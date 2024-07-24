import RequestState from './RequestState.js';

export const schedulingTiles = (tile, tiles3D) => {
	determineFrustumSet(tile, tiles3D);
	markUsedSetLeaves(tile, tiles3D);
	skipTraversal(tile, tiles3D);
	toggleTiles(tile, tiles3D);

	tiles3D.$tilesLoader.lruCache.scheduleUnload();
};

const determineFrustumSet = (tile, tiles3D) => {
	const stats = tiles3D.stats;
	const frameCount = tiles3D.frameCount;
	const errorTarget = tiles3D.errorTarget;
	const maxDepth = tiles3D.maxDepth;
	const loadSiblings = tiles3D.loadSiblings;
	const lruCache = tiles3D.$tilesLoader.lruCache;
	const stopAtEmptyTiles = tiles3D.stopAtEmptyTiles;

	_resetFrameState(tile, frameCount);

	// Early out if this tile is not within view.
	const inFrustum = _tileInView(tile, tiles3D);
	if (inFrustum === false) {
		return false;
	}

	tile.__used = true;
	lruCache.markUsed(tile);

	tile.__inFrustum = true;
	stats.inFrustum++;

	// Early out if this tile has less error than we're targeting but don't stop
	// at an external tile set.
	if ((stopAtEmptyTiles || !tile.__contentEmpty) && !tile.__externalTileSet) {
		// compute the _error and __distanceFromCamera fields
		_calculateError(tile, tiles3D);

		const error = tile.__error;
		if (error <= errorTarget) {
			return true;
		}

		// Early out if we've reached the maximum allowed depth.
		if (maxDepth > 0 && tile.__depth + 1 >= maxDepth) {
			return true;
		}
	}

	// Traverse children and see if any children are in view.
	let anyChildrenUsed = false;
	const children = tile.children;
	for (let i = 0, l = children.length; i < l; i++) {
		const c = children[i];
		const r = determineFrustumSet(c, tiles3D);
		anyChildrenUsed = anyChildrenUsed || r;
	}

	// If there are children within view and we are loading siblings then mark
	// all sibling tiles as used, as well.
	if (anyChildrenUsed && loadSiblings) {
		for (let i = 0, l = children.length; i < l; i++) {
			const c = children[i];
			_recursivelyMarkUsed(c, frameCount, lruCache);
		}
	}

	return true;
};

// Traverse and mark the tiles that are at the leaf nodes of the "used" tree.
const markUsedSetLeaves = (tile, tiles3D) => {
	const stats = tiles3D.stats;
	const frameCount = tiles3D.frameCount;

	if (!_isUsedThisFrame(tile, frameCount)) {
		return;
	}

	stats.used++;

	// This tile is a leaf if none of the children had been used.
	const children = tile.children;
	let anyChildrenUsed = false;
	for (let i = 0, l = children.length; i < l; i++) {
		const c = children[i];
		anyChildrenUsed = anyChildrenUsed || _isUsedThisFrame(c, frameCount);
	}

	if (!anyChildrenUsed) {
		// TODO: This isn't necessarily right because it's possible that a parent tile is considered in the
		// frustum while the child tiles are not, making them unused. If all children have loaded and were properly
		// considered to be in the used set then we shouldn't set ourselves to a leaf here.
		tile.__isLeaf = true;
	} else {
		let childrenWereVisible = false;
		let allChildrenLoaded = true;
		for (let i = 0, l = children.length; i < l; i++) {
			const c = children[i];
			markUsedSetLeaves(c, tiles3D);
			childrenWereVisible = childrenWereVisible || c.__wasSetVisible || c.__childrenWereVisible;

			if (_isUsedThisFrame(c, frameCount)) {
				const childLoaded =
                    c.__allChildrenLoaded ||
                    (!c.__contentEmpty && _isDownloadFinished(c.__loadingState)) ||
                    (c.__externalTileSet && c.__loadingState === RequestState.FAILED);
				allChildrenLoaded = allChildrenLoaded && childLoaded;
			}
		}
		tile.__childrenWereVisible = childrenWereVisible;
		tile.__allChildrenLoaded = allChildrenLoaded;
	}
};

// Skip past tiles we consider unrenderable because they are outside the error threshold.
const skipTraversal = (tile, tiles3D) => {
	const stats = tiles3D.stats;
	const frameCount = tiles3D.frameCount;

	if (!_isUsedThisFrame(tile, frameCount)) {
		return;
	}

	const parent = tile.parent;
	const parentDepthToParent = parent ? parent.__depthFromRenderedParent : -1;
	tile.__depthFromRenderedParent = parentDepthToParent;

	// Request the tile contents or mark it as visible if we've found a leaf.
	const lruCache = tiles3D.$tilesLoader.lruCache;
	if (tile.__isLeaf) {
		tile.__depthFromRenderedParent++;

		if (tile.__loadingState === RequestState.LOADED) {
			if (tile.__inFrustum) {
				tile.__visible = true;
				stats.visible++;
			}
			tile.__active = true;
			stats.active++;
		} else if (!lruCache.isFull() && (!tile.__contentEmpty || tile.__externalTileSet)) {
			tiles3D.$tilesLoader.requestTileContents(tile, tiles3D);
		}

		return;
	}

	const errorRequirement = (tiles3D.errorTarget + 1) * tiles3D.errorThreshold;
	const meetsSSE = tile.__error <= errorRequirement;
	const includeTile = meetsSSE || tile.refine === 'ADD';
	const hasModel = !tile.__contentEmpty;
	const hasContent = hasModel || tile.__externalTileSet;
	const loadedContent = _isDownloadFinished(tile.__loadingState) && hasContent;
	const childrenWereVisible = tile.__childrenWereVisible;
	const children = tile.children;
	const allChildrenHaveContent = tile.__allChildrenLoaded;

	// Increment the relative depth of the node to the nearest rendered parent if it has content
	// and is being rendered.
	if (includeTile && hasModel) {
		tile.__depthFromRenderedParent++;
	}

	// If we've met the SSE requirements and we can load content then fire a fetch.
	if (includeTile && !loadedContent && !lruCache.isFull() && hasContent) {
		tiles3D.$tilesLoader.requestTileContents(tile, tiles3D);
	}

	// Only mark this tile as visible if it meets the screen space error requirements, has loaded content, not
	// all children have loaded yet, and if no children were visible last frame. We want to keep children visible
	// that _were_ visible to avoid a pop in level of detail as the camera moves around and parent / sibling tiles
	// load in.

	// Skip the tile entirely if there's no content to load
	if (
		(meetsSSE && !allChildrenHaveContent && !childrenWereVisible && loadedContent)
        || (tile.refine === 'ADD' && loadedContent)
	) {
		if (tile.__inFrustum) {
			tile.__visible = true;
			stats.visible++;
		}
		tile.__active = true;
		stats.active++;
	}

	// If we're additive then don't stop the traversal here because it doesn't matter whether the children load in
	// at the same rate.
	if (tile.refine !== 'ADD' && meetsSSE && !allChildrenHaveContent && loadedContent) {
		// load the child content if we've found that we've been loaded so we can move down to the next tile
		// layer when the data has loaded.
		for (let i = 0, l = children.length; i < l; i++) {
			const c = children[i];
			if (_isUsedThisFrame(c, frameCount) && !lruCache.isFull()) {
				c.__depthFromRenderedParent = tile.__depthFromRenderedParent + 1;
				_recursivelyLoadTiles(c, c.__depthFromRenderedParent, tiles3D);
			}
		}
	} else {
		for (let i = 0, l = children.length; i < l; i++) {
			const c = children[i];
			if (_isUsedThisFrame(c, frameCount)) {
				skipTraversal(c, tiles3D);
			}
		}
	}
};

// Final traverse to toggle tile visibility.
const toggleTiles = (tile, tiles3D) => {
	const frameCount = tiles3D.frameCount;

	const isUsed = _isUsedThisFrame(tile, frameCount);

	if (isUsed || tile.__usedLastFrame) {
		let setActive = false;
		let setVisible = false;
		if (isUsed) {
			// enable visibility if active due to shadows
			setActive = tile.__active;
			if (tiles3D.displayActiveTiles) {
				setVisible = tile.__active || tile.__visible;
			} else {
				setVisible = tile.__visible;
			}
		}

		// If the active or visible state changed then call the functions.
		if (!tile.__contentEmpty && tile.__loadingState === RequestState.LOADED) {
			if (tile.__wasSetActive !== setActive) {
				tiles3D.$setTileActive(tile, setActive);
			}

			if (tile.__wasSetVisible !== setVisible) {
				tiles3D.$setTileVisible(tile, setVisible);
			}
		}
		tile.__wasSetActive = setActive;
		tile.__wasSetVisible = setVisible;
		tile.__usedLastFrame = isUsed;

		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			const c = children[i];
			toggleTiles(c, tiles3D);
		}
	}
};

// Resets the frame frame information for the given tile
const _resetFrameState = (tile, frameCount) => {
	if (tile.__lastFrameVisited !== frameCount) {
		tile.__lastFrameVisited = frameCount;
		tile.__used = false;
		tile.__inFrustum = false;
		tile.__isLeaf = false;
		tile.__visible = false;
		tile.__active = false;
		tile.__error = Infinity;
		tile.__distanceFromCamera = Infinity;
		tile.__childrenWereVisible = false;
		tile.__allChildrenLoaded = false;
	}
};

const _tileInView = (tile, tiles3D) => {
	const cameraInfos = tiles3D.$cameras.getInfos();

	const cached = tile.cached;
	const boundingVolume = cached.boundingVolume;
	const inFrustum = cached.inFrustum;

	let inView = false;

	for (let i = 0, l = cameraInfos.length; i < l; i++) {
		// Track which camera frustums this tile is in so we can use it
		// to ignore the error calculations for cameras that can't see it
		const frustum = cameraInfos[i].frustum;

		if (boundingVolume.intersectsFrustum(frustum)) {
			inView = true;
			inFrustum[i] = true;
		} else {
			inFrustum[i] = false;
		}
	}

	return inView;
};

const _calculateError = (tile, tiles3D) => {
	const cameraInfos = tiles3D.$cameras.getInfos();

	const cached = tile.cached;
	const inFrustum = cached.inFrustum;
	const boundingVolume = cached.boundingVolume;

	let maxError = -Infinity;
	let minDistance = Infinity;

	for (let i = 0, l = cameraInfos.length; i < l; i++) {
		if (!inFrustum[i]) {
			continue;
		}

		// transform camera position into local frame of the tile bounding box
		const info = cameraInfos[i];
		const invScale = info.invScale;

		let error;
		if (info.isOrthographic) {
			const pixelSize = info.pixelSize;
			error = cached.geometricError / (pixelSize * invScale);
		} else {
			const distance = boundingVolume.distanceToPoint(info.position);
			const scaledDistance = distance * invScale;
			const sseDenominator = info.sseDenominator;
			error = cached.geometricError / (scaledDistance * sseDenominator);

			minDistance = Math.min(minDistance, scaledDistance);
		}

		maxError = Math.max(maxError, error);
	}

	tile.__distanceFromCamera = minDistance;
	tile.__error = maxError;
};

// Recursively mark tiles used down to the next tile with content
const _recursivelyMarkUsed = (tile, frameCount, lruCache) => {
	_resetFrameState(tile, frameCount);

	tile.__used = true;
	lruCache.markUsed(tile);
	if (tile.__contentEmpty) {
		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			_recursivelyMarkUsed(children[i], frameCount, lruCache);
		}
	}
};

// Checks whether this tile was last used on the given frame.
const _isUsedThisFrame = (tile, frameCount) => {
	return tile.__lastFrameVisited === frameCount && tile.__used;
};

function _isDownloadFinished(value) {
	return value === RequestState.LOADED || value === RequestState.FAILED;
}

const _recursivelyLoadTiles = (tile, depthFromRenderedParent, tiles3D) => {
	// Try to load any external tile set children if the external tile set has loaded.
	const doTraverse =
		tile.__contentEmpty && (
			!tile.__externalTileSet ||
			_isDownloadFinished(tile.__loadingState)
		);
	if (doTraverse) {
		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			// don't increment depth to rendered parent here because we should treat
			// the next layer of rendered children as just a single depth away for the
			// sake of sorting.
			const child = children[i];
			child.__depthFromRenderedParent = depthFromRenderedParent;
			_recursivelyLoadTiles(child, depthFromRenderedParent, tiles3D);
		}
	} else {
		tiles3D.$tilesLoader.requestTileContents(tile, tiles3D);
	}
};