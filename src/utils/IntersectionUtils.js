import { Matrix4, Ray, Vector3 } from 't3d';

const _localRay = new Ray();
const _vec = new Vector3();
const _hitArray = [];
const _mat = new Matrix4();

export const distanceSort = (a, b) => {
	return a.distance - b.distance;
};

const intersectTileScene = (tile, ray, renderer, intersects) => {
	const { scene } = tile.cached;

	const lengthBefore = intersects.length;
	const didRaycast = renderer.invokeOnePlugin(plugin => plugin.raycastTile && plugin.raycastTile(tile, scene, ray, intersects));
	if (!didRaycast) {
		scene.traverse(c => {
			// We set the default raycast function to empty so t3d.js doesn't automatically cast against it
			Object.getPrototypeOf(c).raycast.call(c, ray, intersects);
		});

		const lengthAfter = intersects.length;

		// add the tile to intersects
		if (lengthAfter > lengthBefore) {
			for (let i = lengthBefore; i < lengthAfter; i++) {
				intersects[i].tile = tile;
			}
		}
	}
};

function intersectTileSceneFirstHit(tile, ray, renderer) {
	intersectTileScene(tile, ray, renderer, _hitArray);
	_hitArray.sort(distanceSort);

	const hit = _hitArray[0] || null;
	_hitArray.length = 0;
	return hit;
}

function isTileInitialized(tile) {
	return '__used' in tile;
}

// Returns the closest hit when traversing the tree
export const raycastTraverseFirstHit = (tile, renderer, ray, localRay = null) => {
	const { activeTiles } = renderer;

	// get the ray in the local group frame
	if (localRay === null) {
		localRay = _localRay;
		localRay.copy(ray).applyMatrix4(_mat.copy(renderer.worldMatrix).inverse());
	}

	// get a set of intersections so we intersect the nearest one first
	const array = [];
	const children = tile.children;
	for (let i = 0, l = children.length; i < l; i++) {
		const child = children[i];

		if (!isTileInitialized(child) || !child.__used) {
			continue;
		}

		// track the tile and hit distance for sorting
		const boundingVolume = child.cached.boundingVolume;
		if (boundingVolume.intersectRay(localRay, _vec) !== null) {
			_vec.applyMatrix4(renderer.worldMatrix);
			array.push({
				distance: _vec.distanceToSquared(ray.origin),
				tile: child
			});
		}
	}

	// sort them by ascending distance
	array.sort(distanceSort);

	// If the root is active make sure we've checked it
	let bestHit = null;
	let bestHitDistSq = Infinity;
	if (activeTiles.has(tile)) {
		const hit = intersectTileSceneFirstHit(tile, ray, renderer);
		if (hit) {
			bestHit = hit;
			bestHitDistSq = hit.distance * hit.distance;
		}
	}

	// traverse until we find the best hit and early out if a tile bounds
	// couldn't possible include a best hit
	for (let i = 0, l = array.length; i < l; i++) {
		const data = array[i];
		const boundingVolumeDistSq = data.distance;
		const tile = data.tile;
		if (boundingVolumeDistSq > bestHitDistSq) {
			break;
		}

		const hit = raycastTraverseFirstHit(tile, renderer, ray, localRay);
		if (hit) {
			const hitDistSq = hit.distance * hit.distance;
			if (hitDistSq < bestHitDistSq) {
				bestHit = hit;
				bestHitDistSq = hitDistSq;
			}
		}
	}

	return bestHit;
};

export const raycastTraverse = (tile, renderer, ray, intersects, localRay = null) => {
	// if the tile has not been asynchronously initialized then there's no point in
	// traversing the tiles to check intersections.
	if (!isTileInitialized(tile)) {
		return;
	}

	const { activeTiles } = renderer;
	const { boundingVolume } = tile.cached;

	// get the ray in the local group frame
	if (localRay === null) {
		localRay = _localRay;
		localRay.copy(ray).applyMatrix4(_mat.copy(renderer.worldMatrix).inverse());
	}

	// exit early if the tile isn't used or the bounding volume is not intersected
	if (!tile.__used || !boundingVolume.intersectsRay(localRay)) {
		return;
	}

	// only intersect the tile geometry if it's active
	if (activeTiles.has(tile)) {
		intersectTileScene(tile, ray, renderer, intersects);
	}

	const children = tile.children;
	for (let i = 0, l = children.length; i < l; i++) {
		raycastTraverse(children[i], renderer, ray, intersects, localRay);
	}
};