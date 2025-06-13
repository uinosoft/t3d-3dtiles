/**
 * Returns the file extension of the path component of a URL
 * @param {string} url
 * @returns {string} null if no extension found
 */
export const getUrlExtension = url => {
	if (!url) {
		return null;
	}

	const filename = url
		.replace(/[a-z]+:\/\/[^/]+/i, '') 	// remove origin
		.replace(/\?.*$/i, '') 				// remove query
		.replace(/.*\//g, ''); 				// remove path

	const lastPeriod = filename.lastIndexOf('.');
	if (lastPeriod === -1) {
		return null;
	}

	return filename.substring(lastPeriod + 1) || null;
};

export const traverseSet = (tile, beforeCb = null, afterCb = null, parent = null, depth = 0) => {
	if (beforeCb && beforeCb(tile, parent, depth)) {
		if (afterCb) {
			afterCb(tile, parent, depth);
		}

		return;
	}

	const children = tile.children;
	for (let i = 0, l = children.length; i < l; i++) {
		traverseSet(children[i], beforeCb, afterCb, tile, depth + 1);
	}

	if (afterCb) {
		afterCb(tile, parent, depth);
	}
};

/**
 * Traverses the ancestry of the tile up to the root tile.
 */
export function traverseAncestors(tile, callback = null) {
	let current = tile;

	while (current) {
		const depth = current.__depth;
		const parent = current.parent;

		if (callback) {
			callback(current, parent, depth);
		}

		current = parent;
	}
}