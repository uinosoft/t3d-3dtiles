/**
 * Returns the file extension of the path component of a URL
 * @param {string} url
 * @returns {string} null if no extension found
 */
export const getUrlExtension = url => {
	let parsedUrl;
	try {
		parsedUrl = new URL(url, 'http://fakehost.com/');
	} catch (_) {
		// Ignore invalid URLs
		return null;
	}

	const filename = parsedUrl.pathname.split('/').pop();
	const dotIndex = filename.lastIndexOf('.');
	if (dotIndex === -1 || dotIndex === filename.length - 1) {
		// Has no extension or has trailing . character
		return null;
	}

	const extension = filename.substring(dotIndex + 1);
	return extension;
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