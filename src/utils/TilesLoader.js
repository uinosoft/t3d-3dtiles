import { Matrix4, Vector3 } from 't3d';
import { TileBoundingVolume } from '../math/TileBoundingVolume.js';
import RequestState from './RequestState.js';
import { LRUCache } from './LRUCache.js';
import { PriorityQueue } from './PriorityQueue.js';
import { traverseSet, getUrlExtension } from './Utils.js';

export class TilesLoader {

	constructor() {
		this.lruCache = new LRUCache({ unloadPriorityCallback: lruPriorityCallback });
		this.downloadQueue = new PriorityQueue({ maxJobs: 4, priorityCallback });
		this.parseQueue = new PriorityQueue({ maxJobs: 1, priorityCallback });
	}

	preprocessTileSet(json, url, parent = null) {
		const version = json.asset.version;
		const [major, minor] = version.split('.').map(v => parseInt(v));
		console.assert(
			major <= 1,
			'TilesLoader: asset.version is expected to be a 1.x or a compatible version.'
		);

		if (major === 1 && minor > 0) {
			console.warn('TilesLoader: tiles versions at 1.1 or higher have limited support. Some new extensions and features may not be supported.');
		}

		// remove the last file path path-segment from the URL including the trailing slash
		let basePath = url.replace(/\/[^/]*$/, '');
		basePath = new URL(basePath, window.location.href).toString();

		// this.preprocessNode(json.root, basePath, parent);
		traverseSet(
			json.root,
			(node, parent) => preprocessTile(node, parent, basePath),
			null,
			parent,
			parent ? parent.__depth : 0
		);
	}

	requestTileContents(tile, tiles3D) {
		// If the tile is already being loaded then don't
		// start it again.
		if (tile.__loadingState !== RequestState.UNLOADED) {
			return;
		}

		const stats = tiles3D.stats;
		const lruCache = this.lruCache;
		const downloadQueue = this.downloadQueue;
		const parseQueue = this.parseQueue;

		const isExternalTileSet = tile.__externalTileSet;

		lruCache.add(tile, t => {
			if (t.__loadingState === RequestState.LOADING) {
				// Stop the load if it's started
				t.__loadAbort.abort();
				t.__loadAbort = null;
			} else if (isExternalTileSet) {
				t.children.length = 0;
			} else {
				tiles3D.$disposeTile(t);
			}

			// Decrement stats
			if (t.__loadingState === RequestState.LOADING) {
				stats.downloading--;
			} else if (t.__loadingState === RequestState.PARSING) {
				stats.parsing--;
			}

			t.__loadingState = RequestState.UNLOADED;
			t.__loadIndex++;

			downloadQueue.remove(t);
			parseQueue.remove(t);
		});

		// Track a new load index so we avoid the condition where this load is stopped and
		// another begins soon after so we don't parse twice.
		tile.__loadIndex++;
		const loadIndex = tile.__loadIndex;
		const controller = new AbortController();
		const signal = controller.signal;

		stats.downloading++;
		tile.__loadAbort = controller;
		tile.__loadingState = RequestState.LOADING;

		const errorCallback = e => {
			// if it has been unloaded then the tile has been disposed
			if (tile.__loadIndex !== loadIndex) {
				return;
			}

			if (e.name !== 'AbortError') {
				downloadQueue.remove(tile);
				parseQueue.remove(tile);

				if (tile.__loadingState === RequestState.PARSING) {
					stats.parsing--;
				} else if (tile.__loadingState === RequestState.LOADING) {
					stats.downloading--;
				}

				stats.failed++;
				tile.__loadingState = RequestState.FAILED;

				// Handle fetch bug for switching examples in index.html.
				// https://stackoverflow.com/questions/12009423/what-does-status-canceled-for-a-resource-mean-in-chrome-developer-tools
				if (e.message !== 'Failed to fetch') {
					console.error(`TilesLoader: Failed to load tile at url "${tile.content.uri}".`);
					console.error(e);
				}
			} else {
				lruCache.remove(tile);
			}
		};

		let uri = tile.content.uri;
		tiles3D.invokeAllPlugins(plugin => uri = plugin.preprocessURL ? plugin.preprocessURL(uri, tile) : uri);

		if (isExternalTileSet) {
			downloadQueue.add(tile, tileCb => {
				// if it has been unloaded then the tile has been disposed
				if (tileCb.__loadIndex !== loadIndex) {
					return Promise.resolve();
				}

				return tiles3D.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(uri, { ...tiles3D.fetchOptions, signal }));
			}).then(res => {
				if (res.ok) {
					return res.json();
				} else {
					throw new Error(`TilesLoader: Failed to load tileset "${uri}" with status ${res.status} : ${res.statusText}`);
				}
			}).then(json => {
				this.preprocessTileSet(json, uri, tile);
				return json;
			}).then(json => {
				// if it has been unloaded then the tile has been disposed
				if (tile.__loadIndex !== loadIndex) {
					return;
				}

				stats.downloading--;
				tile.__loadAbort = null;
				tile.__loadingState = RequestState.LOADED;

				tile.children.push(json.root);

				tiles3D.dispatchEvent({
					type: 'load-tile-set',
					tileSet: json,
					url: uri
				});
			}).catch(errorCallback);
		} else {
			downloadQueue.add(tile, tileCb => {
				if (tileCb.__loadIndex !== loadIndex) {
					return Promise.resolve();
				}

				return tiles3D.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(uri, { ...tiles3D.fetchOptions, signal }));
			}).then(res => {
				if (tile.__loadIndex !== loadIndex) {
					return;
				}

				if (res.ok) {
					const extension = getUrlExtension(res.url);
					if (extension === 'gltf') {
						return res.json();
					} else {
						return res.arrayBuffer();
					}
				} else {
					throw new Error(`Failed to load model with error code ${res.status}`);
				}
			}).then(buffer => {
				// if it has been unloaded then the tile has been disposed
				if (tile.__loadIndex !== loadIndex) {
					return;
				}

				stats.downloading--;
				stats.parsing++;
				tile.__loadAbort = null;
				tile.__loadingState = RequestState.PARSING;

				return parseQueue.add(tile, parseTile => {
					// if it has been unloaded then the tile has been disposed
					if (parseTile.__loadIndex !== loadIndex) {
						return Promise.resolve();
					}

					const uri = parseTile.content.uri;
					const extension = getUrlExtension(uri);

					return tiles3D.$parseTile(buffer, parseTile, extension);
				});
			}).then(() => {
				// if it has been unloaded then the tile has been disposed
				if (tile.__loadIndex !== loadIndex) {
					return;
				}

				stats.parsing--;
				tile.__loadingState = RequestState.LOADED;

				if (tile.__wasSetVisible) {
					tiles3D.$setTileVisible(tile, true);
				}

				if (tile.__wasSetActive) {
					tiles3D.$setTileActive(tile, true);
				}
			}).catch(errorCallback);
		}
	}

}

/**
 * Function for sorting the evicted LRU items. We should evict the shallowest depth first.
 * @param {Tile} tile
 * @returns number
 */
const lruPriorityCallback = (map, tile) => 1 / (tile.__depthFromRenderedParent + 1);

/**
 * Function for provided to sort all tiles for prioritizing loading.
 *
 * @param {Tile} a
 * @param {Tile} b
 * @returns number
 */
const priorityCallback = (a, b) => {
	if (a.__depth !== b.__depth) {
		// load shallower tiles first
		return a.__depth > b.__depth ? -1 : 1;
	} else if (a.__inFrustum !== b.__inFrustum) {
		// load tiles that are in the frustum at the current depth
		return a.__inFrustum ? 1 : -1;
	} else if (a.__used !== b.__used) {
		// load tiles that have been used
		return a.__used ? 1 : -1;
	} else if (a.__error !== b.__error) {
		// load the tile with the higher error
		return a.__error > b.__error ? 1 : -1;
	} else if (a.__distanceFromCamera !== b.__distanceFromCamera) {
		// and finally visible tiles which have equal error (ex: if geometricError === 0)
		// should prioritize based on distance.
		return a.__distanceFromCamera > b.__distanceFromCamera ? -1 : 1;
	}

	return 0;
};

const preprocessTile = (tile, parentTile, tileSetDir) => {
	if (tile.contents) {
		// TODO: multiple contents (1.1) are not supported yet
		tile.content = tile.contents[0];
	}

	if (tile.content) {
		// Fix old file formats
		if (!('uri' in tile.content) && 'url' in tile.content) {
			tile.content.uri = tile.content.url;
			delete tile.content.url;
		}

		if (tile.content.uri) {
			// tile content uri has to be interpreted relative to the tileset.json
			// tile.content.uri = new URL( tile.content.uri, tileSetDir + '/' ).toString();
			tile.content.uri = `${tileSetDir}/${tile.content.uri}`;
		}

		// NOTE: fix for some cases where tile provide the bounding volume
		// but volumes are not present.
		if (
			tile.content.boundingVolume &&
			!(
				'box' in tile.content.boundingVolume ||
				'sphere' in tile.content.boundingVolume ||
				'region' in tile.content.boundingVolume
			)
		) {
			delete tile.content.boundingVolume;
		}
	}

	tile.parent = parentTile;
	tile.children = tile.children || [];

	const uri = tile.content && tile.content.uri;
	if (uri) {
		// "content" should only indicate loadable meshes, not external tile sets
		const extension = getUrlExtension(tile.content.uri);
		const isExternalTileSet = Boolean(extension && extension.toLowerCase() === 'json');
		tile.__externalTileSet = isExternalTileSet;
		tile.__contentEmpty = isExternalTileSet;
	} else {
		tile.__externalTileSet = false;
		tile.__contentEmpty = true;
	}

	// Expected to be set during calculateError()
	tile.__distanceFromCamera = Infinity;
	tile.__error = Infinity;

	tile.__inFrustum = false;
	tile.__isLeaf = false;

	tile.__usedLastFrame = false;
	tile.__used = false;

	tile.__wasSetVisible = false;
	tile.__visible = false;
	tile.__childrenWereVisible = false;
	tile.__allChildrenLoaded = false;

	tile.__wasSetActive = false;
	tile.__active = false;

	tile.__loadingState = RequestState.UNLOADED;
	tile.__loadIndex = 0;

	tile.__loadAbort = null;

	tile.__depthFromRenderedParent = -1;
	if (parentTile === null) {
		tile.__depth = 0;
		tile.refine = tile.refine || 'REPLACE';
	} else {
		tile.__depth = parentTile.__depth + 1;
		tile.refine = tile.refine || parentTile.refine;
	}

	//

	const transform = new Matrix4();
	if (tile.transform) {
		transform.fromArray(tile.transform);
	}
	if (parentTile) {
		transform.premultiply(parentTile.cached.transform);
	}
	const transformInverse = (new Matrix4()).copy(transform).inverse();

	const transformScale = _vec3_1.setFromMatrixScale(transform);
	const uniformScale = Math.max(transformScale.x, transformScale.y, transformScale.z);
	let geometricError = tile.geometricError * uniformScale;

	const boundingVolume = new TileBoundingVolume();
	if ('box' in tile.boundingVolume) {
		boundingVolume.setOBBData(tile.boundingVolume.box, transform);
	}
	if ('sphere' in tile.boundingVolume) {
		boundingVolume.setSphereData(tile.boundingVolume.sphere, transform);
	}
	if ('region' in tile.boundingVolume) {
		boundingVolume.setRegionData(...tile.boundingVolume.region);
		geometricError = tile.geometricError;
	}

	tile.cached = {
		loadIndex: 0,
		transform,
		transformInverse,

		geometricError, // geometric error applied tile transform scale

		boundingVolume,

		active: false,
		inFrustum: [],

		scene: null,
		geometry: null,
		material: null,

		featureTable: null,
		batchTable: null
	};
};

const _vec3_1 = new Vector3();