import RequestState from './RequestState.js';
import { LRUCache } from './LRUCache.js';
import { PriorityQueue } from './PriorityQueue.js';
import { getUrlExtension } from './Utils.js';

export class TilesLoader {

	constructor() {
		this.lruCache = new LRUCache({ unloadPriorityCallback: lruPriorityCallback });
		this.downloadQueue = new PriorityQueue({ maxJobs: 4, priorityCallback });
		this.parseQueue = new PriorityQueue({ maxJobs: 1, priorityCallback });
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

		const isExternalTileSet = tile.__hasUnrenderableContent;

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

		let uri = new URL(tile.content.uri, tile.__basePath + '/').toString();
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
				tiles3D.preprocessTileSet(json, uri, tile);
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
					tiles3D.invokeOnePlugin(plugin => plugin.setTileVisible && plugin.setTileVisible(tile, true));
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