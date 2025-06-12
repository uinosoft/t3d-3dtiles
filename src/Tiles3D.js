import { EventDispatcher, LoadingManager, Matrix4, Object3D } from 't3d';
import { TileBoundingVolume } from './math/TileBoundingVolume.js';
import { traverseSet, getUrlExtension } from './utils/Utils.js';
import { raycastTraverse, raycastTraverseFirstHit, distanceSort } from './utils/IntersectionUtils.js';
import { CameraList } from './utils/CameraList.js';
import { LRUCache } from './utils/LRUCache.js';
import { PriorityQueue } from './utils/PriorityQueue.js';
import { ModelLoader } from './utils/ModelLoader.js';
import { markUsedTiles, markUsedSetLeaves, markVisibleTiles, toggleTiles } from './utils/TilesScheduler.js';
import { UNLOADED, LOADING, PARSING, LOADED, FAILED } from './constants.js';
import { WGS84_ELLIPSOID } from './math/GeoConstants.js';

const _updateBeforeEvent = { type: 'update-before' };
const _updateAfterEvent = { type: 'update-after' };

const PLUGIN_REGISTERED = Symbol('PLUGIN_REGISTERED');

const INITIAL_FRUSTUM_CULLED = Symbol('INITIAL_FRUSTUM_CULLED');

const tempMat = new Matrix4();
const viewErrorTarget = {
	inView: false,
	error: Infinity
};

const matrixEquals = (matrixA, matrixB, epsilon = Number.EPSILON) => {
	const te = matrixA.elements;
	const me = matrixB.elements;

	for (let i = 0; i < 16; i++) {
		if (Math.abs(te[i] - me[i]) > epsilon) {
			return false;
		}
	}

	return true;
};

// priority queue sort function that takes two tiles to compare. Returning 1 means
// "tile a" is loaded first.
const priorityCallback = (a, b) => {
	if (a.__depthFromRenderedParent !== b.__depthFromRenderedParent) {
		// load shallower tiles first using "depth from rendered parent" to help
		// even out depth disparities caused by non-content parent tiles
		return a.__depthFromRenderedParent > b.__depthFromRenderedParent ? -1 : 1;
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

// lru cache unload callback that takes two tiles to compare. Returning 1 means "tile a"
// is unloaded first.
const lruPriorityCallback = (a, b) => {
	if (a.__depthFromRenderedParent !== b.__depthFromRenderedParent) {
		// dispose of deeper tiles first
		return a.__depthFromRenderedParent > b.__depthFromRenderedParent ? 1 : -1;
	} else if (a.__loadingState !== b.__loadingState) {
		// dispose of tiles that are earlier along in the loading process first
		return a.__loadingState > b.__loadingState ? -1 : 1;
	} else if (a.__lastFrameVisited !== b.__lastFrameVisited) {
		// dispose of least recent tiles first
		return a.__lastFrameVisited > b.__lastFrameVisited ? -1 : 1;
	} else if (a.__hasUnrenderableContent !== b.__hasUnrenderableContent) {
		// dispose of external tile sets last
		return a.__hasUnrenderableContent ? -1 : 1;
	} else if (a.__error !== b.__error) {
		// unload the tile with lower error
		return a.__error > b.__error ? -1 : 1;
	}

	return 0;
};

export class Tiles3D extends Object3D {

	get root() {
		const rootTileSet = this.rootTileSet;
		return rootTileSet ? rootTileSet.root : null;
	}

	constructor(url, manager = new LoadingManager()) {
		super();

		this.ellipsoid = WGS84_ELLIPSOID.clone();

		// options

		this.fetchOptions = {};
		this.errorTarget = 6.0;
		this.errorThreshold = Infinity;
		this.loadSiblings = true;
		this.displayActiveTiles = false;
		this.maxDepth = Infinity;
		this.stopAtEmptyTiles = true;

		this.preprocessURL = null;
		manager.setURLModifier(url => {
			if (this.preprocessURL) {
				return this.preprocessURL(url);
			} else {
				return url;
			}
		});
		this.manager = manager;

		// stats

		this.stats = {
			parsing: 0,
			downloading: 0,
			failed: 0,
			inFrustum: 0,
			used: 0,
			active: 0,
			visible: 0
		};

		this.frameCount = 0;

		this.activeTiles = new Set();
		this.visibleTiles = new Set();
		this.usedSet = new Set();

		this.rootLoadingState = UNLOADED;

		// internals

		this.rootURL = url;
		this.rootTileSet = null;

		this._autoDisableRendererCulling = true;

		this.plugins = [];

		const lruCache = new LRUCache();
		lruCache.unloadPriorityCallback = lruPriorityCallback;

		const downloadQueue = new PriorityQueue();
		downloadQueue.maxJobs = 10;
		downloadQueue.priorityCallback = priorityCallback;

		const parseQueue = new PriorityQueue();
		parseQueue.maxJobs = 1;
		parseQueue.priorityCallback = priorityCallback;

		this.lruCache = lruCache;
		this.downloadQueue = downloadQueue;
		this.parseQueue = parseQueue;

		this.$cameras = new CameraList();
		this.$modelLoader = new ModelLoader(manager);
		this.$events = new EventDispatcher();

		this.lruCache.computeMemoryUsageCallback = tile => tile.cached.bytesUsed ?? null;
	}

	// Plugins
	registerPlugin(plugin) {
		if (plugin[PLUGIN_REGISTERED] === true) {
			throw new Error('Tiles3D: A plugin can only be registered to a single tile set');
		}

		// insert the plugin based on the priority registered on the plugin
		const plugins = this.plugins;
		const priority = plugin.priority || 0;
		let insertionPoint = plugins.length;
		for (let i = 0; i < plugins.length; i++) {
			const otherPriority = plugins[i].priority || 0;
			if (otherPriority > priority) {
				insertionPoint = i;
				break;
			}
		}

		plugins.splice(insertionPoint, 0, plugin);
		plugin[PLUGIN_REGISTERED] = true;
		if (plugin.init) {
			plugin.init(this);
		}
	}

	unregisterPlugin(plugin) {
		const plugins = this.plugins;
		if (typeof plugin === 'string') {
			plugin = this.getPluginByName(name);
		}

		if (plugins.includes(plugin)) {
			const index = plugins.indexOf(plugin);
			plugins.splice(index, 1);
			if (plugin.dispose) {
				plugin.dispose();
			}

			return true;
		}

		return false;
	}

	getPluginByName(name) {
		return this.plugins.find(p => p.name === name) || null;
	}

	traverse(beforecb, aftercb, ensureFullyProcessed = true) {
		if (!this.root) return;

		traverseSet(this.root, (tile, ...args) => {
			if (ensureFullyProcessed) {
				this.ensureChildrenArePreprocessed(tile, true);
			}

			return beforecb ? beforecb(tile, ...args) : false;
		}, aftercb);
	}

	markTileUsed(tile) {
		// save the tile in a separate "used set" so we can mark it as unused
		// before the next tile set traversal
		this.usedSet.add(tile);
		this.lruCache.markUsed(tile);
	}

	// Public API
	update() {
		const { lruCache, usedSet, stats, root } = this;

		if (this.rootLoadingState === UNLOADED) {
			this.rootLoadingState = LOADING;
			this.invokeOnePlugin(plugin => plugin.loadRootTileSet && plugin.loadRootTileSet())
				.then(root => {
					let processedUrl = this.rootURL;
					if (processedUrl !== null) {
						this.invokeAllPlugins(plugin => processedUrl = plugin.preprocessURL ? plugin.preprocessURL(processedUrl, null) : processedUrl);
					}

					this.rootLoadingState = LOADED;
					this.rootTileSet = root;

					this.dispatchEvent({
						type: 'load-tile-set',
						tileSet: root,
						url: processedUrl
					});
				})
				.catch(error => {
					this.rootLoadingState = FAILED;
					console.error(error);

					this.rootTileSet = null;
					this.dispatchEvent({
						type: 'load-error',
						tile: null,
						error,
						url: this.rootURL
					});
				});
		}

		if (!root) {
			return;
		}

		this.dispatchEvent(_updateBeforeEvent);

		this.$cameras.updateInfos(this.worldMatrix);

		stats.inFrustum = 0;
		stats.used = 0;
		stats.active = 0;
		stats.visible = 0;
		this.frameCount++;

		usedSet.forEach(tile => lruCache.markUnused(tile));
		usedSet.clear();

		markUsedTiles(root, this);
		markUsedSetLeaves(root, this);
		markVisibleTiles(root, this);
		toggleTiles(root, this);

		this.lruCache.scheduleUnload();

		this.dispatchEvent(_updateAfterEvent);
	}

	resetFailedTiles() {
		// reset the root tile if it's finished but never loaded
		if (this.rootLoadingState === FAILED) {
			this.rootLoadingState = UNLOADED;
		}

		const stats = this.stats;
		if (stats.failed === 0) {
			return;
		}

		this.traverse(tile => {
			if (tile.__loadingState === FAILED) {
				tile.__loadingState = UNLOADED;
			}
		}, null, false);

		stats.failed = 0;
	}

	dispose() {
		// dispose of all the plugins
		const plugins = [...this.plugins];
		plugins.forEach(plugin => {
			this.unregisterPlugin(plugin);
		});

		const lruCache = this.lruCache;

		// Make sure we've collected all children before disposing of the internal tilesets to avoid
		// dangling children that we inadvertantly skip when deleting the nested tileset.
		const toRemove = [];
		this.traverse(t => {
			toRemove.push(t);
			return false;
		}, null, false);
		for (let i = 0, l = toRemove.length; i < l; i++) {
			lruCache.remove(toRemove[i]);
		}

		this.stats = {
			parsing: 0,
			downloading: 0,
			failed: 0,
			inFrustum: 0,
			used: 0,
			active: 0,
			visible: 0
		};
		this.frameCount = 0;
	}

	dispatchEvent(event) {
		this.$events.dispatchEvent(event);
	}

	fetchData(url, options) {
		return fetch(url, options);
	}

	parseTile(buffer, tile, extension) {
		return this.$modelLoader.loadTileContent(buffer, tile, extension, this)
			.then(scene => {
				scene.traverse(c => {
					c[INITIAL_FRUSTUM_CULLED] = c.frustumCulled; // store initial value
					c.frustumCulled = c.frustumCulled && !this._autoDisableRendererCulling;
				});

				this.dispatchEvent({
					type: 'load-model',
					scene,
					tile
				});
			});
	}

	disposeTile(tile) {
		// This could get called before the tile has finished downloading
		const cached = tile.cached;
		if (cached.scene) {
			const materials = cached.materials;
			const geometry = cached.geometry;
			const textures = cached.textures;
			const parent = cached.scene.parent;

			for (let i = 0, l = geometry.length; i < l; i++) {
				geometry[i].dispose();
			}

			for (let i = 0, l = materials.length; i < l; i++) {
				materials[i].dispose();
			}

			for (let i = 0, l = textures.length; i < l; i++) {
				const texture = textures[i];

				if (texture.image instanceof ImageBitmap) {
					texture.image.close();
				}

				texture.dispose();
			}

			if (parent) {
				parent.remove(cached.scene);
			}

			this.dispatchEvent({
				type: 'dispose-model',
				scene: cached.scene,
				tile
			});

			cached.scene = null;
			cached.materials = null;
			cached.textures = null;
			cached.geometry = null;
			cached.metadata = null;
		}

		tile._loadIndex++;
	}

	preprocessNode(tile, tileSetDir, parentTile = null) {
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

		if (tile.content?.uri) {
			// "content" should only indicate loadable meshes, not external tile sets
			const extension = getUrlExtension(tile.content.uri);

			tile.__hasContent = true;
			tile.__hasUnrenderableContent = Boolean(extension && /json$/.test(extension));
			tile.__hasRenderableContent = !tile.__hasUnrenderableContent;
		} else {
			tile.__hasContent = false;
			tile.__hasUnrenderableContent = false;
			tile.__hasRenderableContent = false;
		}

		// tracker for determining if all the children have been asynchronously
		// processed and are ready to be traversed
		tile.__childrenProcessed = 0;
		if (parentTile) {
			parentTile.__childrenProcessed++;
		}

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

		tile.__loadingState = UNLOADED;

		if (parentTile === null) {
			tile.__depth = 0;
			tile.__depthFromRenderedParent = (tile.__hasRenderableContent ? 1 : 0);
			tile.refine = tile.refine || 'REPLACE';
		} else {
			// increment the "depth from parent" when we encounter a new tile with content
			tile.__depth = parentTile.__depth + 1;
			tile.__depthFromRenderedParent = parentTile.__depthFromRenderedParent + (tile.__hasRenderableContent ? 1 : 0);
			tile.refine = tile.refine || parentTile.refine;
		}

		tile.__basePath = tileSetDir;

		tile.__lastFrameVisited = -1;

		tile.__loadIndex = 0; // TODO remove this
		tile.__loadAbort = null; // TODO remove this

		this.invokeAllPlugins(plugin => {
			plugin !== this && plugin.preprocessNode && plugin.preprocessNode(tile, tileSetDir, parentTile);
		});

		// cached

		const transform = new Matrix4();
		if (tile.transform) {
			transform.fromArray(tile.transform);
		}

		if (parentTile) {
			transform.premultiply(parentTile.cached.transform);
		}

		const transformInverse = new Matrix4().copy(transform).inverse();
		const boundingVolume = new TileBoundingVolume();
		if ('sphere' in tile.boundingVolume) {
			boundingVolume.setSphereData(tile.boundingVolume.sphere, transform);
		}
		if ('box' in tile.boundingVolume) {
			boundingVolume.setOBBData(tile.boundingVolume.box, transform);
		}
		if ('region' in tile.boundingVolume) {
			boundingVolume.setRegionData(this.ellipsoid, ...tile.boundingVolume.region);
		}

		tile.cached = {
			transform,
			transformInverse,

			active: false,

			boundingVolume,

			metadata: null,
			scene: null,
			geometry: null,
			materials: null,
			textures: null,

			// TODO remove this

			inFrustum: [],

			featureTable: null,
			batchTable: null
		};
	}

	setTileActive(tile, active) {
		active ? this.activeTiles.add(tile) : this.activeTiles.delete(tile);
	}

	setTileVisible(tile, visible) {
		const scene = tile.cached.scene;

		if (visible) {
			if (scene) {
				this.add(scene);
				scene.updateMatrix(true);
			}
		} else {
			if (scene) {
				this.remove(scene);
			}
		}

		visible ? this.visibleTiles.add(tile) : this.visibleTiles.delete(tile);

		this.dispatchEvent({
			type: 'tile-visibility-change',
			scene,
			tile,
			visible
		});
	}

	calculateTileViewError(tile, target) {
		// retrieve whether the tile is visible, screen space error, and distance to camera
		// set "inView", "error", "distance"

		const cached = tile.cached;
		const cameraInfo = this.$cameras.getInfos();
		const boundingVolume = cached.boundingVolume;

		let inView = false;
		let inViewError = -Infinity;
		let inViewDistance = Infinity;
		let maxError = -Infinity;
		let minDistance = Infinity;

		for (let i = 0, l = cameraInfo.length; i < l; i++) {
			// calculate the camera error
			const info = cameraInfo[i];
			let error;
			let distance;
			if (info.isOrthographic) {
				const pixelSize = info.pixelSize;
				error = tile.geometricError / pixelSize;
				distance = Infinity;
			} else {
				const sseDenominator = info.sseDenominator;
				distance = boundingVolume.distanceToPoint(info.position);
				error = tile.geometricError / (distance * sseDenominator);
			}

			// Track which camera frustums this tile is in so we can use it
			// to ignore the error calculations for cameras that can't see it
			const frustum = cameraInfo[i].frustum;
			if (boundingVolume.intersectsFrustum(frustum)) {
				inView = true;
				inViewError = Math.max(inViewError, error);
				inViewDistance = Math.min(inViewDistance, distance);
			}

			maxError = Math.max(maxError, error);
			minDistance = Math.min(minDistance, distance);
		}

		// check the plugin visibility
		this.invokeAllPlugins(plugin => {
			if (plugin !== this && plugin.calculateTileViewError) {
				plugin.calculateTileViewError(tile, viewErrorTarget);
				if (viewErrorTarget.inView) {
					inView = true;
					inViewError = Math.max(inViewError, viewErrorTarget.error);
				}

				maxError = Math.max(maxError, viewErrorTarget.error);
			}
		});

		// If the tiles are out of view then use the global distance and error calculated
		if (inView) {
			target.inView = true;
			target.error = inViewError;
			target.distanceToCamera = inViewDistance;
		} else {
			target.inView = false;
			target.error = maxError;
			target.distanceToCamera = minDistance;
		}
	}

	ensureChildrenArePreprocessed(tile, immediate = false) {
		// TODO
	}

	// Private Functions
	preprocessTileSet(json, url, parent = null) {
		const version = json.asset.version;
		const [major, minor] = version.split('.').map(v => parseInt(v));
		console.assert(
			major <= 1,
			'Tiles3D: asset.version is expected to be a 1.x or a compatible version.'
		);

		if (major === 1 && minor > 0) {
			console.warn('Tiles3D: tiles versions at 1.1 or higher have limited support. Some new extensions and features may not be supported.');
		}

		// remove the last file path path-segment from the URL including the trailing slash
		let basePath = url.replace(/\/[^/]*$/, '');
		basePath = new URL(basePath, window.location.href).toString();

		traverseSet(
			json.root,
			(node, parent) => this.preprocessNode(node, basePath, parent),
			null,
			parent,
			parent ? parent.__depth : 0
		);
	}

	loadRootTileSet() {
		// transform the url
		let processedUrl = this.rootURL;
		this.invokeAllPlugins(plugin => processedUrl = plugin.preprocessURL ? plugin.preprocessURL(processedUrl, null) : processedUrl);

		// load the tile set root
		const pr = this
			.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(processedUrl, this.fetchOptions))
			.then(res => {
				if (res.ok) {
					return res.json();
				} else {
					throw new Error(`Tiles3D: Failed to load tileset "${processedUrl}" with status ${res.status} : ${res.statusText}`);
				}
			})
			.then(root => {
				const { extensions = {} } = root;

				// update the ellipsoid based on the extension
				if ('3DTILES_ellipsoid' in extensions) {
					const ext = extensions['3DTILES_ellipsoid'];
					const { ellipsoid } = this;
					ellipsoid.name = ext.body;
					if (ext.radii) {
						ellipsoid.radius.set(...ext.radii);
					} else {
						ellipsoid.radius.set(1, 1, 1);
					}
				}

				this.preprocessTileSet(root, processedUrl);

				return root;
			});

		return pr;
	}

	requestTileContents(tile) {
		// If the tile is already being loaded then don't
		// start it again.
		if (tile.__loadingState !== UNLOADED) {
			return;
		}

		const stats = this.stats;
		const lruCache = this.lruCache;
		const downloadQueue = this.downloadQueue;
		const parseQueue = this.parseQueue;

		const isExternalTileSet = tile.__hasUnrenderableContent;

		lruCache.add(tile, t => {
			if (t.__loadingState === LOADING) {
				// Stop the load if it's started
				t.__loadAbort.abort();
				t.__loadAbort = null;
			} else if (isExternalTileSet) {
				t.children.length = 0;
			} else {
				this.disposeTile(t);
			}

			// Decrement stats
			if (t.__loadingState === LOADING) {
				stats.downloading--;
			} else if (t.__loadingState === PARSING) {
				stats.parsing--;
			}

			t.__loadingState = UNLOADED;
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
		tile.__loadingState = LOADING;

		const errorCallback = e => {
			// if it has been unloaded then the tile has been disposed
			if (tile.__loadIndex !== loadIndex) {
				return;
			}

			if (e.name !== 'AbortError') {
				downloadQueue.remove(tile);
				parseQueue.remove(tile);

				if (tile.__loadingState === PARSING) {
					stats.parsing--;
				} else if (tile.__loadingState === LOADING) {
					stats.downloading--;
				}

				stats.failed++;
				tile.__loadingState = FAILED;

				// Handle fetch bug for switching examples in index.html.
				// https://stackoverflow.com/questions/12009423/what-does-status-canceled-for-a-resource-mean-in-chrome-developer-tools
				if (e.message !== 'Failed to fetch') {
					console.error(`Tiles3D: Failed to load tile at url "${tile.content.uri}".`);
					console.error(e);
				}
			} else {
				lruCache.remove(tile);
			}
		};

		let uri = new URL(tile.content.uri, tile.__basePath + '/').toString();
		this.invokeAllPlugins(plugin => uri = plugin.preprocessURL ? plugin.preprocessURL(uri, tile) : uri);

		if (isExternalTileSet) {
			downloadQueue.add(tile, tileCb => {
				// if it has been unloaded then the tile has been disposed
				if (tileCb.__loadIndex !== loadIndex) {
					return Promise.resolve();
				}

				return this.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(uri, { ...this.fetchOptions, signal }));
			}).then(res => {
				if (res.ok) {
					return res.json();
				} else {
					throw new Error(`Tiles3D: Failed to load tileset "${uri}" with status ${res.status} : ${res.statusText}`);
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
				tile.__loadingState = LOADED;

				tile.children.push(json.root);

				this.dispatchEvent({
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

				return this.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(uri, { ...this.fetchOptions, signal }));
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
				tile.__loadingState = PARSING;

				return parseQueue.add(tile, parseTile => {
					// if it has been unloaded then the tile has been disposed
					if (parseTile.__loadIndex !== loadIndex) {
						return Promise.resolve();
					}

					const uri = parseTile.content.uri;
					const extension = getUrlExtension(uri);

					return this.parseTile(buffer, parseTile, extension);
				});
			}).then(() => {
				// if it has been unloaded then the tile has been disposed
				if (tile.__loadIndex !== loadIndex) {
					return;
				}

				stats.parsing--;
				tile.__loadingState = LOADED;

				if (tile.__wasSetVisible) {
					this.invokeOnePlugin(plugin => plugin.setTileVisible && plugin.setTileVisible(tile, true));
				}

				if (tile.__wasSetActive) {
					this.setTileActive(tile, true);
				}
			}).catch(errorCallback);
		}
	}

	getAttributions(target = []) {
		this.invokeAllPlugins(plugin => plugin !== this && plugin.getAttributions && plugin.getAttributions(target));
		return target;
	}

	invokeOnePlugin(func) {
		const plugins = [...this.plugins, this];
		for (let i = 0; i < plugins.length; i++) {
			const result = func(plugins[i]);
			if (result) {
				return result;
			}
		}

		return null;
	}

	invokeAllPlugins(func) {
		const plugins = [...this.plugins, this];
		const pending = [];
		for (let i = 0; i < plugins.length; i++) {
			const result = func(plugins[i]);
			if (result) {
				pending.push(result);
			}
		}

		return pending.length === 0 ? null : Promise.all(pending);
	}

	//
	get autoDisableRendererCulling() {
		return this._autoDisableRendererCulling;
	}

	set autoDisableRendererCulling(value) {
		if (this._autoDisableRendererCulling !== value) {
			super._autoDisableRendererCulling = value;
			this.traverse(tile => {
				const scene = tile.cached.scene;
				if (scene) {
					scene.traverse(c => {
						c.frustumCulled = c[INITIAL_FRUSTUM_CULLED] && !value;
					});
				}
			});
		}
	}

	setDRACOLoader(dracoLoader) {
		this.$modelLoader.setDRACOLoader(dracoLoader);
	}

	setKTX2Loader(ktx2Loader) {
		this.$modelLoader.setKTX2Loader(ktx2Loader);
	}

	addEventListener(type, listener, thisObject) {
		this.$events.addEventListener(type, listener, thisObject);
	}

	removeEventListener(type, listener, thisObject) {
		this.$events.removeEventListener(type, listener, thisObject);
	}

	addCamera(camera) {
		const success = this.$cameras.add(camera);
		if (success) {
			this.dispatchEvent({ type: 'add-camera', camera });
		}
		return success;
	}

	removeCamera(camera) {
		const success = this.$cameras.remove(camera);
		if (success) {
			this.dispatchEvent({ type: 'delete-camera', camera });
		}
		return success;
	}

	resize(width, height) {
		this.$cameras.setResolution(width, height);
	}

	raycast(ray, intersects) {
		if (!this.root) {
			return null;
		}

		raycastTraverse(this.root, this, ray, intersects);

		intersects.sort(distanceSort);
	}

	raycastFirst(ray) {
		if (!this.root) {
			return null;
		}

		return raycastTraverseFirstHit(this.root, this, ray);
	}

	getBoundingBox(box) {
		if (!this.root) {
			return false;
		}

		const boundingVolume = this.root.cached.boundingVolume;

		if (boundingVolume) {
			boundingVolume.getBoundingBox(box);
			return true;
		} else {
			return false;
		}
	}

	getOrientedBoundingBox(targetBox, targetMatrix) {
		if (!this.root) {
			return false;
		}

		const boundingVolume = this.root.cached.boundingVolume;

		if (boundingVolume) {
			boundingVolume.getOrientedBoundingBox(targetBox, targetMatrix);
			return true;
		} else {
			return false;
		}
	}

	getBoundingSphere(sphere) {
		if (!this.root) {
			return false;
		}

		const boundingVolume = this.root.cached.boundingVolume;

		if (boundingVolume) {
			boundingVolume.getBoundingSphere(sphere);
			return true;
		} else {
			return false;
		}
	}

	forEachLoadedModel(callback) {
		this.traverse(tile => {
			const scene = tile.cached && tile.cached.scene;
			if (scene) {
				callback(scene, tile);
			}
		}, null, false);
	}

	// override Object3D methods
	updateMatrix(force) {
		if (this.matrixAutoUpdate || this.matrixNeedsUpdate) {
			this.matrix.transform(this.position, this.scale, this.quaternion);

			this.matrixNeedsUpdate = false;
			this.worldMatrixNeedsUpdate = true;
		}

		if (this.worldMatrixNeedsUpdate || force) {
			if (this.parent === null) {
				tempMat.copy(this.matrix);
			} else {
				tempMat.multiplyMatrices(this.parent.worldMatrix, this.matrix);
			}
			this.worldMatrixNeedsUpdate = false;

			if (!matrixEquals(tempMat, this.worldMatrix)) {
				this.worldMatrix.copy(tempMat);

				// update children
				// the children will not have to change unless the parent group has updated
				const children = this.children;
				for (let i = 0, l = children.length; i < l; i++) {
					children[i].updateMatrix();
				}
			}
		}
	}

}