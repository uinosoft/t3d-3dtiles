import { EventDispatcher, LoadingManager, Matrix4, Object3D } from 't3d';
import { TileBoundingVolume } from './math/TileBoundingVolume.js';
import { traverseSet, getUrlExtension } from './utils/Utils.js';
import { raycastTraverse, raycastTraverseFirstHit, distanceSort } from './utils/IntersectionUtils.js';
import RequestState from './utils/RequestState.js';
import { CameraList } from './utils/CameraList.js';
import { TilesLoader } from './utils/TilesLoader.js';
import { ModelLoader } from './utils/ModelLoader.js';
import { schedulingTiles } from './utils/TilesScheduler.js';
import { WGS84_ELLIPSOID } from './math/GeoConstants.js';

export class Tiles3D extends Object3D {

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

		// internals

		this.rootURL = url;
		this._rootTileSet = null;

		this._autoDisableRendererCulling = true;

		this.plugins = [];

		this.$cameras = new CameraList();
		this.$tilesLoader = new TilesLoader();
		this.$modelLoader = new ModelLoader(manager);
		this.$events = new EventDispatcher();
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

	fetchData(url, options) {
		return fetch(url, options);
	}

	get rootTileSet() {
		const rootTileSet = this._rootTileSet;

		if (!rootTileSet) {
			this._rootTileSet = this.invokeOnePlugin(plugin => plugin.loadRootTileSet && plugin.loadRootTileSet())
				.then(root => {
					let processedUrl = this.rootURL;
					if (processedUrl !== null) {
						this.invokeAllPlugins(plugin => processedUrl = plugin.preprocessURL ? plugin.preprocessURL(processedUrl, null) : processedUrl);
					}

					this._rootTileSet = root;

					this.dispatchEvent({
						type: 'load-tile-set',
						tileSet: root,
						url: processedUrl
					});
				})
				.catch(err => {
					console.error(err);
					this._rootTileSet = err;
				});
			return null;
		} else if (rootTileSet instanceof Promise || rootTileSet instanceof Error) {
			return null;
		}

		return rootTileSet;
	}

	get root() {
		const rootTileSet = this.rootTileSet;
		return rootTileSet ? rootTileSet.root : null;
	}

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

	dispatchEvent(event) {
		this.$events.dispatchEvent(event);
	}

	markTileUsed(tile) {
		this.$tilesLoader.lruCache.markUsed(tile);
	}

	update() {
		const rootTile = this.root;
		if (!rootTile) return;

		this.dispatchEvent(_updateBeforeEvent);

		const { stats } = this;
		stats.inFrustum = 0;
		stats.used = 0;
		stats.active = 0;
		stats.visible = 0;

		this.frameCount++;

		this.$cameras.updateInfos(this.worldMatrix);

		schedulingTiles(rootTile, this);

		this.dispatchEvent(_updateAfterEvent);
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
		}, null);
	}

	traverse(beforecb, aftercb) {
		const rootTile = this.root;
		if (!rootTile) return;
		traverseSet(rootTile, beforecb, aftercb);
	}

	resetFailedTiles() {
		const stats = this.stats;
		if (stats.failed === 0) {
			return;
		}

		this.traverse(tile => {
			if (tile.__loadingState === RequestState.FAILED) {
				tile.__loadingState = RequestState.UNLOADED;
			}
		});

		stats.failed = 0;
	}

	dispose() {
		// dispose of all the plugins
		this.plugins.forEach(plugin => {
			this.unregisterPlugin(plugin);
		});

		const lruCache = this.$tilesLoader.lruCache;
		this.traverse(tile => {
			lruCache.remove(tile);
		});

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

		tile.__loadingState = RequestState.UNLOADED;

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

	setTileVisible(tile, visible) {
		const scene = tile.cached.scene;
		if (!scene) {
			return;
		}
		const visibleTiles = this.visibleTiles;
		if (visible) {
			this.add(scene);
			visibleTiles.add(tile);
			scene.updateMatrix(true); // TODO: remove this?
		} else {
			this.remove(scene);
			visibleTiles.delete(tile);
		}

		this.dispatchEvent({
			type: 'tile-visibility-change',
			scene,
			tile,
			visible
		});
	}

	$setTileActive(tile, active) {
		const activeTiles = this.activeTiles;
		if (active) {
			activeTiles.add(tile);
		} else {
			activeTiles.delete(tile);
		}
	}

	$disposeTile(tile) {
		// This could get called before the tile has finished downloading
		const cached = tile.cached;
		if (cached.scene) {
			const materials = cached.materials;
			const geometry = cached.geometry;
			const textures = cached.textures;

			for (let i = 0, l = geometry.length; i < l; i++) {
				geometry[i].dispose();
			}

			for (let i = 0, l = materials.length; i < l; i++) {
				materials[i].dispose();
			}

			for (let i = 0, l = textures.length; i < l; i++) {
				const texture = textures[i];
				texture.dispose();
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
		}

		tile._loadIndex++;
	}

}

const _updateBeforeEvent = { type: 'update-before' };
const _updateAfterEvent = { type: 'update-after' };

const PLUGIN_REGISTERED = Symbol('PLUGIN_REGISTERED');

const INITIAL_FRUSTUM_CULLED = Symbol('INITIAL_FRUSTUM_CULLED');

const tempMat = new Matrix4();

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