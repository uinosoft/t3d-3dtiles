import { TilesRendererBase } from '../core/renderer/tiles/TilesRendererBase.js';
import { B3DMLoader } from '../loaders/B3DMLoader.js';
import { I3DMLoader } from '../loaders/I3DMLoader.js';
import { PNTSLoader } from '../loaders/PNTSLoader.js';
import { CMPTLoader } from '../loaders/CMPTLoader.js';
import { TilesGroup } from './TilesGroup.js';
import { EventDispatcher, LoadingManager, Matrix4, Vector3 } from 't3d';
import { raycastTraverse, raycastTraverseFirstHit } from './raycastTraverse.js';
import { readMagicBytes } from '../core/renderer/utilities/readMagicBytes.js';
import { TileBoundingVolume } from '../math/TileBoundingVolume.js';
import { CameraList } from './CameraList.js';
import { estimateBytesUsed } from './utilities.js';
import { WGS84_ELLIPSOID } from '../math/GeoConstants.js';
import { TileGLTFLoader } from '../loaders/TileGLTFLoader.js';

const INITIAL_FRUSTUM_CULLED = Symbol('INITIAL_FRUSTUM_CULLED');
const viewErrorTarget = {
	inView: false,
	error: Infinity
};

const X_AXIS = new Vector3(1, 0, 0);
const Y_AXIS = new Vector3(0, 1, 0);

function updateFrustumCulled(object, toInitialValue) {
	object.traverse(c => {
		c.frustumCulled = c[INITIAL_FRUSTUM_CULLED] && toInitialValue;
	});
}

const _updateBeforeEvent = { type: 'update-before' };
const _updateAfterEvent = { type: 'update-after' };

export class TilesRenderer extends TilesRendererBase {

	get autoDisableRendererCulling() {
		return this._autoDisableRendererCulling;
	}

	set autoDisableRendererCulling(value) {
		if (this._autoDisableRendererCulling !== value) {
			this._autoDisableRendererCulling = value;
			this.forEachLoadedModel(scene => {
				updateFrustumCulled(scene, !value);
			});
		}
	}

	get optimizeRaycast() {
		return this._optimizeRaycast;
	}

	set optimizeRaycast(v) {
		console.warn('TilesRenderer: The "optimizeRaycast" option has been deprecated.');
		this._optimizeRaycast = v;
	}

	constructor(...args) {
		super(...args);
		this.group = new TilesGroup(this);
		this.ellipsoid = WGS84_ELLIPSOID.clone();
		this.$cameras = new CameraList();
		this._optimizeRaycast = true;
		this._upRotationMatrix = new Matrix4();
		this._bytesUsed = new WeakMap();

		// flag indicating whether frustum culling should be disabled
		this._autoDisableRendererCulling = true;

		const manager = new LoadingManager();
		manager.setURLModifier(url => {
			if (this.preprocessURL) {
				return this.preprocessURL(url);
			} else {
				return url;
			}
		});
		this.manager = manager;

		// saved for event dispatcher functions
		this._listeners = {};

		// loaders
		const b3dmLoader = new B3DMLoader(manager);
		const i3dmLoader = new I3DMLoader(manager);
		const pntsLoader = new PNTSLoader(manager);
		const cmptLoader = new CMPTLoader(manager);
		const gltfLoader = new TileGLTFLoader(manager);
		this._loaders = new Map([
			['b3dm', b3dmLoader],
			['i3dm', i3dmLoader],
			['pnts', pntsLoader],
			['cmpt', cmptLoader],
			['gltf', gltfLoader]
		]);
	}

	addEventListener(...args) {
		EventDispatcher.prototype.addEventListener.call(this, ...args);
	}

	hasEventListener(...args) {
		EventDispatcher.prototype.hasEventListener.call(this, ...args);
	}

	removeEventListener(...args) {
		EventDispatcher.prototype.removeEventListener.call(this, ...args);
	}

	dispatchEvent(...args) {
		EventDispatcher.prototype.dispatchEvent.call(this, ...args);
	}

	// Public API
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

	getBoundingSphere(target) {
		if (!this.root) {
			return false;
		}

		const boundingVolume = this.root.cached.boundingVolume;

		if (boundingVolume) {
			boundingVolume.getBoundingSphere(target);
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

	raycast(ray, intersects) {
		if (!this.root) {
			return null;
		}

		raycastTraverse(this, this.root, ray, intersects);
	}

	raycastFirst(ray) {
		if (!this.root) {
			return null;
		}

		return raycastTraverseFirstHit(this, this.root, ray);
	}

	addCamera(camera) {
		const success = this.$cameras.add(camera);
		if (success) {
			this.dispatchEvent({ type: 'add-camera', camera });
		}
		return success;
	}

	setResolution(width, height) {
		this.$cameras.setResolution(width, height);
	}

	removeCamera(camera) {
		const success = this.$cameras.remove(camera);
		if (success) {
			this.dispatchEvent({ type: 'delete-camera', camera });
		}
		return success;
	}

	/* Overriden */
	loadRootTileSet(...args) {
		return super.loadRootTileSet(...args)
			.then(root => {
				// cache the gltf tile set rotation matrix
				const { asset, extensions = {} } = root;
				const upAxis = asset && asset.gltfUpAxis || 'y';
				switch (upAxis.toLowerCase()) {
					case 'x':
						this._upRotationMatrix.makeRotationAxis(Y_AXIS, -Math.PI / 2);
						break;

					case 'y':
						this._upRotationMatrix.makeRotationAxis(X_AXIS, Math.PI / 2);
						break;
					default:
						this._upRotationMatrix.identity();
						break;
				}

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

				return root;
			});
	}

	update() {
		// check if the plugins that can block the tile updates require it
		let needsUpdate = null;
		this.invokeAllPlugins(plugin => {
			if (plugin.doTilesNeedUpdate) {
				const res = plugin.doTilesNeedUpdate();
				if (needsUpdate === null) {
					needsUpdate = res;
				} else {
					needsUpdate = Boolean(needsUpdate || res);
				}
			}
		});

		if (needsUpdate === false) {
			this.dispatchEvent({ type: 'update-before' });
			this.dispatchEvent({ type: 'update-after' });
			return;
		}

		this.dispatchEvent(_updateBeforeEvent);

		this.$cameras.updateInfos(this.group);

		super.update();

		this.dispatchEvent(_updateAfterEvent);

		// check for cameras _after_ base update so we can enable pre-loading the root tile set
		const cameras = this.$cameras._cameras;
		if (cameras.length === 0 && this.root) {
			let found = false;
			this.invokeAllPlugins(plugin => found = found || Boolean(plugin !== this && plugin.calculateTileViewError));
			if (found === false) {
				console.warn('TilesRenderer: no cameras defined. Cannot update 3d tiles.');
			}
		}
	}

	preprocessNode(tile, tileSetDir, parentTile = null) {
		// TODO: multiple contents (1.1) are not supported yet
		if (tile.contents) {
			tile.content = tile.contents[0];
		}

		super.preprocessNode(tile, tileSetDir, parentTile);

		const transform = new Matrix4();
		if (tile.transform) {
			transform.fromArray(tile.transform);
		}

		if (parentTile) {
			transform.premultiply(parentTile.cached.transform);
		}

		const transformInverse = new Matrix4().copy(transform).invert();
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

			featureTable: null,
			batchTable: null
		};
	}

	async parseTile(buffer, tile, extension, uri, abortSignal) {
		const cached = tile.cached;
		const uriSplits = uri.split(/[\\\/]/g); // eslint-disable-line no-useless-escape
		uriSplits.pop();
		const workingPath = uriSplits.join('/');
		const fetchOptions = this.fetchOptions;

		let promise = null;

		const cachedTransform = cached.transform;
		const upRotationMatrix = this._upRotationMatrix;
		const fileType = (readMagicBytes(buffer) || extension).toLowerCase();

		switch (fileType) {
			case 'b3dm': {
				promise = this._loaders.get('b3dm').load(uri, {
					fetchOptions,
					path: workingPath,
					buffer,
					adjustmentTransform: upRotationMatrix.clone()
				});
				break;
			}
			case 'pnts': {
				promise = this._loaders.get('pnts').load(uri, {
					fetchOptions,
					path: workingPath,
					buffer
				});
				break;
			}
			case 'i3dm': {
				promise = this._loaders.get('i3dm').load(uri, {
					fetchOptions,
					path: workingPath,
					buffer,
					adjustmentTransform: upRotationMatrix.clone()
				});
				break;
			}
			case 'cmpt': {
				promise = this._loaders.get('i3dm').load(uri, {
					fetchOptions,
					path: workingPath,
					buffer,
					adjustmentTransform: upRotationMatrix.clone()
				});
				break;
			}
			case 'gltf':
			case 'glb': {
				promise = this._loaders.get('gltf').load(uri, {
					fetchOptions,
					path: workingPath,
					buffer
				}).then(result => {
					// apply the local up-axis correction rotation
					// GLTFLoader seems to never set a transformation on the root scene object so
					// any transformations applied to it can be assumed to be applied after load
					// (such as applying RTC_CENTER) meaning they should happen _after_ the z-up
					// rotation fix which is why "multiply" happens here.
					const { root: scene } = result;
					scene.matrix
						.multiply(upRotationMatrix)
						.decompose(scene.position, scene.quaternion, scene.scale);
					return result;
				});
				break;
			}
			default: {
				promise = this.invokeOnePlugin(plugin => plugin.parseToMesh && plugin.parseToMesh(buffer, tile, extension, uri, abortSignal));
				break;
			}
		}

		// wait for the tile to load
		const result = await promise;
		if (result === null) {
			throw new Error(`TilesRenderer: Content type "${fileType}" not supported.`);
		}

		// get the scene data
		let scene;
		let metadata;
		if (result.isObject3D) {
			scene = result;
			metadata = null;
		} else {
			scene = result.root;
			metadata = result;
		}

		// ensure the matrix is up to date in case the scene has a transform applied
		scene.updateMatrix();
		scene.matrix.premultiply(cachedTransform);
		scene.matrix.decompose(scene.position, scene.quaternion, scene.scale);

		// wait for extra processing by plugins if needed
		await this.invokeAllPlugins(plugin => {
			return plugin.processTileModel && plugin.processTileModel(scene, tile);
		});

		// frustum culling
		scene.traverse(c => {
			c[INITIAL_FRUSTUM_CULLED] = c.frustumCulled;
		});
		updateFrustumCulled(scene, !this.autoDisableRendererCulling);

		// collect all original geometries, materials, etc to be disposed of later
		const materials = [];
		const geometry = [];
		const textures = [];
		scene.traverse(c => {
			if (c.geometry) {
				geometry.push(c.geometry);
			}

			if (c.material) {
				const material = c.material;
				materials.push(c.material);

				for (const key in material) {
					const value = material[key];
					if (value && value.isTexture) {
						textures.push(value);
					}
				}
			}
		});

		// exit early if a new request has already started
		if (abortSignal.aborted) {
			// dispose of any image bitmaps that have been opened.
			// TODO: share this code with the "disposeTile" code below, possibly allow for the tiles
			// renderer base to trigger a disposal of unneeded data
			for (let i = 0, l = textures.length; i < l; i++) {
				const texture = textures[i];

				if (texture.image instanceof ImageBitmap) {
					texture.image.close();
				}

				texture.dispose();
			}

			return;
		}

		scene.batchTable = result.batchTable;
		scene.featureTable = result.featureTable;

		cached.materials = materials;
		cached.geometry = geometry;
		cached.textures = textures;
		cached.scene = scene;
		cached.metadata = metadata;
	}

	disposeTile(tile) {
		super.disposeTile(tile);

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
	}

	setTileVisible(tile, visible) {
		const scene = tile.cached.scene;
		const group = this.group;

		if (visible) {
			if (scene) {
				group.add(scene);
				scene.updateMatrix(true);
			}
		} else {
			if (scene) {
				group.remove(scene);
			}
		}

		super.setTileVisible(tile, visible);

		this.dispatchEvent({
			type: 'tile-visibility-change',
			scene,
			tile,
			visible
		});
	}

	calculateBytesUsed(tile, scene) {
		const bytesUsed = this._bytesUsed;
		if (!bytesUsed.has(tile) && scene) {
			bytesUsed.set(tile, estimateBytesUsed(scene));
		}

		return bytesUsed.get(tile) ?? null;
	}

	calculateTileViewError(tile, target) {
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
				// avoid dividing 0 by 0 which can result in NaN. If the distance to the tile is
				// 0 then the error should be infinity.
				const sseDenominator = info.sseDenominator;
				distance = boundingVolume.distanceToPoint(info.position);
				error = distance === 0 ? Infinity : tile.geometricError / (distance * sseDenominator);
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

	dispose() {
		super.dispose();
		this.group.removeFromParent();
	}

	setDRACOLoader(dracoLoader) {
		this._loaders.get('b3dm').setDRACOLoader(dracoLoader);
		this._loaders.get('i3dm').setDRACOLoader(dracoLoader);
		this._loaders.get('cmpt').setDRACOLoader(dracoLoader);
		this._loaders.get('gltf').setDRACOLoader(dracoLoader);
	}

	setKTX2Loader(ktx2Loader) {
		this._loaders.get('b3dm').setKTX2Loader(ktx2Loader);
		this._loaders.get('i3dm').setKTX2Loader(ktx2Loader);
		this._loaders.get('cmpt').setKTX2Loader(ktx2Loader);
		this._loaders.get('gltf').setKTX2Loader(ktx2Loader);
	}

}