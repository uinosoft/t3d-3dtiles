
import { Matrix4, Vector3 } from 't3d';
import { B3DMLoader } from '../loaders/B3DMLoader.js';
import { I3DMLoader } from '../loaders/I3DMLoader.js';
import { PNTSLoader } from '../loaders/PNTSLoader.js';
import { CMPTLoader } from '../loaders/CMPTLoader.js';
import { TileGLTFLoader } from '../loaders/TileGLTFLoader.js';

export class ModelLoader {

	constructor(manager) {
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
			['gltf', gltfLoader],
			['glb', gltfLoader]
		]);
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

	loadTileContent(buffer, tile, extension, tiles3D) {
		tile._loadIndex = tile._loadIndex || 0;
		tile._loadIndex++;

		const uri = tile.content.uri;
		const uriSplits = uri.split(/[\\\/]/g); // eslint-disable-line no-useless-escape
		uriSplits.pop();
		const workingPath = uriSplits.join('/');
		const fetchOptions = tiles3D.fetchOptions;

		const loadIndex = tile._loadIndex;
		let promise = null;

		const upAxis = tiles3D.rootTileSet.asset && tiles3D.rootTileSet.asset.gltfUpAxis || 'y';
		const cached = tile.cached;
		const cachedTransform = cached.transform;

		switch (upAxis.toLowerCase()) {
			case 'x':
				mat4_1.makeRotationAxis(Y_AXIS, -Math.PI / 2);
				break;
			case 'y':
				mat4_1.makeRotationAxis(X_AXIS, Math.PI / 2);
				break;
			case 'z':
				mat4_1.identity();
				break;
		}

		const loader = this._loaders.get(extension);

		if (loader) {
			const config = {
				fetchOptions,
				path: workingPath,
				buffer
			};

			if (extension === 'b3dm' || extension === 'i3dm' || extension === 'cmpt') {
				config.adjustmentTransform = mat4_1.clone();
			}

			promise = loader.load(uri, config);
		} else {
			console.warn(`TilesRenderer: Content type "${extension}" not supported.`);
			promise = Promise.resolve(null);
		}

		return promise.then(resource => {
			const scene = resource.root;

			if (tile._loadIndex !== loadIndex || !scene) {
				return;
			}

			// ensure the matrix is up to date in case the scene has a transform applied
			scene.updateMatrix();

			// apply the local up-axis correction rotation
			// GLTFLoader seems to never set a transformation on the root scene object so
			// any transformations applied to it can be assumed to be applied after load
			// (such as applying RTC_CENTER) meaning they should happen _after_ the z-up
			// rotation fix which is why "multiply" happens here.
			if (extension === 'glb' || extension === 'gltf') {
				scene.matrix.multiply(mat4_1);
			}

			scene.matrix.premultiply(cachedTransform);
			scene.matrix.decompose(scene.position, scene.quaternion, scene.scale);

			cached.scene = scene;
			cached.featureTable = resource.featureTable;
			cached.batchTable = resource.batchTable;

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

			cached.materials = materials;
			cached.geometry = geometry;
			cached.textures = textures;

			return scene;
		});
	}

}

const mat4_1 = new Matrix4();

const X_AXIS = new Vector3(1, 0, 0);
const Y_AXIS = new Vector3(0, 1, 0);