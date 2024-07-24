import { Color3, Object3D } from 't3d';
import { Box3Helper } from 't3d/addons/objects/Box3Helper.js';
import { SphereHelper } from 't3d/addons/objects/SphereHelper.js';

export class Tiles3DHelper extends Object3D {

	constructor(tiles3D) {
		super();

		const boxGroup = new Object3D();
		boxGroup.name = 'Tiles3DHelper.boxGroup';
		this.add(boxGroup);

		const sphereGroup = new Object3D();
		sphereGroup.name = 'Tiles3DHelper.sphereGroup';
		this.add(sphereGroup);

		const helperMap = new Map();

		function getHelperCache(tile) {
			let helperCache = helperMap.get(tile);
			if (!helperCache) {
				helperCache = {
					box: null,
					sphere: null
				};
				helperMap.set(tile, helperCache);
			}
			return helperCache;
		}

		tiles3D.addEventListener('TileVisibilityChanged', ({ scene, tile, visible }) => {
			const helperCache = getHelperCache(tile);

			if (visible) {
				if (helperCache.box !== null) {
					boxGroup.add(helperCache.box);
					helperCache.box.updateMatrix(true);
				}

				if (helperCache.sphere !== null) {
					sphereGroup.add(helperCache.sphere);
					helperCache.sphere.updateMatrix(true);
				}
			} else {
				if (helperCache.box !== null) {
					boxGroup.remove(helperCache.box);
				}

				if (helperCache.sphere !== null) {
					sphereGroup.remove(helperCache.sphere);
				}
			}
		});

		tiles3D.addEventListener('TileLoaded', ({ scene, tile }) => {
			const helperCache = getHelperCache(tile);
			const boundingVolume = tile.cached.boundingVolume;

			if (boundingVolume.obb) {
				// Create debug bounding box
				// In some cases the bounding box may have a scale of 0 in one dimension resulting
				// in the NaNs in an extracted rotation so we disable matrix updates instead.
				const boxWrapper = new Object3D();
				boxWrapper.name = 'Tiles3DHelper.boxWrapper';
				boxWrapper.matrix.copy(boundingVolume.obb._originBoxTransform);
				boxWrapper.matrixAutoUpdate = false;
				boxWrapper.matrixNeedsUpdate = false;

				const boxHelper = new Box3Helper(boundingVolume.obb._originBox, getIndexedRandomColor(tile.__depth));
				boxHelper.name = 'Box3Helper';
				boxHelper.raycast = emptyRaycast;
				boxWrapper.add(boxHelper);

				helperCache.box = boxWrapper;

				if (tiles3D.visibleTiles.has(tile)) {
					boxGroup.add(boxWrapper);
					boxWrapper.updateMatrix(true);
				}
			}

			if (boundingVolume.sphere) {
				const sphereHelper = new SphereHelper(boundingVolume.sphere, getIndexedRandomColor(tile.__depth));
				sphereHelper.name = 'SphereHelper';
				sphereHelper.raycast = emptyRaycast;

				helperCache.sphere = sphereHelper;

				if (tiles3D.visibleTiles.has(tile)) {
					sphereGroup.add(sphereHelper);
					sphereHelper.updateMatrix(true);
				}
			}
		});

		tiles3D.addEventListener('TileDisposed', ({ scene, tile }) => {
			const helperCache = getHelperCache(tile);

			if (helperCache.box !== null) {
				helperCache.box.children[0].geometry.dispose();
				helperCache.box.children[0].material.dispose();
				delete helperCache.box;
			}

			if (helperCache.sphere !== null) {
				helperCache.sphere.geometry.dispose();
				helperCache.sphere.material.dispose();
				delete helperCache.sphere;
			}
		});

		this.boxGroup = boxGroup;
		this.sphereGroup = sphereGroup;
	}

}

function emptyRaycast() {}

const _color3_1 = new Color3();

const _colors = {};

// Return a consistant random color for an index
const getIndexedRandomColor = index => {
	if (!_colors[index]) {
		const h = Math.random();
		const s = 0.5 + Math.random() * 0.5;
		const l = 0.375 + Math.random() * 0.25;

		_colors[index] = _color3_1.setHSL(h, s, l).getHex();
	}
	return _colors[index];
};