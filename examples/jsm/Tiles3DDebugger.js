import { PBRMaterial, PointsMaterial, SHADING_TYPE, DRAW_MODE, Sphere } from 't3d';
import { InstancedPBRMaterial, InstancedBasicMaterial } from 't3d-3dtiles';
import { Tiles3DHelper } from './helpers/Tiles3DHelper.js';

export class Tiles3DDebugger {

	constructor(tiles3D) {
		this.helper = new Tiles3DHelper(tiles3D);
		tiles3D.add(this.helper);

		this.colorMode = NONE;

		this.maxDebugDepth = -1;
		this.maxDebugDistance = -1;
		this.maxDebugError = -1;

		this.customColorCallback = null;

		this.getDebugColor = (value, target) => {
			target.setRGB(value, value, value);
		};

		let extremeDebugDepth = -1;
		let extremeDebugError = -1;

		tiles3D.addEventListener('TileSetLoaded', ({ json, url }) => {
			// only update for the root tileset

			if (url !== tiles3D.rootURL) return;

			// initialize the extreme values of the hierarchy

			let maxDepth = -1;
			tiles3D.traverse(tile => {
				maxDepth = Math.max(maxDepth, tile.__depth);
			});

			let maxError = -1;
			tiles3D.traverse(tile => {
				maxError = Math.max(maxError, tile.cached.geometricError);
			});

			extremeDebugDepth = maxDepth;
			extremeDebugError = maxError;
		});

		this.update = () => {
			if (!tiles3D.root) return;

			// get max values to use for materials
			let maxDepth = -1;
			if (this.maxDebugDepth === -1) {
				maxDepth = extremeDebugDepth;
			} else {
				maxDepth = this.maxDebugDepth;
			}

			let maxError = -1;
			if (this.maxDebugError === -1) {
				maxError = extremeDebugError;
			} else {
				maxError = this.maxDebugError;
			}

			let maxDistance = -1;
			if (this.maxDebugDistance === -1) {
				const boundingVolume = tiles3D.root.cached.boundingVolume;
				const rootSphere = boundingVolume.getBoundingSphere(_sphere_1);
				maxDistance = rootSphere.radius;
			} else {
				maxDistance = this.maxDebugDistance;
			}

			const errorTarget = tiles3D.errorTarget;
			const colorMode = this.colorMode;
			const visibleTiles = tiles3D.visibleTiles;

			visibleTiles.forEach(tile => {
				const scene = tile.cached.scene;

				// create a random color per-tile
				let h, s, l;
				if (colorMode === RANDOM_COLOR) {
					h = Math.random();
					s = 0.5 + Math.random() * 0.5;
					l = 0.375 + Math.random() * 0.25;
				}

				scene.traverse(c => {
					if (!c[ORIGINAL_MATERIAL]) {
						c[ORIGINAL_MATERIAL] = c.material;
					}

					if (colorMode === RANDOM_NODE_COLOR) {
						h = Math.random();
						s = 0.5 + Math.random() * 0.5;
						l = 0.375 + Math.random() * 0.25;
					}

					if (c.material) {
						const originalMaterial = c[ORIGINAL_MATERIAL];

						if (colorMode === NONE && c.material !== originalMaterial) {
							c.material.dispose();
							c.material = c[ORIGINAL_MATERIAL];
						} else if (colorMode !== NONE && c.material === originalMaterial) {
							if (c.material.drawMode === DRAW_MODE.POINTS) {
								const pointsMaterial = new PointsMaterial();
								pointsMaterial.size = originalMaterial.size;
								pointsMaterial.sizeAttenuation = originalMaterial.sizeAttenuation;
								c.material = pointsMaterial;
							} else {
								if (c.material.isInstancedPBRMaterial) {
									c.material = new InstancedPBRMaterial();
								} else if (c.material.isInstancedBasicMaterial) {
									c.material = new InstancedBasicMaterial();
								} else {
									c.material = new PBRMaterial();
								}

								c.material.shading = SHADING_TYPE.FLAT_SHADING;
							}
						}

						if (colorMode !== RANDOM_COLOR) {
							delete c.material[HAS_RANDOM_COLOR];
						}

						if (colorMode !== RANDOM_NODE_COLOR) {
							delete c.material[HAS_RANDOM_NODE_COLOR];
						}

						switch (colorMode) {
							case DEPTH: {
								const val = tile.__depth / maxDepth;
								this.getDebugColor(val, c.material.diffuse);
								break;
							}
							case RELATIVE_DEPTH: {
								const val = tile.__depthFromRenderedParent / maxDepth;
								this.getDebugColor(val, c.material.diffuse);
								break;
							}
							case SCREEN_ERROR: {
								const val = tile.__error / errorTarget;
								if (val > 1.0) {
									c.material.diffuse.setRGB(1.0, 0.0, 0.0);
								} else {
									this.getDebugColor(val, c.material.diffuse);
								}
								break;
							}
							case GEOMETRIC_ERROR: {
								const val = Math.min(tile.cached.geometricError / maxError, 1);
								this.getDebugColor(val, c.material.diffuse);
								break;
							}
							case DISTANCE: {
								// We don't update the distance if the geometric error is 0.0 so
								// it will always be black.
								const val = Math.min(tile.__distanceFromCamera / maxDistance, 1);
								this.getDebugColor(val, c.material.diffuse);
								break;
							}
							case IS_LEAF: {
								if (!tile.children || tile.children.length === 0) {
									this.getDebugColor(1.0, c.material.diffuse);
								} else {
									this.getDebugColor(0.0, c.material.diffuse);
								}
								break;
							}
							case RANDOM_NODE_COLOR: {
								if (!c.material[HAS_RANDOM_NODE_COLOR]) {
									c.material.diffuse.setHSL(h, s, l);
									c.material[HAS_RANDOM_NODE_COLOR] = true;
								}
								break;
							}
							case RANDOM_COLOR: {
								if (!c.material[HAS_RANDOM_COLOR]) {
									c.material.diffuse.setHSL(h, s, l);
									c.material[HAS_RANDOM_COLOR] = true;
								}
								break;
							}
							case CUSTOM_COLOR: {
								if (this.customColorCallback) {
									this.customColorCallback(tile, c);
								} else {
									console.warn('Tiles3DDebugger: customColorCallback not defined');
								}
								break;
							}
						}
					}
				});
			});
		};
	}

}

const _sphere_1 = new Sphere();

const ORIGINAL_MATERIAL = Symbol('ORIGINAL_MATERIAL');
const HAS_RANDOM_COLOR = Symbol('HAS_RANDOM_COLOR');
const HAS_RANDOM_NODE_COLOR = Symbol('HAS_RANDOM_NODE_COLOR');

export const NONE = 0;
export const SCREEN_ERROR = 1;
export const GEOMETRIC_ERROR = 2;
export const DISTANCE = 3;
export const DEPTH = 4;
export const RELATIVE_DEPTH = 5;
export const IS_LEAF = 6;
export const RANDOM_COLOR = 7;
export const RANDOM_NODE_COLOR = 8;
export const CUSTOM_COLOR = 9;