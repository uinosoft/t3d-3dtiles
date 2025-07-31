export { ImplicitTilingPlugin } from './core/plugins/ImplicitTilingPlugin.js';

export * from './tiles/TilesRenderer.js';

export * from './math/OBB.js';

export * from './loaders/B3DMLoader.js';
export * from './loaders/CMPTLoader.js';
export * from './loaders/TileGLTFLoader.js';
export * from './loaders/I3DMLoader.js';
export * from './loaders/PNTSLoader.js';

export { GlobeControls } from './controls/GlobeControls.js';
export { EnvironmentControls } from './controls/EnvironmentControls.js';

export { TilesFadePlugin } from './plugins/fade/TilesFadePlugin.js';
export { CesiumIonAuthPlugin } from './plugins/CesiumIonAuthPlugin.js';
export { DebugTilesPlugin } from './plugins/DebugTilesPlugin.js';
export { ReorientationPlugin } from './plugins/ReorientationPlugin.js';
export { QuantizedMeshPlugin } from './plugins/QuantizedMeshPlugin.js';
export { XYZTilesPlugin, TMSTilesPlugin } from './plugins/images/EPSGTilesPlugin.js';

export { LoadParser as DebugLoadParser } from './loaders/parsers/LoadParser.js';

// Exporting some prototype methods for Matrix4, Vector3, Camera and Object3D classes

import { Vector3, Quaternion, Triangle, Object3D, MathUtils, Camera } from 't3d';

const _quaternion = new Quaternion();

Vector3.prototype.applyEuler = function(euler) {
	return this.applyQuaternion(_quaternion.setFromEuler(euler));
};

Vector3.prototype.applyAxisAngle = function(axis, angle) {
	return this.applyQuaternion(_quaternion.setFromAxisAngle(axis, angle));
};

Triangle.prototype.setFromAttributeAndIndices = function(attribute, i0, i1, i2) {
	const array = attribute.buffer.array;
	const itemSize = attribute.size;
	const offset = attribute.offset;
	this.a.fromArray(array, i0 * itemSize + offset);
	this.b.fromArray(array, i1 * itemSize + offset);
	this.c.fromArray(array, i2 * itemSize + offset);
};

Object3D.prototype.removeFromParent = function() {
	const parent = this.parent;

	if (parent !== null) {
		parent.remove(this);
	}

	return this;
};

MathUtils.DEG2RAD = Math.PI / 180;
MathUtils.RAD2DEG = 180 / Math.PI;

const _setOrtho = Camera.prototype.setOrtho;
Camera.prototype.setOrtho = function(left, right, bottom, top, near, far) {
	this.left = left;
	this.right = right;
	this.bottom = bottom;
	this.top = top;
	this.near = near;
	this.far = far;

	this.zoom = 1;

	this.isPerspectiveCamera = false;
	this.isOrthographicCamera = true;

	_setOrtho.call(this, left, right, bottom, top, near, far);
};

const _setPerspective = Camera.prototype.setPerspective;
Camera.prototype.setPerspective = function(fov, aspect, near, far) {
	this.fov = fov;
	this.aspect = aspect;
	this.near = near;
	this.far = far;

	this.isPerspectiveCamera = true;
	this.isOrthographicCamera = false;

	_setPerspective.call(this, fov, aspect, near, far);
};

Camera.prototype.updateProjectionMatrix = function() {
	if (this.isOrthographicCamera) {
		this.setOrtho(
			this.left, this.right, this.bottom, this.top, this.near, this.far
		);
	} else if (this.isPerspectiveCamera) {
		this.setPerspective(
			this.fov, this.aspect, this.near, this.far
		);
	}

	return this;
};