import { Vector2, Vector3, Matrix4 } from 't3d';
import { FastFrustum } from '../math/FastFrustum.js';

const _mat4_1 = new Matrix4();
const _vec3_1 = new Vector3();

export class CameraList {

	constructor() {
		this._cameras = [];
		this._resolution = new Vector2();
		this._infos = [];
	}

	add(camera) {
		const cameras = this._cameras;

		if (cameras.indexOf(camera) === -1) {
			cameras.push(camera);
			return true;
		}

		return false;
	}

	remove(camera) {
		const cameras = this._cameras;
		const index = cameras.indexOf(camera);

		if (index !== -1) {
			cameras.splice(index, 1);
			return true;
		}

		return false;
	}

	setResolution(width, height) {
		this._resolution.set(width, height);
	}

	updateInfos(group) {
		const cameras = this._cameras;
		const infos = this._infos;
		const resolution = this._resolution;

		// automatically scale the array of infos to match the cameras
		while (infos.length > cameras.length) {
			infos.pop();
		}

		while (infos.length < cameras.length) {
			infos.push({
				frustum: new FastFrustum(), // in origin space
				isOrthographic: false,
				sseDenominator: -1, // used if isOrthographic is false
				position: new Vector3(), // in origin space
				invScale: -1,
				pixelSize: 0 // used if isOrthographic is true
			});
		}

		// extract scale of group container
		_vec3_1.setFromMatrixScale(group.worldMatrixInverse);
		if (Math.abs(Math.max(_vec3_1.x - _vec3_1.y, _vec3_1.x - _vec3_1.z)) > 1e-6) {
			console.warn('TilesRenderer : Non uniform scale used for tile which may cause issues when calculating screen space error.');
		}

		// store the camera cameraInfo in the 3d tiles root frame
		for (let i = 0, l = infos.length; i < l; i++) {
			const camera = cameras[i];
			const info = infos[i];
			const frustum = info.frustum;
			const position = info.position;

			const cameraResolutionX = resolution.x * (camera.rect.z - camera.rect.x);
			const cameraResolutionY = resolution.y * (camera.rect.w - camera.rect.y);

			if (cameraResolutionX === 0 || cameraResolutionY === 0) {
				console.warn('CameraList.updateInfos(): Resolution for camera error calculation is not set.');
			}

			// Read the calculated projection matrix directly to support custom Camera implementations
			const projection = camera.projectionMatrix.elements;

			// The last element of the projection matrix is 1 for orthographic, 0 for perspective
			info.isOrthographic = projection[15] === 1;

			if (info.isOrthographic) {
				// the view width and height are used to populate matrix elements 0 and 5.
				const w = 2 / projection[0];
				const h = 2 / projection[5];
				info.pixelSize = Math.max(h / cameraResolutionY, w / cameraResolutionX);
			} else {
				// the vertical FOV is used to populate matrix element 5.
				info.sseDenominator = (2 / projection[5]) / cameraResolutionY;
			}

			// get frustum in origin space
			_mat4_1.copy(group.worldMatrix).premultiply(camera.projectionViewMatrix);

			frustum.setFromMatrix(_mat4_1);
			frustum.updateCache();

			// get camera position in origin space
			position.setFromMatrixPosition(camera.worldMatrix).applyMatrix4(group.worldMatrixInverse);
		}
	}

	getInfos() {
		return this._infos;
	}

}