import { Vector2, Vector3, Matrix4 } from 't3d';
import { FastFrustum } from '../math/FastFrustum.js';

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

	updateInfos(originMatrix) {
		const cameras = this._cameras;
		const cameraCount = cameras.length;
		const infos = this._infos;
		const resolution = this._resolution;

		if (cameraCount === 0) {
			console.warn('CameraList.updateInfos(): No camera added.');
			return;
		}

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

		// get inverse scale of origin matrix

		_mat4_1.copy(originMatrix).inverse();

		const invScaleX = _vec3_1.setFromMatrixColumn(_mat4_1, 0).getLength();
		const invScaleY = _vec3_1.setFromMatrixColumn(_mat4_1, 1).getLength();
		const invScaleZ = _vec3_1.setFromMatrixColumn(_mat4_1, 2).getLength();

		if (Math.abs(Math.max(invScaleX - invScaleY, invScaleX - invScaleZ)) > 1e-6) {
			console.warn('CameraList.updateInfos(): Non uniform scale used for tile which may cause issues when calculating screen space error.');
		}

		const invScale = invScaleX;
		const invOriginMatrix = _mat4_1;

		// update the camera infos

		for (let i = 0, l = infos.length; i < l; i++) {
			const camera = cameras[i];
			const info = infos[i];

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

			info.invScale = invScale;

			// get frustum in origin space
			_mat4_2.copy(originMatrix).premultiply(camera.projectionViewMatrix);
			info.frustum.setFromMatrix(_mat4_2);
			info.frustum.updateCache();

			// get camera position in origin space
			info.position.setFromMatrixPosition(camera.worldMatrix).applyMatrix4(invOriginMatrix);
		}
	}

	getInfos() {
		return this._infos;
	}

}

const _mat4_1 = new Matrix4();
const _mat4_2 = new Matrix4();
const _vec3_1 = new Vector3();