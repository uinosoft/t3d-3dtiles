import { Box3, Matrix4, Plane, Ray, Vector3 } from 't3d';
import { OBB } from './OBB.js';

export class TileOBB extends OBB {

	constructor() {
		super();

		// cache obb points and planes
		// to speed up intersection test with frustum
		this._points = new Array(8).fill().map(() => new Vector3());
		this._planes = new Array(6).fill().map(() => new Plane());

		// cache obb origin box and transform matrix4x4
		// to speed up intersection test with ray and error calculation
		this._originBox = new Box3();
		this._originBoxTransform = new Matrix4();
		this._originBoxTransformInverse = new Matrix4();
	}

	updateCache() {
		this.getPoints(this._points);
		this.getPlanes(this._planes);

		this.toBoundingBoxWithTransform(
			this._originBox, this._originBoxTransform
		);
		this._originBoxTransformInverse.copy(this._originBoxTransform).inverse();
	}

	containsPoint(point) {
		_vec3_1.copy(point).applyMatrix4(this._originBoxTransformInverse);
		return this.box.containsPoint(_vec3_1);
	}

	intersectsRay(ray) {
		_ray_1.copy(ray).applyMatrix4(this._originBoxTransformInverse);
		return _ray_1.intersectsBox(this._originBox);
	}

	intersectRay(ray, target) {
		_ray_1.copy(ray).applyMatrix4(this._originBoxTransformInverse);

		if (_ray_1.intersectBox(this._originBox, target)) {
			return target.applyMatrix4(this._originBoxTransform);
		}

		return null;
	}

	// optimized intersection test with frustum
	intersectsFrustum(frustum) {
		for (let i = 0; i < 6; i++) {
			const plane = frustum.planes[i];
			let maxDistance = -Infinity;
			for (let j = 0; j < 8; j++) {
				const v = this._points[j];
				const dist = plane.distanceToPoint(v);
				maxDistance = maxDistance < dist ? dist : maxDistance;
			}

			if (maxDistance < 0) {
				return false;
			}
		}

		// do the opposite check using the obb planes to avoid false positives
		for (let i = 0; i < 6; i++) {
			const plane = this._planes[i];
			let maxDistance = -Infinity;
			for (let j = 0; j < 8; j++) {
				const v = frustum.points[j];
				const dist = plane.distanceToPoint(v);
				maxDistance = maxDistance < dist ? dist : maxDistance;
			}

			if (maxDistance < 0) {
				return false;
			}
		}

		return true;
	}

	distanceToPoint(point) {
		// originBoxTransformInverse has no scale,
		// so we don't need to scale the distance
		_vec3_1.copy(point).applyMatrix4(this._originBoxTransformInverse);
		return this._originBox.distanceToPoint(_vec3_1);
	}

	getBoundingSphere(target) {
		return this.box.getBoundingSphere(target);
	}

	getBoundingBox(target) {
		return target.setFromPoints(this._points);
	}

}

const _ray_1 = new Ray();
const _vec3_1 = new Vector3();