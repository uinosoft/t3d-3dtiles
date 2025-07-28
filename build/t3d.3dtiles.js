// t3d-3dtiles
(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports, require('t3d')) :
	typeof define === 'function' && define.amd ? define(['exports', 't3d'], factory) :
	(global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.t3d = global.t3d || {}, global.t3d));
})(this, (function (exports, t3d) { 'use strict';

	/**
	 * An oriented bounding box.
	 */
	class OBB {
		/**
				* Create a new OBB.
				* @param {Box3} [box] - The axis-aligned bounding box.
				* @param {Matrix3} [rotation] - The rotation matrix.
				*/
		constructor(box = new t3d.Box3(), rotation = new t3d.Matrix3()) {
			this.box = box;
			this.rotation = rotation;
		}

		/**
				* Set the OBB from a center point and axes with half sizes.
				* @param {Vector3} center - The center of the OBB.
				* @param {Vector3} axisX - The X axis with half size.
				* @param {Vector3} axisY - The Y axis with half size.
				* @param {Vector3} axisZ - The Z axis with half size.
				* @return {OBB} A reference to this OBB.
				*/
		setFromCenterAndAxes(center, axisX, axisY, axisZ) {
			_vec3_1$4.copy(axisX);
			_vec3_2$1.copy(axisY);
			_vec3_3$1.copy(axisZ);
			const scaleX = _vec3_1$4.getLength();
			const scaleY = _vec3_2$1.getLength();
			const scaleZ = _vec3_3$1.getLength();
			_vec3_1$4.normalize();
			_vec3_2$1.normalize();
			_vec3_3$1.normalize();

			// handle the case where the box has a dimension of 0 in one axis
			if (scaleX === 0) {
				_vec3_1$4.crossVectors(_vec3_2$1, _vec3_3$1);
			}
			if (scaleY === 0) {
				_vec3_2$1.crossVectors(_vec3_1$4, _vec3_3$1);
			}
			if (scaleZ === 0) {
				_vec3_3$1.crossVectors(_vec3_1$4, _vec3_2$1);
			}
			this.rotation.set(_vec3_1$4.x, _vec3_2$1.x, _vec3_3$1.x, _vec3_1$4.y, _vec3_2$1.y, _vec3_3$1.y, _vec3_1$4.z, _vec3_2$1.z, _vec3_3$1.z);
			const halfSize = _vec3_1$4.set(scaleX, scaleY, scaleZ);
			this.box.min.copy(center).sub(halfSize);
			this.box.max.copy(center).add(halfSize);
			return this;
		}

		/**
				* Transforms this OBB with the supplied matrix.
				* @param {Matrix4} matrix - The transformation matrix.
				* @return {OBB} A reference to this OBB.
				*/
		applyMatrix4(matrix) {
			const e = matrix.elements;
			let sx = _vec3_1$4.set(e[0], e[1], e[2]).getLength();
			const sy = _vec3_1$4.set(e[4], e[5], e[6]).getLength();
			const sz = _vec3_1$4.set(e[8], e[9], e[10]).getLength();
			const det = matrix.determinant();
			if (det < 0) sx = -sx;
			_mat3_1$2.setFromMatrix4(matrix);
			const invSX = 1 / sx;
			const invSY = 1 / sy;
			const invSZ = 1 / sz;
			_mat3_1$2.elements[0] *= invSX;
			_mat3_1$2.elements[1] *= invSX;
			_mat3_1$2.elements[2] *= invSX;
			_mat3_1$2.elements[3] *= invSY;
			_mat3_1$2.elements[4] *= invSY;
			_mat3_1$2.elements[5] *= invSY;
			_mat3_1$2.elements[6] *= invSZ;
			_mat3_1$2.elements[7] *= invSZ;
			_mat3_1$2.elements[8] *= invSZ;
			this.rotation.multiply(_mat3_1$2);
			const center = this.box.getCenter(_vec3_1$4);
			const halfSize = this.box.getSize(_vec3_2$1).multiplyScalar(0.5);
			halfSize.x *= sx;
			halfSize.y *= sy;
			halfSize.z *= sz;

			// https://github.com/mrdoob/three.js/issues/21753
			center.applyMatrix4(matrix);
			this.box.min.copy(center).sub(halfSize);
			this.box.max.copy(center).add(halfSize);
			return this;
		}

		/**
				* Get the 8 corner points of the OBB.
				* @param {Vector3[]} points - The array to store the points.
				* @return {Vector3[]} The array of points.
				*/
		getPoints(points) {
			const center = this.box.getCenter(_vec3_1$4);
			const min = _vec3_2$1.subVectors(this.box.min, center);
			const max = _vec3_3$1.subVectors(this.box.max, center);
			let index = 0;
			for (let x = -1; x <= 1; x += 2) {
				for (let y = -1; y <= 1; y += 2) {
					for (let z = -1; z <= 1; z += 2) {
						points[index].set(x < 0 ? min.x : max.x, y < 0 ? min.y : max.y, z < 0 ? min.z : max.z).applyMatrix3(this.rotation).add(center);
						index++;
					}
				}
			}
			return points;
		}

		/**
				* Get the 6 planes of the OBB.
				* @param {Plane[]} planes - The array to store the planes.
				* @return {Plane[]} The array of planes.
				*/
		getPlanes(planes) {
			const center = this.box.getCenter(_vec3_1$4);
			const worldMin = _vec3_2$1.subVectors(this.box.min, center).applyMatrix3(this.rotation).add(center);
			const worldMax = _vec3_3$1.subVectors(this.box.max, center).applyMatrix3(this.rotation).add(center);
			_vec3_1$4.set(0, 0, 1).applyMatrix3(this.rotation).normalize();
			planes[0].setFromNormalAndCoplanarPoint(_vec3_1$4, worldMin);
			planes[1].setFromNormalAndCoplanarPoint(_vec3_1$4, worldMax);
			planes[1].normal.negate();
			planes[1].constant *= -1;
			_vec3_1$4.set(0, 1, 0).applyMatrix3(this.rotation).normalize();
			planes[2].setFromNormalAndCoplanarPoint(_vec3_1$4, worldMin);
			planes[3].setFromNormalAndCoplanarPoint(_vec3_1$4, worldMax);
			planes[3].normal.negate();
			planes[3].constant *= -1;
			_vec3_1$4.set(1, 0, 0).applyMatrix3(this.rotation).normalize();
			planes[4].setFromNormalAndCoplanarPoint(_vec3_1$4, worldMin);
			planes[5].setFromNormalAndCoplanarPoint(_vec3_1$4, worldMax);
			planes[5].normal.negate();
			planes[5].constant *= -1;
			return planes;
		}
		containsPoint(point) {
			const obbStruct = getOBBStruct(this, a);
			const v = _vec3_1$4.subVectors(point, obbStruct.c);

			// project _vec3_4 onto each axis and check if these points lie inside the OBB

			return Math.abs(v.dot(obbStruct.u[0])) <= obbStruct.e[0] && Math.abs(v.dot(obbStruct.u[1])) <= obbStruct.e[1] && Math.abs(v.dot(obbStruct.u[2])) <= obbStruct.e[2];
		}

		/**
		 * Reference: Closest Point on OBB to Point in Real-Time Collision Detection
		 * by Christer Ericson (chapter 5.1.4)
		 */
		clampPoint(point, result) {
			const obbStruct = getOBBStruct(this, a);
			const v = _vec3_1$4.subVectors(point, obbStruct.c);

			// start at the center position of the OBB

			result.copy(obbStruct.c);

			// project the target onto the OBB axes and walk towards that point

			const x = Math.max(Math.min(v.dot(obbStruct.u[0]), obbStruct.e[0]), -obbStruct.e[0]);
			result.add(obbStruct.u[0].multiplyScalar(x));
			const y = Math.max(Math.min(v.dot(obbStruct.u[1]), obbStruct.e[1]), -obbStruct.e[1]);
			result.add(obbStruct.u[1].multiplyScalar(y));
			const z = Math.max(Math.min(v.dot(obbStruct.u[2]), obbStruct.e[2]), -obbStruct.e[2]);
			result.add(obbStruct.u[2].multiplyScalar(z));
			return result;
		}
		intersectsSphere(sphere) {
			// find the point on the OBB closest to the sphere center
			this.clampPoint(sphere.center, closestPoint);
			// if that point is inside the sphere, the OBB and sphere intersect
			return closestPoint.distanceToSquared(sphere.center) <= sphere.radius * sphere.radius;
		}
		intersectsOBB(obb, epsilon = Number.EPSILON) {
			// prepare data structures (the code uses the same nomenclature like the reference)

			getOBBStruct(this, a);
			getOBBStruct(obb, b);

			// compute rotation matrix expressing b in a's coordinate frame

			for (let i = 0; i < 3; i++) {
				for (let j = 0; j < 3; j++) {
					R[i][j] = a.u[i].dot(b.u[j]);
				}
			}

			// compute translation vector

			const v1 = _vec3_1$4.subVectors(b.c, a.c);

			// bring translation into a's coordinate frame

			t[0] = v1.dot(a.u[0]);
			t[1] = v1.dot(a.u[1]);
			t[2] = v1.dot(a.u[2]);

			// compute common subexpressions. Add in an epsilon term to
			// counteract arithmetic errors when two edges are parallel and
			// their cross product is (near) null

			for (let i = 0; i < 3; i++) {
				for (let j = 0; j < 3; j++) {
					AbsR[i][j] = Math.abs(R[i][j]) + epsilon;
				}
			}
			let ra, rb;

			// test axes L = A0, L = A1, L = A2

			for (let i = 0; i < 3; i++) {
				ra = a.e[i];
				rb = b.e[0] * AbsR[i][0] + b.e[1] * AbsR[i][1] + b.e[2] * AbsR[i][2];
				if (Math.abs(t[i]) > ra + rb) return false;
			}

			// test axes L = B0, L = B1, L = B2

			for (let i = 0; i < 3; i++) {
				ra = a.e[0] * AbsR[0][i] + a.e[1] * AbsR[1][i] + a.e[2] * AbsR[2][i];
				rb = b.e[i];
				if (Math.abs(t[0] * R[0][i] + t[1] * R[1][i] + t[2] * R[2][i]) > ra + rb) return false;
			}

			// test axis L = A0 x B0

			ra = a.e[1] * AbsR[2][0] + a.e[2] * AbsR[1][0];
			rb = b.e[1] * AbsR[0][2] + b.e[2] * AbsR[0][1];
			if (Math.abs(t[2] * R[1][0] - t[1] * R[2][0]) > ra + rb) return false;

			// test axis L = A0 x B1

			ra = a.e[1] * AbsR[2][1] + a.e[2] * AbsR[1][1];
			rb = b.e[0] * AbsR[0][2] + b.e[2] * AbsR[0][0];
			if (Math.abs(t[2] * R[1][1] - t[1] * R[2][1]) > ra + rb) return false;

			// test axis L = A0 x B2

			ra = a.e[1] * AbsR[2][2] + a.e[2] * AbsR[1][2];
			rb = b.e[0] * AbsR[0][1] + b.e[1] * AbsR[0][0];
			if (Math.abs(t[2] * R[1][2] - t[1] * R[2][2]) > ra + rb) return false;

			// test axis L = A1 x B0

			ra = a.e[0] * AbsR[2][0] + a.e[2] * AbsR[0][0];
			rb = b.e[1] * AbsR[1][2] + b.e[2] * AbsR[1][1];
			if (Math.abs(t[0] * R[2][0] - t[2] * R[0][0]) > ra + rb) return false;

			// test axis L = A1 x B1

			ra = a.e[0] * AbsR[2][1] + a.e[2] * AbsR[0][1];
			rb = b.e[0] * AbsR[1][2] + b.e[2] * AbsR[1][0];
			if (Math.abs(t[0] * R[2][1] - t[2] * R[0][1]) > ra + rb) return false;

			// test axis L = A1 x B2

			ra = a.e[0] * AbsR[2][2] + a.e[2] * AbsR[0][2];
			rb = b.e[0] * AbsR[1][1] + b.e[1] * AbsR[1][0];
			if (Math.abs(t[0] * R[2][2] - t[2] * R[0][2]) > ra + rb) return false;

			// test axis L = A2 x B0

			ra = a.e[0] * AbsR[1][0] + a.e[1] * AbsR[0][0];
			rb = b.e[1] * AbsR[2][2] + b.e[2] * AbsR[2][1];
			if (Math.abs(t[1] * R[0][0] - t[0] * R[1][0]) > ra + rb) return false;

			// test axis L = A2 x B1

			ra = a.e[0] * AbsR[1][1] + a.e[1] * AbsR[0][1];
			rb = b.e[0] * AbsR[2][2] + b.e[2] * AbsR[2][0];
			if (Math.abs(t[1] * R[0][1] - t[0] * R[1][1]) > ra + rb) return false;

			// test axis L = A2 x B2

			ra = a.e[0] * AbsR[1][2] + a.e[1] * AbsR[0][2];
			rb = b.e[0] * AbsR[2][1] + b.e[1] * AbsR[2][0];
			if (Math.abs(t[1] * R[0][2] - t[0] * R[1][2]) > ra + rb) return false;

			// since no separating axis is found, the OBBs must be intersecting

			return true;
		}

		/**
				* To AABB with matrix4 transform.
				* This method can ensure the AABB center is origin,
				* and transform only contains rotation and offset.
				* @param {Box3} box3 - The target box.
				* @param {Matrix4} transform - The transformation matrix.
				*/
		toBoundingBoxWithTransform(box3, transform) {
			const center = this.box.getCenter(_vec3_1$4);
			box3.min.copy(this.box.min).sub(center);
			box3.max.copy(this.box.max).sub(center);
			const e = this.rotation.elements;
			transform.set(e[0], e[3], e[6], center.x, e[1], e[4], e[7], center.y, e[2], e[5], e[8], center.z, 0, 0, 0, 1);
		}
	}
	const closestPoint = new t3d.Vector3();
	const _vec3_1$4 = new t3d.Vector3();
	const _vec3_2$1 = new t3d.Vector3();
	const _vec3_3$1 = new t3d.Vector3();
	const _mat3_1$2 = new t3d.Matrix3();
	const R = [[], [], []];
	const AbsR = [[], [], []];
	const t = [];
	const _halfSize = new t3d.Vector3();
	const a = {
		c: new t3d.Vector3(),
		// center
		u: [new t3d.Vector3(), new t3d.Vector3(), new t3d.Vector3()],
		// basis vectors
		e: [] // half width
	};
	const b = {
		c: new t3d.Vector3(),
		// center
		u: [new t3d.Vector3(), new t3d.Vector3(), new t3d.Vector3()],
		// basis vectors
		e: [] // half width
	};
	function getOBBStruct(obb, target) {
		obb.box.getCenter(target.c);
		const rotationElements = obb.rotation.elements;
		target.u[0].fromArray(rotationElements, 0);
		target.u[1].fromArray(rotationElements, 3);
		target.u[2].fromArray(rotationElements, 6);
		obb.box.getSize(_halfSize).multiplyScalar(0.5).toArray(target.e);
		return target;
	}

	class TileOBB extends OBB {
		constructor() {
			super();

			// cache obb points and planes
			// to speed up intersection test with frustum
			this._points = new Array(8).fill().map(() => new t3d.Vector3());
			this._planes = new Array(6).fill().map(() => new t3d.Plane());

			// cache obb origin box and transform matrix4x4
			// to speed up intersection test with ray and error calculation
			this._originBox = new t3d.Box3();
			this._originBoxTransform = new t3d.Matrix4();
			this._originBoxTransformInverse = new t3d.Matrix4();
		}
		updateCache() {
			this.getPoints(this._points);
			this.getPlanes(this._planes);
			this.toBoundingBoxWithTransform(this._originBox, this._originBoxTransform);
			this._originBoxTransformInverse.copy(this._originBoxTransform).inverse();
		}
		containsPoint(point) {
			_vec3_1$3.copy(point).applyMatrix4(this._originBoxTransformInverse);
			return this.box.containsPoint(_vec3_1$3);
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
			_vec3_1$3.copy(point).applyMatrix4(this._originBoxTransformInverse);
			return this._originBox.distanceToPoint(_vec3_1$3);
		}
		getBoundingSphere(target) {
			return this.box.getBoundingSphere(target);
		}
		getBoundingBox(target) {
			return target.setFromPoints(this._points);
		}
	}
	const _ray_1 = new t3d.Ray();
	const _vec3_1$3 = new t3d.Vector3();

	// Cesium / 3D tiles Spheroid:
	// - Up is Z at 90 degrees latitude
	// - 0, 0 latitude, longitude is X axis
	//			Z
	//			|
	//			|
	//			.----- Y
	//		 /
	//	 X

	// t3d.js Spherical Coordinates
	// - Up is Y at 90 degrees latitude
	// - 0, 0 latitude, longitude is Z
	//			Y
	//			|
	//			|
	//			.----- X
	//		 /
	//	 Z

	function swapToGeoFrame(target) {
		const {
			x,
			y,
			z
		} = target;
		target.x = z;
		target.y = x;
		target.z = y;
	}
	function latitudeToSphericalPhi(latitude) {
		return -latitude + Math.PI / 2;
	}

	class Ellipsoid {
		constructor(radius = new t3d.Vector3(1, 1, 1)) {
			this.name = '';
			this.radius = radius;
		}
		intersectRay(ray, target) {
			_matrix$1.makeScale(...this.radius.toArray([])).invert();
			_sphere$2.center.set(0, 0, 0);
			_sphere$2.radius = 1;
			_ray$3.copy(ray).applyMatrix4(_matrix$1);
			if (_ray$3.intersectSphere(_sphere$2, target)) {
				_matrix$1.makeScale(...this.radius.toArray([]));
				target.applyMatrix4(_matrix$1);
				return target;
			} else {
				return null;
			}
		}

		// returns a frame with Z indicating altitude
		// Y pointing north
		// X pointing east
		getEastNorthUpFrame(lat, lon, target) {
			this.getEastNorthUpAxes(lat, lon, _vecX, _vecY, _vecZ, _pos$4);
			return target.makeBasis(_vecX, _vecY, _vecZ).setPosition(_pos$4);
		}
		getEastNorthUpAxes(lat, lon, vecEast, vecNorth, vecUp, point = _pos$4) {
			this.getCartographicToPosition(lat, lon, 0, point);
			this.getCartographicToNormal(lat, lon, vecUp); // up
			vecEast.set(-point.y, point.x, 0).normalize(); // east
			vecNorth.crossVectors(vecUp, vecEast).normalize(); // north
		}
		getRotationMatrixFromAzElRoll(lat, lon, az, el, roll, target, frame = ENU_FRAME) {
			this.getEastNorthUpFrame(lat, lon, _matrix$1);
			_euler.set(el, roll, -az, 'ZXY');
			target.makeRotationFromEuler(_euler).premultiply(_matrix$1).setPosition(0, 0, 0);

			// Add in the orientation adjustment for objects and cameras so "forward" and "up" are oriented
			// correctly
			if (frame === CAMERA_FRAME) {
				_euler.set(Math.PI / 2, 0, 0, 'XYZ');
				_matrix2.makeRotationFromEuler(_euler);
				target.multiply(_matrix2);
			} else if (frame === OBJECT_FRAME) {
				_euler.set(-Math.PI / 2, 0, Math.PI, 'XYZ');
				_matrix2.makeRotationFromEuler(_euler);
				target.multiply(_matrix2);
			}
			return target;
		}
		getCartographicToPosition(lat, lon, height, target) {
			// From Cesium function Ellipsoid.cartographicToCartesian
			// https://github.com/CesiumGS/cesium/blob/665ec32e813d5d6fe906ec3e87187f6c38ed5e49/packages/engine/Source/Core/Ellipsoid.js#L396
			this.getCartographicToNormal(lat, lon, _norm$2);
			const radius = this.radius;
			_vec$6.copy(_norm$2);
			_vec$6.x *= radius.x ** 2;
			_vec$6.y *= radius.y ** 2;
			_vec$6.z *= radius.z ** 2;
			const gamma = Math.sqrt(_norm$2.dot(_vec$6));
			_vec$6.multiplyScalar(1 / gamma);
			return target.copy(_vec$6).addScaledVector(_norm$2, height);
		}
		getPositionToCartographic(pos, target) {
			// From Cesium function Ellipsoid.cartesianToCartographic
			// https://github.com/CesiumGS/cesium/blob/665ec32e813d5d6fe906ec3e87187f6c38ed5e49/packages/engine/Source/Core/Ellipsoid.js#L463
			this.getPositionToSurfacePoint(pos, _vec$6);
			this.getPositionToNormal(pos, _norm$2);
			const heightDelta = _vec2$1.subVectors(pos, _vec$6);
			target.lon = Math.atan2(_norm$2.y, _norm$2.x);
			target.lat = Math.asin(_norm$2.z);
			target.height = Math.sign(heightDelta.dot(pos)) * heightDelta.getLength();
			return target;
		}
		getCartographicToNormal(lat, lon, target) {
			_spherical.set(1, latitudeToSphericalPhi(lat), lon);
			target.setFromSpherical(_spherical).normalize();

			// swap frame from the t3d.js frame to the geo coord frame
			swapToGeoFrame(target);
			return target;
		}
		getPositionToNormal(pos, target) {
			const radius = this.radius;
			target.copy(pos);
			target.x /= radius.x ** 2;
			target.y /= radius.y ** 2;
			target.z /= radius.z ** 2;
			target.normalize();
			return target;
		}
		getPositionToSurfacePoint(pos, target) {
			// From Cesium function Ellipsoid.scaleToGeodeticSurface
			// https://github.com/CesiumGS/cesium/blob/d11b746e5809ac115fcff65b7b0c6bdfe81dcf1c/packages/engine/Source/Core/scaleToGeodeticSurface.js#L25
			const radius = this.radius;
			const invRadiusSqX = 1 / radius.x ** 2;
			const invRadiusSqY = 1 / radius.y ** 2;
			const invRadiusSqZ = 1 / radius.z ** 2;
			const x2 = pos.x * pos.x * invRadiusSqX;
			const y2 = pos.y * pos.y * invRadiusSqY;
			const z2 = pos.z * pos.z * invRadiusSqZ;

			// Compute the squared ellipsoid norm.
			const squaredNorm = x2 + y2 + z2;
			const ratio = Math.sqrt(1.0 / squaredNorm);

			// As an initial approximation, assume that the radial intersection is the projection point.
			const intersection = _vec$6.copy(pos).multiplyScalar(ratio);
			if (squaredNorm < CENTER_EPS) {
				return !isFinite(ratio) ? null : target.copy(intersection);
			}

			// Use the gradient at the intersection point in place of the true unit normal.
			// The difference in magnitude will be absorbed in the multiplier.
			const gradient = _vec2$1.set(intersection.x * invRadiusSqX * 2.0, intersection.y * invRadiusSqY * 2.0, intersection.z * invRadiusSqZ * 2.0);

			// Compute the initial guess at the normal vector multiplier, lambda.
			let lambda = (1.0 - ratio) * pos.getLength() / (0.5 * gradient.getLength());
			let correction = 0.0;
			let func, denominator;
			let xMultiplier, yMultiplier, zMultiplier;
			let xMultiplier2, yMultiplier2, zMultiplier2;
			let xMultiplier3, yMultiplier3, zMultiplier3;
			do {
				lambda -= correction;
				xMultiplier = 1.0 / (1.0 + lambda * invRadiusSqX);
				yMultiplier = 1.0 / (1.0 + lambda * invRadiusSqY);
				zMultiplier = 1.0 / (1.0 + lambda * invRadiusSqZ);
				xMultiplier2 = xMultiplier * xMultiplier;
				yMultiplier2 = yMultiplier * yMultiplier;
				zMultiplier2 = zMultiplier * zMultiplier;
				xMultiplier3 = xMultiplier2 * xMultiplier;
				yMultiplier3 = yMultiplier2 * yMultiplier;
				zMultiplier3 = zMultiplier2 * zMultiplier;
				func = x2 * xMultiplier2 + y2 * yMultiplier2 + z2 * zMultiplier2 - 1.0;

				// "denominator" here refers to the use of this expression in the velocity and acceleration
				// computations in the sections to follow.
				denominator = x2 * xMultiplier3 * invRadiusSqX + y2 * yMultiplier3 * invRadiusSqY + z2 * zMultiplier3 * invRadiusSqZ;
				const derivative = -2 * denominator;
				correction = func / derivative;
			} while (Math.abs(func) > EPSILON12);
			return target.set(pos.x * xMultiplier, pos.y * yMultiplier, pos.z * zMultiplier);
		}
		calculateHorizonDistance(latitude, elevation) {
			// from https://aty.sdsu.edu/explain/atmos_refr/horizon.html
			// OG = sqrt ( 2 R h + h2 ) .
			const effectiveRadius = this.calculateEffectiveRadius(latitude);
			return Math.sqrt(2 * effectiveRadius * elevation + elevation ** 2);
		}
		calculateEffectiveRadius(latitude) {
			// This radius represents the distance from the center of the ellipsoid to the surface along the normal at the given latitude.
			// from https://en.wikipedia.org/wiki/Earth_radius#Prime_vertical
			// N = a / sqrt(1 - e^2 * sin^2(phi))
			const semiMajorAxis = this.radius.x;
			const semiMinorAxis = this.radius.z;
			const eSquared = 1 - semiMinorAxis ** 2 / semiMajorAxis ** 2;
			const phi = latitude * t3d.MathUtils.DEG2RAD;
			const sinPhiSquared = Math.sin(phi) ** 2;
			const N = semiMajorAxis / Math.sqrt(1 - eSquared * sinPhiSquared);
			return N;
		}
		getPositionElevation(pos) {
			// logic from "getPositionToCartographic"
			this.getPositionToSurfacePoint(pos, _vec$6);
			const heightDelta = _vec2$1.subVectors(pos, _vec$6);
			return Math.sign(heightDelta.dot(pos)) * heightDelta.getLength();
		}
		copy(source) {
			this.radius.copy(source.radius);
			return this;
		}
		clone() {
			return new this.constructor().copy(this);
		}
	}
	const _spherical = new t3d.Spherical();
	const _norm$2 = new t3d.Vector3();
	const _vec$6 = new t3d.Vector3();
	const _vec2$1 = new t3d.Vector3();
	const _matrix$1 = new t3d.Matrix4();
	const _matrix2 = new t3d.Matrix4();
	const _sphere$2 = new t3d.Sphere();
	const _euler = new t3d.Euler();
	const _vecX = new t3d.Vector3();
	const _vecY = new t3d.Vector3();
	const _vecZ = new t3d.Vector3();
	const _pos$4 = new t3d.Vector3();
	const _ray$3 = new t3d.Ray();
	const EPSILON12 = 1e-12;
	const CENTER_EPS = 0.1;
	const ENU_FRAME = 0;
	const CAMERA_FRAME = 1;
	const OBJECT_FRAME = 2;

	class EllipsoidRegion extends Ellipsoid {
		constructor(radius = new t3d.Vector3(1, 1, 1), latRange = new t3d.Vector2(-HALF_PI, HALF_PI), lonRange = new t3d.Vector2(0, 2 * PI), heightRange = new t3d.Vector2(0, 1)) {
			super(radius);
			this.latRange = latRange;
			this.lonRange = lonRange;
			this.heightRange = heightRange;
		}

		// refer to https://github.com/CesiumGS/cesium/blob/1.119/packages/engine/Source/Core/OrientedBoundingBox.js#L343
		// refer to https://github.com/NASA-AMMOS/3DTilesRendererJS/blob/c46d59c674e9ac1e652b9e9b65849bf02a645a6a/src/three/math/EllipsoidRegion.js#L120
		getOrientedBoundingBox(target) {
			resetPool();
			const {
				latRange,
				lonRange
			} = this;
			const latRangeValue = latRange.y - latRange.x;
			if (latRangeValue < PI / 2) {
				// get the midway point for the region
				const midLat = t3d.MathUtils.mapLinear(0.5, 0, 1, latRange.x, latRange.y);
				const midLon = t3d.MathUtils.mapLinear(0.5, 0, 1, lonRange.x, lonRange.y);

				// get the frame matrix for the box - works well for smaller regions
				this.getCartographicToNormal(midLat, midLon, _orthoZ);
				_orthoY.set(0, 0, 1);
				_orthoX.crossVectors(_orthoY, _orthoZ);
				_orthoY.crossVectors(_orthoX, _orthoZ);
			} else {
				_orthoX.set(1, 0, 0);
				_orthoY.set(0, 1, 0);
				_orthoZ.set(0, 0, 1);
			}
			target.rotation.set(_orthoX.x, _orthoY.x, _orthoZ.x, _orthoX.y, _orthoY.y, _orthoZ.y, _orthoX.z, _orthoY.z, _orthoZ.z);

			// transform the points into the local frame
			_invMatrix$1.setFromMatrix3(target.rotation).inverse();
			const points = this._getPoints(true);

			// get the center of the region
			target.box.makeEmpty();
			for (let i = 0, l = points.length; i < l; i++) {
				_center$1.copy(points[i]).applyMatrix4(_invMatrix$1);
				target.box.expandByPoint(_center$1);
			}
			target.box.getCenter(_center$1);
			_center$1.applyMatrix3(target.rotation);
			for (let i = 0, l = points.length; i < l; i++) {
				points[i].sub(_center$1).applyMatrix4(_invMatrix$1).add(_center$1);
			}
			target.box.makeEmpty();
			target.box.setFromPoints(points);
		}
		getBoundingSphere(target) {
			resetPool();
			const points = this._getPoints(true);
			target.makeEmpty();
			target.setFromPoints(points);
		}
		_getPoints(usePool = false) {
			const {
				latRange,
				lonRange,
				heightRange
			} = this;
			const midLat = t3d.MathUtils.mapLinear(0.5, 0, 1, latRange.x, latRange.y);
			const midLon = t3d.MathUtils.mapLinear(0.5, 0, 1, lonRange.x, lonRange.y);
			const lonOffset = Math.floor(lonRange.x / HALF_PI) * HALF_PI;
			const latlon = [[-PI / 2, 0], [PI / 2, 0], [0, lonOffset], [0, lonOffset + PI / 2], [0, lonOffset + PI], [0, lonOffset + 3 * PI / 2], [latRange.x, lonRange.y], [latRange.y, lonRange.y], [latRange.x, lonRange.x], [latRange.y, lonRange.x], [0, lonRange.x], [0, lonRange.y], [midLat, midLon], [latRange.x, midLon], [latRange.y, midLon], [midLat, lonRange.x], [midLat, lonRange.y]];
			const target = [];
			const total = latlon.length;
			for (let z = 0; z <= 1; z++) {
				const height = t3d.MathUtils.mapLinear(z, 0, 1, heightRange.x, heightRange.y);
				for (let i = 0, l = total; i < l; i++) {
					const [lat, lon] = latlon[i];
					if (lat >= latRange.x && lat <= latRange.y && lon >= lonRange.x && lon <= lonRange.y) {
						const v = getVector(usePool);
						target.push(v);
						this.getCartographicToPosition(lat, lon, height, v);
					}
				}
			}
			return target;
		}
	}
	const _orthoX = new t3d.Vector3();
	const _orthoY = new t3d.Vector3();
	const _orthoZ = new t3d.Vector3();
	const _center$1 = new t3d.Vector3();
	const _invMatrix$1 = new t3d.Matrix4();
	const PI = Math.PI;
	const HALF_PI = PI / 2;
	let _poolIndex = 0;
	const _pointsPool = [];
	function getVector(usePool = false) {
		if (!usePool) {
			return new t3d.Vector3();
		}
		if (!_pointsPool[_poolIndex]) {
			_pointsPool[_poolIndex] = new t3d.Vector3();
		}
		_poolIndex++;
		return _pointsPool[_poolIndex - 1];
	}
	function resetPool() {
		_poolIndex = 0;
	}

	class TileBoundingVolume {
		constructor() {
			this.sphere = null;
			this.obb = null;
			this.region = null;
		}
		setOBBData(data, transform) {
			const obb = new TileOBB();
			obb.setFromCenterAndAxes(_vec3_4.set(data[0], data[1], data[2]), _vec3_1$2.set(data[3], data[4], data[5]), _vec3_2.set(data[6], data[7], data[8]), _vec3_3.set(data[9], data[10], data[11])).applyMatrix4(transform);
			obb.updateCache();
			this.obb = obb;
		}
		setSphereData(data, transform) {
			const sphere = new t3d.Sphere();
			sphere.center.set(data[0], data[1], data[2]);
			sphere.radius = data[3];
			sphere.applyMatrix4(transform);
			this.sphere = sphere;
		}
		setRegionData(ellipsoid, west, south, east, north, minHeight, maxHeight) {
			const region = new EllipsoidRegion(ellipsoid.radius.clone(), new t3d.Vector2(south, north), new t3d.Vector2(west, east), new t3d.Vector2(minHeight, maxHeight));
			this.region = region;
			const obb = new TileOBB();
			region.getOrientedBoundingBox(obb);
			obb.updateCache();
			this.obb = obb;

			// const sphere = new Sphere();
			// obb.getBoundingSphere(sphere);
			// this.sphere = sphere;
		}
		intersectsRay(ray) {
			const sphere = this.sphere;
			const obb = this.obb;

			// Early out if we don't hit this tile sphere
			if (sphere && !ray.intersectsSphere(sphere)) {
				return false;
			}

			// Early out if we don't this this tile box
			if (obb && !obb.intersectsRay(ray)) {
				return false;
			}
			return true;
		}
		intersectRay(ray, target) {
			const sphere = this.sphere;
			const obb = this.obb;
			let sphereDistSq = -Infinity;
			let obbDistSq = -Infinity;
			if (sphere) {
				if (ray.intersectSphere(sphere, _vec3_1$2)) {
					sphereDistSq = sphere.containsPoint(ray.origin) ? 0 : ray.origin.distanceToSquared(_vec3_1$2);
				}
			}
			if (obb) {
				if (obb.intersectRay(ray, _vec3_1$2)) {
					obbDistSq = obb.containsPoint(ray.origin) ? 0 : ray.origin.distanceToSquared(_vec3_1$2);
				}
			}

			// if we didn't hit anything then exit
			const furthestDist = Math.max(sphereDistSq, obbDistSq);
			if (furthestDist === -Infinity) {
				return null;
			}

			// get the furthest hit point if needed
			return ray.at(Math.sqrt(furthestDist), target);
		}
		distanceToPoint(point) {
			const sphere = this.sphere;
			const obb = this.obb;
			let sphereDistance = -Infinity;
			let obbDistance = -Infinity;
			if (sphere) {
				// Sphere#distanceToPoint is negative inside the sphere, whereas Box3#distanceToPoint is
				// zero inside the box. Clipping the distance to a minimum of zero ensures that both
				// types of bounding volume behave the same way.
				sphereDistance = Math.max(sphere.distanceToPoint(point), 0);
			}
			if (obb) {
				obbDistance = obb.distanceToPoint(point);
			}

			// return the further distance of the two volumes
			return sphereDistance > obbDistance ? sphereDistance : obbDistance;
		}
		intersectsFrustum(frustum) {
			const sphere = this.sphere;
			if (sphere && !frustum.intersectsSphere(sphere)) {
				return false;
			}
			const obb = this.obb;
			if (obb && !obb.intersectsFrustum(frustum)) {
				return false;
			}

			// if we don't have a sphere or obb then just say we did intersect
			return Boolean(sphere || obb);
		}
		getOrientedBoundingBox(targetBox, targetMatrix) {
			if (this.obb) {
				targetBox.copy(this.obb._originBox);
				targetMatrix.copy(this.obb._originBoxTransform);
			} else {
				this.getBoundingBox(targetBox);
				targetMatrix.identity();
			}
		}
		getBoundingBox(target) {
			if (this.sphere) {
				return this.sphere.getBoundingBox(target);
			} else {
				return this.obb.getBoundingBox(target);
			}
		}
		getBoundingSphere(target) {
			if (this.sphere) {
				return target.copy(this.sphere);
			} else {
				return this.obb.getBoundingSphere(target);
			}
		}
	}
	const _vec3_1$2 = new t3d.Vector3();
	const _vec3_2 = new t3d.Vector3();
	const _vec3_3 = new t3d.Vector3();
	const _vec3_4 = new t3d.Vector3();

	/**
	 * Returns the file extension of the path component of a URL
	 * @param {string} url
	 * @returns {string} null if no extension found
	 */
	function getUrlExtension(url) {
		if (!url) {
			return null;
		}
		const filename = url.replace(/[a-z]+:\/\/[^/]+/i, '') // remove origin
		.replace(/\?.*$/i, '') // remove query
		.replace(/.*\//g, ''); // remove path

		const lastPeriod = filename.lastIndexOf('.');
		if (lastPeriod === -1) {
			return null;
		}
		return filename.substring(lastPeriod + 1) || null;
	}

	const _localRay = new t3d.Ray();
	const _vec$5 = new t3d.Vector3();
	const _hitArray = [];
	const _mat = new t3d.Matrix4();
	function distanceSort(a, b) {
		return a.distance - b.distance;
	}
	function intersectTileScene(tile, ray, renderer, intersects) {
		const {
			scene
		} = tile.cached;
		const lengthBefore = intersects.length;
		const didRaycast = renderer.invokeOnePlugin(plugin => plugin.raycastTile && plugin.raycastTile(tile, scene, ray, intersects));
		if (!didRaycast) {
			scene.traverse(c => {
				// We set the default raycast function to empty so t3d.js doesn't automatically cast against it
				Object.getPrototypeOf(c).raycast.call(c, ray, intersects);
			});
			const lengthAfter = intersects.length;

			// add the tile to intersects
			if (lengthAfter > lengthBefore) {
				for (let i = lengthBefore; i < lengthAfter; i++) {
					intersects[i].tile = tile;
				}
			}
		}
	}
	function intersectTileSceneFirstHit(tile, ray, renderer) {
		intersectTileScene(tile, ray, renderer, _hitArray);
		_hitArray.sort(distanceSort);
		const hit = _hitArray[0] || null;
		_hitArray.length = 0;
		return hit;
	}
	function isTileInitialized(tile) {
		return '__used' in tile;
	}

	// Returns the closest hit when traversing the tree
	function raycastTraverseFirstHit(renderer, tile, ray, localRay = null) {
		const {
			activeTiles
		} = renderer;

		// get the ray in the local group frame
		if (localRay === null) {
			localRay = _localRay;
			localRay.copy(ray).applyMatrix4(_mat.copy(renderer.worldMatrix).inverse());
		}

		// get a set of intersections so we intersect the nearest one first
		const array = [];
		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			const child = children[i];
			if (!isTileInitialized(child) || !child.__used) {
				continue;
			}

			// track the tile and hit distance for sorting
			const boundingVolume = child.cached.boundingVolume;
			if (boundingVolume.intersectRay(localRay, _vec$5) !== null) {
				_vec$5.applyMatrix4(renderer.worldMatrix);
				array.push({
					distance: _vec$5.distanceToSquared(ray.origin),
					tile: child
				});
			}
		}

		// sort them by ascending distance
		array.sort(distanceSort);

		// If the root is active make sure we've checked it
		let bestHit = null;
		let bestHitDistSq = Infinity;
		if (activeTiles.has(tile)) {
			const hit = intersectTileSceneFirstHit(tile, ray, renderer);
			if (hit) {
				bestHit = hit;
				bestHitDistSq = hit.distance * hit.distance;
			}
		}

		// traverse until we find the best hit and early out if a tile bounds
		// couldn't possible include a best hit
		for (let i = 0, l = array.length; i < l; i++) {
			const data = array[i];
			const boundingVolumeDistSq = data.distance;
			const tile = data.tile;
			if (boundingVolumeDistSq > bestHitDistSq) {
				break;
			}
			const hit = raycastTraverseFirstHit(renderer, tile, ray, localRay);
			if (hit) {
				const hitDistSq = hit.distance * hit.distance;
				if (hitDistSq < bestHitDistSq) {
					bestHit = hit;
					bestHitDistSq = hitDistSq;
				}
			}
		}
		return bestHit;
	}
	function raycastTraverse(renderer, tile, ray, intersects, localRay = null) {
		// if the tile has not been asynchronously initialized then there's no point in
		// traversing the tiles to check intersections.
		if (!isTileInitialized(tile)) {
			return;
		}
		const {
			activeTiles
		} = renderer;
		const {
			boundingVolume
		} = tile.cached;

		// get the ray in the local group frame
		if (localRay === null) {
			localRay = _localRay;
			localRay.copy(ray).applyMatrix4(_mat.copy(renderer.worldMatrix).inverse());
		}

		// exit early if the tile isn't used or the bounding volume is not intersected
		if (!tile.__used || !boundingVolume.intersectsRay(localRay)) {
			return;
		}

		// only intersect the tile geometry if it's active
		if (activeTiles.has(tile)) {
			intersectTileScene(tile, ray, renderer, intersects);
		}
		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			raycastTraverse(renderer, children[i], ray, intersects, localRay);
		}
	}

	class FastFrustum extends t3d.Frustum {
		constructor() {
			super();
			this.points = new Array(8).fill().map(() => new t3d.Vector3());
		}
		updateCache() {
			const {
				planes,
				points
			} = this;
			const planeIntersections = [[planes[0], planes[3], planes[4]],
			// Near top left
			[planes[1], planes[3], planes[4]],
			// Near top right
			[planes[0], planes[2], planes[4]],
			// Near bottom left
			[planes[1], planes[2], planes[4]],
			// Near bottom right
			[planes[0], planes[3], planes[5]],
			// Far top left
			[planes[1], planes[3], planes[5]],
			// Far top right
			[planes[0], planes[2], planes[5]],
			// Far bottom left
			[planes[1], planes[2], planes[5]] // Far bottom right
			];
			planeIntersections.forEach((planes, index) => {
				findIntersectionPoint(planes[0], planes[1], planes[2], points[index]);
			});
		}
	}
	const _mat3_1$1 = new t3d.Matrix3();

	// Solve a system of equations to find the point where the three planes intersect
	function findIntersectionPoint(plane1, plane2, plane3, target) {
		// Create the matrix A using the normals of the planes as rows
		const A = _mat3_1$1.set(plane1.normal.x, plane1.normal.y, plane1.normal.z, plane2.normal.x, plane2.normal.y, plane2.normal.z, plane3.normal.x, plane3.normal.y, plane3.normal.z);

		// Create the vector B using the constants of the planes
		target.set(-plane1.constant, -plane2.constant, -plane3.constant);

		// Solve for X by applying the inverse matrix to B
		target.applyMatrix3(A.inverse());
		return target;
	}

	class CameraList {
		constructor() {
			this._cameras = [];
			this._resolution = new t3d.Vector2();
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
					frustum: new FastFrustum(),
					// in origin space
					isOrthographic: false,
					sseDenominator: -1,
					// used if isOrthographic is false
					position: new t3d.Vector3(),
					// in origin space
					invScale: -1,
					pixelSize: 0 // used if isOrthographic is true
				});
			}

			// get inverse scale of origin matrix

			_mat4_1.copy(originMatrix).inverse();
			const invScaleX = _vec3_1$1.setFromMatrixColumn(_mat4_1, 0).getLength();
			const invScaleY = _vec3_1$1.setFromMatrixColumn(_mat4_1, 1).getLength();
			const invScaleZ = _vec3_1$1.setFromMatrixColumn(_mat4_1, 2).getLength();
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
					info.sseDenominator = 2 / projection[5] / cameraResolutionY;
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
	const _mat4_1 = new t3d.Matrix4();
	const _mat4_2 = new t3d.Matrix4();
	const _vec3_1$1 = new t3d.Vector3();

	const GIGABYTE_BYTES = 2 ** 30;
	class LRUCache {
		get unloadPriorityCallback() {
			return this._unloadPriorityCallback;
		}
		set unloadPriorityCallback(cb) {
			if (cb.length === 1) {
				console.warn('LRUCache: "unloadPriorityCallback" function has been changed to take two arguments.');
				this._unloadPriorityCallback = (a, b) => {
					const valA = cb(a);
					const valB = cb(b);
					if (valA < valB) return -1;
					if (valA > valB) return 1;
					return 0;
				};
			} else {
				this._unloadPriorityCallback = cb;
			}
		}
		constructor() {
			// options
			this.minSize = 6000;
			this.maxSize = 8000;
			this.minBytesSize = 0.3 * GIGABYTE_BYTES;
			this.maxBytesSize = 0.4 * GIGABYTE_BYTES;
			this.unloadPercent = 0.05;
			this.autoMarkUnused = true;

			// "itemSet" doubles as both the list of the full set of items currently
			// stored in the cache (keys) as well as a map to the time the item was last
			// used so it can be sorted appropriately.
			this.itemSet = new Map();
			this.itemList = [];
			this.usedSet = new Set();
			this.callbacks = new Map();
			this.unloadingHandle = -1;
			this.cachedBytes = 0;
			this.bytesMap = new Map();
			this.loadedSet = new Set();
			this._unloadPriorityCallback = null;
			this.computeMemoryUsageCallback = () => null;
			const itemSet = this.itemSet;
			this.defaultPriorityCallback = item => itemSet.get(item);
		}

		// Returns whether or not the cache has reached the maximum size
		isFull() {
			return this.itemSet.size >= this.maxSize || this.cachedBytes >= this.maxBytesSize;
		}
		getMemoryUsage(item) {
			return this.bytesMap.get(item) ?? null;
		}
		add(item, removeCb) {
			const itemSet = this.itemSet;
			if (itemSet.has(item)) {
				return false;
			}
			if (this.isFull()) {
				return false;
			}
			const usedSet = this.usedSet;
			const itemList = this.itemList;
			const callbacks = this.callbacks;
			const bytesMap = this.bytesMap;
			itemList.push(item);
			usedSet.add(item);
			itemSet.set(item, Date.now());
			callbacks.set(item, removeCb);

			// computeMemoryUsageCallback can return "null" if memory usage is not known, yet
			const bytes = this.computeMemoryUsageCallback(item);
			this.cachedBytes += bytes || 0;
			bytesMap.set(item, bytes);
			return true;
		}
		has(item) {
			return this.itemSet.has(item);
		}
		remove(item) {
			const usedSet = this.usedSet;
			const itemSet = this.itemSet;
			const itemList = this.itemList;
			const bytesMap = this.bytesMap;
			const callbacks = this.callbacks;
			const loadedSet = this.loadedSet;
			if (itemSet.has(item)) {
				this.cachedBytes -= bytesMap.get(item) || 0;
				bytesMap.delete(item);
				callbacks.get(item)(item);
				const index = itemList.indexOf(item);
				itemList.splice(index, 1);
				usedSet.delete(item);
				itemSet.delete(item);
				callbacks.delete(item);
				loadedSet.delete(item);
				return true;
			}
			return false;
		}

		// Marks whether tiles in the cache have been completely loaded or not. Tiles that have not been completely
		// loaded are subject to being disposed early if the cache is full above its max size limits, even if they
		// are marked as used.
		setLoaded(item, value) {
			const {
				itemSet,
				loadedSet
			} = this;
			if (itemSet.has(item)) {
				if (value === true) {
					loadedSet.add(item);
				} else {
					loadedSet.delete(item);
				}
			}
		}
		updateMemoryUsage(item) {
			const itemSet = this.itemSet;
			const bytesMap = this.bytesMap;
			if (!itemSet.has(item)) {
				return;
			}
			this.cachedBytes -= bytesMap.get(item) || 0;
			const bytes = this.computeMemoryUsageCallback(item);
			bytesMap.set(item, bytes);
			this.cachedBytes += bytes;
		}
		markUsed(item) {
			const itemSet = this.itemSet;
			const usedSet = this.usedSet;
			if (itemSet.has(item) && !usedSet.has(item)) {
				itemSet.set(item, Date.now());
				usedSet.add(item);
			}
		}
		markUnused(item) {
			this.usedSet.delete(item);
		}
		markAllUnused() {
			this.usedSet.clear();
		}

		// TODO: this should be renamed because it's not necessarily unloading all unused content
		// Maybe call it "cleanup" or "unloadToMinSize"
		unloadUnusedContent() {
			const {
				unloadPercent,
				minSize,
				maxSize,
				itemList,
				itemSet,
				usedSet,
				loadedSet,
				callbacks,
				bytesMap,
				minBytesSize,
				maxBytesSize
			} = this;
			const unused = itemList.length - usedSet.size;
			const unloaded = itemList.length - loadedSet.size;
			const excessNodes = Math.max(Math.min(itemList.length - minSize, unused), 0);
			const excessBytes = this.cachedBytes - minBytesSize;
			const unloadPriorityCallback = this.unloadPriorityCallback || this.defaultPriorityCallback;
			let needsRerun = false;
			const hasNodesToUnload = excessNodes > 0 && unused > 0 || unloaded && itemList.length > maxSize;
			const hasBytesToUnload = unused && this.cachedBytes > minBytesSize || unloaded && this.cachedBytes > maxBytesSize;
			if (hasBytesToUnload || hasNodesToUnload) {
				// used items should be at the end of the array, "unloaded" items in the middle of the array
				itemList.sort((a, b) => {
					const usedA = usedSet.has(a);
					const usedB = usedSet.has(b);
					if (usedA === usedB) {
						const loadedA = loadedSet.has(a);
						const loadedB = loadedSet.has(b);
						if (loadedA === loadedB) {
							// Use the sort function otherwise
							// higher priority should be further to the left
							return -unloadPriorityCallback(a, b);
						} else {
							return loadedA ? 1 : -1;
						}
					} else {
						// If one is used and the other is not move the used one towards the end of the array
						return usedA ? 1 : -1;
					}
				});

				// address corner cases where the minSize might be zero or smaller than maxSize - minSize,
				// which would result in a very small or no items being unloaded.
				const maxUnload = Math.max(minSize * unloadPercent, excessNodes * unloadPercent);
				const nodesToUnload = Math.ceil(Math.min(maxUnload, unused, excessNodes));
				const maxBytesUnload = Math.max(unloadPercent * excessBytes, unloadPercent * minBytesSize);
				const bytesToUnload = Math.min(maxBytesUnload, excessBytes);
				let removedNodes = 0;
				let removedBytes = 0;

				// evict up to the max node or bytes size, keeping one more item over the max bytes limit
				// so the "full" function behaves correctly.
				while (this.cachedBytes - removedBytes > maxBytesSize || itemList.length - removedNodes > maxSize) {
					const item = itemList[removedNodes];
					const bytes = bytesMap.get(item) || 0;
					if (usedSet.has(item) && loadedSet.has(item) || this.cachedBytes - removedBytes - bytes < maxBytesSize && itemList.length - removedNodes <= maxSize) {
						break;
					}
					removedBytes += bytes;
					removedNodes++;
				}

				// evict up to the min node or bytes size, keeping one more item over the min bytes limit
				// so we're meeting it
				while (removedBytes < bytesToUnload || removedNodes < nodesToUnload) {
					const item = itemList[removedNodes];
					const bytes = bytesMap.get(item) || 0;
					if (usedSet.has(item) || this.cachedBytes - removedBytes - bytes < minBytesSize && removedNodes >= nodesToUnload) {
						break;
					}
					removedBytes += bytes;
					removedNodes++;
				}

				// remove the nodes
				itemList.splice(0, removedNodes).forEach(item => {
					this.cachedBytes -= bytesMap.get(item) || 0;
					callbacks.get(item)(item);
					bytesMap.delete(item);
					itemSet.delete(item);
					callbacks.delete(item);
					loadedSet.delete(item);
					usedSet.delete(item);
				});

				// if we didn't remove enough nodes or we still have excess bytes and there are nodes to removed
				// then we want to fire another round of unloading
				needsRerun = removedNodes < excessNodes || removedBytes < excessBytes && removedNodes < unused;
				needsRerun = needsRerun && removedNodes > 0;
			}
			if (needsRerun) {
				this.unloadingHandle = requestAnimationFrame(() => this.scheduleUnload());
			}
		}
		scheduleUnload() {
			cancelAnimationFrame(this.unloadingHandle);
			if (!this.scheduled) {
				this.scheduled = true;
				queueMicrotask(() => {
					this.scheduled = false;
					this.unloadUnusedContent();
				});
			}
		}
	}

	class PriorityQueue {
		// returns whether tasks are queued or actively running
		get running() {
			return this.items.length !== 0 || this.currJobs !== 0;
		}
		constructor() {
			// options
			this.maxJobs = 6;
			this.items = [];
			this.callbacks = new Map();
			this.currJobs = 0;
			this.scheduled = false;
			this.autoUpdate = true;
			this.priorityCallback = () => {
				throw new Error('PriorityQueue: PriorityCallback function not defined.');
			};

			// Customizable scheduling callback. Default using requestAnimationFrame()
			this.schedulingCallback = func => {
				requestAnimationFrame(func);
			};
			this._runjobs = () => {
				this.scheduled = false;
				this.tryRunJobs();
			};
		}
		sort() {
			const priorityCallback = this.priorityCallback;
			const items = this.items;
			items.sort(priorityCallback);
		}
		has(item) {
			return this.callbacks.has(item);
		}
		add(item, callback) {
			const data = {
				callback,
				reject: null,
				resolve: null,
				promise: null
			};
			data.promise = new Promise((resolve, reject) => {
				const items = this.items;
				const callbacks = this.callbacks;
				data.resolve = resolve;
				data.reject = reject;
				items.push(item);
				callbacks.set(item, data);
				if (this.autoUpdate) {
					this.scheduleJobRun();
				}
			});
			return data.promise;
		}
		remove(item) {
			const items = this.items;
			const callbacks = this.callbacks;
			const index = items.indexOf(item);
			if (index !== -1) {
				// reject the promise to ensure there are no dangling promises - add a
				// catch here to handle the case where the promise was never used anywhere
				// else.
				const info = callbacks.get(item);
				info.promise.catch(() => {});
				info.reject(new Error('PriorityQueue: Item removed.'));
				items.splice(index, 1);
				callbacks.delete(item);
			}
		}
		tryRunJobs() {
			this.sort();
			const items = this.items;
			const callbacks = this.callbacks;
			const maxJobs = this.maxJobs;
			let iterated = 0;
			const completedCallback = () => {
				this.currJobs--;
				if (this.autoUpdate) {
					this.scheduleJobRun();
				}
			};
			while (maxJobs > this.currJobs && items.length > 0 && iterated < maxJobs) {
				this.currJobs++;
				iterated++;
				const item = items.pop();
				const {
					callback,
					resolve,
					reject
				} = callbacks.get(item);
				callbacks.delete(item);
				let result;
				try {
					result = callback(item);
				} catch (err) {
					reject(err);
					completedCallback();
				}
				if (result instanceof Promise) {
					result.then(resolve).catch(reject).finally(completedCallback);
				} else {
					resolve(result);
					completedCallback();
				}
			}
		}
		scheduleJobRun() {
			if (!this.scheduled) {
				this.schedulingCallback(this._runjobs);
				this.scheduled = true;
			}
		}
	}

	// FAILED is negative so lru cache priority sorting will unload it first
	const FAILED = -1;
	const UNLOADED = 0;
	const LOADING = 1;
	const PARSING = 2;
	const LOADED = 3;

	// https://en.wikipedia.org/wiki/World_Geodetic_System
	// https://en.wikipedia.org/wiki/Flattening
	const WGS84_RADIUS = 6378137;
	const WGS84_HEIGHT = 6356752.314245179;

	const viewErrorTarget$1 = {
		inView: false,
		error: Infinity,
		distance: Infinity
	};
	function isDownloadFinished(value) {
		return value === LOADED || value === FAILED;
	}

	// Checks whether this tile was last used on the given frame.
	function isUsedThisFrame(tile, frameCount) {
		return tile.__lastFrameVisited === frameCount && tile.__used;
	}
	function areChildrenProcessed(tile) {
		return tile.__childrenProcessed === tile.children.length;
	}

	// Resets the frame frame information for the given tile
	function resetFrameState(tile, renderer) {
		if (tile.__lastFrameVisited !== renderer.frameCount) {
			tile.__lastFrameVisited = renderer.frameCount;
			tile.__used = false;
			tile.__inFrustum = false;
			tile.__isLeaf = false;
			tile.__visible = false;
			tile.__active = false;
			tile.__error = Infinity;
			tile.__distanceFromCamera = Infinity;
			tile.__childrenWereVisible = false;
			tile.__allChildrenLoaded = false;

			// update tile frustum and error state
			renderer.calculateTileViewError(tile, viewErrorTarget$1);
			tile.__inFrustum = viewErrorTarget$1.inView;
			tile.__error = viewErrorTarget$1.error;
			tile.__distanceFromCamera = viewErrorTarget$1.distance;
		}
	}

	// Recursively mark tiles used down to the next tile with content
	function recursivelyMarkUsed(tile, renderer) {
		renderer.ensureChildrenArePreprocessed(tile);
		resetFrameState(tile, renderer);
		markUsed(tile, renderer);

		// don't traverse if the children have not been processed, yet
		if (!tile.__hasRenderableContent && areChildrenProcessed(tile)) {
			const children = tile.children;
			for (let i = 0, l = children.length; i < l; i++) {
				recursivelyMarkUsed(children[i], renderer);
			}
		}
	}

	// Recursively traverses to the next tiles with unloaded renderable content to load them
	function recursivelyLoadNextRenderableTiles(tile, renderer) {
		renderer.ensureChildrenArePreprocessed(tile);

		// exit the recursion if the tile hasn't been used this frame
		if (isUsedThisFrame(tile, renderer.frameCount)) {
			// queue this tile to download content
			if (tile.__hasContent && tile.__loadingState === UNLOADED && !renderer.lruCache.isFull()) {
				renderer.queueTileForDownload(tile);
			}
			if (areChildrenProcessed(tile)) {
				// queue any used child tiles
				const children = tile.children;
				for (let i = 0, l = children.length; i < l; i++) {
					recursivelyLoadNextRenderableTiles(children[i], renderer);
				}
			}
		}
	}

	// Mark a tile as being used by current view
	function markUsed(tile, renderer) {
		if (tile.__used) {
			return;
		}
		tile.__used = true;
		renderer.markTileUsed(tile);
		renderer.stats.used++;
		if (tile.__inFrustum === true) {
			renderer.stats.inFrustum++;
		}
	}

	// Returns whether the tile can be traversed to the next layer of children by checking the tile metrics
	function canTraverse(tile, renderer) {
		// If we've met the error requirements then don't load further
		if (tile.__error <= renderer.errorTarget) {
			return false;
		}

		// Early out if we've reached the maximum allowed depth.
		if (renderer.maxDepth > 0 && tile.__depth + 1 >= renderer.maxDepth) {
			return false;
		}

		// Early out if the children haven't been processed, yet
		if (!areChildrenProcessed(tile)) {
			return false;
		}
		return true;
	}

	// Helper function for traversing a tile set. If `beforeCb` returns `true` then the
	// traversal will end early.
	function traverseSet(tile, beforeCb = null, afterCb = null) {
		const stack = [];

		// A stack-based, depth-first traversal, storing
		// triplets (tile, parent, depth) in the stack array.

		stack.push(tile);
		stack.push(null);
		stack.push(0);
		while (stack.length > 0) {
			const depth = stack.pop();
			const parent = stack.pop();
			const tile = stack.pop();
			if (beforeCb && beforeCb(tile, parent, depth)) {
				if (afterCb) {
					afterCb(tile, parent, depth);
				}
				return;
			}
			const children = tile.children;

			// Children might be undefined if the tile has not been preprocessed yet
			if (children) {
				for (let i = children.length - 1; i >= 0; i--) {
					stack.push(children[i]);
					stack.push(tile);
					stack.push(depth + 1);
				}
			}
			if (afterCb) {
				afterCb(tile, parent, depth);
			}
		}
	}
	function markUsedTiles(tile, renderer) {
		// determine frustum set is run first so we can ensure the preprocessing of all the necessary
		// child tiles has happened here.
		renderer.ensureChildrenArePreprocessed(tile);
		resetFrameState(tile, renderer);
		if (!tile.__inFrustum) {
			return;
		}
		if (!canTraverse(tile, renderer)) {
			markUsed(tile, renderer);
			return;
		}

		// Traverse children and see if any children are in view.
		let anyChildrenUsed = false;
		let anyChildrenInFrustum = false;
		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			const c = children[i];
			markUsedTiles(c, renderer);
			anyChildrenUsed = anyChildrenUsed || isUsedThisFrame(c, renderer.frameCount);
			anyChildrenInFrustum = anyChildrenInFrustum || c.__inFrustum;
		}

		// Disabled for now because this will cause otherwise unused children to be added to the lru cache
		// if none of the children are in the frustum then this tile shouldn't be displayed.
		// Otherwise this can cause load oscillation as parents are traversed and loaded and then determined
		// to not be used because children aren't visible. See #1165.
		// if (tile.refine === 'REPLACE' && !anyChildrenInFrustum && children.length !== 0 && !tile.__hasUnrenderableContent) {
		// 	// TODO: we're not checking tiles with unrenderable content here since external tile sets might look like they're in the frustum,
		// 	// load the children, then the children indicate that it's not visible, causing it to be unloaded. Then it will be loaded again.
		// 	// The impact when including external tile set roots in the check is more significant but can't be used unless we keep external tile
		// 	// sets around even when they're not needed. See issue #741.

		// 	// TODO: what if we mark the tile as not in the frustum but we _do_ mark it as used? Then we can stop frustum traversal and at least
		// 	// prevent tiles from rendering unless they're needed.
		// 	console.log('FAILED');
		// 	tile.__inFrustum = false;
		// 	return;
		// }

		markUsed(tile, renderer);

		// If this is a tile that needs children loaded to refine then recursively load child
		// tiles until error is met
		if (anyChildrenUsed && tile.refine === 'REPLACE') {
			for (let i = 0, l = children.length; i < l; i++) {
				const c = children[i];
				recursivelyMarkUsed(c, renderer);
			}
		}
	}

	// Traverse and mark the tiles that are at the leaf nodes of the "used" tree.
	function markUsedSetLeaves(tile, renderer) {
		const frameCount = renderer.frameCount;
		if (!isUsedThisFrame(tile, frameCount)) {
			return;
		}

		// This tile is a leaf if none of the children had been used.
		const children = tile.children;
		let anyChildrenUsed = false;
		for (let i = 0, l = children.length; i < l; i++) {
			const c = children[i];
			anyChildrenUsed = anyChildrenUsed || isUsedThisFrame(c, frameCount);
		}
		if (!anyChildrenUsed) {
			tile.__isLeaf = true;
		} else {
			let childrenWereVisible = false;
			let allChildrenLoaded = true;
			for (let i = 0, l = children.length; i < l; i++) {
				const c = children[i];
				markUsedSetLeaves(c, renderer);
				childrenWereVisible = childrenWereVisible || c.__wasSetVisible || c.__childrenWereVisible;
				if (isUsedThisFrame(c, frameCount)) {
					// consider a child to be loaded if
					// - the children's children have been loaded
					// - the tile content has loaded
					// - the tile is completely empty - ie has no children and no content
					// - the child tile set has tried to load but failed
					const childLoaded = c.__allChildrenLoaded || c.__hasRenderableContent && isDownloadFinished(c.__loadingState) || c.__hasUnrenderableContent && c.__loadingState === FAILED;
					allChildrenLoaded = allChildrenLoaded && childLoaded;
				}
			}
			tile.__childrenWereVisible = childrenWereVisible;
			tile.__allChildrenLoaded = allChildrenLoaded;
		}
	}

	// Skip past tiles we consider unrenderable because they are outside the error threshold.
	function markVisibleTiles(tile, renderer) {
		const stats = renderer.stats;
		if (!isUsedThisFrame(tile, renderer.frameCount)) {
			return;
		}

		// Request the tile contents or mark it as visible if we've found a leaf.
		const lruCache = renderer.lruCache;
		if (tile.__isLeaf) {
			if (tile.__loadingState === LOADED) {
				if (tile.__inFrustum) {
					tile.__visible = true;
					stats.visible++;
				}
				tile.__active = true;
				stats.active++;
			} else if (!lruCache.isFull() && tile.__hasContent) {
				renderer.queueTileForDownload(tile);
			}
			return;
		}
		const children = tile.children;
		const hasContent = tile.__hasContent;
		const loadedContent = isDownloadFinished(tile.__loadingState) && hasContent;
		const errorRequirement = (renderer.errorTarget + 1) * renderer.errorThreshold;
		const meetsSSE = tile.__error <= errorRequirement;
		const childrenWereVisible = tile.__childrenWereVisible;

		// NOTE: We can "trickle" root tiles in by enabling these lines.
		// Don't wait for all children tiles to load if this tile set has empty tiles at the root
		// const emptyRootTile = tile.__depthFromRenderedParent === 0;
		// const allChildrenLoaded = tile.__allChildrenLoaded || emptyRootTile;

		// If we've met the SSE requirements and we can load content then fire a fetch.
		const allChildrenLoaded = tile.__allChildrenLoaded;
		const includeTile = meetsSSE || tile.refine === 'ADD';
		if (includeTile && !loadedContent && !lruCache.isFull() && hasContent) {
			renderer.queueTileForDownload(tile);
		}

		// Only mark this tile as visible if it meets the screen space error requirements, has loaded content, not
		// all children have loaded yet, and if no children were visible last frame. We want to keep children visible
		// that _were_ visible to avoid a pop in level of detail as the camera moves around and parent / sibling tiles
		// load in.

		// Skip the tile entirely if there's no content to load
		if (meetsSSE && !allChildrenLoaded && !childrenWereVisible && loadedContent || tile.refine === 'ADD' && loadedContent) {
			if (tile.__inFrustum) {
				tile.__visible = true;
				stats.visible++;
			}
			tile.__active = true;
			stats.active++;
		}

		// If we're additive then don't stop the traversal here because it doesn't matter whether the children load in
		// at the same rate.
		if (tile.refine === 'REPLACE' && meetsSSE && !allChildrenLoaded) {
			// load the child content if we've found that we've been loaded so we can move down to the next tile
			// layer when the data has loaded.
			for (let i = 0, l = children.length; i < l; i++) {
				const c = children[i];
				if (isUsedThisFrame(c, renderer.frameCount)) {
					recursivelyLoadNextRenderableTiles(c, renderer);
				}
			}
		} else {
			for (let i = 0, l = children.length; i < l; i++) {
				markVisibleTiles(children[i], renderer);
			}
		}
	}

	// Final traverse to toggle tile visibility.
	const toggleTiles = (tile, renderer) => {
		const isUsed = isUsedThisFrame(tile, renderer.frameCount);
		if (isUsed || tile.__usedLastFrame) {
			let setActive = false;
			let setVisible = false;
			if (isUsed) {
				// enable visibility if active due to shadows
				setActive = tile.__active;
				if (renderer.displayActiveTiles) {
					setVisible = tile.__active || tile.__visible;
				} else {
					setVisible = tile.__visible;
				}
			} else {
				// if the tile was used last frame but not this one then there's potential for the tile
				// to not have been visited during the traversal, meaning it hasn't been reset and has
				// stale values. This ensures the values are not stale.
				resetFrameState(tile, renderer);
			}

			// If the active or visible state changed then call the functions.
			if (tile.__hasRenderableContent && tile.__loadingState === LOADED) {
				if (tile.__wasSetActive !== setActive) {
					renderer.invokeOnePlugin(plugin => plugin.setTileActive && plugin.setTileActive(tile, setActive));
				}
				if (tile.__wasSetVisible !== setVisible) {
					renderer.invokeOnePlugin(plugin => plugin.setTileVisible && plugin.setTileVisible(tile, setVisible));
				}
			}
			tile.__wasSetActive = setActive;
			tile.__wasSetVisible = setVisible;
			tile.__usedLastFrame = isUsed;
			const children = tile.children;
			for (let i = 0, l = children.length; i < l; i++) {
				const c = children[i];
				toggleTiles(c, renderer);
			}
		}
	};

	/**
	 * Traverses the ancestry of the tile up to the root tile.
	 */
	function traverseAncestors(tile, callback = null) {
		let current = tile;
		while (current) {
			const depth = current.__depth;
			const parent = current.parent;
			if (callback) {
				callback(current, parent, depth);
			}
			current = parent;
		}
	}

	const WGS84_ELLIPSOID = new Ellipsoid(new t3d.Vector3(WGS84_RADIUS, WGS84_RADIUS, WGS84_HEIGHT));
	WGS84_ELLIPSOID.name = 'WGS84 Earth';

	// function that rate limits the amount of time a function can be called to once
	// per frame, initially queuing a new call for the next frame.
	function throttle(callback) {
		let handle = null;
		return () => {
			if (handle === null) {
				handle = requestAnimationFrame(() => {
					handle = null;
					callback();
				});
			}
		};
	}

	class ImageBitmapLoader extends t3d.Loader {
		constructor(manager) {
			super(manager);
			if (typeof createImageBitmap === 'undefined') {
				console.warn('ImageBitmapLoader: createImageBitmap() not supported.');
			}
			if (typeof fetch === 'undefined') {
				console.warn('ImageBitmapLoader: fetch() not supported.');
			}
			this.options = {
				premultiplyAlpha: 'none'
			};
		}
		setOptions(options) {
			this.options = options;
			return this;
		}
		load(url, onLoad, _onProgress, onError) {
			if (url === undefined) url = '';
			if (this.path !== undefined) url = this.path + url;
			url = this.manager.resolveURL(url);
			const scope = this;
			const fetchOptions = {};
			fetchOptions.credentials = this.crossOrigin === 'anonymous' ? 'same-origin' : 'include';
			fetchOptions.headers = this.requestHeader;
			fetch(url, fetchOptions).then(function (res) {
				return res.blob();
			}).then(function (blob) {
				return createImageBitmap(blob, Object.assign(scope.options, {
					colorSpaceConversion: 'none'
				}));
			}).then(function (imageBitmap) {
				if (onLoad) onLoad(imageBitmap);
				scope.manager.itemEnd(url);
			}).catch(function (e) {
				if (onError) onError(e);
				scope.manager.itemError(url);
				scope.manager.itemEnd(url);
			});
			scope.manager.itemStart(url);
		}
	}

	const _vec4_1 = new t3d.Vector4();
	class GLTFUtils {
		constructor() {}
		static extractUrlBase(url) {
			const parts = url.split('/');
			parts.pop();
			return (parts.length < 1 ? '.' : parts.join('/')) + '/';
		}

		// url: aa.bin ;	path:example/resource/model/		 (for example)example/resource/model/aa.bin
		static resolveURL(url, path) {
			// Invalid URL
			if (typeof url !== 'string' || url === '') return '';

			// Absolute URL http://,https://,//
			if (/^(https?:)?\/\//i.test(url)) return url;

			// Data URI
			if (/^data:/i.test(url)) return url;

			// Blob URL
			if (/^blob:/i.test(url)) return url;

			// Relative URL
			return path + url;
		}
		static decodeText(array) {
			if (typeof TextDecoder !== 'undefined') {
				return new TextDecoder().decode(array);
			}

			// Avoid the String.fromCharCode.apply(null, array) shortcut, which
			// throws a "maximum call stack size exceeded" error for large arrays.

			let s = '';
			for (let i = 0, il = array.length; i < il; i++) {
				// Implicitly assumes little-endian.
				s += String.fromCharCode(array[i]);
			}
			try {
				// merges multi-byte utf-8 characters.

				return decodeURIComponent(escape(s));
			} catch (e) {
				// see #16358
				return s;
			}
		}
		static parseGLB(glb) {
			const UINT32_LENGTH = 4;
			const GLB_HEADER_MAGIC = 0x46546C67; // 'glTF'
			const GLB_HEADER_LENGTH = 12;
			const GLB_CHUNK_TYPES = {
				JSON: 0x4E4F534A,
				BIN: 0x004E4942
			};
			const dataView = new DataView(glb);
			const header = {
				magic: dataView.getUint32(0, true),
				version: dataView.getUint32(UINT32_LENGTH, true),
				length: dataView.getUint32(2 * UINT32_LENGTH, true)
			};
			if (header.magic !== GLB_HEADER_MAGIC) {
				console.error('Invalid glb magic number. Expected 0x46546C67, found 0x' + header.magic.toString(16));
				return null;
			} else if (header.version < 2.0) {
				console.error('GLTFLoader: Legacy binary file detected.');
			}
			let chunkLength = dataView.getUint32(GLB_HEADER_LENGTH, true);
			let chunkType = dataView.getUint32(GLB_HEADER_LENGTH + UINT32_LENGTH, true);
			if (chunkType !== GLB_CHUNK_TYPES.JSON) {
				console.error('Invalid glb chunk type. Expected 0x4E4F534A, found 0x' + chunkType.toString(16));
				return null;
			}
			const glTFData = new Uint8Array(glb, GLB_HEADER_LENGTH + 2 * UINT32_LENGTH, chunkLength);
			const gltf = JSON.parse(GLTFUtils.decodeText(glTFData));
			const buffers = [];
			let byteOffset = GLB_HEADER_LENGTH + 2 * UINT32_LENGTH + chunkLength;
			while (byteOffset < header.length) {
				chunkLength = dataView.getUint32(byteOffset, true);
				chunkType = dataView.getUint32(byteOffset + UINT32_LENGTH, true);
				if (chunkType !== GLB_CHUNK_TYPES.BIN) {
					console.error('Invalid glb chunk type. Expected 0x004E4942, found 0x' + chunkType.toString(16));
					return null;
				}
				const currentOffset = byteOffset + 2 * UINT32_LENGTH;
				const buffer = glb.slice(currentOffset, currentOffset + chunkLength);
				buffers.push(buffer);
				byteOffset += chunkLength + 2 * UINT32_LENGTH;
			}
			return {
				gltf,
				buffers
			};
		}

		// Reference: github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_mesh_quantization#encoding-quantized-data
		static getNormalizedComponentScale(constructor) {
			if (constructor === Int8Array) {
				return 1 / 127;
			} else if (constructor === Uint8Array) {
				return 1 / 255;
			} else if (constructor === Int16Array) {
				return 1 / 32767;
			} else if (constructor === Uint16Array) {
				return 1 / 65535;
			} else {
				throw new Error('Unsupported normalized accessor component type.');
			}
		}
		static normalizeSkinWeights(skinWeight) {
			const offset = skinWeight.offset;
			const buffer = skinWeight.buffer;
			const stride = buffer.stride;
			for (let i = 0, l = buffer.count; i < l; i++) {
				_vec4_1.fromArray(buffer.array, i * stride + offset);
				const scale = 1.0 / _vec4_1.getManhattanLength();
				if (scale !== Infinity) {
					_vec4_1.multiplyScalar(scale);
				} else {
					_vec4_1.set(1, 0, 0, 0); // do something reasonable
				}
				_vec4_1.toArray(buffer.array, i * stride + offset);
			}
		}
	}

	let IndexParser$1 = class IndexParser {
		static parse(context, loader) {
			const {
				url
			} = context;
			return loader.loadFile(url, 'arraybuffer').then(data => {
				const magic = GLTFUtils.decodeText(new Uint8Array(data, 0, 4));
				if (magic === 'glTF') {
					const glbData = GLTFUtils.parseGLB(data);
					context.gltf = glbData.gltf;
					context.buffers = glbData.buffers;
				} else {
					const gltfString = GLTFUtils.decodeText(new Uint8Array(data));
					context.gltf = JSON.parse(gltfString);
				}
			});
		}
	};

	// Marks the special nodes/meshes in json for efficient parse.
	class ReferenceParser {
		static parse(context, loader) {
			const {
				gltf,
				path
			} = context;
			const {
				nodes = [],
				skins = [],
				meshes = [],
				buffers,
				images
			} = gltf;

			// Nothing in the node definition indicates whether it is a Bone or an
			// Object3D. Use the skins' joint references to mark bones.
			skins.forEach(skin => {
				const {
					joints = []
				} = skin;
				joints.forEach(joint => {
					nodes[joint].isBone = true;
				});
			});

			// Nothing in the mesh definition indicates whether it is
			// a SkinnedMesh or Mesh. Use the node's mesh reference
			// to mark SkinnedMesh if node has skin.
			nodes.forEach(node => {
				if (node.mesh !== undefined) {
					if (node.skin !== undefined) {
						meshes[node.mesh].isSkinned = true;
					}
				}
			});

			// setup loading list for detail load progress
			if (loader.detailLoadProgress) {
				const loadItems = new Set();
				if (buffers) {
					buffers.forEach(buffer => {
						if (!buffer.uri) {
							// glb or other
							return;
						}
						const bufferUrl = GLTFUtils.resolveURL(buffer.uri, path);
						loadItems.add(bufferUrl);
					});
				}
				if (images) {
					images.forEach((image, index) => {
						const {
							uri,
							bufferView: bufferViewIndex
						} = image;
						let imageUrl = uri;
						if (bufferViewIndex !== undefined) {
							imageUrl = 'blob<' + index + '>'; // fake url for blob image
						}
						imageUrl = GLTFUtils.resolveURL(imageUrl, path);
						loadItems.add(imageUrl);
					});
				}
				loadItems.forEach(item => loader.manager.itemStart(item));
				context.loadItems = loadItems;
			}
		}
	}

	class Validator {
		static parse(context) {
			const {
				gltf: {
					asset: {
						version
					}
				}
			} = context;
			const gltfVersion = Number(version);
			if (!(gltfVersion >= 2 && gltfVersion < 3)) {
				throw 'Only support gltf 2.x.';
			}
		}
	}

	class BufferParser {
		static parse(context, loader) {
			const {
				gltf,
				loadItems
			} = context;
			if (context.buffers !== null) {
				// buffers have been parsed
				return null;
			} else {
				return Promise.all(gltf.buffers.map(buffer => {
					const bufferUrl = GLTFUtils.resolveURL(buffer.uri, context.path);
					if (loader.detailLoadProgress) {
						loadItems.delete(bufferUrl);
					}
					const promise = loader.loadFile(bufferUrl, 'arraybuffer').then(buffer => {
						if (loader.detailLoadProgress) {
							loader.manager.itemEnd(bufferUrl);
						}
						return buffer;
					});
					if (loader.detailLoadProgress) {
						promise.catch(() => loader.manager.itemEnd(bufferUrl));
					}
					return promise;
				})).then(buffers => {
					context.buffers = buffers;
				});
			}
		}
	}

	class BufferViewParser {
		static parse(context, loader) {
			const {
				buffers,
				gltf
			} = context;
			if (!gltf.bufferViews) return;
			const meshoptExt = loader.extensions.get('EXT_meshopt_compression');
			return Promise.all(gltf.bufferViews.map(bufferView => {
				const {
					buffer,
					byteOffset = 0,
					byteLength = 0
				} = bufferView;
				if (bufferView.extensions) {
					const {
						EXT_meshopt_compression
					} = bufferView.extensions;
					if (EXT_meshopt_compression && meshoptExt) {
						return meshoptExt.loadBufferView(EXT_meshopt_compression, buffers, loader.getMeshoptDecoder());
					}
				}
				const arrayBuffer = buffers[buffer];
				return arrayBuffer.slice(byteOffset, byteOffset + byteLength);
			})).then(bufferViews => {
				context.bufferViews = bufferViews;
			});
		}
	}

	class ImageParser {
		static parse(context, loader) {
			const {
				gltf,
				bufferViews,
				path,
				loadItems
			} = context;
			if (!gltf.images) return;
			const basisuExt = loader.extensions.get('KHR_texture_basisu');
			return Promise.all(gltf.images.map((params, index) => {
				const {
					uri,
					bufferView: bufferViewIndex,
					mimeType,
					name: imageName
				} = params;
				let isObjectURL = false;
				let sourceUrl = uri || '';
				if (bufferViewIndex !== undefined) {
					const bufferViewData = bufferViews[bufferViewIndex];
					const blob = new Blob([bufferViewData], {
						type: mimeType
					});
					sourceUrl = URL.createObjectURL(blob);
					isObjectURL = true;
				}
				const imageUrl = GLTFUtils.resolveURL(sourceUrl, path);
				if (loader.detailLoadProgress) {
					loadItems.delete(imageUrl);
				}
				let promise;
				if (mimeType && mimeType.includes('ktx2') && basisuExt) {
					promise = basisuExt.loadTextureData(imageUrl, loader.getKTX2Loader()).then(transcodeResult => {
						if (loader.detailLoadProgress) {
							if (isObjectURL) {
								loader.manager.itemEnd(GLTFUtils.resolveURL('blob<' + index + '>', path));
							} else {
								loader.manager.itemEnd(imageUrl);
							}
						}
						return transcodeResult;
					});
				} else {
					const param = {
						loader,
						imageUrl,
						imageName,
						isObjectURL,
						sourceUrl,
						index,
						path
					};
					if (mimeType && (mimeType.includes('avif') || mimeType.includes('webp'))) {
						promise = detectSupport(mimeType).then(isSupported => {
							if (isSupported) return loadImage(param);
							throw new Error('GLTFLoader: WebP or AVIF required by asset but unsupported.');
						});
					} else {
						return loadImage(param);
					}
				}
				if (loader.detailLoadProgress) {
					promise.catch(() => loader.manager.itemEnd(imageUrl));
				}
				return promise;
			})).then(images => {
				context.images = images;
			});
		}
	}
	function detectSupport(mimeType) {
		const isSupported = new Promise(resolve => {
			// Lossy test image.
			const image = new Image();
			if (mimeType.includes('avif')) {
				image.src = 'data:image/avif;base64,AAAAIGZ0eXBhdmlmAAAAAGF2aWZtaWYxbWlhZk1BMUIAAADybWV0YQAAAAAAAAAoaGRscgAAAAAAAAAAcGljdAAAAAAAAAAAAAAAAGxpYmF2aWYAAAAADnBpdG0AAAAAAAEAAAAeaWxvYwAAAABEAAABAAEAAAABAAABGgAAABcAAAAoaWluZgAAAAAAAQAAABppbmZlAgAAAAABAABhdjAxQ29sb3IAAAAAamlwcnAAAABLaXBjbwAAABRpc3BlAAAAAAAAAAEAAAABAAAAEHBpeGkAAAAAAwgICAAAAAxhdjFDgQAMAAAAABNjb2xybmNseAACAAIABoAAAAAXaXBtYQAAAAAAAAABAAEEAQKDBAAAAB9tZGF0EgAKCBgABogQEDQgMgkQAAAAB8dSLfI=';
			} else {
				image.src = 'data:image/webp;base64,UklGRiIAAABXRUJQVlA4IBYAAAAwAQCdASoBAAEADsD+JaQAA3AAAAAA';
			}
			image.onload = () => {
				resolve(image.height === 1);
			};
		});
		return isSupported;
	}
	function loadImage(param) {
		const {
			loader,
			imageUrl,
			imageName,
			isObjectURL,
			sourceUrl,
			index,
			path
		} = param;
		const promise = loader.loadImage(imageUrl).then(image => {
			image.__name = imageName;
			if (isObjectURL === true) {
				URL.revokeObjectURL(sourceUrl);
			}
			if (loader.detailLoadProgress) {
				if (isObjectURL) {
					loader.manager.itemEnd(GLTFUtils.resolveURL('blob<' + index + '>', path));
				} else {
					loader.manager.itemEnd(imageUrl);
				}
			}
			return image;
		});
		return promise;
	}

	const ATTRIBUTES = {
		POSITION: 'a_Position',
		NORMAL: 'a_Normal',
		TANGENT: 'a_Tangent',
		TEXCOORD_0: 'a_Uv',
		TEXCOORD_1: 'a_Uv2',
		TEXCOORD_2: 'a_Uv3',
		TEXCOORD_3: 'a_Uv4',
		TEXCOORD_4: 'a_Uv5',
		TEXCOORD_5: 'a_Uv6',
		TEXCOORD_6: 'a_Uv7',
		TEXCOORD_7: 'a_Uv8',
		COLOR_0: 'a_Color',
		WEIGHTS_0: 'skinWeight',
		JOINTS_0: 'skinIndex',
		TEXCOORD0: 'a_Uv',
		// deprecated
		TEXCOORD: 'a_Uv',
		// deprecated
		COLOR0: 'a_Color',
		// deprecated
		COLOR: 'a_Color',
		// deprecated
		WEIGHT: 'skinWeight',
		// deprecated
		JOINT: 'skinIndex' // deprecated
	};
	const ALPHA_MODES = {
		MASK: 'MASK',
		BLEND: 'BLEND'
	};
	const ACCESSOR_TYPE_SIZES = {
		'SCALAR': 1,
		'VEC2': 2,
		'VEC3': 3,
		'VEC4': 4,
		'MAT2': 4,
		'MAT3': 9,
		'MAT4': 16
	};
	const ACCESSOR_COMPONENT_TYPES = {
		5120: Int8Array,
		5121: Uint8Array,
		5122: Int16Array,
		5123: Uint16Array,
		5125: Uint32Array,
		5126: Float32Array
	};
	const WEBGL_FILTERS = {
		9728: t3d.TEXTURE_FILTER.NEAREST,
		9729: t3d.TEXTURE_FILTER.LINEAR,
		9984: t3d.TEXTURE_FILTER.NEAREST_MIPMAP_NEAREST,
		9985: t3d.TEXTURE_FILTER.LINEAR_MIPMAP_NEAREST,
		9986: t3d.TEXTURE_FILTER.NEAREST_MIPMAP_LINEAR,
		9987: t3d.TEXTURE_FILTER.LINEAR_MIPMAP_LINEAR
	};
	const WEBGL_WRAPPINGS = {
		33071: t3d.TEXTURE_WRAP.CLAMP_TO_EDGE,
		33648: t3d.TEXTURE_WRAP.MIRRORED_REPEAT,
		10497: t3d.TEXTURE_WRAP.REPEAT
	};
	const WEBGL_DRAW_MODES = {
		POINTS: 0,
		LINES: 1,
		LINE_LOOP: 2,
		LINE_STRIP: 3,
		TRIANGLE_STRIP: 5,
		TRIANGLE_FAN: 6
	};

	class TextureParser {
		static parse(context) {
			const {
				gltf,
				images
			} = context;
			if (!gltf.textures) return;
			const textureCache = new Map();
			return Promise.all(gltf.textures.map((params, index) => {
				const {
					sampler,
					source = 0,
					name: textureName
				} = params;
				let sourceIndex = source,
					isTextureData = false;
				if (params.extensions) {
					const {
						KHR_texture_basisu
					} = params.extensions;
					if (KHR_texture_basisu) {
						sourceIndex = KHR_texture_basisu.source;
						isTextureData = true;
					} else if (Object.values(params.extensions).length && Object.values(params.extensions)[0].hasOwnProperty('source')) {
						sourceIndex = Object.values(params.extensions)[0].source;
					} else {
						console.warn('GLTFLoader: unknown texture extension');
					}
				}
				const cacheKey = sourceIndex + ':' + sampler;
				if (textureCache.has(cacheKey)) {
					return textureCache.get(cacheKey);
				}
				const texture = new t3d.Texture2D();
				if (isTextureData) {
					const {
						image,
						mipmaps,
						type,
						format,
						minFilter,
						magFilter,
						generateMipmaps,
						encoding,
						premultiplyAlpha
					} = images[sourceIndex];
					texture.image = image;
					texture.mipmaps = mipmaps;
					texture.type = type;
					texture.format = format;
					texture.minFilter = minFilter;
					texture.magFilter = magFilter;
					texture.generateMipmaps = generateMipmaps;
					texture.encoding = encoding;
					texture.premultiplyAlpha = premultiplyAlpha;
				} else {
					texture.image = images[sourceIndex];
				}
				texture.version++;
				texture.name = textureName || texture.image.__name || `texture_${index}`;
				texture.flipY = false;
				const {
					mimeType,
					uri
				} = gltf.images[sourceIndex];
				texture.userData.mimeType = mimeType || getImageURIMimeType(uri);
				const samplers = gltf.samplers || {};
				parseSampler(texture, samplers[sampler]);
				textureCache.set(cacheKey, texture);
				return texture;
			})).then(textures => {
				context.textures = textures;
				textureCache.clear();
			});
		}
	}
	function parseSampler(texture, sampler = {}) {
		const {
			magFilter,
			minFilter,
			wrapS,
			wrapT
		} = sampler;
		texture.magFilter = WEBGL_FILTERS[magFilter] || t3d.TEXTURE_FILTER.LINEAR;
		texture.minFilter = WEBGL_FILTERS[minFilter] || t3d.TEXTURE_FILTER.LINEAR_MIPMAP_LINEAR;
		texture.wrapS = WEBGL_WRAPPINGS[wrapS] || t3d.TEXTURE_WRAP.REPEAT;
		texture.wrapT = WEBGL_WRAPPINGS[wrapT] || t3d.TEXTURE_WRAP.REPEAT;
	}

	// only for jpeg, png, webp
	// because other should get mimeType from glTF image object
	function getImageURIMimeType(uri) {
		if (uri.startsWith('data:image/')) {
			// early return for data URIs
			if (uri.startsWith('data:image/jpeg')) return 'image/jpeg';
			if (uri.startsWith('data:image/webp')) return 'image/webp';
			return 'image/png';
		} else {
			if (uri.search(/\.jpe?g($|\?)/i) > 0) return 'image/jpeg';
			if (uri.search(/\.webp($|\?)/i) > 0) return 'image/webp';
			return 'image/png';
		}
	}

	let MaterialParser$1 = class MaterialParser {
		static parse(context, loader) {
			const {
				gltf,
				textures
			} = context;
			if (!gltf.materials) return;
			const transformExt = loader.extensions.get('KHR_texture_transform');
			const materials = [];
			for (let i = 0; i < gltf.materials.length; i++) {
				const {
					extensions = {},
					pbrMetallicRoughness,
					normalTexture,
					occlusionTexture,
					emissiveTexture,
					emissiveFactor,
					alphaMode,
					alphaCutoff,
					doubleSided,
					name = ''
				} = gltf.materials[i];
				let material = null;
				const materialExtNames = loader.autoParseConfig.materials;

				// TODO: refactor invoke method
				for (let j = 0; j < materialExtNames.length; j++) {
					const extName = materialExtNames[j];
					const extParams = extensions[extName];
					const ext = loader.extensions.get(extName);
					if (extParams && ext && ext.getMaterial) {
						material = ext.getMaterial();
						break;
					}
				}
				material = material || new t3d.PBRMaterial();
				material.name = name;
				for (let j = 0; j < materialExtNames.length; j++) {
					const extName = materialExtNames[j];
					const extParams = extensions[extName];
					const ext = loader.extensions.get(extName);
					if (extParams && ext && ext.parseParams) {
						ext.parseParams(material, extParams, textures, transformExt);
					}
				}
				const {
					KHR_materials_unlit,
					KHR_materials_pbrSpecularGlossiness
				} = extensions;
				if (pbrMetallicRoughness) {
					const {
						baseColorFactor,
						baseColorTexture,
						metallicFactor,
						roughnessFactor,
						metallicRoughnessTexture
					} = pbrMetallicRoughness;
					if (Array.isArray(baseColorFactor)) {
						material.diffuse.fromArray(baseColorFactor);
						material.opacity = baseColorFactor[3] !== undefined ? baseColorFactor[3] : 1;
					}
					if (baseColorTexture) {
						material.diffuseMap = textures[baseColorTexture.index];
						material.diffuseMapCoord = baseColorTexture.texCoord || 0;
						if (material.diffuseMap) {
							material.diffuseMap.encoding = t3d.TEXEL_ENCODING_TYPE.SRGB;
							transformExt && transformExt.handleMaterialMap(material, 'diffuseMap', baseColorTexture);
						}
					}
					if (!KHR_materials_unlit && !KHR_materials_pbrSpecularGlossiness) {
						material.metalness = metallicFactor !== undefined ? metallicFactor : 1;
						material.roughness = roughnessFactor !== undefined ? roughnessFactor : 1;
						if (metallicRoughnessTexture) {
							material.metalnessMap = textures[metallicRoughnessTexture.index];
							material.roughnessMap = textures[metallicRoughnessTexture.index];
							// metallicRoughnessTexture transform not supported yet
						}
					}
				}
				if (emissiveFactor) {
					material.emissive.fromArray(emissiveFactor);
				}
				if (emissiveTexture) {
					material.emissiveMap = textures[emissiveTexture.index];
					material.emissiveMapCoord = emissiveTexture.texCoord || 0;
					if (material.emissiveMap) {
						material.emissiveMap.encoding = t3d.TEXEL_ENCODING_TYPE.SRGB;
						transformExt && transformExt.handleMaterialMap(material, 'emissiveMap', emissiveTexture);
					}
				}
				if (occlusionTexture) {
					material.aoMap = textures[occlusionTexture.index];
					material.aoMapCoord = occlusionTexture.texCoord || 0;
					if (occlusionTexture.strength !== undefined) {
						material.aoMapIntensity = occlusionTexture.strength;
					}
					if (material.aoMap) {
						transformExt && transformExt.handleMaterialMap(material, 'aoMap', occlusionTexture);
					}
				}
				if (!KHR_materials_unlit) {
					if (normalTexture) {
						material.normalMap = textures[normalTexture.index];
						material.normalScale.set(1, -1);
						if (normalTexture.scale !== undefined) {
							// fix flip y for normal map
							// https://github.com/mrdoob/three.js/issues/11438#issuecomment-507003995
							material.normalScale.set(normalTexture.scale, -normalTexture.scale);
						}

						// normal map transform not supported yet
					}
				}
				material.side = doubleSided === true ? t3d.DRAW_SIDE.DOUBLE : t3d.DRAW_SIDE.FRONT;
				if (alphaMode === ALPHA_MODES.BLEND) {
					material.transparent = true;
				} else {
					material.transparent = false;
					if (alphaMode === ALPHA_MODES.MASK) {
						material.alphaTest = alphaCutoff !== undefined ? alphaCutoff : 0.5;
					}
				}
				materials[i] = material;
			}
			context.materials = materials;
		}
	};

	class AccessorParser {
		static parse(context) {
			const {
				bufferViews,
				gltf
			} = context;
			if (!gltf.accessors) return;
			const interleavedBufferCache = new Map();
			const accessors = gltf.accessors.map(accessor => {
				const {
					bufferView: bufferViewIndex,
					type,
					componentType,
					count,
					byteOffset = 0,
					normalized = false,
					sparse
				} = accessor;
				if (bufferViewIndex === undefined && sparse === undefined) {
					// Ignore empty accessors, which may be used to declare runtime
					// information about attributes coming from another source (e.g. Draco compression extension).
					return null;
				}

				// Get buffer view infos
				const bufferView = bufferViewIndex !== undefined ? bufferViews[bufferViewIndex] : null;
				const byteStride = bufferViewIndex !== undefined ? gltf.bufferViews[bufferViewIndex].byteStride : undefined;

				// Get accessor infos
				const itemSize = ACCESSOR_TYPE_SIZES[type];
				const TypedArray = ACCESSOR_COMPONENT_TYPES[componentType];
				const elementBytes = TypedArray.BYTES_PER_ELEMENT;
				const itemBytes = elementBytes * itemSize; // For VEC3: itemSize is 3, elementBytes is 4, itemBytes is 12.

				let array, attribute;
				if (byteStride && byteStride !== itemBytes) {
					// The buffer is interleaved
					// Each "slice" of the buffer, as defined by 'count' elements of 'byteStride' bytes, gets its own InterleavedBuffer
					// This makes sure that IBA.count reflects accessor.count properly
					const ibSlice = Math.floor(byteOffset / byteStride);
					const ibCacheKey = 'Buffer:' + bufferViewIndex + ':' + componentType + ':' + ibSlice + ':' + count;
					let ib = interleavedBufferCache.get(ibCacheKey);
					if (!ib) {
						// Use the full buffer if it's interleaved.
						array = new TypedArray(bufferView, ibSlice * byteStride, count * byteStride / elementBytes);

						// Integer parameters to IB/IBA are in array elements, not bytes.
						ib = new t3d.Buffer(array, byteStride / elementBytes);
						interleavedBufferCache.set(ibCacheKey, ib);
					}
					attribute = new t3d.Attribute(ib, itemSize, byteOffset % byteStride / elementBytes, normalized);
				} else {
					if (bufferView === null) {
						array = new TypedArray(count * itemSize);
					} else {
						array = new TypedArray(bufferView, byteOffset, count * itemSize);
					}
					attribute = new t3d.Attribute(new t3d.Buffer(array, itemSize), itemSize, 0, normalized);
				}

				// https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#sparse-accessors
				if (sparse) {
					const itemSizeIndices = ACCESSOR_TYPE_SIZES.SCALAR;
					const TypedArrayIndices = ACCESSOR_COMPONENT_TYPES[sparse.indices.componentType];
					const byteOffsetIndices = sparse.indices.byteOffset || 0;
					const byteOffsetValues = sparse.values.byteOffset || 0;
					const sparseIndices = new TypedArrayIndices(bufferViews[sparse.indices.bufferView], byteOffsetIndices, sparse.count * itemSizeIndices);
					const sparseValues = new TypedArray(bufferViews[sparse.values.bufferView], byteOffsetValues, sparse.count * itemSize);
					if (bufferView !== null) {
						// Avoid modifying the original ArrayBuffer, if the bufferView wasn't initialized with zeroes.
						attribute = new t3d.Attribute(attribute.buffer.clone(), attribute.size, attribute.offset, attribute.normalized);
					}
					const buffer = attribute.buffer;
					for (let i = 0, il = sparseIndices.length; i < il; i++) {
						const index = sparseIndices[i];
						buffer.array[index * attribute.size] = sparseValues[i * itemSize];
						if (itemSize >= 2) buffer.array[index * attribute.size + 1] = sparseValues[i * itemSize + 1];
						if (itemSize >= 3) buffer.array[index * attribute.size + 2] = sparseValues[i * itemSize + 2];
						if (itemSize >= 4) buffer.array[index * attribute.size + 3] = sparseValues[i * itemSize + 3];
						if (itemSize >= 5) throw new Error('Unsupported itemSize in sparse Attribute.');
					}
				}
				return attribute;
			});
			interleavedBufferCache.clear();
			context.accessors = accessors;
		}
	}

	let PrimitiveParser$1 = class PrimitiveParser {
		static parse(context, loader) {
			const {
				gltf,
				accessors,
				materials,
				bufferViews
			} = context;
			if (!gltf.meshes) return;
			const dracoExt = loader.extensions.get('KHR_draco_mesh_compression');
			const materialCache = new Map();
			const geometryPromiseCache = new Map();
			const meshPromises = [];
			for (let i = 0; i < gltf.meshes.length; i++) {
				const gltfMesh = gltf.meshes[i];
				const primitivePromises = [];
				for (let j = 0; j < gltfMesh.primitives.length; j++) {
					const gltfPrimitive = gltfMesh.primitives[j];
					const {
						extensions = {},
						mode,
						material
					} = gltfPrimitive;
					const {
						KHR_draco_mesh_compression
					} = extensions;
					let geometryPromise;
					const geometryKey = createGeometryKey$1(gltfPrimitive);
					if (geometryPromiseCache.has(geometryKey)) {
						geometryPromise = geometryPromiseCache.get(geometryKey);
					} else {
						if (KHR_draco_mesh_compression && dracoExt) {
							geometryPromise = dracoExt.getGeometry(KHR_draco_mesh_compression, bufferViews, gltfPrimitive.attributes, gltf.accessors, loader.getDRACOLoader());
						} else {
							geometryPromise = Promise.resolve(new t3d.Geometry());
						}
						geometryPromise = geometryPromise.then(geometry => {
							parseGeometryFromGLTFPrimitive$1(geometry, gltfPrimitive, gltf, accessors);
							return geometry;
						});
						geometryPromiseCache.set(geometryKey, geometryPromise);
					}
					const primitivePromise = geometryPromise.then(geometry => {
						const primitive = {
							mode,
							geometry,
							material: material === undefined ? new t3d.PBRMaterial() : materials[material],
							weights: Object.keys(geometry.morphAttributes).length > 0 && gltfMesh.weights ? gltfMesh.weights.slice(0) : undefined,
							skinned: gltfMesh.isSkinned
						};
						assignFinalMaterial$1(primitive, materialCache);
						return primitive;
					});
					primitivePromises.push(primitivePromise);
				}
				meshPromises.push(Promise.all(primitivePromises));
			}
			materialCache.clear();
			geometryPromiseCache.clear();
			return Promise.all(meshPromises).then(primitives => {
				context.primitives = primitives;
			});
		}
	};
	function parseGeometryFromGLTFPrimitive$1(geometry, gltfPrimitive, gltf, accessors) {
		const {
			attributes,
			indices,
			targets
		} = gltfPrimitive;

		// set attributes

		for (const attributeSemantic in attributes) {
			const accessorIdx = attributes[attributeSemantic];
			const attributeName = ATTRIBUTES[attributeSemantic] === undefined ? attributeSemantic : ATTRIBUTES[attributeSemantic];
			// Skip attributes already provided by e.g. Draco extension.
			if (attributeName in geometry.attributes) continue;
			geometry.addAttribute(attributeName, accessors[accessorIdx]);
		}

		// set index

		if (indices !== undefined && !geometry.index) {
			geometry.setIndex(accessors[indices]);
		}

		// compute bounds

		const {
			boundingBox,
			boundingSphere
		} = geometry;
		if (attributes.POSITION !== undefined) {
			const accessorIdx = attributes.POSITION;
			const accessor = gltf.accessors[accessorIdx];
			if (accessor.min && accessor.max) {
				boundingBox.min.fromArray(accessor.min);
				boundingBox.max.fromArray(accessor.max);
				if (accessor.normalized) {
					const boxScale = GLTFUtils.getNormalizedComponentScale(ACCESSOR_COMPONENT_TYPES[accessor.componentType]);
					boundingBox.min.multiplyScalar(boxScale);
					boundingBox.max.multiplyScalar(boxScale);
				}
			} else {
				geometry.computeBoundingBox();
			}
		} else {
			geometry.computeBoundingBox();
		}
		boundingBox.getCenter(boundingSphere.center);
		boundingSphere.radius = boundingBox.min.distanceTo(boundingBox.max) / 2;

		// set morph targets

		if (targets) {
			let hasMorphPosition = false;
			let hasMorphNormal = false;
			for (let i = 0, il = targets.length; i < il; i++) {
				const target = targets[i];
				if (target.POSITION !== undefined) hasMorphPosition = true;
				if (target.NORMAL !== undefined) hasMorphNormal = true;
				if (hasMorphPosition && hasMorphNormal) break;
			}
			if (hasMorphPosition || hasMorphNormal) {
				const morphPositions = [];
				const morphNormals = [];
				for (let i = 0, il = targets.length; i < il; i++) {
					const target = targets[i];
					if (hasMorphPosition) {
						morphPositions.push(target.POSITION !== undefined ? accessors[target.POSITION] : geometry.attributes[ATTRIBUTES.POSITION]);
					}
					if (hasMorphNormal) {
						morphNormals.push(target.NORMAL !== undefined ? accessors[target.NORMAL] : geometry.attributes[ATTRIBUTES.NORMAL]);
					}
				}
				if (hasMorphPosition) {
					geometry.morphAttributes.position = morphPositions;
				}
				if (hasMorphNormal) {
					geometry.morphAttributes.normal = morphNormals;
				}
			}
		}
		return geometry;
	}
	function assignFinalMaterial$1(primitive, materialCache) {
		let {
			geometry,
			material,
			skinned,
			mode
		} = primitive;

		// If the material will be modified later on, clone it now.
		const useVertexTangents = geometry.attributes[ATTRIBUTES.TANGENT] !== undefined;
		const useVertexColors = geometry.attributes[ATTRIBUTES.COLOR_0] !== undefined;
		const useFlatShading = geometry.attributes[ATTRIBUTES.NORMAL] === undefined;
		const useSkinning = skinned;
		if (mode === WEBGL_DRAW_MODES.POINTS) {
			const cacheKey = 'PointsMaterial:' + material.id;
			let pointsMaterial = materialCache.get(cacheKey);
			if (!pointsMaterial) {
				pointsMaterial = new t3d.PointsMaterial();
				t3d.Material.prototype.copy.call(pointsMaterial, material);
				pointsMaterial.diffuse.copy(material.diffuse);
				pointsMaterial.diffuseMap = material.map;
				pointsMaterial.drawMode = mode;
				pointsMaterial.acceptLight = false; // PointsMaterial doesn't support lights yet
				materialCache.set(cacheKey, pointsMaterial);
			}
			material = pointsMaterial;
		} else if (mode === WEBGL_DRAW_MODES.LINES || mode === WEBGL_DRAW_MODES.LINE_STRIP || mode === WEBGL_DRAW_MODES.LINE_LOOP) {
			const cacheKey = 'BasicMaterial:' + material.id;
			let basicMaterial = materialCache.get(cacheKey);
			if (!basicMaterial) {
				basicMaterial = new t3d.BasicMaterial();
				basicMaterial.envMap = undefined; // force close env map
				basicMaterial.diffuse.copy(material.diffuse);
				basicMaterial.diffuseMap = material.diffuseMap;
				basicMaterial.drawMode = mode;
				materialCache.set(cacheKey, basicMaterial);
			}
			material = basicMaterial;
		} else if (mode === WEBGL_DRAW_MODES.TRIANGLE_STRIP) {
			// TODO
			console.warn('TRIANGLE_STRIP will be removed later.');
			material.drawMode = WEBGL_DRAW_MODES.TRIANGLE_STRIP;
		} else if (mode === WEBGL_DRAW_MODES.TRIANGLE_FAN) {
			// TODO
			console.warn('TRIANGLE_FAN will be removed later.');
			material.drawMode = WEBGL_DRAW_MODES.TRIANGLE_FAN;
		}
		if (useVertexTangents || useVertexColors || useFlatShading || useSkinning) {
			let cacheKey = 'ClonedMaterial:' + material.id + ':';
			if (useVertexTangents) cacheKey += 'vertex-tangents:';
			if (useVertexColors) {
				if (geometry.attributes[ATTRIBUTES.COLOR_0].size === 3) {
					cacheKey += 'vertex-colors-rgb:';
				} else if (geometry.attributes[ATTRIBUTES.COLOR_0].size === 4) {
					cacheKey += 'vertex-colors-rgba:';
				}
			}
			if (useFlatShading) cacheKey += 'flat-shading:';
			let cachedMaterial = materialCache.get(cacheKey);
			if (!cachedMaterial) {
				cachedMaterial = material.clone();
				if (useVertexTangents) {
					cachedMaterial.vertexTangents = true;

					// revert flip y fix for tangents
					// https://github.com/mrdoob/three.js/issues/11438#issuecomment-507003995
					if (cachedMaterial.normalMap) {
						cachedMaterial.normalScale.y *= -1;
					}
				}
				if (useVertexColors) {
					if (geometry.attributes[ATTRIBUTES.COLOR_0].size === 3) {
						cachedMaterial.vertexColors = t3d.VERTEX_COLOR.RGB;
					} else if (geometry.attributes[ATTRIBUTES.COLOR_0].size === 4) {
						cachedMaterial.vertexColors = t3d.VERTEX_COLOR.RGBA;
					} else {
						console.warn('Illegal vertex color size: ' + geometry.attributes[ATTRIBUTES.COLOR_0].size);
					}
				}
				if (useFlatShading) {
					cachedMaterial.shading = t3d.SHADING_TYPE.FLAT_SHADING;
				}
			}
			material = cachedMaterial;
		}
		primitive.material = material;
	}
	function createGeometryKey$1(primitive) {
		const dracoExtension = primitive.extensions && primitive.extensions.KHR_draco_mesh_compression;
		let geometryKey;
		if (dracoExtension) {
			geometryKey = 'draco:' + dracoExtension.bufferView + ':' + dracoExtension.indices + ':' + createAttributesKey$1(dracoExtension.attributes);
		} else {
			geometryKey = primitive.indices + ':' + createAttributesKey$1(primitive.attributes) + ':' + primitive.mode;
		}
		if (primitive.targets) {
			for (let i = 0, il = primitive.targets.length; i < il; i++) {
				geometryKey += ':' + createAttributesKey$1(primitive.targets[i]);
			}
		}
		return geometryKey;
	}
	function createAttributesKey$1(attributes) {
		let attributesKey = '';
		const keys = Object.keys(attributes).sort();
		for (let i = 0, il = keys.length; i < il; i++) {
			attributesKey += keys[i] + ':' + attributes[keys[i]] + ';';
		}
		return attributesKey;
	}

	class NodeParser {
		static parse(context, loader) {
			const {
				gltf: {
					nodes: gltfNodes,
					cameras: gltfCameras,
					extensions: gltfExtensions
				}
			} = context;
			if (!gltfNodes) return;
			const lightsExt = loader.extensions.get('KHR_lights_punctual');
			const instancingExt = loader.extensions.get('EXT_mesh_gpu_instancing');
			const cameras = [];
			const lights = [];
			const nodes = gltfNodes.map(gltfNode => {
				const {
					matrix,
					translation,
					rotation,
					scale,
					camera: cameraID,
					mesh: meshID,
					extensions = {}
				} = gltfNode;
				const {
					KHR_lights_punctual,
					EXT_mesh_gpu_instancing
				} = extensions;
				let node = null;
				if (gltfNode.isBone) {
					// .isBone isn't in glTF spec. Marked in IndexParser
					node = new t3d.Bone();
				} else if (meshID !== undefined) {
					if (EXT_mesh_gpu_instancing && instancingExt) {
						node = instancingExt.getInstancedMesh(context, gltfNode);
					} else {
						node = createMesh(context, gltfNode);
					}
				} else if (cameraID !== undefined) {
					node = createCamera(gltfCameras[cameraID]);
					cameras.push(node);
				} else if (KHR_lights_punctual && lightsExt) {
					const lightIndex = KHR_lights_punctual.light;
					const gltfLights = gltfExtensions.KHR_lights_punctual.lights;
					node = lightsExt.getLight(gltfLights[lightIndex]);
					lights.push(node);
				} else {
					node = new t3d.Object3D();
				}
				node.name = gltfNode.name || '';
				if (!!node.name && node.children.length > 0) {
					for (let i = 0; i < node.children.length; i++) {
						node.children[i].name = node.name + '_' + i;
					}
				}
				if (matrix !== undefined) {
					node.matrix.fromArray(matrix);
					node.matrix.decompose(node.position, node.quaternion, node.scale);
				} else {
					if (translation !== undefined) {
						node.position.fromArray(translation);
					}
					if (rotation !== undefined) {
						node.quaternion.fromArray(rotation);
					}
					if (scale !== undefined) {
						node.scale.fromArray(scale);
					}
				}
				return node;
			});
			context.nodes = nodes;
			context.cameras = cameras;
			context.lights = lights;
		}
	}
	function createCamera(cameraDef) {
		const {
			orthographic,
			perspective,
			type
		} = cameraDef;
		const camera = new t3d.Camera();
		if (type == 'perspective') {
			const {
				aspectRatio,
				yfov,
				zfar,
				znear
			} = perspective;
			camera.setPerspective(yfov, aspectRatio || 1, znear || 1, zfar || 2e6);
		} else if (type == 'orthographic') {
			const {
				xmag,
				ymag,
				zfar,
				znear
			} = orthographic;
			// https:// github.com/KhronosGroup/glTF/issues/1663
			camera.setOrtho(-xmag, xmag, -ymag, ymag, znear || 1, zfar || 2e6);
		}
		return camera;
	}
	function createMesh(context, gltfNode) {
		const {
			primitives
		} = context;
		const {
			mesh: meshID,
			skin: skinID
		} = gltfNode;
		const meshes = primitives[meshID].map(primitive => {
			const {
				geometry,
				material,
				weights
			} = primitive;
			let mesh;
			if (skinID !== undefined) {
				mesh = new t3d.SkinnedMesh(geometry, material);
				if (geometry.attributes.skinWeight && !geometry.attributes.skinWeight.normalized) {
					GLTFUtils.normalizeSkinWeights(geometry.attributes.skinWeight);
				}
			} else {
				mesh = new t3d.Mesh(geometry, material);
				if (weights) {
					mesh.morphTargetInfluences = weights.slice();
				}
			}
			return mesh;
		});
		if (meshes.length > 1) {
			const parent = new t3d.Object3D();
			meshes.forEach(mesh => parent.add(mesh));
			return parent;
		} else {
			return meshes[0];
		}
	}

	class SkinParser {
		static parse(context) {
			const {
				gltf,
				accessors,
				nodes
			} = context;
			const gltfSkins = gltf.skins;
			if (!gltfSkins) return;
			const skins = gltfSkins.map(skin => {
				const {
					inverseBindMatrices,
					joints
				} = skin;
				const attribute = accessors[inverseBindMatrices];
				const bones = [];
				const boneInverses = [];
				joints.forEach((jointId, index) => {
					const jointNode = nodes[jointId];
					if (jointNode) {
						bones.push(jointNode);
						const boneInverse = new t3d.Matrix4();
						if (attribute) {
							boneInverse.fromArray(attribute.buffer.array, index * 16);
						}
						boneInverses.push(boneInverse);
					} else {
						console.warn('Joint ' + jointId + ' could not be found.');
					}
				});
				return new t3d.Skeleton(bones, boneInverses);
			});
			context.skins = skins;

			// Bind all skined meshes
			nodes.forEach((node, index) => {
				const {
					skin: skinID
				} = gltf.nodes[index];
				if (skinID !== undefined) {
					node.traverse(function (mesh) {
						if (!mesh.isSkinnedMesh) return;
						mesh.bind(skins[skinID], mesh.worldMatrix); // TODO need updateMatrix ?
					});
				}
			});
		}
	}

	class SceneParser {
		static parse(context) {
			const {
				gltf,
				nodes
			} = context;
			const roots = gltf.scenes.map(sceneDef => {
				const {
					name: sceneName = '',
					nodes: nodeIds = []
				} = sceneDef;
				const group = new t3d.Object3D();
				group.name = sceneName;
				for (let i = 0; i < nodeIds.length; i++) {
					buildNodeHierachy(nodeIds[i], group, gltf.nodes, nodes);
				}
				return group;
			});
			context.roots = roots;
			context.root = roots[gltf.scene || 0];
		}
	}
	function buildNodeHierachy(nodeId, parentNode, gltfNodes, nodes) {
		const node = nodes[nodeId];
		const nodeDef = gltfNodes[nodeId];
		parentNode.add(node);
		if (nodeDef.children) {
			const children = nodeDef.children;
			for (let i = 0, il = children.length; i < il; i++) {
				const child = children[i];
				buildNodeHierachy(child, node, gltfNodes, nodes);
			}
		}
	}

	class AnimationParser {
		static parse(context, loader) {
			const {
				gltf,
				nodes,
				accessors
			} = context;
			const {
				animations
			} = gltf;
			if (!animations) return;
			const pointerExt = loader.extensions.get('KHR_animation_pointer');
			const animationClips = animations.map((gltfAnimation, index) => {
				const {
					channels,
					samplers,
					name = `animation_${index}`
				} = gltfAnimation;
				const trackInfos = [];
				let duration = 0;
				for (let i = 0; i < channels.length; i++) {
					const gltfChannel = channels[i];
					const gltfSampler = samplers[gltfChannel.sampler];
					if (!gltfSampler) continue;
					const targetDef = gltfChannel.target;
					const inputAccessor = accessors[gltfSampler.input];
					const input = new inputAccessor.buffer.array.constructor(inputAccessor.buffer.array);
					const outputAccessor = accessors[gltfSampler.output];
					const output = new Float32Array(outputAccessor.buffer.array);
					if (outputAccessor.normalized) {
						const scale = GLTFUtils.getNormalizedComponentScale(outputAccessor.buffer.array.constructor);
						for (let j = 0, jl = output.length; j < jl; j++) {
							output[j] *= scale;
						}
					}
					duration = Math.max(duration, input[input.length - 1]);
					if (pointerExt && targetDef.extensions && targetDef.extensions['KHR_animation_pointer']) {
						pointerExt.getTrackInfos(context, targetDef.extensions['KHR_animation_pointer'], input, output, gltfSampler.interpolation, trackInfos);
					} else {
						const target = nodes[targetDef.node !== undefined ? targetDef.node : targetDef.id]; // Note: targetDef.id is deprecated.

						if (!target) continue;
						let TypedKeyframeTrack, propertyPath;
						if (targetDef.path === 'rotation') {
							TypedKeyframeTrack = t3d.QuaternionKeyframeTrack;
							propertyPath = 'quaternion';
						} else if (targetDef.path === 'weights') {
							TypedKeyframeTrack = t3d.NumberKeyframeTrack;
							propertyPath = 'morphTargetInfluences';
						} else if (targetDef.path === 'translation') {
							TypedKeyframeTrack = t3d.VectorKeyframeTrack;
							propertyPath = 'position';
						} else if (targetDef.path === 'scale') {
							TypedKeyframeTrack = t3d.VectorKeyframeTrack;
							propertyPath = 'scale';
						} else {
							continue;
						}
						trackInfos.push({
							TypedKeyframeTrack,
							target,
							propertyPath,
							times: input,
							values: output,
							interpolation: gltfSampler.interpolation
						});
					}
				}
				const tracks = [];
				trackInfos.forEach(trackInfo => {
					const {
						TypedKeyframeTrack,
						target,
						propertyPath,
						times,
						values,
						interpolation
					} = trackInfo;
					const interpolant = getInterpolant(interpolation, TypedKeyframeTrack === t3d.QuaternionKeyframeTrack);
					if (propertyPath === 'morphTargetInfluences') {
						// node may be a Object3D (glTF mesh with several primitives) or a Mesh.
						target.traverse(object => {
							if (object.isMesh && object.morphTargetInfluences) {
								const track = new TypedKeyframeTrack(object, propertyPath, times, values, interpolant);
								tracks.push(track);
							}
						});
					} else {
						const track = new TypedKeyframeTrack(target, propertyPath, times, values, interpolant);
						tracks.push(track);
					}
				});
				return new t3d.KeyframeClip(name, tracks, duration);
			});
			context.animations = animationClips;
		}
	}
	function getInterpolant(type, quaternion) {
		switch (type) {
			case 'STEP':
				return t3d.StepInterpolant;
			case 'CUBICSPLINE':
				return quaternion ? t3d.QuaternionCubicSplineInterpolant : t3d.CubicSplineInterpolant;
			case 'LINEAR':
			default:
				return quaternion ? t3d.QuaternionLinearInterpolant : t3d.LinearInterpolant;
		}
	}

	let resourceId = 0;
	class GLTFResource {
		constructor() {
			this.id = ++resourceId;
			this.url = ''; // url string
			this.path = ''; // path string
			this.options = null; // load options
			this.gltf = null; // gltf json after IndexParser
			this.loadItems = null; // String[] after IndexParser. Store all urls that need to load.
			this.buffers = null; // ArrayBuffer[] after BufferParser
			this.bufferViews = null; // ArrayBuffer[] after BufferViewParser
			this.images = null; // Image[] after ImageParser
			this.textures = null; // Texture2D[] after TextureParser
			this.materials = null; // Material[] after MaterialParser
			this.accessors = null; // Attribute[] after AccessorParser
			this.primitives = null; // { mode, geometry, material, weights, skinned }[] after PrimitiveParser
			this.nodes = null; // Object3D[] after NodeParser
			this.cameras = null; // Camera[] after NodeParser
			this.lights = null; // Light[] after NodeParser
			this.skins = null; // Skeleton[] after SkinParser
			this.root = null; // root after SceneParser
			this.roots = null; // root[] after SceneParser
			this.animations = null; // KeyframeClip[] after AnimationParser
		}
	}

	/**
	 * meshopt BufferView Compression Extension
	 *
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Vendor/EXT_meshopt_compression
	 */
	class EXT_meshopt_compression {
		static loadBufferView(extensionDef, buffers, meshoptDecoder) {
			const buffer = buffers[extensionDef.buffer];
			if (!meshoptDecoder || !meshoptDecoder.supported) {
				throw new Error('GLTFLoader: setMeshoptDecoder must be called before loading compressed files.');
			}
			const byteOffset = extensionDef.byteOffset || 0;
			const byteLength = extensionDef.byteLength || 0;
			const count = extensionDef.count;
			const stride = extensionDef.byteStride;
			const source = new Uint8Array(buffer, byteOffset, byteLength);
			if (meshoptDecoder.decodeGltfBufferAsync) {
				return meshoptDecoder.decodeGltfBufferAsync(count, stride, source, extensionDef.mode, extensionDef.filter).then(res => res.buffer);
			} else {
				// Support for MeshoptDecoder 0.18 or earlier, without decodeGltfBufferAsync
				return meshoptDecoder.ready.then(() => {
					const result = new ArrayBuffer(count * stride);
					meshoptDecoder.decodeGltfBuffer(new Uint8Array(result), count, stride, source, extensionDef.mode, extensionDef.filter);
					return result;
				});
			}
		}
	}

	/**
	 * KHR_draco_mesh_compression extension
	 * https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_draco_mesh_compression
	 */
	class KHR_draco_mesh_compression {
		static getGeometry(params, bufferViews, attributes, accessors, dracoLoader) {
			const {
				bufferView: bufferViewIndex,
				attributes: gltfAttributeMap
			} = params;
			if (!dracoLoader) {
				throw new Error('GLTFLoader: No DRACOLoader instance provided.');
			}
			const attributeMap = {};
			for (const attributeSemantic in gltfAttributeMap) {
				const attributeName = ATTRIBUTES[attributeSemantic] === undefined ? attributeSemantic : ATTRIBUTES[attributeSemantic];
				attributeMap[attributeName] = gltfAttributeMap[attributeSemantic];
			}
			const attributeNormalizedMap = {};
			const attributeTypeMap = {};
			for (const attributeNameItem in attributes) {
				const attributeName = ATTRIBUTES[attributeNameItem] || attributeNameItem.toLowerCase();
				if (gltfAttributeMap[attributeNameItem] !== undefined) {
					const accessorDef = accessors[attributes[attributeNameItem]];
					const componentType = ACCESSOR_COMPONENT_TYPES[accessorDef.componentType];
					attributeTypeMap[attributeName] = componentType.name;
					attributeNormalizedMap[attributeName] = accessorDef.normalized === true;
				}
			}
			const bufferView = bufferViews[bufferViewIndex];
			return new Promise(function (resolve) {
				dracoLoader.decodeDracoFile(bufferView, function (geometry) {
					for (const attributeName in geometry.attributes) {
						const attribute = geometry.attributes[attributeName];
						const normalized = attributeNormalizedMap[attributeName];
						if (normalized !== undefined) attribute.normalized = normalized;
					}
					resolve(geometry);
				}, attributeMap, attributeTypeMap);
			});
		}
	}

	/**
	 * KHR_lights_punctual extension
	 * https://github.com/KhronosGroup/glTF/blob/master/extensions/2.0/Khronos/KHR_lights_punctual/README.md
	 */
	class KHR_lights_punctual {
		static getLight(params) {
			const {
				color,
				intensity = 1,
				type,
				range,
				spot
			} = params;
			let lightNode;
			if (type === 'directional') {
				lightNode = new t3d.DirectionalLight();
			} else if (type === 'point') {
				lightNode = new t3d.PointLight();
				if (range !== undefined) {
					lightNode.distance = range;
				}

				// https://github.com/KhronosGroup/glTF/blob/master/extensions/2.0/Khronos/KHR_lights_punctual/README.md#range-property
				// lightNode.decay = 2;
			} else if (type === 'spot') {
				lightNode = new t3d.SpotLight();
				if (range !== undefined) {
					lightNode.distance = range;
				}

				// https://github.com/KhronosGroup/glTF/blob/master/extensions/2.0/Khronos/KHR_lights_punctual/README.md#range-property
				// lightNode.decay = 2;

				if (spot) {
					const {
						innerConeAngle = 0,
						outerConeAngle = Math.PI / 4
					} = spot;
					lightNode.angle = outerConeAngle;
					lightNode.penumbra = 1.0 - innerConeAngle / outerConeAngle;
				}
			} else {
				throw new Error('Unexpected light type: ' + type);
			}
			if (color) {
				lightNode.color.fromArray(color);
			}
			lightNode.intensity = intensity;
			return lightNode;
		}
	}

	/**
	 * KHR_materials_clearcoat extension
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_materials_clearcoat
	 */
	class KHR_materials_clearcoat {
		static getMaterial() {
			return new t3d.PBRMaterial();
		}
		static parseParams(material, extension, textures) {
			if (material.constructor !== t3d.PBRMaterial) return;
			const {
				clearcoatFactor,
				clearcoatTexture,
				clearcoatRoughnessFactor,
				clearcoatRoughnessTexture,
				clearcoatNormalTexture
			} = extension;
			if (clearcoatFactor) {
				material.clearcoat = clearcoatFactor;
			}
			if (clearcoatTexture) {
				material.clearcoatMap = textures[clearcoatTexture.index];
				// material does not yet support the transform of clearcoatMap.
				// parseTextureTransform(material, 'clearcoatMap', clearcoatTexture.extensions);
			}
			if (clearcoatRoughnessFactor) {
				material.clearcoatRoughness = clearcoatRoughnessFactor;
			}
			if (clearcoatRoughnessTexture) {
				material.clearcoatRoughnessMap = textures[clearcoatRoughnessTexture.index];
				// material does not yet support the transform of clearcoatRoughnessMap.
				// parseTextureTransform(material, 'clearcoatRoughnessMap', clearcoatRoughnessTexture.extensions);
			}
			if (clearcoatNormalTexture) {
				material.clearcoatNormalMap = textures[clearcoatNormalTexture.index];
				// material does not yet support the transform of clearcoatNormalMap.
				// parseTextureTransform(material, 'clearcoatNormalMap', clearcoatNormalTexture.extensions);
				if (clearcoatNormalTexture.scale) {
					const scale = clearcoatNormalTexture.scale;
					material.clearcoatNormalScale = new t3d.Vector2(scale, scale);
				}
			}
		}
	}

	/**
	 * KHR_materials_pbrSpecularGlossiness extension
	 * https://github.com/KhronosGroup/glTF/blob/main/extensions/2.0/Archived/KHR_materials_pbrSpecularGlossiness/README.md
	 */
	class KHR_materials_pbrSpecularGlossiness {
		static getMaterial() {
			return new t3d.PBR2Material();
		}
		static parseParams(material, params, textures, transformExt) {
			if (material.constructor !== t3d.PBR2Material) return;
			const {
				diffuseFactor,
				diffuseTexture,
				specularFactor,
				glossinessFactor,
				specularGlossinessTexture
			} = params;
			if (Array.isArray(diffuseFactor)) {
				material.diffuse.fromArray(diffuseFactor);
				material.opacity = diffuseFactor[3] || 1;
			}
			if (diffuseTexture) {
				material.diffuseMap = textures[diffuseTexture.index];
				material.diffuseMapCoord = diffuseTexture.texCoord || 0;
				if (material.diffuseMap) {
					material.diffuseMap.encoding = t3d.TEXEL_ENCODING_TYPE.SRGB;
					transformExt && transformExt.handleMaterialMap(material, 'diffuseMap', diffuseTexture);
				}
			}
			material.glossiness = glossinessFactor !== undefined ? glossinessFactor : 1.0;
			if (Array.isArray(specularFactor)) {
				material.specular.fromArray(specularFactor);
			}
			if (specularGlossinessTexture) {
				material.glossinessMap = textures[specularGlossinessTexture.index];
				material.specularMap = textures[specularGlossinessTexture.index];
				// specularGlossinessTexture transform not supported yet
			}
		}
	}

	/**
	 * KHR_materials_unlit extension
	 * https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_materials_unlit
	 */
	class KHR_materials_unlit {
		static getMaterial() {
			return new t3d.BasicMaterial();
		}
	}

	/**
	 * BasisU Texture Extension
	 *
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_texture_basisu
	 */
	class KHR_texture_basisu {
		static loadTextureData(url, ktx2Loader) {
			return new Promise((resolve, reject) => {
				ktx2Loader.load(url, resolve, undefined, reject);
			});
		}
	}

	/**
	 * KHR_texture_transform extension
	 * https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_texture_transform
	 */
	class KHR_texture_transform {
		static handleMaterialMap(material, mapType, textureDef) {
			if (!textureDef.extensions) return;
			const extDef = textureDef.extensions.KHR_texture_transform;
			if (!extDef) return;

			// If texCoord is present, it overrides the texture's texCoord
			if (extDef.texCoord !== undefined) {
				material[mapType + 'Coord'] = extDef.texCoord;
			}
			const transform = material[mapType + 'Transform'];
			if (!transform) return;
			if (extDef.offset !== undefined) {
				transform.offset.fromArray(extDef.offset);
			}
			if (extDef.rotation !== undefined) {
				transform.rotation = extDef.rotation;
			}
			if (extDef.scale !== undefined) {
				transform.scale.fromArray(extDef.scale);
			}
			transform.updateMatrix();
		}
	}

	class KHR_animation_pointer {
		static getTrackInfos(context, extensionDef, input, output, interpolation, trackInfos) {
			const {
				pointer
			} = extensionDef;
			const segments = pointer.replace(/^\//, '').split('/');
			const type = segments[0];
			const index = parseInt(segments[1]);
			const property = segments[segments.length - 1];
			const searchArray = context[type];
			if (!searchArray) return;
			const target = searchArray[index];
			if (!target) return;
			let TypedKeyframeTrack, propertyPath, TypedKeyframeTrack2, propertyPath2;
			if (property === 'rotation') {
				TypedKeyframeTrack = t3d.QuaternionKeyframeTrack;
				propertyPath = 'quaternion';
			} else if (property === 'weights') {
				TypedKeyframeTrack = t3d.NumberKeyframeTrack;
				propertyPath = 'morphTargetInfluences';
			} else if (property === 'translation') {
				TypedKeyframeTrack = t3d.VectorKeyframeTrack;
				propertyPath = 'position';
			} else if (property === 'scale') {
				TypedKeyframeTrack = t3d.VectorKeyframeTrack;
				propertyPath = 'scale';
			} else if (property === 'baseColorFactor') {
				TypedKeyframeTrack = t3d.ColorKeyframeTrack;
				propertyPath = 'diffuse';
				TypedKeyframeTrack2 = t3d.NumberKeyframeTrack;
				propertyPath2 = 'opacity';
			} else if (property === 'metallicFactor') {
				TypedKeyframeTrack = t3d.NumberKeyframeTrack;
				propertyPath = 'metalness';
			} else if (property === 'roughnessFactor') {
				TypedKeyframeTrack = t3d.NumberKeyframeTrack;
				propertyPath = 'roughness';
			} else if (property === 'emissiveFactor') {
				TypedKeyframeTrack = t3d.VectorKeyframeTrack;
				propertyPath = 'emissive';
			} else if (segments[segments.length - 2] === 'KHR_texture_transform') {
				TypedKeyframeTrack = t3d.VectorKeyframeTrack;
				const textureProperty = segments[segments.length - 4];
				if (textureProperty === 'baseColorTexture') {
					propertyPath = 'diffuseMapTransform.' + property;
				} else if (textureProperty === 'emissiveTexture') {
					propertyPath = 'emissiveMapTransform.' + property;
				} else {
					return;
				}
			} else {
				return;
			}
			if (property === 'baseColorFactor') {
				// Separate the alpha channel from the color
				const color3Output = new Float32Array(output.length / 4 * 3);
				const alphaOutput = new Float32Array(output.length / 4);
				for (let i = 0; i < output.length / 4; i++) {
					color3Output[i * 3] = output[i * 4];
					color3Output[i * 3 + 1] = output[i * 4 + 1];
					color3Output[i * 3 + 2] = output[i * 4 + 2];
					alphaOutput[i] = output[i * 4 + 3];
				}
				trackInfos.push({
					TypedKeyframeTrack,
					target,
					propertyPath,
					times: input,
					values: color3Output,
					interpolation
				});
				trackInfos.push({
					TypedKeyframeTrack: TypedKeyframeTrack2,
					target,
					propertyPath: propertyPath2,
					times: input,
					values: alphaOutput,
					interpolation
				});
			} else {
				trackInfos.push({
					TypedKeyframeTrack,
					target,
					propertyPath,
					times: input,
					values: output,
					interpolation
				});
			}
		}
	}

	const DefaultParsePipeline = [IndexParser$1, ReferenceParser, Validator, BufferParser, BufferViewParser, ImageParser, TextureParser, MaterialParser$1, AccessorParser, PrimitiveParser$1, NodeParser, SkinParser, SceneParser, AnimationParser];
	const DefaultExtensions = new Map([['EXT_meshopt_compression', EXT_meshopt_compression], ['KHR_draco_mesh_compression', KHR_draco_mesh_compression], ['KHR_lights_punctual', KHR_lights_punctual], ['KHR_materials_clearcoat', KHR_materials_clearcoat], ['KHR_materials_pbrSpecularGlossiness', KHR_materials_pbrSpecularGlossiness], ['KHR_materials_unlit', KHR_materials_unlit], ['KHR_mesh_quantization', {}],
	// This is supported by default
	['KHR_texture_basisu', KHR_texture_basisu], ['KHR_texture_transform', KHR_texture_transform], ['KHR_animation_pointer', KHR_animation_pointer]]);
	class GLTFLoader {
		constructor(manager = t3d.DefaultLoadingManager, parsers = DefaultParsePipeline, extensions = DefaultExtensions) {
			this.manager = manager;

			// If ture, loading manager will dispatch progress for every buffer and image.
			// otherwise, loading manager will only dispatch progress for the whole gltf resource.
			this.detailLoadProgress = true;

			// If set false, need add Promise.catch to catch errors.
			this.autoLogError = true;
			this.extensions = new Map(extensions);

			// Indicate which extensions can be parsed in a uniform way.
			this.autoParseConfig = {
				materials: ['KHR_materials_clearcoat', 'KHR_materials_pbrSpecularGlossiness', 'KHR_materials_unlit', 'KHR_materials_transmission', 'KHR_materials_ior', 'KHR_materials_volume', 'KHR_materials_dispersion']
			};
			this._parsers = parsers.slice(0);
			this._dracoLoader = null;
			this._meshoptDecoder = null;
			this._ktx2Loader = null;
			this._fileLoader = new t3d.FileLoader();
			const userAgent = navigator.userAgent;
			const isSafari = /^((?!chrome|android).)*safari/i.test(userAgent) === true;
			const safariMatch = userAgent.match(/Version\/(\d+)/);
			const safariVersion = isSafari && safariMatch ? parseInt(safariMatch[1], 10) : -1;
			const isFirefox = userAgent.indexOf('Firefox') > -1;
			const firefoxVersion = isFirefox ? userAgent.match(/Firefox\/([0-9]+)\./)[1] : -1;
			if (typeof createImageBitmap === 'undefined' || isSafari && safariVersion < 17 || isFirefox && firefoxVersion < 98) {
				this._imageLoader = new t3d.ImageLoader();
			} else {
				this._imageLoader = new ImageBitmapLoader();
			}
		}
		load(url, options = {}) {
			this.manager.itemStart(url);
			return new Promise((resolve, reject) => {
				const resource = new GLTFResource();
				resource.url = url;
				resource.path = GLTFUtils.extractUrlBase(url);
				resource.options = options;
				this._parse(resource).then(resolve).then(() => this.manager.itemEnd(url)).catch(e => {
					if (this.autoLogError) {
						console.error(e);
					}
					if (this.detailLoadProgress && resource.loadItems) {
						resource.loadItems.forEach(item => {
							this.manager.itemEnd(item);
						});
					}
					this.manager.itemError(url);
					this.manager.itemEnd(url);
					reject(`Error loading glTF model from ${url} .`);
				});
			});
		}
		_parse(context) {
			let lastParser;
			return new Promise((resolve, reject) => {
				this._parsers.forEach(parser => {
					if (lastParser) {
						lastParser = lastParser.then(() => parser.parse(context, this));
					} else {
						lastParser = parser.parse(context, this);
					}
				});
				if (lastParser) {
					lastParser.then(() => resolve(context)).catch(reject);
				} else {
					resolve(context);
				}
			});
		}
		setDRACOLoader(dracoLoader) {
			this._dracoLoader = dracoLoader;
			return this;
		}
		getDRACOLoader() {
			return this._dracoLoader;
		}
		setMeshoptDecoder(meshoptDecoder) {
			this._meshoptDecoder = meshoptDecoder;
			return this;
		}
		getMeshoptDecoder() {
			return this._meshoptDecoder;
		}
		setKTX2Loader(ktx2Loader) {
			this._ktx2Loader = ktx2Loader;
			return this;
		}
		getKTX2Loader() {
			return this._ktx2Loader;
		}
		loadFile(url, type = 'json') {
			this._fileLoader.setResponseType(type);
			return new Promise((resolve, reject) => {
				url = this.manager.resolveURL(url);
				this._fileLoader.load(url, resolve, undefined, reject);
			});
		}
		loadImage(url) {
			return new Promise((resolve, reject) => {
				url = this.manager.resolveURL(url);
				this._imageLoader.load(url, resolve, undefined, reject);
			});
		}
		insertParser(parser, index) {
			this._parsers.splice(index, 0, parser);
		}
		replaceParser(parser, index) {
			this._parsers.splice(index, 1, parser);
		}
	}

	/**
	 * This parser is used to parse the header of a 3D Tiles resource.
	 * For 'b3dm', 'i3dm', 'pnts' and 'cmpt' formats.
	 */
	class HeaderParser {
		static parse(context, loader) {
			const buffer = context.options.buffer;

			// TODO: this should be able to take a uint8array with an offset and length
			const dataView = new DataView(buffer);

			// 32-byte header for i3dm and 28-byte header for the others.

			// 4 bytes for the magic, can be 'b3dm', 'i3dm', 'pnts', 'cmpt' for now.
			// 'vctr' is not supported yet.

			const magic = String.fromCharCode(dataView.getUint8(0)) + String.fromCharCode(dataView.getUint8(1)) + String.fromCharCode(dataView.getUint8(2)) + String.fromCharCode(dataView.getUint8(3));
			const urlExtension = getUrlExtension(context.url);
			if (magic !== urlExtension) {
				throw `Not a ${urlExtension} type resource, with url ${context.url}!`;
			}

			// 4 bytes for the version number.

			const version = dataView.getUint32(4, true);
			if (version !== 1) {
				throw `${urlExtension} version must be 1, with url ${context.url}!`;
			}

			// 4 bytes for the byte length of the entire tile content.

			const byteLength = dataView.getUint32(8, true);
			if (byteLength !== buffer.byteLength) {
				throw `${urlExtension} data byte length check failed, with url ${context.url}!`;
			}

			// output the header information to the context and return
			// if the tile content is cmpt.

			if (urlExtension === 'cmpt') {
				// 4 bytes for the tiles length
				const tilesLength = dataView.getUint32(12, true);
				context.header = {
					magic,
					version,
					byteLength,
					tilesLength
				};
				return;
			}

			// 4 bytes for the byte length of the feature table JSON.

			const featureTableJSONByteLength = dataView.getUint32(12, true);

			// 4 bytes for the byte length of the feature table binary.

			const featureTableBinaryByteLength = dataView.getUint32(16, true);

			// 4 bytes for the byte length of the batch table JSON.

			const batchTableJSONByteLength = dataView.getUint32(20, true);

			// 4 bytes for the byte length of the batch table binary.

			const batchTableBinaryByteLength = dataView.getUint32(24, true);

			// 4 bytes for the gltf format if the tile content format is i3dm.

			let gltfFormat = null;
			if (urlExtension === 'i3dm') {
				gltfFormat = dataView.getUint32(28, true);
			}

			// output the header information to the context.

			context.header = {
				magic,
				version,
				byteLength,
				featureTableJSONByteLength,
				featureTableBinaryByteLength,
				batchTableJSONByteLength,
				batchTableBinaryByteLength,
				gltfFormat
			};
			return;
		}
	}

	class FeatureTable {
		constructor(buffer, start, headerLength, binLength) {
			this.buffer = buffer;
			this.binOffset = start + headerLength;
			this.binLength = binLength;
			let header = null;
			if (headerLength !== 0) {
				const headerData = new Uint8Array(buffer, start, headerLength);
				header = JSON.parse(GLTFUtils.decodeText(headerData));
			} else {
				header = {};
			}
			this.header = header;
		}
		getKeys() {
			return Object.keys(this.header);
		}
		getData(key, count, defaultComponentType = null, defaultType = null) {
			const header = this.header;
			if (!(key in header)) {
				return null;
			}
			const feature = header[key];
			if (!(feature instanceof Object)) {
				return feature;
			} else if (Array.isArray(feature)) {
				return feature;
			} else {
				const {
					buffer,
					binOffset,
					binLength
				} = this;
				const byteOffset = feature.byteOffset || 0;
				const featureType = feature.type || defaultType;
				const featureComponentType = feature.componentType || defaultComponentType;
				if ('type' in feature && defaultType && feature.type !== defaultType) {
					throw new Error('FeatureTable: Specified type does not match expected type.');
				}
				let stride;
				switch (featureType) {
					case 'SCALAR':
						stride = 1;
						break;
					case 'VEC2':
						stride = 2;
						break;
					case 'VEC3':
						stride = 3;
						break;
					case 'VEC4':
						stride = 4;
						break;
					default:
						throw new Error(`FeatureTable: Feature type not provided for "${key}".`);
				}
				let data;
				const arrayStart = binOffset + byteOffset;
				const arrayLength = count * stride;
				switch (featureComponentType) {
					case 'BYTE':
						data = new Int8Array(buffer, arrayStart, arrayLength);
						break;
					case 'UNSIGNED_BYTE':
						data = new Uint8Array(buffer, arrayStart, arrayLength);
						break;
					case 'SHORT':
						data = new Int16Array(buffer, arrayStart, arrayLength);
						break;
					case 'UNSIGNED_SHORT':
						data = new Uint16Array(buffer, arrayStart, arrayLength);
						break;
					case 'INT':
						data = new Int32Array(buffer, arrayStart, arrayLength);
						break;
					case 'UNSIGNED_INT':
						data = new Uint32Array(buffer, arrayStart, arrayLength);
						break;
					case 'FLOAT':
						data = new Float32Array(buffer, arrayStart, arrayLength);
						break;
					case 'DOUBLE':
						data = new Float64Array(buffer, arrayStart, arrayLength);
						break;
					default:
						throw new Error(`FeatureTable: Feature component type not provided for "${key}".`);
				}
				const dataEnd = arrayStart + arrayLength * data.BYTES_PER_ELEMENT;
				if (dataEnd > binOffset + binLength) {
					throw new Error('FeatureTable: Feature data read outside binary body length.');
				}
				return data;
			}
		}
	}

	class BatchTable extends FeatureTable {
		constructor(buffer, batchSize, start, headerLength, binLength) {
			super(buffer, start, headerLength, binLength);
			this.batchSize = batchSize;
		}
		getData(key, componentType = null, type = null) {
			return super.getData(key, this.batchSize, componentType, type);
		}
	}

	/**
	 * This parser is used to parse the feature table and batch table of a 3D Tiles resource.
	 * For 'b3dm', 'i3dm' and 'pnts' formats.
	 */
	class TableParser {
		static parse(context, loader) {
			const {
				header,
				options
			} = context;
			const buffer = options.buffer;
			const featureTableStart = header.magic === 'i3dm' ? 32 : 28;
			const featureTableEnd = featureTableStart + header.featureTableJSONByteLength + header.featureTableBinaryByteLength;
			const batchTableStart = featureTableEnd;
			const batchTableEnd = batchTableStart + header.batchTableJSONByteLength + header.batchTableBinaryByteLength;

			// parse the feature table

			const featureTableBuffer = buffer.slice(featureTableStart, featureTableEnd);
			const featureTable = new FeatureTable(featureTableBuffer, 0, header.featureTableJSONByteLength, header.featureTableBinaryByteLength);

			// parse the batch table

			let batchSize;
			if (header.magic === 'b3dm') {
				batchSize = featureTable.getData('BATCH_LENGTH');
			} else if (header.magic === 'i3dm') {
				batchSize = featureTable.getData('INSTANCES_LENGTH');
			} else if (header.magic === 'pnts') {
				batchSize = featureTable.getData('BATCH_LENGTH') || featureTable.getData('POINTS_LENGTH');
			} else {
				throw `Unrecognized magic: ${header.magic}!`;
			}
			const batchTableBuffer = buffer.slice(batchTableStart, batchTableEnd);
			const batchTable = new BatchTable(batchTableBuffer, batchSize, 0, header.batchTableJSONByteLength, header.batchTableBinaryByteLength);

			// output the tables to the context

			context.featureTable = featureTable;
			context.batchTable = batchTable;
			context.batchTableEnd = batchTableEnd;
		}
	}

	class B3DMParser {
		static parse(context, loader) {
			const glbBytes = new Uint8Array(context.options.buffer, context.batchTableEnd, context.header.byteLength - context.batchTableEnd);
			const glbData = GLTFUtils.parseGLB(glbBytes.slice().buffer);
			context.gltf = glbData.gltf;
			context.buffers = glbData.buffers;
		}
	}

	class B3DMRootParser {
		static parse(context, loader) {
			const {
				root,
				featureTable,
				options
			} = context;

			// fix rtc center

			const rtcCenter = featureTable.getData('RTC_CENTER');
			if (rtcCenter) {
				root.position.x += rtcCenter[0];
				root.position.y += rtcCenter[1];
				root.position.z += rtcCenter[2];
			}
			if (options.adjustmentTransform) {
				root.matrix.transform(root.position, root.scale, root.quaternion);
				root.matrix.multiply(options.adjustmentTransform);
				root.matrix.decompose(root.position, root.quaternion, root.scale);
			}
		}
	}

	/**
	 * KHR_techniques_webgl extension
	 * https://github.com/KhronosGroup/glTF/blob/main/extensions/2.0/Archived/KHR_techniques_webgl/README.md
	 * This extension has been archived, so we only provide a basic implementation.
	 */
	class KHR_techniques_webgl {
		static getMaterial() {
			return new t3d.PBRMaterial();
		}
		static parseParams(material, extension, textures) {
			const {
				values
			} = extension;
			const {
				u_diffuse
			} = values;
			if (u_diffuse) {
				material.diffuseMap = textures[u_diffuse.index];
				material.diffuseMapCoord = u_diffuse.texCoord || 0;
			}
		}
	}

	/**
	 * B3DMLoader is a loader for the B3DM format.
	 */
	class B3DMLoader extends GLTFLoader {
		constructor(manager) {
			super(manager, [HeaderParser,
			// insert HeaderParser
			TableParser,
			// insert TableParser
			B3DMParser,
			// insert B3DMParser
			ReferenceParser, Validator, BufferParser, BufferViewParser, ImageParser, TextureParser, MaterialParser$1, AccessorParser, PrimitiveParser$1, NodeParser, SkinParser, SceneParser, AnimationParser, B3DMRootParser // insert B3DMRootParser
			]);
			this.extensions.set('KHR_techniques_webgl', KHR_techniques_webgl);
			this.autoParseConfig.materials.push('KHR_techniques_webgl');
		}
	}

	class I3DMParser {
		static parse(context, loader) {
			const bodyBytes = new Uint8Array(context.options.buffer, context.batchTableEnd, context.header.byteLength - context.batchTableEnd);
			let promise = null;
			if (context.header.gltfFormat === 1) {
				promise = Promise.resolve(bodyBytes);
			} else {
				const externalUri = GLTFUtils.resolveURL(GLTFUtils.decodeText(bodyBytes), context.path);
				promise = loader.loadFile(externalUri, 'arraybuffer').then(buffer => new Uint8Array(buffer));
			}
			return promise.then(glbBytes => {
				const glbData = GLTFUtils.parseGLB(glbBytes.slice().buffer);
				context.gltf = glbData.gltf;
				context.buffers = glbData.buffers;
			});
		}
	}

	const instancing_pars_vert = `
#ifdef USE_INSTANCING
	attribute mat4 instanceMatrix;
	uniform mat4 instanceOffset;
#endif
`;
	const instancing_position_vert = `
#ifdef USE_INSTANCING
	mat4 instancingMatrix = inverseMat4(instanceOffset) * instanceMatrix * instanceOffset;
	transformed = (instancingMatrix * vec4(transformed, 1.0)).xyz;
#endif
`;
	const instancing_normal_vert = `
#ifdef USE_INSTANCING
	mat4 instancingNormalMatrix = transposeMat4(inverseMat4(instancingMatrix));

	objectNormal = (instancingNormalMatrix * vec4(objectNormal, 0.0)).xyz;

	#ifdef USE_TANGENT
		objectTangent = (instancingNormalMatrix * vec4(objectTangent, 0.0)).xyz;
	#endif
#endif
`;

	// InstancedPBRMaterial

	let pbr_vert = t3d.ShaderLib.pbr_vert;
	pbr_vert = pbr_vert.replace('#include <logdepthbuf_pars_vert>', `
#include <logdepthbuf_pars_vert>
${instancing_pars_vert}
`);
	pbr_vert = pbr_vert.replace('#include <pvm_vert>', `
${instancing_position_vert}
#include <pvm_vert>
`);
	pbr_vert = pbr_vert.replace('#include <normal_vert>', `
${instancing_normal_vert}
#include <normal_vert>
`);
	class InstancedPBRMaterial extends t3d.PBRMaterial {
		constructor(sourceMaterial) {
			super();
			this.type = t3d.MATERIAL_TYPE.SHADER;
			if (sourceMaterial) {
				this.copy(sourceMaterial);
			}
			this.shaderName = 'InstancedPBR';
			this.vertexShader = pbr_vert;
			this.fragmentShader = t3d.ShaderLib.pbr_frag;
			this.defines.USE_INSTANCING = true;
			this.uniforms.instanceOffset = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);
		}
	}

	// InstancedBasicMaterial

	let basic_vert = t3d.ShaderLib.basic_vert;
	basic_vert = basic_vert.replace('#include <logdepthbuf_pars_vert>', `
	#include <logdepthbuf_pars_vert>
	${instancing_pars_vert}
`);
	basic_vert = basic_vert.replace('#include <pvm_vert>', `
	${instancing_position_vert}
	#include <pvm_vert>
`);
	class InstancedBasicMaterial extends t3d.BasicMaterial {
		constructor(sourceMaterial) {
			super();
			this.type = t3d.MATERIAL_TYPE.SHADER;
			if (sourceMaterial) {
				this.copy(sourceMaterial);
			}
			this.shaderName = 'InstancedBasic';
			this.vertexShader = basic_vert;
			this.fragmentShader = t3d.ShaderLib.basic_frag;
			this.defines.USE_INSTANCING = true;
			this.uniforms.instanceOffset = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);
		}
	}

	// InstancedDepthMaterial

	let depth_vert = t3d.ShaderLib.depth_vert;
	depth_vert = depth_vert.replace('#include <logdepthbuf_pars_vert>', `
	#include <logdepthbuf_pars_vert>
	${instancing_pars_vert}
`);
	depth_vert = depth_vert.replace('#include <pvm_vert>', `
	${instancing_position_vert}
	#include <pvm_vert>
`);

	class MaterialParser {
		static parse(context, loader) {
			const {
				gltf,
				textures
			} = context;
			if (!gltf.materials) return;
			const transformExt = loader.extensions.get('KHR_texture_transform');
			const materials = [];
			for (let i = 0; i < gltf.materials.length; i++) {
				const {
					extensions = {},
					pbrMetallicRoughness,
					normalTexture,
					occlusionTexture,
					emissiveTexture,
					emissiveFactor,
					alphaMode,
					alphaCutoff,
					doubleSided,
					name = ''
				} = gltf.materials[i];
				let material = null;
				const materialExtNames = loader.autoParseConfig.materials;

				// TODO: refactor invoke method
				for (let j = 0; j < materialExtNames.length; j++) {
					const extName = materialExtNames[j];
					const extParams = extensions[extName];
					const ext = loader.extensions.get(extName);
					if (extParams && ext && ext.getMaterial) {
						material = ext.getMaterial();
						break;
					}
				}
				material = material || new InstancedPBRMaterial(); // @parser-modification - instanced materials
				material.name = name;
				for (let j = 0; j < materialExtNames.length; j++) {
					const extName = materialExtNames[j];
					const extParams = extensions[extName];
					const ext = loader.extensions.get(extName);
					if (extParams && ext && ext.parseParams) {
						ext.parseParams(material, extParams, textures, transformExt);
					}
				}
				const {
					KHR_materials_unlit,
					KHR_materials_pbrSpecularGlossiness
				} = extensions;
				if (pbrMetallicRoughness) {
					const {
						baseColorFactor,
						baseColorTexture,
						metallicFactor,
						roughnessFactor,
						metallicRoughnessTexture
					} = pbrMetallicRoughness;
					if (Array.isArray(baseColorFactor)) {
						material.diffuse.fromArray(baseColorFactor);
						material.opacity = baseColorFactor[3] !== undefined ? baseColorFactor[3] : 1;
					}
					if (baseColorTexture) {
						material.diffuseMap = textures[baseColorTexture.index];
						material.diffuseMapCoord = baseColorTexture.texCoord || 0;
						if (material.diffuseMap) {
							material.diffuseMap.encoding = t3d.TEXEL_ENCODING_TYPE.SRGB;
							transformExt && transformExt.handleMaterialMap(material, 'diffuseMap', baseColorTexture);
						}
					}
					if (!KHR_materials_unlit && !KHR_materials_pbrSpecularGlossiness) {
						material.metalness = metallicFactor !== undefined ? metallicFactor : 1;
						material.roughness = roughnessFactor !== undefined ? roughnessFactor : 1;
						if (metallicRoughnessTexture) {
							material.metalnessMap = textures[metallicRoughnessTexture.index];
							material.roughnessMap = textures[metallicRoughnessTexture.index];
							// metallicRoughnessTexture transform not supported yet
						}
					}
				}
				if (emissiveFactor) {
					material.emissive.fromArray(emissiveFactor);
				}
				if (emissiveTexture) {
					material.emissiveMap = textures[emissiveTexture.index];
					material.emissiveMapCoord = emissiveTexture.texCoord || 0;
					if (material.emissiveMap) {
						material.emissiveMap.encoding = t3d.TEXEL_ENCODING_TYPE.SRGB;
						transformExt && transformExt.handleMaterialMap(material, 'emissiveMap', emissiveTexture);
					}
				}
				if (occlusionTexture) {
					material.aoMap = textures[occlusionTexture.index];
					material.aoMapCoord = occlusionTexture.texCoord || 0;
					if (occlusionTexture.strength !== undefined) {
						material.aoMapIntensity = occlusionTexture.strength;
					}
					if (material.aoMap) {
						transformExt && transformExt.handleMaterialMap(material, 'aoMap', occlusionTexture);
					}
				}
				if (!KHR_materials_unlit) {
					if (normalTexture) {
						material.normalMap = textures[normalTexture.index];
						material.normalScale.set(1, -1);
						if (normalTexture.scale !== undefined) {
							// fix flip y for normal map
							// https://github.com/mrdoob/three.js/issues/11438#issuecomment-507003995
							material.normalScale.set(normalTexture.scale, -normalTexture.scale);
						}

						// normal map transform not supported yet
					}
				}
				material.side = doubleSided === true ? t3d.DRAW_SIDE.DOUBLE : t3d.DRAW_SIDE.FRONT;
				if (alphaMode === ALPHA_MODES.BLEND) {
					material.transparent = true;
				} else {
					material.transparent = false;
					if (alphaMode === ALPHA_MODES.MASK) {
						material.alphaTest = alphaCutoff !== undefined ? alphaCutoff : 0.5;
					}
				}
				materials[i] = material;
			}
			context.materials = materials;
		}
	}

	class PrimitiveParser {
		static parse(context, loader) {
			const {
				gltf,
				accessors,
				materials,
				bufferViews
			} = context;
			if (!gltf.meshes) return;
			const dracoExt = loader.extensions.get('KHR_draco_mesh_compression');
			const materialCache = new Map();
			const geometryPromiseCache = new Map();
			const meshPromises = [];
			for (let i = 0; i < gltf.meshes.length; i++) {
				const gltfMesh = gltf.meshes[i];
				const primitivePromises = [];
				for (let j = 0; j < gltfMesh.primitives.length; j++) {
					const gltfPrimitive = gltfMesh.primitives[j];
					const {
						extensions = {},
						mode,
						material
					} = gltfPrimitive;
					const {
						KHR_draco_mesh_compression
					} = extensions;
					let geometryPromise;
					const geometryKey = createGeometryKey(gltfPrimitive);
					if (geometryPromiseCache.has(geometryKey)) {
						geometryPromise = geometryPromiseCache.get(geometryKey);
					} else {
						if (KHR_draco_mesh_compression && dracoExt) {
							geometryPromise = dracoExt.getGeometry(KHR_draco_mesh_compression, bufferViews, gltfPrimitive.attributes, gltf.accessors, loader.getDRACOLoader());
						} else {
							geometryPromise = Promise.resolve(new t3d.Geometry());
						}
						geometryPromise = geometryPromise.then(geometry => {
							parseGeometryFromGLTFPrimitive(geometry, gltfPrimitive, gltf, accessors);
							return geometry;
						});
						geometryPromiseCache.set(geometryKey, geometryPromise);
					}
					const primitivePromise = geometryPromise.then(geometry => {
						const primitive = {
							mode,
							geometry,
							material: material === undefined ? new InstancedPBRMaterial() : materials[material],
							// @parser-modification - instanced materials
							weights: Object.keys(geometry.morphAttributes).length > 0 && gltfMesh.weights ? gltfMesh.weights.slice(0) : undefined,
							skinned: gltfMesh.isSkinned
						};
						assignFinalMaterial(primitive, materialCache);
						return primitive;
					});
					primitivePromises.push(primitivePromise);
				}
				meshPromises.push(Promise.all(primitivePromises));
			}
			materialCache.clear();
			geometryPromiseCache.clear();
			return Promise.all(meshPromises).then(primitives => {
				context.primitives = primitives;
			});
		}
	}
	function parseGeometryFromGLTFPrimitive(geometry, gltfPrimitive, gltf, accessors) {
		const {
			attributes,
			indices,
			targets
		} = gltfPrimitive;

		// set attributes

		for (const attributeSemantic in attributes) {
			const accessorIdx = attributes[attributeSemantic];
			const attributeName = ATTRIBUTES[attributeSemantic] === undefined ? attributeSemantic : ATTRIBUTES[attributeSemantic];
			// Skip attributes already provided by e.g. Draco extension.
			if (attributeName in geometry.attributes) continue;
			geometry.addAttribute(attributeName, accessors[accessorIdx]);
		}

		// set index

		if (indices !== undefined && !geometry.index) {
			geometry.setIndex(accessors[indices]);
		}

		// compute bounds

		const {
			boundingBox,
			boundingSphere
		} = geometry;
		if (attributes.POSITION !== undefined) {
			const accessorIdx = attributes.POSITION;
			const accessor = gltf.accessors[accessorIdx];
			if (accessor.min && accessor.max) {
				boundingBox.min.fromArray(accessor.min);
				boundingBox.max.fromArray(accessor.max);
				if (accessor.normalized) {
					const boxScale = GLTFUtils.getNormalizedComponentScale(ACCESSOR_COMPONENT_TYPES[accessor.componentType]);
					boundingBox.min.multiplyScalar(boxScale);
					boundingBox.max.multiplyScalar(boxScale);
				}
			} else {
				geometry.computeBoundingBox();
			}
		} else {
			geometry.computeBoundingBox();
		}
		boundingBox.getCenter(boundingSphere.center);
		boundingSphere.radius = boundingBox.min.distanceTo(boundingBox.max) / 2;

		// set morph targets

		if (targets) {
			let hasMorphPosition = false;
			let hasMorphNormal = false;
			for (let i = 0, il = targets.length; i < il; i++) {
				const target = targets[i];
				if (target.POSITION !== undefined) hasMorphPosition = true;
				if (target.NORMAL !== undefined) hasMorphNormal = true;
				if (hasMorphPosition && hasMorphNormal) break;
			}
			if (hasMorphPosition || hasMorphNormal) {
				const morphPositions = [];
				const morphNormals = [];
				for (let i = 0, il = targets.length; i < il; i++) {
					const target = targets[i];
					if (hasMorphPosition) {
						morphPositions.push(target.POSITION !== undefined ? accessors[target.POSITION] : geometry.attributes[ATTRIBUTES.POSITION]);
					}
					if (hasMorphNormal) {
						morphNormals.push(target.NORMAL !== undefined ? accessors[target.NORMAL] : geometry.attributes[ATTRIBUTES.NORMAL]);
					}
				}
				if (hasMorphPosition) {
					geometry.morphAttributes.position = morphPositions;
				}
				if (hasMorphNormal) {
					geometry.morphAttributes.normal = morphNormals;
				}
			}
		}
		return geometry;
	}
	function assignFinalMaterial(primitive, materialCache) {
		let {
			geometry,
			material,
			skinned,
			mode
		} = primitive;

		// If the material will be modified later on, clone it now.
		const useVertexTangents = geometry.attributes[ATTRIBUTES.TANGENT] !== undefined;
		const useVertexColors = geometry.attributes[ATTRIBUTES.COLOR_0] !== undefined;
		const useFlatShading = geometry.attributes[ATTRIBUTES.NORMAL] === undefined;
		const useSkinning = skinned;
		if (mode === WEBGL_DRAW_MODES.POINTS) {
			const cacheKey = 'PointsMaterial:' + material.id;
			let pointsMaterial = materialCache.get(cacheKey);
			if (!pointsMaterial) {
				pointsMaterial = new t3d.PointsMaterial();
				t3d.Material.prototype.copy.call(pointsMaterial, material);
				pointsMaterial.diffuse.copy(material.diffuse);
				pointsMaterial.diffuseMap = material.map;
				pointsMaterial.drawMode = mode;
				pointsMaterial.acceptLight = false; // PointsMaterial doesn't support lights yet
				materialCache.set(cacheKey, pointsMaterial);
			}
			material = pointsMaterial;
		} else if (mode === WEBGL_DRAW_MODES.LINES || mode === WEBGL_DRAW_MODES.LINE_STRIP || mode === WEBGL_DRAW_MODES.LINE_LOOP) {
			const cacheKey = 'BasicMaterial:' + material.id;
			let basicMaterial = materialCache.get(cacheKey);
			if (!basicMaterial) {
				basicMaterial = new t3d.BasicMaterial();
				basicMaterial.envMap = undefined; // force close env map
				basicMaterial.diffuse.copy(material.diffuse);
				basicMaterial.diffuseMap = material.diffuseMap;
				basicMaterial.drawMode = mode;
				materialCache.set(cacheKey, basicMaterial);
			}
			material = basicMaterial;
		} else if (mode === WEBGL_DRAW_MODES.TRIANGLE_STRIP) {
			// TODO
			console.warn('TRIANGLE_STRIP will be removed later.');
			material.drawMode = WEBGL_DRAW_MODES.TRIANGLE_STRIP;
		} else if (mode === WEBGL_DRAW_MODES.TRIANGLE_FAN) {
			// TODO
			console.warn('TRIANGLE_FAN will be removed later.');
			material.drawMode = WEBGL_DRAW_MODES.TRIANGLE_FAN;
		}
		if (useVertexTangents || useVertexColors || useFlatShading || useSkinning) {
			let cacheKey = 'ClonedMaterial:' + material.id + ':';
			if (useVertexTangents) cacheKey += 'vertex-tangents:';
			if (useVertexColors) {
				if (geometry.attributes[ATTRIBUTES.COLOR_0].size === 3) {
					cacheKey += 'vertex-colors-rgb:';
				} else if (geometry.attributes[ATTRIBUTES.COLOR_0].size === 4) {
					cacheKey += 'vertex-colors-rgba:';
				}
			}
			if (useFlatShading) cacheKey += 'flat-shading:';
			let cachedMaterial = materialCache.get(cacheKey);
			if (!cachedMaterial) {
				cachedMaterial = material.clone();
				if (useVertexTangents) {
					cachedMaterial.vertexTangents = true;

					// revert flip y fix for tangents
					// https://github.com/mrdoob/three.js/issues/11438#issuecomment-507003995
					if (cachedMaterial.normalMap) {
						cachedMaterial.normalScale.y *= -1;
					}
				}
				if (useVertexColors) {
					if (geometry.attributes[ATTRIBUTES.COLOR_0].size === 3) {
						cachedMaterial.vertexColors = t3d.VERTEX_COLOR.RGB;
					} else if (geometry.attributes[ATTRIBUTES.COLOR_0].size === 4) {
						cachedMaterial.vertexColors = t3d.VERTEX_COLOR.RGBA;
					} else {
						console.warn('Illegal vertex color size: ' + geometry.attributes[ATTRIBUTES.COLOR_0].size);
					}
				}
				if (useFlatShading) {
					cachedMaterial.shading = t3d.SHADING_TYPE.FLAT_SHADING;
				}
			}
			material = cachedMaterial;
		}
		primitive.material = material;
	}
	function createGeometryKey(primitive) {
		const dracoExtension = primitive.extensions && primitive.extensions.KHR_draco_mesh_compression;
		let geometryKey;
		if (dracoExtension) {
			geometryKey = 'draco:' + dracoExtension.bufferView + ':' + dracoExtension.indices + ':' + createAttributesKey(dracoExtension.attributes);
		} else {
			geometryKey = primitive.indices + ':' + createAttributesKey(primitive.attributes) + ':' + primitive.mode;
		}
		if (primitive.targets) {
			for (let i = 0, il = primitive.targets.length; i < il; i++) {
				geometryKey += ':' + createAttributesKey(primitive.targets[i]);
			}
		}
		return geometryKey;
	}
	function createAttributesKey(attributes) {
		let attributesKey = '';
		const keys = Object.keys(attributes).sort();
		for (let i = 0, il = keys.length; i < il; i++) {
			attributesKey += keys[i] + ':' + attributes[keys[i]] + ';';
		}
		return attributesKey;
	}

	class I3DMRootParser {
		static parse(context, loader) {
			const {
				featureTable,
				root,
				options
			} = context;
			const INSTANCES_LENGTH = featureTable.getData('INSTANCES_LENGTH');
			const POSITION = featureTable.getData('POSITION', INSTANCES_LENGTH, 'FLOAT', 'VEC3');
			const NORMAL_UP = featureTable.getData('NORMAL_UP', INSTANCES_LENGTH, 'FLOAT', 'VEC3');
			const NORMAL_RIGHT = featureTable.getData('NORMAL_RIGHT', INSTANCES_LENGTH, 'FLOAT', 'VEC3');
			const SCALE = featureTable.getData('SCALE', INSTANCES_LENGTH, 'FLOAT', 'SCALAR');
			const SCALE_NON_UNIFORM = featureTable.getData('SCALE_NON_UNIFORM', INSTANCES_LENGTH, 'FLOAT', 'VEC3');

			// check unsupported features

			[
			// Global Properties
			'QUANTIZED_VOLUME_OFFSET', 'QUANTIZED_VOLUME_SCALE', 'EAST_NORTH_UP',
			// Per-instance Properties
			'POSITION_QUANTIZED', 'NORMAL_UP_OCT32P', 'NORMAL_RIGHT_OCT32P'].forEach(feature => {
				if (feature in featureTable.header) {
					console.warn(`I3DMLoader: Unsupported FeatureTable feature "${feature}" detected.`);
				}
			});

			// set instance matrix for all geometries

			const averageVector = new t3d.Vector3();
			for (let i = 0; i < INSTANCES_LENGTH; i++) {
				averageVector.x += POSITION[i * 3 + 0] / INSTANCES_LENGTH;
				averageVector.y += POSITION[i * 3 + 1] / INSTANCES_LENGTH;
				averageVector.z += POSITION[i * 3 + 2] / INSTANCES_LENGTH;
			}
			const instances = [];
			root.traverse(child => {
				if (child.isMesh) {
					const {
						geometry
					} = child;
					geometry.instanceCount = INSTANCES_LENGTH;
					const instanceMatrix = new t3d.Attribute(new t3d.Buffer(new Float32Array(INSTANCES_LENGTH * 16), 16), 16);
					instanceMatrix.divisor = 1;
					geometry.addAttribute('instanceMatrix', instanceMatrix);

					// Center the instance around an average point to avoid jitter at large scales.
					// Transform the average vector by matrix world so we can account for any existing
					// transforms of the instanced mesh.
					child.updateMatrix(true);
					child.position.copy(averageVector).applyMatrix4(child.worldMatrix);
					instances.push(child);
				}
			});
			for (let i = 0; i < INSTANCES_LENGTH; i++) {
				// position
				tempPos.fromArray(POSITION, i * 3).sub(averageVector);

				// rotation
				if (NORMAL_UP) {
					tempUp.fromArray(NORMAL_UP, i * 3);
					tempRight.fromArray(NORMAL_RIGHT, i * 3);
					tempFwd.crossVectors(tempRight, tempUp).normalize();
					tempMat$1.set(tempRight.x, tempUp.x, tempFwd.x, 0, tempRight.y, tempUp.y, tempFwd.y, 0, tempRight.z, tempUp.z, tempFwd.z, 0, 0, 0, 0, 1);
					tempQuat.setFromRotationMatrix(tempMat$1);
				} else {
					tempQuat.set(0, 0, 0, 1);
				}

				// scale
				if (SCALE) {
					tempSca.set(SCALE[i], SCALE[i], SCALE[i]);
				} else if (SCALE_NON_UNIFORM) {
					tempSca.fromArray(SCALE_NON_UNIFORM, i * 3);
				} else {
					tempSca.set(1, 1, 1);
				}

				// TODO instance matrix should be applied to model root
				tempMat$1.transform(tempPos, tempSca, tempQuat).multiply(options.adjustmentTransform);
				for (let j = 0, l = instances.length; j < l; j++) {
					const {
						geometry
					} = instances[j];
					const instanceArray = geometry.getAttribute('instanceMatrix').buffer.array;
					tempMat$1.toArray(instanceArray, i * 16);
					geometry.version++;
				}
			}

			// fix rtc center

			const rtcCenter = featureTable.getData('RTC_CENTER');
			if (rtcCenter) {
				root.position.x += rtcCenter[0];
				root.position.y += rtcCenter[1];
				root.position.z += rtcCenter[2];
			}
		}
	}
	const tempFwd = new t3d.Vector3();
	const tempUp = new t3d.Vector3();
	const tempRight = new t3d.Vector3();
	const tempPos = new t3d.Vector3();
	const tempQuat = new t3d.Quaternion();
	const tempSca = new t3d.Vector3();
	const tempMat$1 = new t3d.Matrix4();

	class KHR_techniques_webgl_i extends KHR_techniques_webgl {
		static getMaterial() {
			return new InstancedPBRMaterial();
		}
	}

	class KHR_materials_unlit_i {
		static getMaterial() {
			return new InstancedBasicMaterial();
		}
	}

	class KHR_materials_pbrSpecularGlossiness_i extends KHR_materials_pbrSpecularGlossiness {
		static getMaterial() {
			const material = new InstancedPBRMaterial(); // fallback to InstancedPBRMaterial
			material.specular = new t3d.Color3(0x111111);
			return material;
		}
	}

	class KHR_materials_clearcoat_i extends KHR_materials_clearcoat {
		static getMaterial() {
			return new InstancedPBRMaterial();
		}
	}

	/**
	 * I3DMLoader is a loader for the I3DM format.
	 */
	class I3DMLoader extends GLTFLoader {
		constructor(manager) {
			super(manager, [HeaderParser,
			// insert HeaderParser
			TableParser,
			// insert TableParser
			I3DMParser,
			// insert I3DMParser
			ReferenceParser, Validator, BufferParser, BufferViewParser, ImageParser, TextureParser, MaterialParser,
			// replace MaterialParser
			AccessorParser, PrimitiveParser,
			// replace PrimitiveParser
			NodeParser, SkinParser, SceneParser, AnimationParser, I3DMRootParser // insert I3DMSetupParser
			]);
			this.extensions.set('KHR_techniques_webgl', KHR_techniques_webgl_i);
			this.extensions.set('KHR_materials_unlit', KHR_materials_unlit_i);
			this.extensions.set('KHR_materials_pbrSpecularGlossiness', KHR_materials_pbrSpecularGlossiness_i);
			this.extensions.set('KHR_materials_clearcoat', KHR_materials_clearcoat_i);
			this.autoParseConfig.materials.push('KHR_techniques_webgl');
		}
	}

	class PNTSRootParser {
		static parse(context, loader) {
			const {
				featureTable
			} = context;
			const POINTS_LENGTH = featureTable.getData('POINTS_LENGTH');
			const POSITION = featureTable.getData('POSITION', POINTS_LENGTH, 'FLOAT', 'VEC3');
			const RGB = featureTable.getData('RGB', POINTS_LENGTH, 'UNSIGNED_BYTE', 'VEC3');
			const RGBA = featureTable.getData('RGBA', POINTS_LENGTH, 'UNSIGNED_BYTE', 'VEC4');

			// check unsupported features

			[
			// Global Properties
			'QUANTIZED_VOLUME_OFFSET', 'QUANTIZED_VOLUME_SCALE', 'CONSTANT_RGBA', 'BATCH_LENGTH',
			// Per-point Properties
			'POSITION_QUANTIZED', 'RGB565', 'NORMAL', 'NORMAL_OCT16P', 'BATCH_ID'].forEach(feature => {
				if (feature in featureTable.header) {
					console.warn(`PNTSLoader: Unsupported FeatureTable feature "${feature}" detected.`);
				}
			});

			// generate root

			const geometry = new t3d.Geometry();
			geometry.addAttribute('a_Position', new t3d.Attribute(new t3d.Buffer(POSITION, 3), 3, 0, true));
			geometry.computeBoundingBox();
			geometry.computeBoundingSphere();
			const material = new t3d.PointsMaterial();
			material.size = 2;
			material.sizeAttenuation = false;
			if (RGB !== null) {
				geometry.addAttribute('a_Color', new t3d.Attribute(new t3d.Buffer(RGB, 3), 3, 0, true));
				material.vertexColors = t3d.VERTEX_COLOR.RGB;
			} else if (RGBA !== null) {
				geometry.addAttribute('a_Color', new t3d.Attribute(new t3d.Buffer(RGBA, 4), 4, 0, true));
				material.vertexColors = t3d.VERTEX_COLOR.RGBA;
			}
			const root = new t3d.Mesh(geometry, material);

			// output mesh to root

			context.root = root;

			// fix rtc center

			const rtcCenter = featureTable.getData('RTC_CENTER');
			if (rtcCenter) {
				root.position.x += rtcCenter[0];
				root.position.y += rtcCenter[1];
				root.position.z += rtcCenter[2];
			}
		}
	}

	/**
	 * PNTSLoader is a loader for the PNTS format.
	 */
	class PNTSLoader extends GLTFLoader {
		constructor(manager) {
			super(manager, [HeaderParser, TableParser, PNTSRootParser]);
		}
	}

	class CMPTParser {
		static parse(context, loader) {
			const buffer = context.options.buffer;
			const tilesLength = context.header.tilesLength;
			const tiles = [];
			let offset = 16;
			for (let i = 0; i < tilesLength; i++) {
				const tileView = new DataView(buffer, offset, 12);
				const tileMagic = String.fromCharCode(tileView.getUint8(0)) + String.fromCharCode(tileView.getUint8(1)) + String.fromCharCode(tileView.getUint8(2)) + String.fromCharCode(tileView.getUint8(3));
				const tileVersion = tileView.getUint32(4, true);
				const byteLength = tileView.getUint32(8, true);
				const tileBuffer = new Uint8Array(buffer, offset, byteLength);
				tiles.push({
					type: tileMagic,
					buffer: tileBuffer,
					version: tileVersion
				});
				offset += byteLength;
			}
			context.tiles = tiles;
		}
	}

	class CMPTRootParser {
		static parse(context, loader) {
			const {
				tiles,
				options,
				path
			} = context;
			const adjustmentTransform = options.adjustmentTransform;
			const promises = [];
			for (const i in tiles) {
				const {
					type,
					buffer
				} = tiles[i];
				const config = {
					fetchOptions: options.fetchOptions,
					path,
					buffer: buffer.slice().buffer
				};
				if (type === 'b3dm' || type === 'i3dm') {
					config.adjustmentTransform = adjustmentTransform;
				}
				const _loader = loader._loaders.get(type);
				if (_loader) {
					promises.push(_loader.load(`${path}/temp.${type}`, config));
				}
			}
			return Promise.all(promises).then(results => {
				const group = new t3d.Object3D();
				results.forEach(result => {
					group.add(result.root);
				});
				return {
					tiles: results,
					root: group
				};
			});
		}
	}

	/**
	 * CMPTLoader is a loader for the CMPT format.
	 */
	class CMPTLoader extends GLTFLoader {
		constructor(manager) {
			super(manager, [HeaderParser, CMPTParser, CMPTRootParser]);
			const b3dmLoader = new B3DMLoader(manager);
			const i3dmLoader = new I3DMLoader(manager);
			const pntsLoader = new PNTSLoader(manager);
			this._loaders = new Map([['b3dm', b3dmLoader], ['i3dm', i3dmLoader], ['pnts', pntsLoader]]);
		}
		setDRACOLoader(dracoLoader) {
			for (const loader of this._loaders.values()) {
				loader.setDRACOLoader(dracoLoader);
			}
			return super.setDRACOLoader(dracoLoader);
		}
		setKTX2Loader(ktx2Loader) {
			for (const loader of this._loaders.values()) {
				loader.setKTX2Loader(ktx2Loader);
			}
			return super.setKTX2Loader(ktx2Loader);
		}
	}

	class IndexParser {
		static parse(context, loader) {
			const {
				url,
				options
			} = context;
			const buffer = options.buffer;
			const _isGLB = isGLB(url);
			if (_isGLB) {
				const glbData = GLTFUtils.parseGLB(buffer);
				context.gltf = glbData.gltf;
				context.buffers = glbData.buffers;
			} else {
				context.gltf = buffer;
			}
		}
	}
	const isGLB = url => {
		return getUrlExtension(url) === 'glb';
	};

	/**
	 * TileGLTFLoader is a gltf loader
	 * that extends the default gltf loader with additional parsers
	 * to support the tile gltf format.
	 */
	class TileGLTFLoader extends GLTFLoader {
		constructor(manager) {
			super(manager);
			this.replaceParser(IndexParser, 0);
		}
	}

	function readMagicBytes(bufferOrDataView) {
		if (bufferOrDataView === null || bufferOrDataView.byteLength < 4) {
			return '';
		}
		let view;
		if (bufferOrDataView instanceof DataView) {
			view = bufferOrDataView;
		} else {
			view = new DataView(bufferOrDataView);
		}
		if (String.fromCharCode(view.getUint8(0)) === '{') {
			return null;
		}
		let magicBytes = '';
		for (let i = 0; i < 4; i++) {
			magicBytes += String.fromCharCode(view.getUint8(i));
		}
		return magicBytes;
	}

	const _updateBeforeEvent = {
		type: 'update-before'
	};
	const _updateAfterEvent = {
		type: 'update-after'
	};
	const _tilesLoadStartEvent = {
		type: 'tiles-load-start'
	};
	const _tilesLoadEndEvent = {
		type: 'tiles-load-end'
	};
	const PLUGIN_REGISTERED = Symbol('PLUGIN_REGISTERED');
	const INITIAL_FRUSTUM_CULLED = Symbol('INITIAL_FRUSTUM_CULLED');
	const tempMat = new t3d.Matrix4();
	const viewErrorTarget = {
		inView: false,
		error: Infinity
	};
	const X_AXIS = new t3d.Vector3(1, 0, 0);
	const Y_AXIS = new t3d.Vector3(0, 1, 0);
	function updateFrustumCulled(object, toInitialValue) {
		object.traverse(c => {
			c.frustumCulled = c[INITIAL_FRUSTUM_CULLED] && toInitialValue;
		});
	}
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
	class Tiles3D extends t3d.Object3D {
		get root() {
			const rootTileSet = this.rootTileSet;
			return rootTileSet ? rootTileSet.root : null;
		}
		get loadProgress() {
			const {
				stats,
				isLoading
			} = this;
			const loading = stats.downloading + stats.parsing;
			const total = stats.inCacheSinceLoad + (isLoading ? 1 : 0);
			return total === 0 ? 1.0 : 1.0 - loading / total;
		}
		constructor(url, manager = new t3d.LoadingManager()) {
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
				inCacheSinceLoad: 0,
				inCache: 0,
				parsing: 0,
				downloading: 0,
				failed: 0,
				inFrustum: 0,
				used: 0,
				active: 0,
				visible: 0
			};
			this.frameCount = 0;

			// callbacks
			this._dispatchNeedsUpdateEvent = throttle(() => {
				this.dispatchEvent({
					type: 'needs-update'
				});
			});
			this.activeTiles = new Set();
			this.visibleTiles = new Set();
			this.usedSet = new Set();
			this.rootLoadingState = UNLOADED;

			// internals

			this.rootURL = url;
			this.rootTileSet = null;
			this._autoDisableRendererCulling = true;
			this.plugins = [];
			this.queuedTiles = [];
			this.cachedSinceLoadComplete = new Set();
			this.isLoading = false;
			const lruCache = new LRUCache();
			lruCache.unloadPriorityCallback = lruPriorityCallback;
			const downloadQueue = new PriorityQueue();
			downloadQueue.maxJobs = 10;
			downloadQueue.priorityCallback = priorityCallback;
			const parseQueue = new PriorityQueue();
			parseQueue.maxJobs = 1;
			parseQueue.priorityCallback = priorityCallback;
			const processNodeQueue = new PriorityQueue();
			processNodeQueue.maxJobs = 25;
			processNodeQueue.priorityCallback = priorityCallback;
			processNodeQueue.log = true;
			this.lruCache = lruCache;
			this.downloadQueue = downloadQueue;
			this.parseQueue = parseQueue;
			this.processNodeQueue = processNodeQueue;
			this.$cameras = new CameraList();
			this.lruCache.computeMemoryUsageCallback = tile => tile.cached.bytesUsed ?? null;
			const b3dmLoader = new B3DMLoader(manager);
			const i3dmLoader = new I3DMLoader(manager);
			const pntsLoader = new PNTSLoader(manager);
			const cmptLoader = new CMPTLoader(manager);
			const gltfLoader = new TileGLTFLoader(manager);
			this._loaders = new Map([['b3dm', b3dmLoader], ['i3dm', i3dmLoader], ['pnts', pntsLoader], ['cmpt', cmptLoader], ['gltf', gltfLoader]]);
			this._upRotationMatrix = new t3d.Matrix4();
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
		queueTileForDownload(tile) {
			if (tile.__loadingState !== UNLOADED) {
				return;
			}
			this.queuedTiles.push(tile);
		}
		markTileUsed(tile) {
			// save the tile in a separate "used set" so we can mark it as unused
			// before the next tile set traversal
			this.usedSet.add(tile);
			this.lruCache.markUsed(tile);
		}

		// Public API
		update() {
			const {
				lruCache,
				usedSet,
				stats,
				root,
				downloadQueue,
				parseQueue,
				processNodeQueue
			} = this;
			if (this.rootLoadingState === UNLOADED) {
				this.rootLoadingState = LOADING;
				this.invokeOnePlugin(plugin => plugin.loadRootTileSet && plugin.loadRootTileSet()).then(root => {
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
				}).catch(error => {
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

			// TODO: This will only sort for one tile set. We may want to store this queue on the
			// LRUCache so multiple tile sets can use it at once
			// start the downloads of the tiles as needed
			const queuedTiles = this.queuedTiles;
			queuedTiles.sort(lruCache.unloadPriorityCallback);
			for (let i = 0, l = queuedTiles.length; i < l && !lruCache.isFull(); i++) {
				this.requestTileContents(queuedTiles[i]);
			}
			queuedTiles.length = 0;
			this.lruCache.scheduleUnload();

			// if all tasks have finished and we've been marked as actively loading then fire the completion event
			const runningTasks = downloadQueue.running || parseQueue.running || processNodeQueue.running;
			if (runningTasks === false && this.isLoading === true) {
				this.cachedSinceLoadComplete.clear();
				stats.inCacheSinceLoad = 0;
				this.dispatchEvent(_tilesLoadEndEvent);
				this.isLoading = false;
			}
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
		dispatchEvent(...args) {
			t3d.EventDispatcher.prototype.dispatchEvent.call(this, ...args);
		}
		fetchData(url, options) {
			return fetch(url, options);
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
				case 'b3dm':
					{
						promise = this._loaders.get('b3dm').load(uri, {
							fetchOptions,
							path: workingPath,
							buffer,
							adjustmentTransform: upRotationMatrix.clone()
						});
						break;
					}
				case 'pnts':
					{
						promise = this._loaders.get('pnts').load(uri, {
							fetchOptions,
							path: workingPath,
							buffer
						});
						break;
					}
				case 'i3dm':
					{
						promise = this._loaders.get('i3dm').load(uri, {
							fetchOptions,
							path: workingPath,
							buffer,
							adjustmentTransform: upRotationMatrix.clone()
						});
						break;
					}
				case 'cmpt':
					{
						promise = this._loaders.get('i3dm').load(uri, {
							fetchOptions,
							path: workingPath,
							buffer,
							adjustmentTransform: upRotationMatrix.clone()
						});
						break;
					}
				case 'gltf':
				case 'glb':
					{
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
							const {
								root: scene
							} = result;
							scene.matrix.multiply(upRotationMatrix).decompose(scene.position, scene.quaternion, scene.scale);
							return result;
						});
						break;
					}
				default:
					{
						promise = this.invokeOnePlugin(plugin => plugin.parseToMesh && plugin.parseToMesh(buffer, tile, extension, uri, abortSignal));
						break;
					}
			}

			// wait for the tile to load
			const result = await promise;
			if (result === null) {
				throw new Error(`Tiles3D: Content type "${fileType}" not supported.`);
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

			// wait for extra processing by plugins if needed
			await this.invokeAllPlugins(plugin => {
				return plugin.processTileModel && plugin.processTileModel(scene, tile);
			});

			// ensure the matrix is up to date in case the scene has a transform applied
			scene.updateMatrix();
			scene.matrix.premultiply(cachedTransform);
			scene.matrix.decompose(scene.position, scene.quaternion, scene.scale);
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
			cached.materials = materials;
			cached.geometry = geometry;
			cached.textures = textures;
			cached.scene = scene;
			cached.metadata = metadata;
			// cached.bytesUsed = estimateBytesUsed(scene);
			cached.featureTable = result.featureTable;
			cached.batchTable = result.batchTable;
		}
		disposeTile(tile) {
			// TODO: are these necessary? Are we disposing tiles when they are currently visible?
			if (tile.__visible) {
				this.invokeOnePlugin(plugin => plugin.setTileVisible && plugin.setTileVisible(tile, false));
				tile.__visible = false;
			}
			if (tile.__active) {
				this.invokeOnePlugin(plugin => plugin.setTileActive && plugin.setTileActive(tile, false));
				tile.__active = false;
			}

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
				if (tile.content.boundingVolume && !('box' in tile.content.boundingVolume || 'sphere' in tile.content.boundingVolume || 'region' in tile.content.boundingVolume)) {
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
				tile.__depthFromRenderedParent = tile.__hasRenderableContent ? 1 : 0;
				tile.refine = tile.refine || 'REPLACE';
			} else {
				// increment the "depth from parent" when we encounter a new tile with content
				tile.__depth = parentTile.__depth + 1;
				tile.__depthFromRenderedParent = parentTile.__depthFromRenderedParent + (tile.__hasRenderableContent ? 1 : 0);
				tile.refine = tile.refine || parentTile.refine;
			}
			tile.__basePath = tileSetDir;
			tile.__lastFrameVisited = -1;
			this.invokeAllPlugins(plugin => {
				plugin !== this && plugin.preprocessNode && plugin.preprocessNode(tile, tileSetDir, parentTile);
			});

			// cached

			const transform = new t3d.Matrix4();
			if (tile.transform) {
				transform.fromArray(tile.transform);
			}
			if (parentTile) {
				transform.premultiply(parentTile.cached.transform);
			}
			const transformInverse = new t3d.Matrix4().copy(transform).inverse();
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
			const children = tile.children;
			for (let i = 0, l = children.length; i < l; i++) {
				const child = children[i];
				if ('__depth' in child) {
					// the child has already been processed
					break;
				} else if (immediate) {
					// process the node immediately and make sure we don't double process it
					this.processNodeQueue.remove(child);
					this.preprocessNode(child, tile.__basePath, tile);
				} else {
					// queue the node for processing if it hasn't been already
					if (!this.processNodeQueue.has(child)) {
						this.processNodeQueue.add(child, child => {
							this.preprocessNode(child, tile.__basePath, tile);
							this._dispatchNeedsUpdateEvent();
						});
					}
				}
			}
		}

		// Private Functions
		preprocessTileSet(json, url, parent = null) {
			const version = json.asset.version;
			const [major, minor] = version.split('.').map(v => parseInt(v));
			console.assert(major <= 1, 'Tiles3D: asset.version is expected to be a 1.x or a compatible version.');
			if (major === 1 && minor > 0) {
				console.warn('Tiles3D: tiles versions at 1.1 or higher have limited support. Some new extensions and features may not be supported.');
			}

			// remove the last file path path-segment from the URL including the trailing slash
			let basePath = url.replace(/\/[^/]*$/, '');
			basePath = new URL(basePath, window.location.href).toString();
			this.preprocessNode(json.root, basePath, parent);
		}
		loadRootTileSet() {
			// transform the url
			let processedUrl = this.rootURL;
			this.invokeAllPlugins(plugin => processedUrl = plugin.preprocessURL ? plugin.preprocessURL(processedUrl, null) : processedUrl);

			// load the tile set root
			const pr = this.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(processedUrl, this.fetchOptions)).then(res => {
				if (res.ok) {
					return res.json();
				} else {
					throw new Error(`Tiles3D: Failed to load tileset "${processedUrl}" with status ${res.status} : ${res.statusText}`);
				}
			}).then(root => {
				// cache the gltf tile set rotation matrix
				const {
					asset,
					extensions = {}
				} = root;
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
					const {
						ellipsoid
					} = this;
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
			let isExternalTileSet = false;
			let externalTileset = null;
			let uri = new URL(tile.content.uri, tile.__basePath + '/').toString();
			this.invokeAllPlugins(plugin => uri = plugin.preprocessURL ? plugin.preprocessURL(uri, tile) : uri);
			const stats = this.stats;
			const lruCache = this.lruCache;
			const downloadQueue = this.downloadQueue;
			const parseQueue = this.parseQueue;
			const extension = getUrlExtension(uri);

			// track an abort controller and pass-through the below conditions if aborted
			const controller = new AbortController();
			const signal = controller.signal;
			const addedSuccessfully = lruCache.add(tile, t => {
				// Stop the load if it's started
				controller.abort();
				if (isExternalTileSet) {
					t.children.length = 0;
					t.__childrenProcessed = 0;
				} else {
					this.invokeAllPlugins(plugin => {
						plugin.disposeTile && plugin.disposeTile(t);
					});
				}

				// Decrement stats
				stats.inCache--;
				if (this.cachedSinceLoadComplete.has(tile)) {
					this.cachedSinceLoadComplete.delete(tile);
					stats.inCacheSinceLoad--;
				}
				if (t.__loadingState === LOADING) {
					stats.downloading--;
				} else if (t.__loadingState === PARSING) {
					stats.parsing--;
				}
				t.__loadingState = UNLOADED;
				downloadQueue.remove(t);
				parseQueue.remove(t);
			});

			// if we couldn't add the tile to the lru cache because it's full then skip
			if (!addedSuccessfully) {
				return;
			}

			// check if this is the beginning of a new set of tiles to load and dispatch and event
			if (!this.isLoading) {
				this.isLoading = true;
				this.dispatchEvent(_tilesLoadStartEvent);
			}
			this.cachedSinceLoadComplete.add(tile);
			stats.inCacheSinceLoad++;
			stats.inCache++;
			stats.downloading++;
			tile.__loadingState = LOADING;
			downloadQueue.add(tile, downloadTile => {
				if (signal.aborted) {
					return Promise.resolve();
				}
				const res = this.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(uri, {
					...this.fetchOptions,
					signal
				}));
				this.dispatchEvent({
					type: 'tile-download-start',
					tile
				});
				return res;
			}).then(res => {
				if (signal.aborted) {
					return;
				}
				if (!(res instanceof Response)) {
					return res;
				} else if (res.ok) {
					return extension === 'json' || extension === 'gltf' ? res.json() : res.arrayBuffer();
				} else {
					throw new Error(`Failed to load model with error code ${res.status}`);
				}
			}).then(content => {
				// if it has been unloaded then the tile has been disposed
				if (signal.aborted) {
					return;
				}
				stats.downloading--;
				stats.parsing++;
				tile.__loadingState = PARSING;
				return parseQueue.add(tile, parseTile => {
					// if it has been unloaded then the tile has been disposed
					if (signal.aborted) {
						return Promise.resolve();
					}
					if (extension === 'json' && content.root) {
						this.preprocessTileSet(content, uri, tile);
						tile.children.push(content.root);
						externalTileset = content;
						isExternalTileSet = true;
						return Promise.resolve();
					} else {
						return this.invokeOnePlugin(plugin => plugin.parseTile && plugin.parseTile(content, parseTile, extension, uri, signal));
					}
				});
			}).then(() => {
				// if it has been unloaded then the tile has been disposed
				if (signal.aborted) {
					return;
				}
				stats.parsing--;
				tile.__loadingState = LOADED;
				lruCache.setLoaded(tile, true);

				// If the memory of the item hasn't been registered yet then that means the memory usage hasn't
				// been accounted for by the cache yet so we need to check if it fits or if we should remove it.
				if (lruCache.getMemoryUsage(tile) === null) {
					if (lruCache.isFull() && lruCache.computeMemoryUsageCallback(tile) > 0) {
						// And if the cache is full due to newly loaded memory then lets discard this tile - it will
						// be loaded again later from the disk cache if needed.
						lruCache.remove(tile);
					} else {
						// Otherwise update the item to the latest known value
						lruCache.updateMemoryUsage(tile);
					}
				}

				// dispatch an event indicating that this model has completed and that a new
				// call to "update" is needed.
				this.dispatchEvent({
					type: 'needs-update'
				});
				this.dispatchEvent({
					type: 'load-content'
				});
				if (isExternalTileSet) {
					this.dispatchEvent({
						type: 'load-tile-set',
						tileSet: externalTileset,
						url: uri
					});
				}
				if (tile.cached.scene) {
					this.dispatchEvent({
						type: 'load-model',
						scene: tile.cached.scene,
						tile
					});
				}
			}).catch(error => {
				// if it has been unloaded then the tile has been disposed
				if (signal.aborted) {
					return;
				}
				if (error.name !== 'AbortError') {
					downloadQueue.remove(tile);
					parseQueue.remove(tile);
					if (tile.__loadingState === PARSING) {
						stats.parsing--;
					} else if (tile.__loadingState === LOADING) {
						stats.downloading--;
					}
					stats.failed++;
					console.error(`Tiles3D: Failed to load tile at url "${tile.content.uri}".`);
					console.error(error);
					tile.__loadingState = FAILED;
					lruCache.setLoaded(tile, true);
					this.dispatchEvent({
						type: 'load-error',
						tile,
						error,
						url: uri
					});
				} else {
					lruCache.remove(tile);
				}
			});
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
				this._autoDisableRendererCulling = value;
				this.forEachLoadedModel(scene => {
					updateFrustumCulled(scene, !value);
				});
			}
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
		addEventListener(...args) {
			t3d.EventDispatcher.prototype.addEventListener.call(this, ...args);
		}
		hasEventListener(...args) {
			t3d.EventDispatcher.prototype.hasEventListener.call(this, ...args);
		}
		removeEventListener(...args) {
			t3d.EventDispatcher.prototype.removeEventListener.call(this, ...args);
		}
		addCamera(camera) {
			const success = this.$cameras.add(camera);
			if (success) {
				this.dispatchEvent({
					type: 'add-camera',
					camera
				});
			}
			return success;
		}
		removeCamera(camera) {
			const success = this.$cameras.remove(camera);
			if (success) {
				this.dispatchEvent({
					type: 'delete-camera',
					camera
				});
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
			raycastTraverse(this, this.root, ray, intersects);
		}
		raycastFirst(ray) {
			if (!this.root) {
				return null;
			}
			return raycastTraverseFirstHit(this, this.root, ray);
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

	class PivotPointMesh extends t3d.Mesh {
		constructor() {
			super(new t3d.PlaneGeometry(0, 0), new PivotMaterial());
			this.renderOrder = Infinity;
		}
	}
	class PivotMaterial extends t3d.ShaderMaterial {
		constructor() {
			super(pivotShader);
			this.depthWrite = false;
			this.depthTest = false;
			this.transparent = true;
		}
	}
	const pivotShader = {
		name: 'PivotPoint',
		uniforms: {
			resolution: [512, 512],
			size: 15,
			thickness: 2,
			opacity: 1
		},
		vertexShader: /* glsl */`
		attribute vec3 a_Position; 
		attribute vec2 a_Uv;

		uniform mat4 u_ProjectionView;
		uniform mat4 u_Model;

		uniform float pixelRatio;
		uniform float size;
		uniform float thickness;
		uniform vec2 resolution;

		varying vec2 v_Uv;

		void main() {
			v_Uv = a_Uv;

			float aspect = resolution.x / resolution.y;
			vec2 offset = a_Uv * 2.0 - vec2(1.0);
			offset.y *= aspect;

			vec4 screenPoint = u_ProjectionView * u_Model * vec4(a_Position, 1.0);
			screenPoint.xy += offset * (size + thickness) * screenPoint.w / resolution.x;

			gl_Position = screenPoint;
		}
	`,
		fragmentShader: /* glsl */`
		uniform float size;
		uniform float thickness;
		uniform float opacity;

		varying vec2 v_Uv;

		void main() {
			float ht = 0.5 * thickness;
			float planeDim = size + thickness;
			float offset = (planeDim - ht - 2.0) / planeDim;
			float texelThickness = ht / planeDim;

			vec2 vec = v_Uv * 2.0 - vec2(1.0);
			float dist = abs(length(vec) - offset);
			float fw = fwidth(dist) * 0.5;
			float a = smoothstep(texelThickness - fw, texelThickness + fw, dist);

			gl_FragColor = vec4(1, 1, 1, opacity * (1.0 - a));
		}
	`
	};

	const _vec$4 = new t3d.Vector2();
	const _vec2 = new t3d.Vector2();
	class PointerTracker {
		constructor() {
			this.domElement = null;
			this.buttons = 0;
			this.pointerType = null;
			this.pointerOrder = [];
			this.previousPositions = {};
			this.pointerPositions = {};
			this.startPositions = {};
			this.pointerSetThisFrame = {};
			this.hoverPosition = new t3d.Vector2();
			this.hoverSet = false;
		}
		reset() {
			this.buttons = 0;
			this.pointerType = null;
			this.pointerOrder = [];
			this.previousPositions = {};
			this.pointerPositions = {};
			this.startPositions = {};
			this.pointerSetThisFrame = {};
			this.hoverPosition = new t3d.Vector2();
			this.hoverSet = false;
		}

		// The pointers can be set multiple times per frame so track whether the pointer has
		// been set this frame or not so we don't overwrite the previous position and lose information
		// about pointer movement
		updateFrame() {
			const {
				previousPositions,
				pointerPositions
			} = this;
			for (const id in pointerPositions) {
				previousPositions[id].copy(pointerPositions[id]);
			}
		}
		setHoverEvent(e) {
			if (e.pointerType === 'mouse' || e.type === 'wheel') {
				this.getAdjustedPointer(e, this.hoverPosition);
				this.hoverSet = true;
			}
		}
		getLatestPoint(target) {
			if (this.pointerType !== null) {
				this.getCenterPoint(target);
				return target;
			} else if (this.hoverSet) {
				target.copy(this.hoverPosition);
				return target;
			} else {
				return null;
			}
		}

		// get the pointer position in the coordinate system of the target element
		getAdjustedPointer(e, target) {
			const domRef = this.domElement ? this.domElement : e.target;
			const rect = domRef.getBoundingClientRect();
			const x = e.clientX - rect.left;
			const y = e.clientY - rect.top;
			target.set(x, y);
		}
		addPointer(e) {
			const id = e.pointerId;
			const position = new t3d.Vector2();
			this.getAdjustedPointer(e, position);
			this.pointerOrder.push(id);
			this.pointerPositions[id] = position;
			this.previousPositions[id] = position.clone();
			this.startPositions[id] = position.clone();
			if (this.getPointerCount() === 1) {
				this.pointerType = e.pointerType;
				this.buttons = e.buttons;
			}
		}
		updatePointer(e) {
			const id = e.pointerId;
			if (!(id in this.pointerPositions)) {
				return false;
			}
			this.getAdjustedPointer(e, this.pointerPositions[id]);
			return true;
		}
		deletePointer(e) {
			const id = e.pointerId;
			const pointerOrder = this.pointerOrder;
			pointerOrder.splice(pointerOrder.indexOf(id), 1);
			delete this.pointerPositions[id];
			delete this.previousPositions[id];
			delete this.startPositions[id];
			if (this.getPointerCount() === 0) {
				this.buttons = 0;
				this.pointerType = null;
			}
		}
		getPointerCount() {
			return this.pointerOrder.length;
		}
		getCenterPoint(target, pointerPositions = this.pointerPositions) {
			const pointerOrder = this.pointerOrder;
			if (this.getPointerCount() === 1 || this.getPointerType() === 'mouse') {
				const id = pointerOrder[0];
				target.copy(pointerPositions[id]);
				return target;
			} else if (this.getPointerCount() === 2) {
				const id0 = this.pointerOrder[0];
				const id1 = this.pointerOrder[1];
				const p0 = pointerPositions[id0];
				const p1 = pointerPositions[id1];
				target.addVectors(p0, p1).multiplyScalar(0.5);
				return target;
			}
			return null;
		}
		getPreviousCenterPoint(target) {
			return this.getCenterPoint(target, this.previousPositions);
		}
		getStartCenterPoint(target) {
			return this.getCenterPoint(target, this.startPositions);
		}
		getMoveDistance() {
			this.getCenterPoint(_vec$4);
			this.getPreviousCenterPoint(_vec2);
			return _vec$4.sub(_vec2).getLength();
		}
		getTouchPointerDistance(pointerPositions = this.pointerPositions) {
			if (this.getPointerCount() <= 1 || this.getPointerType() === 'mouse') {
				return 0;
			}
			const {
				pointerOrder
			} = this;
			const id0 = pointerOrder[0];
			const id1 = pointerOrder[1];
			const p0 = pointerPositions[id0];
			const p1 = pointerPositions[id1];
			return p0.distanceTo(p1);
		}
		getPreviousTouchPointerDistance() {
			return this.getTouchPointerDistance(this.previousPositions);
		}
		getStartTouchPointerDistance() {
			return this.getTouchPointerDistance(this.startPositions);
		}
		getPointerType() {
			return this.pointerType;
		}
		isPointerTouch() {
			return this.getPointerType() === 'touch';
		}
		getPointerButtons() {
			return this.buttons;
		}
		isLeftClicked() {
			return Boolean(this.buttons & 1);
		}
		isRightClicked() {
			return Boolean(this.buttons & 2);
		}
	}

	const _matrix = new t3d.Matrix4();
	const _ray$2 = new t3d.Ray();
	const _vec$3 = new t3d.Vector3();

	// helper function for constructing a matrix for rotating around a point
	function makeRotateAroundPoint(point, quat, target) {
		target.makeTranslation(-point.x, -point.y, -point.z);
		_matrix.makeRotationFromQuaternion(quat);
		target.premultiply(_matrix);
		_matrix.makeTranslation(point.x, point.y, point.z);
		target.premultiply(_matrix);
		return target;
	}

	// get the three.js pointer coords from an event
	function mouseToCoords(clientX, clientY, element, target) {
		target.x = (clientX - element.offsetLeft) / element.clientWidth * 2 - 1;
		target.y = -((clientY - element.offsetTop) / element.clientHeight) * 2 + 1;
		if (target.isVector3) {
			target.z = 0;
		}
	}

	// Returns an estimate of the closest point on the ellipsoid to the ray. Returns
	// the surface intersection if they collide.
	function closestRayEllipsoidSurfacePointEstimate(ray, ellipsoid, target) {
		if (ellipsoid.intersectRay(ray, target)) {
			return target;
		} else {
			_matrix.makeScale(...ellipsoid.radius).invert();
			_ray$2.copy(ray).applyMatrix4(_matrix);
			_vec$3.set(0, 0, 0);
			_ray$2.closestPointToPoint(_vec$3, target).normalize();
			_matrix.makeScale(...ellipsoid.radius);
			return target.applyMatrix4(_matrix);
		}
	}

	// find the closest ray on the horizon when the ray passes above the sphere
	function closestRaySpherePointFromRotation(ray, radius, target) {
		const hypotenuse = ray.origin.getLength();

		// angle inside the sphere
		const theta = Math.acos(radius / hypotenuse);

		// the direction to the camera
		target.copy(ray.origin).multiplyScalar(-1).normalize();

		// get the normal of the plane the ray and origin lie in
		const rotationVec = _vec$3.crossVectors(target, ray.direction).normalize();

		// rotate the camera direction by angle and scale it to the surface
		target.multiplyScalar(-1).applyAxisAngle(rotationVec, -theta).normalize().multiplyScalar(radius);
	}

	// custom version of set raycaster from camera that relies on the underlying matrices
	// so the ray origin is position at the camera near clip.
	function setRaycasterFromCamera(raycaster, coords, camera) {
		const ray = raycaster instanceof t3d.Ray ? raycaster : raycaster.ray;
		const {
			origin,
			direction
		} = ray;

		// get the origin and direction of the frustum ray
		origin.set(coords.x, coords.y, -1).unproject(camera);
		direction.set(coords.x, coords.y, 1).unproject(camera).sub(origin);

		// normalize the ray direction
		direction.normalize();
	}

	const NONE$1 = 0;
	const DRAG = 1;
	const ROTATE = 2;
	const ZOOM = 3;
	const WAITING = 4;
	const DRAG_PLANE_THRESHOLD = 0.05;
	const DRAG_UP_THRESHOLD = 0.025;
	const _rotMatrix$1 = /* @__PURE__ */new t3d.Matrix4();
	const _delta = /* @__PURE__ */new t3d.Vector3();
	const _vec$2 = /* @__PURE__ */new t3d.Vector3();
	const _forward$1 = /* @__PURE__ */new t3d.Vector3();
	const _right$1 = /* @__PURE__ */new t3d.Vector3();
	const _rotationAxis = /* @__PURE__ */new t3d.Vector3();
	const _quaternion$2 = /* @__PURE__ */new t3d.Quaternion();
	const _plane = /* @__PURE__ */new t3d.Plane();
	const _localUp = /* @__PURE__ */new t3d.Vector3();
	const _mouseBefore = /* @__PURE__ */new t3d.Vector3();
	const _mouseAfter = /* @__PURE__ */new t3d.Vector3();
	const _identityQuat = /* @__PURE__ */new t3d.Quaternion();
	const _ray$1 = /* @__PURE__ */new t3d.Ray();
	const _zoomPointPointer = /* @__PURE__ */new t3d.Vector2();
	const _pointer$1 = /* @__PURE__ */new t3d.Vector2();
	const _prevPointer = /* @__PURE__ */new t3d.Vector2();
	const _deltaPointer = /* @__PURE__ */new t3d.Vector2();
	const _centerPoint = /* @__PURE__ */new t3d.Vector2();
	const _startCenterPoint = /* @__PURE__ */new t3d.Vector2();
	const _changeEvent = {
		type: 'change'
	};
	const _startEvent = {
		type: 'start'
	};
	const _endEvent = {
		type: 'end'
	};
	class EnvironmentControls extends t3d.EventDispatcher {
		get enabled() {
			return this._enabled;
		}
		set enabled(v) {
			if (v !== this.enabled) {
				this._enabled = v;
				this.resetState();
				this.pointerTracker.reset();
				if (!this.enabled) {
					this.dragInertia.set(0, 0, 0);
					this.rotationInertia.set(0, 0);
				}
			}
		}
		constructor(scene = null, camera = null, domElement = null, tilesRenderer = null) {
			super();
			this.isEnvironmentControls = true;
			this.domElement = null;
			this.camera = null;
			this.scene = null;
			this.tilesRenderer = null;

			// settings
			this._enabled = true;
			this.cameraRadius = 5;
			this.rotationSpeed = 1;
			this.minAltitude = 0;
			this.maxAltitude = 0.45 * Math.PI;
			this.minDistance = 10;
			this.maxDistance = Infinity;
			this.minZoom = 0;
			this.maxZoom = Infinity;
			this.zoomSpeed = 1;
			this.adjustHeight = true;
			this.enableDamping = false;
			this.dampingFactor = 0.15;
			this.fallbackPlane = new t3d.Plane(new t3d.Vector3(0, 1, 0), 0);
			this.useFallbackPlane = true;

			// settings for GlobeControls
			this.reorientOnDrag = true;
			this.scaleZoomOrientationAtEdges = false;

			// internal state
			this.state = NONE$1;
			this.pointerTracker = new PointerTracker();
			this.needsUpdate = false;
			this.actionHeightOffset = 0;
			this.pivotPoint = new t3d.Vector3();

			// used for zoom
			this.zoomDirectionSet = false;
			this.zoomPointSet = false;
			this.zoomDirection = new t3d.Vector3();
			this.zoomPoint = new t3d.Vector3();
			this.zoomDelta = 0;

			// fields used for inertia
			this.rotationInertiaPivot = new t3d.Vector3();
			this.rotationInertia = new t3d.Vector2();
			this.dragInertia = new t3d.Vector3();
			this.inertiaTargetDistance = Infinity; // track the distance from the camera that we want to use to calculate the inertia end threshold
			this.inertiaStableFrames = 0; // the number of frames that the camera has not moved while the user is interacting

			// circular pivot mesh
			this.pivotMesh = new PivotPointMesh();
			this.pivotMesh.raycast = () => {};
			this.pivotMesh.scale.setScalar(0.25);

			// raycaster
			this.raycaster = new t3d.Raycaster();
			this.up = new t3d.Vector3(0, 1, 0);
			this._detachCallback = null;
			this._upInitialized = false;
			this._lastUsedState = NONE$1;
			this._zoomPointWasSet = false;

			// always update the zoom target point in case the tiles are changing
			this._tilesOnChangeCallback = () => this.zoomPointSet = false;

			// init
			if (domElement) this.attach(domElement);
			if (camera) this.setCamera(camera);
			if (scene) this.setScene(scene);
			if (tilesRenderer) this.setTilesRenderer(tilesRenderer);
		}
		setScene(scene) {
			this.scene = scene;
		}
		setCamera(camera) {
			this.camera = camera;
			this._upInitialized = false;
			this.zoomDirectionSet = false;
			this.zoomPointSet = false;
			this.needsUpdate = true;
			this.raycaster.camera = camera;
			this.resetState();
		}
		setTilesRenderer(tilesRenderer) {
			// TODO: what if a scene has multiple tile sets?
			if (this.tilesRenderer) {
				this.tilesRenderer.removeEventListener('tile-visibility-change', this._tilesOnChangeCallback);
			}
			this.tilesRenderer = tilesRenderer;
			if (this.tilesRenderer !== null) {
				this.tilesRenderer.addEventListener('tile-visibility-change', this._tilesOnChangeCallback);
				if (this.scene === null) {
					this.setScene(this.tilesRenderer);
				}
			}
		}
		attach(domElement) {
			if (this.domElement) {
				throw new Error('EnvironmentControls: Controls already attached to element');
			}

			// set the touch action to none so the browser does not
			// drag the page to refresh or scroll
			this.domElement = domElement;
			this.pointerTracker.domElement = domElement;
			domElement.style.touchAction = 'none';
			const contextMenuCallback = e => {
				// exit early if the controls are disabled
				if (!this.enabled) {
					return;
				}
				e.preventDefault();
			};
			const pointerdownCallback = e => {
				// exit early if the controls are disabled
				if (!this.enabled) {
					return;
				}
				e.preventDefault();
				const {
					camera,
					raycaster,
					domElement,
					up,
					pivotMesh,
					pointerTracker,
					scene,
					pivotPoint,
					enabled
				} = this;

				// init the pointer
				pointerTracker.addPointer(e);
				this.needsUpdate = true;

				// handle cases where we need to capture the pointer or
				// reset state when we have too many pointers
				if (pointerTracker.isPointerTouch()) {
					pivotMesh.visible = false;
					if (pointerTracker.getPointerCount() === 0) {
						domElement.setPointerCapture(e.pointerId);
					} else if (pointerTracker.getPointerCount() > 2) {
						this.resetState();
						return;
					}
				}

				// the "pointer" for zooming and rotating should be based on the center point
				pointerTracker.getCenterPoint(_pointer$1);
				mouseToCoords(_pointer$1.x, _pointer$1.y, domElement, _pointer$1);
				setRaycasterFromCamera(raycaster, _pointer$1, camera);

				// prevent the drag distance from getting too severe by limiting the drag point
				// to a reasonable angle and reasonable distance with the drag plane
				const dot = Math.abs(raycaster.ray.direction.dot(up));
				if (dot < DRAG_PLANE_THRESHOLD || dot < DRAG_UP_THRESHOLD) {
					return;
				}

				// find the hit point
				const hit = this._raycast(raycaster);
				if (hit) {
					// if two fingers, right click, or shift click are being used then we trigger
					// a rotation action to begin
					if (pointerTracker.getPointerCount() === 2 || pointerTracker.isRightClicked() || pointerTracker.isLeftClicked() && e.shiftKey) {
						this.setState(pointerTracker.isPointerTouch() ? WAITING : ROTATE);
						pivotPoint.copy(hit.point);
						pivotMesh.position.copy(hit.point);
						pivotMesh.visible = pointerTracker.isPointerTouch() ? false : enabled;
						pivotMesh.updateMatrix();
						scene.add(pivotMesh);
					} else if (pointerTracker.isLeftClicked()) {
						// if the clicked point is coming from below the plane then don't perform the drag
						this.setState(DRAG);
						pivotPoint.copy(hit.point);
						pivotMesh.position.copy(hit.point);
						pivotMesh.updateMatrix();
						scene.add(pivotMesh);
					}
				}
			};
			let _pointerMoveQueued = false;
			const pointermoveCallback = e => {
				// exit early if the controls are disabled
				if (!this.enabled) {
					return;
				}
				e.preventDefault();
				const {
					pivotMesh,
					enabled
				} = this;

				// whenever the pointer moves we need to re-derive the zoom direction and point
				this.zoomDirectionSet = false;
				this.zoomPointSet = false;
				if (this.state !== NONE$1) {
					this.needsUpdate = true;
				}
				const {
					pointerTracker
				} = this;
				pointerTracker.setHoverEvent(e);
				if (!pointerTracker.updatePointer(e)) {
					return;
				}
				if (pointerTracker.isPointerTouch() && pointerTracker.getPointerCount() === 2) {
					// We queue this event to ensure that all pointers have been updated
					if (!_pointerMoveQueued) {
						_pointerMoveQueued = true;
						queueMicrotask(() => {
							_pointerMoveQueued = false;

							// adjust the pointer position to be the center point
							pointerTracker.getCenterPoint(_centerPoint);

							// detect zoom transition
							const startDist = pointerTracker.getStartTouchPointerDistance();
							const pointerDist = pointerTracker.getTouchPointerDistance();
							const separateDelta = pointerDist - startDist;
							if (this.state === NONE$1 || this.state === WAITING) {
								// check which direction was moved in first - if the pointers are pinching then
								// it's a zoom. But if they move in parallel it's a rotation
								pointerTracker.getCenterPoint(_centerPoint);
								pointerTracker.getStartCenterPoint(_startCenterPoint);

								// adjust the drag requirement by the dpr
								const dragThreshold = 2.0 * window.devicePixelRatio;
								const parallelDelta = _centerPoint.distanceTo(_startCenterPoint);
								if (Math.abs(separateDelta) > dragThreshold || parallelDelta > dragThreshold) {
									if (Math.abs(separateDelta) > parallelDelta) {
										this.setState(ZOOM);
										this.zoomDirectionSet = false;
									} else {
										this.setState(ROTATE);
									}
								}
							}
							if (this.state === ZOOM) {
								const previousDist = pointerTracker.getPreviousTouchPointerDistance();
								this.zoomDelta += pointerDist - previousDist;
								pivotMesh.visible = false;
							} else if (this.state === ROTATE) {
								pivotMesh.visible = enabled;
							}
						});
					}
				}

				// TODO: we have the potential to fire change multiple times per frame - should we debounce?
				this.dispatchEvent(_changeEvent);
			};
			const pointerupCallback = e => {
				// exit early if the controls are disabled
				if (!this.enabled) {
					return;
				}
				const {
					pointerTracker
				} = this;
				pointerTracker.deletePointer(e);
				if (pointerTracker.getPointerType() === 'touch' && pointerTracker.getPointerCount() === 0) {
					domElement.releasePointerCapture(e.pointerId);
				}
				this.resetState();
				this.needsUpdate = true;
			};
			const wheelCallback = e => {
				// exit early if the controls are disabled
				if (!this.enabled) {
					return;
				}
				e.preventDefault();
				const {
					pointerTracker
				} = this;
				pointerTracker.setHoverEvent(e);
				pointerTracker.updatePointer(e);

				// TODO: do we need events here?
				this.dispatchEvent(_startEvent);
				let delta;
				switch (e.deltaMode) {
					case 2:
						// Pages
						delta = e.deltaY * 800;
						break;
					case 1:
						// Lines
						delta = e.deltaY * 40;
						break;
					case 0:
						// Pixels
						delta = e.deltaY;
						break;
				}

				// use LOG to scale the scroll delta and hopefully normalize them across platforms
				const deltaSign = Math.sign(delta);
				const normalizedDelta = Math.abs(delta);
				this.zoomDelta -= 0.25 * deltaSign * normalizedDelta;
				this.needsUpdate = true;
				this._lastUsedState = ZOOM;
				this.dispatchEvent(_endEvent);
			};
			const pointerenterCallback = e => {
				// exit early if the controls are disabled
				if (!this.enabled) {
					return;
				}
				const {
					pointerTracker
				} = this;
				if (e.buttons !== pointerTracker.getPointerButtons()) {
					pointerTracker.deletePointer(e);
					this.resetState();
				}
			};
			domElement.addEventListener('contextmenu', contextMenuCallback);
			domElement.addEventListener('pointerdown', pointerdownCallback);
			domElement.addEventListener('pointermove', pointermoveCallback);
			domElement.addEventListener('pointerup', pointerupCallback);
			domElement.addEventListener('wheel', wheelCallback, {
				passive: false
			});
			domElement.addEventListener('pointerenter', pointerenterCallback);
			this._detachCallback = () => {
				domElement.removeEventListener('contextmenu', contextMenuCallback);
				domElement.removeEventListener('pointerdown', pointerdownCallback);
				domElement.removeEventListener('pointermove', pointermoveCallback);
				domElement.removeEventListener('pointerup', pointerupCallback);
				domElement.removeEventListener('wheel', wheelCallback);
				domElement.removeEventListener('pointerenter', pointerenterCallback);
			};
		}

		// override-able functions for retrieving the up direction at a point
		getUpDirection(point, target) {
			target.copy(this.up);
		}
		getCameraUpDirection(target) {
			this.getUpDirection(this.camera.position, target);
		}

		// returns the active / last used pivot point for the scene
		getPivotPoint(target) {
			let result = null;

			// get the last interacted point as the focus
			if (this._lastUsedState === ZOOM) {
				if (this._zoomPointWasSet) {
					result = target.copy(this.zoomPoint);
				}
			} else if (this._lastUsedState === ROTATE || this._lastUsedState === DRAG) {
				result = target.copy(this.pivotPoint);
			}

			// If the last used point is outside the camera view then skip it
			const {
				camera,
				raycaster
			} = this;
			if (result !== null) {
				_vec$2.copy(result).project(camera);
				if (_vec$2.x < -1 || _vec$2.x > 1 || _vec$2.y < -1 || _vec$2.y > 1) {
					result = null;
				}
			}

			// default to the raycast hit if we have not result or the hit is closer to the camera
			// set a ray in the local ellipsoid frame
			setRaycasterFromCamera(raycaster, {
				x: 0,
				y: 0
			}, camera);
			const hit = this._raycast(raycaster);
			if (hit) {
				if (result === null || hit.distance < result.distanceTo(raycaster.ray.origin)) {
					result = target.copy(hit.point);
				}
			}
			return result;
		}
		detach() {
			this.domElement = null;
			if (this._detachCallback) {
				this._detachCallback();
				this._detachCallback = null;
				this.pointerTracker.reset();
			}
		}
		resetState() {
			if (this.state !== NONE$1) {
				this.dispatchEvent(_endEvent);
			}
			this.state = NONE$1;
			this.pivotMesh.removeFromParent();
			this.pivotMesh.visible = this.enabled;
			this.actionHeightOffset = 0;
		}
		setState(state = this.state, fireEvent = true) {
			if (this.state === state) {
				return;
			}
			if (this.state === NONE$1 && fireEvent) {
				this.dispatchEvent(_startEvent);
			}
			this.pivotMesh.visible = this.enabled;
			this.dragInertia.set(0, 0, 0);
			this.rotationInertia.set(0, 0);
			this.inertiaStableFrames = 0;
			this.state = state;
			if (state !== NONE$1 && state !== WAITING) {
				this._lastUsedState = state;
			}
		}
		update(deltaTime = 64 / 1000) {
			if (!this.enabled || !this.camera || deltaTime === 0) {
				return;
			}
			const {
				camera,
				cameraRadius,
				pivotPoint,
				up,
				state,
				adjustHeight
			} = this;
			camera.updateMatrix();

			// set the "up" vector immediately so it's available in the following functions
			this.getCameraUpDirection(_localUp);
			if (!this._upInitialized) {
				this._upInitialized = true;
				this.up.copy(_localUp);
			}

			// update the actions
			const inertiaNeedsUpdate = this._inertiaNeedsUpdate();
			if (this.needsUpdate || inertiaNeedsUpdate) {
				const zoomDelta = this.zoomDelta;
				this._updateZoom();
				this._updatePosition(deltaTime);
				this._updateRotation(deltaTime);
				if (state === DRAG || state === ROTATE) {
					_forward$1.set(0, 0, -1).transformDirection(camera.worldMatrix);
					this.inertiaTargetDistance = _vec$2.copy(this.pivotPoint).sub(camera.position).dot(_forward$1);
				} else if (state === NONE$1) {
					this._updateInertia(deltaTime);
				}
				if (state !== NONE$1 || zoomDelta !== 0 || inertiaNeedsUpdate) {
					this.dispatchEvent(_changeEvent);
				}
				this.needsUpdate = false;
			}

			// update the up direction based on where the camera moved to
			// if using an orthographic camera then rotate around drag pivot
			// reuse the "hit" information since it can be slow to perform multiple hits
			// TODO Orthographic Camera
			const hit = camera.isOrthographicCamera ? null : adjustHeight && this._getPointBelowCamera() || null;
			const rotationPoint = camera.isOrthographicCamera ? pivotPoint : hit && hit.point || null;
			this.getCameraUpDirection(_localUp);
			this._setFrame(_localUp, rotationPoint);

			// when dragging the camera and drag point may be moved
			// to accommodate terrain so we try to move it back down
			// to the original point.
			if ((this.state === DRAG || this.state === ROTATE) && this.actionHeightOffset !== 0) {
				const {
					actionHeightOffset
				} = this;
				camera.position.addScaledVector(up, -actionHeightOffset);
				pivotPoint.addScaledVector(up, -actionHeightOffset);

				// adjust the height
				if (hit) {
					hit.distance -= actionHeightOffset;
				}
			}
			this.actionHeightOffset = 0;
			if (hit) {
				const dist = hit.distance;
				if (dist < cameraRadius) {
					const delta = cameraRadius - dist;
					camera.position.addScaledVector(up, delta);
					pivotPoint.addScaledVector(up, delta);
					this.actionHeightOffset = delta;
				}
			}
			this.pointerTracker.updateFrame();
		}

		// updates the camera to position it based on the constraints of the controls
		adjustCamera(camera) {
			const {
				adjustHeight,
				cameraRadius
			} = this;
			if (camera.isPerspectiveCamera) {
				// adjust the camera height
				this.getUpDirection(camera.position, _localUp);
				const hit = adjustHeight && this._getPointBelowCamera(camera.position, _localUp) || null;
				if (hit) {
					const dist = hit.distance;
					if (dist < cameraRadius) {
						camera.position.addScaledVector(_localUp, cameraRadius - dist);
					}
				}
			}
		}
		dispose() {
			this.detach();
		}

		// private
		_updateInertia(deltaTime) {
			// update the damping of momentum variables
			const {
				rotationInertia,
				pivotPoint,
				dragInertia,
				enableDamping,
				dampingFactor,
				camera,
				cameraRadius,
				minDistance,
				inertiaTargetDistance
			} = this;
			if (!this.enableDamping || this.inertiaStableFrames > 1) {
				dragInertia.set(0, 0, 0);
				rotationInertia.set(0, 0, 0);
				return;
			}

			// Based on Freya Holmer's frame-rate independent lerp function
			const factor = Math.pow(2, -deltaTime / dampingFactor);
			const stableDistance = Math.max(camera.near, cameraRadius, minDistance, inertiaTargetDistance);
			const resolution = 2 * 1e3;
			const pixelWidth = 2 / resolution;
			const pixelThreshold = 0.25 * pixelWidth;

			// scale the residual rotation motion
			if (rotationInertia.getLengthSquared() > 0) {
				// calculate two screen points at 1 pixel apart in our notional resolution so we can stop when the delta is ~ 1 pixel
				// projected into world space
				setRaycasterFromCamera(_ray$1, _vec$2.set(0, 0, -1), camera);
				_ray$1.applyMatrix4(camera.viewMatrix);
				_ray$1.direction.normalize();
				_ray$1.recast(-_ray$1.direction.dot(_ray$1.origin)).at(stableDistance / _ray$1.direction.z, _vec$2);
				_vec$2.applyMatrix4(camera.worldMatrix);
				setRaycasterFromCamera(_ray$1, _delta.set(pixelThreshold, pixelThreshold, -1), camera);
				_ray$1.applyMatrix4(camera.viewMatrix);
				_ray$1.direction.normalize();
				_ray$1.recast(-_ray$1.direction.dot(_ray$1.origin)).at(stableDistance / _ray$1.direction.z, _delta);
				_delta.applyMatrix4(camera.worldMatrix);

				// get implied angle
				_vec$2.sub(pivotPoint).normalize();
				_delta.sub(pivotPoint).normalize();

				// calculate the rotation threshold
				const threshold = _vec$2.angleTo(_delta) / deltaTime;
				rotationInertia.multiplyScalar(factor);
				if (rotationInertia.getLengthSquared() < threshold ** 2 || !enableDamping) {
					rotationInertia.set(0, 0);
				}
			}

			// scale the residual translation motion
			if (dragInertia.getLengthSquared() > 0) {
				// calculate two screen points at 1 pixel apart in our notional resolution so we can stop when the delta is ~ 1 pixel
				// projected into world space
				setRaycasterFromCamera(_ray$1, _vec$2.set(0, 0, -1), camera);
				_ray$1.applyMatrix4(camera.viewMatrix);
				_ray$1.direction.normalize();
				_ray$1.recast(-_ray$1.direction.dot(_ray$1.origin)).at(stableDistance / _ray$1.direction.z, _vec$2);
				_vec$2.applyMatrix4(camera.worldMatrix);
				setRaycasterFromCamera(_ray$1, _delta.set(pixelThreshold, pixelThreshold, -1), camera);
				_ray$1.applyMatrix4(camera.viewMatrix);
				_ray$1.direction.normalize();
				_ray$1.recast(-_ray$1.direction.dot(_ray$1.origin)).at(stableDistance / _ray$1.direction.z, _delta);
				_delta.applyMatrix4(camera.worldMatrix);

				// calculate movement threshold
				const threshold = _vec$2.distanceTo(_delta) / deltaTime;
				dragInertia.multiplyScalar(factor);
				if (dragInertia.getLengthSquared() < threshold ** 2 || !enableDamping) {
					dragInertia.set(0, 0, 0);
				}
			}

			// apply the inertia changes
			if (rotationInertia.getLengthSquared() > 0) {
				this._applyRotation(rotationInertia.x * deltaTime, rotationInertia.y * deltaTime, pivotPoint);
			}
			if (dragInertia.getLengthSquared() > 0) {
				camera.position.addScaledVector(dragInertia, deltaTime);
				camera.updateMatrix();
			}
		}
		_inertiaNeedsUpdate() {
			const {
				rotationInertia,
				dragInertia
			} = this;
			return rotationInertia.getLengthSquared() !== 0 || dragInertia.getLengthSquared() !== 0;
		}
		_updateZoom() {
			const {
				zoomPoint,
				zoomDirection,
				camera,
				minDistance,
				maxDistance,
				pointerTracker,
				domElement,
				minZoom,
				maxZoom,
				zoomSpeed,
				state
			} = this;
			let scale = this.zoomDelta;
			this.zoomDelta = 0;

			// get the latest hover / touch point
			if (!pointerTracker.getLatestPoint(_pointer$1) || scale === 0 && state !== ZOOM) {
				return;
			}

			// reset momentum
			this.rotationInertia.set(0, 0);
			this.dragInertia.set(0, 0, 0);
			if (camera.isOrthographicCamera) {
				// update the zoom direction
				this._updateZoomDirection();

				// zoom straight into the globe if we haven't hit anything
				const zoomIntoPoint = this.zoomPointSet || this._updateZoomPoint();

				// get the mouse position before zoom
				_mouseBefore.unproject(camera);

				// zoom the camera
				const normalizedDelta = Math.pow(0.95, Math.abs(scale * 0.05));
				let scaleFactor = scale > 0 ? 1 / Math.abs(normalizedDelta) : normalizedDelta;
				scaleFactor *= zoomSpeed;
				if (scaleFactor > 1) {
					if (maxZoom < camera.zoom * scaleFactor) {
						scaleFactor = 1;
					}
				} else {
					if (minZoom > camera.zoom * scaleFactor) {
						scaleFactor = 1;
					}
				}

				// TODO
				// camera.zoom *= scaleFactor;
				// camera.updateProjectionMatrix();

				// adjust the surface point to be in the same position if the globe is hovered over
				if (zoomIntoPoint) {
					// get the mouse position after zoom
					mouseToCoords(_pointer$1.x, _pointer$1.y, domElement, _mouseAfter);
					_mouseAfter.unproject(camera);

					// shift the camera on the near plane so the mouse is in the same spot
					camera.position.sub(_mouseAfter).add(_mouseBefore);
					camera.updateMatrix();
				}
			} else {
				// initialize the zoom direction
				this._updateZoomDirection();

				// track the zoom direction we're going to use
				const finalZoomDirection = _vec$2.copy(zoomDirection);
				if (this.zoomPointSet || this._updateZoomPoint()) {
					const dist = zoomPoint.distanceTo(camera.position);

					// scale the distance based on how far there is to move
					if (scale < 0) {
						const remainingDistance = Math.min(0, dist - maxDistance);
						scale = scale * dist * zoomSpeed * 0.0025;
						scale = Math.max(scale, remainingDistance);
					} else {
						const remainingDistance = Math.max(0, dist - minDistance);
						scale = scale * Math.max(dist - minDistance, 0) * zoomSpeed * 0.0025;
						scale = Math.min(scale, remainingDistance);
					}
					camera.position.addScaledVector(zoomDirection, scale);
					camera.updateMatrix();
				} else {
					// if we're zooming into nothing then use the distance from the ground to scale movement
					const hit = this._getPointBelowCamera();
					if (hit) {
						const dist = hit.distance;
						finalZoomDirection.set(0, 0, -1).transformDirection(camera.worldMatrix);
						camera.position.addScaledVector(finalZoomDirection, scale * dist * 0.01);
						camera.updateMatrix();
					}
				}
			}
		}
		_updateZoomDirection() {
			if (this.zoomDirectionSet) {
				return;
			}
			const {
				domElement,
				raycaster,
				camera,
				zoomDirection,
				pointerTracker
			} = this;
			pointerTracker.getLatestPoint(_pointer$1);
			mouseToCoords(_pointer$1.x, _pointer$1.y, domElement, _mouseBefore);
			setRaycasterFromCamera(raycaster, _mouseBefore, camera);
			zoomDirection.copy(raycaster.ray.direction).normalize();
			this.zoomDirectionSet = true;
		}

		// update the point being zoomed in to based on the zoom direction
		_updateZoomPoint() {
			const {
				camera,
				zoomDirectionSet,
				zoomDirection,
				raycaster,
				zoomPoint,
				pointerTracker,
				domElement
			} = this;
			this._zoomPointWasSet = false;
			if (!zoomDirectionSet) {
				return false;
			}

			// If using an orthographic camera we have to account for the mouse position when picking the point
			if (camera.isOrthographicCamera && pointerTracker.getLatestPoint(_zoomPointPointer)) {
				mouseToCoords(_zoomPointPointer.x, _zoomPointPointer.y, domElement, _zoomPointPointer);
				setRaycasterFromCamera(raycaster, _zoomPointPointer, camera);
			} else {
				raycaster.ray.origin.copy(camera.position);
				raycaster.ray.direction.copy(zoomDirection);
				raycaster.near = 0;
				raycaster.far = Infinity;
			}

			// get the hit point
			const hit = this._raycast(raycaster);
			if (hit) {
				zoomPoint.copy(hit.point);
				this.zoomPointSet = true;
				this._zoomPointWasSet = true;
				return true;
			}
			return false;
		}

		// returns the point below the camera
		_getPointBelowCamera(point = this.camera.position, up = this.up) {
			const {
				raycaster
			} = this;
			raycaster.ray.direction.copy(up).multiplyScalar(-1);
			raycaster.ray.origin.copy(point).addScaledVector(up, 1e5);
			raycaster.near = 0;
			raycaster.far = Infinity;
			const hit = this._raycast(raycaster);
			if (hit) {
				hit.distance -= 1e5;
			}
			return hit;
		}

		// update the drag action
		_updatePosition(deltaTime) {
			const {
				raycaster,
				camera,
				pivotPoint,
				up,
				pointerTracker,
				domElement,
				state,
				dragInertia
			} = this;
			if (state === DRAG) {
				// get the pointer and plane
				pointerTracker.getCenterPoint(_pointer$1);
				mouseToCoords(_pointer$1.x, _pointer$1.y, domElement, _pointer$1);
				_plane.setFromNormalAndCoplanarPoint(up, pivotPoint);
				setRaycasterFromCamera(raycaster, _pointer$1, camera);

				// prevent the drag distance from getting too severe by limiting the drag point
				// to a reasonable angle with the drag plane
				if (Math.abs(raycaster.ray.direction.dot(up)) < DRAG_PLANE_THRESHOLD) {
					// rotate the pointer direction down to the correct angle for horizontal dragging
					const angle = Math.acos(DRAG_PLANE_THRESHOLD);
					_rotationAxis.crossVectors(raycaster.ray.direction, up).normalize();
					raycaster.ray.direction.copy(up).applyAxisAngle(_rotationAxis, angle).multiplyScalar(-1);
				}

				// TODO: dragging causes the camera to rise because we're getting "pushed" up by lower resolution tiles and
				// don't lower back down. We should maintain a target height above tiles where possible
				// prevent the drag from inverting

				// if we drag to a point that's near the edge of the earth then we want to prevent it
				// from wrapping around and causing unexpected rotations
				this.getUpDirection(pivotPoint, _localUp);
				if (Math.abs(raycaster.ray.direction.dot(_localUp)) < DRAG_UP_THRESHOLD) {
					const angle = Math.acos(DRAG_UP_THRESHOLD);
					_rotationAxis.crossVectors(raycaster.ray.direction, _localUp).normalize();
					raycaster.ray.direction.copy(_localUp).applyAxisAngle(_rotationAxis, angle).multiplyScalar(-1);
				}

				// find the point on the plane that we should drag to
				if (raycaster.ray.intersectPlane(_plane, _vec$2)) {
					_delta.subVectors(pivotPoint, _vec$2);
					camera.position.add(_delta);
					camera.updateMatrix();

					// update the drag inertia
					_delta.multiplyScalar(1 / deltaTime);
					if (pointerTracker.getMoveDistance() / deltaTime < 2 * window.devicePixelRatio) {
						this.inertiaStableFrames++;
					} else {
						dragInertia.copy(_delta);
						this.inertiaStableFrames = 0;
					}
				}
			}
		}
		_updateRotation(deltaTime) {
			const {
				pivotPoint,
				pointerTracker,
				domElement,
				state,
				rotationInertia
			} = this;
			if (state === ROTATE) {
				// get the rotation motion and divide out the container height to normalize for element size
				pointerTracker.getCenterPoint(_pointer$1);
				pointerTracker.getPreviousCenterPoint(_prevPointer);
				_deltaPointer.subVectors(_pointer$1, _prevPointer).multiplyScalar(2 * Math.PI / domElement.clientHeight);
				this._applyRotation(_deltaPointer.x, _deltaPointer.y, pivotPoint);

				// update rotation inertia
				_deltaPointer.multiplyScalar(1 / deltaTime);
				if (pointerTracker.getMoveDistance() / deltaTime < 2 * window.devicePixelRatio) {
					this.inertiaStableFrames++;
				} else {
					rotationInertia.copy(_deltaPointer);
					this.inertiaStableFrames = 0;
				}
			}
		}
		_applyRotation(x, y, pivotPoint) {
			if (x === 0 && y === 0) {
				return;
			}
			const {
				camera,
				minAltitude,
				maxAltitude,
				rotationSpeed
			} = this;
			const azimuth = -x * rotationSpeed;
			let altitude = y * rotationSpeed;

			// calculate current angles and clamp
			_forward$1.set(0, 0, 1).transformDirection(camera.worldMatrix);
			this.getUpDirection(pivotPoint, _localUp);

			// get the signed angle relative to the top down view
			_vec$2.crossVectors(_localUp, _forward$1).normalize();
			_right$1.set(1, 0, 0).transformDirection(camera.worldMatrix).normalize();
			const sign = Math.sign(_vec$2.dot(_right$1));
			const angle = sign * _localUp.angleTo(_forward$1);

			// clamp the rotation to be within the provided limits
			// clamp to 0 here, as well, so we don't "pop" to the the value range
			if (altitude > 0) {
				altitude = Math.min(angle - minAltitude - 1e-2, altitude);
				altitude = Math.max(0, altitude);
			} else {
				altitude = Math.max(angle - maxAltitude, altitude);
				altitude = Math.min(0, altitude);
			}

			// rotate around the up axis
			_quaternion$2.setFromAxisAngle(_localUp, azimuth);
			makeRotateAroundPoint(pivotPoint, _quaternion$2, _rotMatrix$1);
			camera.worldMatrix.premultiply(_rotMatrix$1);

			// get a rotation axis for altitude and rotate
			_rotationAxis.set(-1, 0, 0).transformDirection(camera.worldMatrix);
			_quaternion$2.setFromAxisAngle(_rotationAxis, altitude);
			makeRotateAroundPoint(pivotPoint, _quaternion$2, _rotMatrix$1);
			camera.worldMatrix.premultiply(_rotMatrix$1);

			// update the transform members
			camera.worldMatrix.decompose(camera.position, camera.quaternion, _vec$2);
		}

		// sets the "up" axis for the current surface of the tile set
		_setFrame(newUp, pivot) {
			const {
				up,
				camera,
				state,
				zoomPoint,
				zoomDirectionSet,
				zoomPointSet,
				reorientOnDrag,
				scaleZoomOrientationAtEdges
			} = this;
			camera.updateMatrix();

			// get the amount needed to rotate
			_quaternion$2.setFromUnitVectors(up, newUp);

			// If we're zooming then reorient around the zoom point
			const action = state;
			if (zoomDirectionSet && (zoomPointSet || this._updateZoomPoint())) {
				this.getUpDirection(zoomPoint, _vec$2);
				if (scaleZoomOrientationAtEdges) {
					let amt = Math.max(_vec$2.dot(up) - 0.6, 0) / 0.4;
					amt = t3d.MathUtils.mapLinear(amt, 0, 0.5, 0, 1);
					amt = Math.min(amt, 1);

					// scale the value if we're using an orthographic camera so
					// GlobeControls works correctly
					if (camera.isOrthographicCamera) {
						amt *= 0.1;
					}
					_quaternion$2.slerpQuaternions(_quaternion$2, _identityQuat, 1.0 - amt);
				}

				// rotates the camera position around the point being zoomed in to
				makeRotateAroundPoint(zoomPoint, _quaternion$2, _rotMatrix$1);
				camera.worldMatrix.premultiply(_rotMatrix$1);
				camera.worldMatrix.decompose(camera.position, camera.quaternion, _vec$2);

				// recompute the zoom direction after updating rotation to align with frame
				this.zoomDirectionSet = false;
				this._updateZoomDirection();
			} else if (action === DRAG && reorientOnDrag) {
				// If we're dragging then reorient around the drag point

				// NOTE: We used to derive the pivot point here by getting the point below the camera
				// but decided to pass it in via "update" to avoid multiple ray casts

				if (pivot) {
					// perform a simple realignment by rotating the camera around the pivot
					makeRotateAroundPoint(pivot, _quaternion$2, _rotMatrix$1);
					camera.worldMatrix.premultiply(_rotMatrix$1);
					camera.worldMatrix.decompose(camera.position, camera.quaternion, _vec$2);
				}
			}
			up.copy(newUp);
			camera.updateMatrix();
		}
		_raycast(raycaster) {
			const {
				scene,
				useFallbackPlane,
				fallbackPlane
			} = this;
			const result = raycaster.intersectObject(scene, true)[0] || null;
			if (result) {
				return result;
			} else if (useFallbackPlane) {
				// if we don't hit any geometry then try to intersect the fallback
				// plane so the camera can still be manipulated
				const plane = fallbackPlane;
				if (raycaster.ray.intersectPlane(plane, _vec$2)) {
					const planeHit = {
						point: _vec$2.clone(),
						distance: raycaster.ray.origin.distanceTo(_vec$2)
					};
					return planeHit;
				}
			}
			return null;
		}
	}

	const _invMatrix = /* @__PURE__ */new t3d.Matrix4();
	const _rotMatrix = /* @__PURE__ */new t3d.Matrix4();
	const _pos$3 = /* @__PURE__ */new t3d.Vector3();
	const _vec$1 = /* @__PURE__ */new t3d.Vector3();
	const _center = /* @__PURE__ */new t3d.Vector3();
	const _forward = /* @__PURE__ */new t3d.Vector3();
	const _right = /* @__PURE__ */new t3d.Vector3();
	const _targetRight = /* @__PURE__ */new t3d.Vector3();
	const _globalUp = /* @__PURE__ */new t3d.Vector3();
	const _quaternion$1 = /* @__PURE__ */new t3d.Quaternion();
	const _zoomPointUp = /* @__PURE__ */new t3d.Vector3();
	const _toCenter = /* @__PURE__ */new t3d.Vector3();
	const _ray = /* @__PURE__ */new t3d.Ray();
	const _ellipsoid = /* @__PURE__ */new Ellipsoid();
	const _pointer = /* @__PURE__ */new t3d.Vector2();
	const _latLon = {};

	// hand picked minimum elevation to tune far plane near surface
	const MIN_ELEVATION = 2550;
	class GlobeControls extends EnvironmentControls {
		get ellipsoid() {
			return this.tilesRenderer ? this.tilesRenderer.ellipsoid : null;
		}
		get tilesGroup() {
			return this.tilesRenderer ? this.tilesRenderer : null;
		}
		constructor(scene = null, camera = null, domElement = null, tilesRenderer = null) {
			// store which mode the drag stats are in
			super(scene, camera, domElement);
			this.isGlobeControls = true;
			this._dragMode = 0;
			this._rotationMode = 0;
			this.maxZoom = 0.01;
			this.nearMargin = 0.25;
			this.farMargin = 0;
			this.useFallbackPlane = false;
			this.reorientOnDrag = false;
			this.globeInertia = new t3d.Quaternion();
			this.globeInertiaFactor = 0;
			this.setTilesRenderer(tilesRenderer);
		}
		setScene(scene) {
			if (scene === null && this.tilesRenderer !== null) {
				super.setScene(this.tilesRenderer);
			} else {
				super.setScene(scene);
			}
		}
		getPivotPoint(target) {
			const {
				camera,
				tilesGroup,
				ellipsoid
			} = this;

			// get camera values
			_forward.set(0, 0, -1).transformDirection(camera.worldMatrix);

			// set a ray in the local ellipsoid frame
			_ray.origin.copy(camera.position);
			_ray.direction.copy(_forward);
			_invMatrix.copy(tilesGroup.worldMatrix).invert();
			_ray.applyMatrix4(_invMatrix);

			// get the estimated closest point
			closestRayEllipsoidSurfacePointEstimate(_ray, ellipsoid, _vec$1);
			_vec$1.applyMatrix4(tilesGroup.worldMatrix);

			// use the closest point if no pivot was provided or it's closer
			if (super.getPivotPoint(target) === null || target.distanceTo(_ray.origin) > _vec$1.distanceTo(_ray.origin)) {
				target.copy(_vec$1);
			}
			return target;
		}

		// get the vector to the center of the provided globe
		getVectorToCenter(target) {
			const {
				tilesGroup,
				camera
			} = this;
			return target.setFromMatrixPosition(tilesGroup.worldMatrix).sub(camera.position);
		}

		// get the distance to the center of the globe
		getDistanceToCenter() {
			return this.getVectorToCenter(_vec$1).getLength();
		}
		getUpDirection(point, target) {
			// get the "up" direction based on the wgs84 ellipsoid
			const {
				tilesGroup,
				ellipsoid
			} = this;
			_invMatrix.copy(tilesGroup.worldMatrix).invert();
			_vec$1.copy(point).applyMatrix4(_invMatrix);
			ellipsoid.getPositionToNormal(_vec$1, target);
			target.transformDirection(tilesGroup.worldMatrix);
		}
		getCameraUpDirection(target) {
			const {
				tilesGroup,
				ellipsoid,
				camera
			} = this;
			if (camera.isOrthographicCamera) {
				this._getVirtualOrthoCameraPosition(_vec$1);
				_invMatrix.copy(tilesGroup.worldMatrix).invert();
				_vec$1.applyMatrix4(_invMatrix);
				ellipsoid.getPositionToNormal(_vec$1, target);
				target.transformDirection(tilesGroup.worldMatrix);
			} else {
				this.getUpDirection(camera.position, target);
			}
		}
		update(deltaTime = 64 / 1000) {
			if (!this.enabled || !this.tilesGroup || !this.camera || deltaTime === 0) {
				return;
			}
			const {
				camera,
				pivotMesh
			} = this;

			// if we're outside the transition threshold then we toggle some reorientation behavior
			// when adjusting the up frame while moving the camera
			if (this._isNearControls()) {
				this.scaleZoomOrientationAtEdges = this.zoomDelta < 0;
			} else {
				if (this.state !== NONE$1 && this._dragMode !== 1 && this._rotationMode !== 1) {
					pivotMesh.visible = false;
				}
				this.scaleZoomOrientationAtEdges = false;
			}

			// fire basic controls update
			super.update(deltaTime);

			// update the camera planes and the ortho camera position
			this.adjustCamera(camera);
		}

		// Updates the passed camera near and far clip planes to encapsulate the ellipsoid from the
		// current position in addition to adjusting the height.
		adjustCamera(camera) {
			super.adjustCamera(camera);
			const {
				tilesGroup,
				ellipsoid,
				nearMargin,
				farMargin
			} = this;
			const maxRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
			if (camera.isPerspectiveCamera) {
				// adjust the clip planes
				const distanceToCenter = _vec$1.setFromMatrixPosition(tilesGroup.worldMatrix).sub(camera.position).getLength();

				// update the projection matrix
				// interpolate from the 25% radius margin around the globe down to the surface
				// so we can avoid z fighting when near value is too far at a high altitude
				const margin = nearMargin * maxRadius;
				const alpha = t3d.MathUtils.clamp((distanceToCenter - maxRadius) / margin, 0, 1);
				const minNear = t3d.MathUtils.lerp(1, 1000, alpha);
				camera.near = Math.max(minNear, distanceToCenter - maxRadius - margin);

				// update the far plane to the horizon distance
				_invMatrix.copy(tilesGroup.worldMatrix).invert();
				_pos$3.copy(camera.position).applyMatrix4(_invMatrix);
				ellipsoid.getPositionToCartographic(_pos$3, _latLon);

				// use a minimum elevation for computing the horizon distance to avoid the far clip
				// plane approaching zero or clipping mountains over the horizon in the distance as
				// the camera goes to or below sea level.
				const elevation = Math.max(ellipsoid.getPositionElevation(_pos$3), MIN_ELEVATION);
				const horizonDistance = ellipsoid.calculateHorizonDistance(_latLon.lat, elevation);
				camera.far = horizonDistance + 0.1 + maxRadius * farMargin;
				camera.updateProjectionMatrix();
			} else {
				this._getVirtualOrthoCameraPosition(camera.position, camera);
				camera.updateMatrix();
				_invMatrix.copy(camera.worldMatrix).invert();
				_vec$1.setFromMatrixPosition(tilesGroup.worldMatrix).applyMatrix4(_invMatrix);
				const distanceToCenter = -_vec$1.z;
				camera.near = distanceToCenter - maxRadius * (1 + nearMargin);
				camera.far = distanceToCenter + 0.1 + maxRadius * farMargin;

				// adjust the position of the ortho camera such that the near value is 0
				camera.position.addScaledVector(_forward, camera.near);
				camera.far -= camera.near;
				camera.near = 0;
				camera.updateProjectionMatrix();
				camera.updateMatrix();
			}
		}

		// resets the "stuck" drag modes
		resetState() {
			super.resetState();
			this._dragMode = 0;
			this._rotationMode = 0;
		}
		_updateInertia(deltaTime) {
			super._updateInertia(deltaTime);
			const {
				globeInertia,
				enableDamping,
				dampingFactor,
				camera,
				cameraRadius,
				minDistance,
				inertiaTargetDistance,
				tilesGroup
			} = this;
			if (!this.enableDamping || this.inertiaStableFrames > 1) {
				this.globeInertiaFactor = 0;
				this.globeInertia.identity();
				return;
			}
			const factor = Math.pow(2, -deltaTime / dampingFactor);
			const stableDistance = Math.max(camera.near, cameraRadius, minDistance, inertiaTargetDistance);
			const resolution = 2 * 1e3;
			const pixelWidth = 2 / resolution;
			const pixelThreshold = 0.25 * pixelWidth;
			_center.setFromMatrixPosition(tilesGroup.worldMatrix);
			if (this.globeInertiaFactor !== 0) {
				// calculate two screen points at 1 pixel apart in our notional resolution so we can stop when the delta is ~ 1 pixel
				// projected into world space
				setRaycasterFromCamera(_ray, _vec$1.set(0, 0, -1), camera);
				_ray.applyMatrix4(camera.viewMatrix);
				_ray.direction.normalize();
				_ray.recast(-_ray.direction.dot(_ray.origin)).at(stableDistance / _ray.direction.z, _vec$1);
				_vec$1.applyMatrix4(camera.worldMatrix);
				setRaycasterFromCamera(_ray, _pos$3.set(pixelThreshold, pixelThreshold, -1), camera);
				_ray.applyMatrix4(camera.viewMatrix);
				_ray.direction.normalize();
				_ray.recast(-_ray.direction.dot(_ray.origin)).at(stableDistance / _ray.direction.z, _pos$3);
				_pos$3.applyMatrix4(camera.worldMatrix);

				// get implied angle
				_vec$1.sub(_center).normalize();
				_pos$3.sub(_center).normalize();
				this.globeInertiaFactor *= factor;
				const threshold = _vec$1.angleTo(_pos$3) / deltaTime;
				const globeAngle = 2 * Math.acos(globeInertia.w) * this.globeInertiaFactor;
				if (globeAngle < threshold || !enableDamping) {
					this.globeInertiaFactor = 0;
					globeInertia.identity();
				}
			}
			if (this.globeInertiaFactor !== 0) {
				// ensure our w component is non-one if the xyz values are
				// non zero to ensure we can animate
				if (globeInertia.w === 1 && (globeInertia.x !== 0 || globeInertia.y !== 0 || globeInertia.z !== 0)) {
					globeInertia.w = Math.min(globeInertia.w, 1 - 1e-9);
				}

				// construct the rotation matrix
				_center.setFromMatrixPosition(tilesGroup.worldMatrix);
				_quaternion$1.identity().slerp(globeInertia, this.globeInertiaFactor * deltaTime);
				makeRotateAroundPoint(_center, _quaternion$1, _rotMatrix);

				// apply the rotation
				camera.worldMatrix.premultiply(_rotMatrix);
				camera.worldMatrix.decompose(camera.position, camera.quaternion, _vec$1);
			}
		}
		_inertiaNeedsUpdate() {
			return super._inertiaNeedsUpdate() || this.globeInertiaFactor !== 0;
		}
		_updatePosition(deltaTime) {
			if (this.state === DRAG) {
				// save the drag mode state so we can update the pivot mesh visuals in "update"
				if (this._dragMode === 0) {
					this._dragMode = this._isNearControls() ? 1 : -1;
				}
				const {
					raycaster,
					camera,
					pivotPoint,
					pointerTracker,
					domElement,
					tilesGroup
				} = this;

				// reuse cache variables
				const pivotDir = _pos$3;
				const newPivotDir = _targetRight;

				// get the pointer and ray
				pointerTracker.getCenterPoint(_pointer);
				mouseToCoords(_pointer.x, _pointer.y, domElement, _pointer);
				setRaycasterFromCamera(raycaster, _pointer, camera);
				_invMatrix.copy(tilesGroup.worldMatrix).invert();

				// transform to ellipsoid frame
				raycaster.ray.applyMatrix4(_invMatrix);

				// construct an ellipsoid that matches a sphere with the radius of the globe so
				// the drag position matches where the initial click was
				const pivotRadius = _vec$1.copy(pivotPoint).applyMatrix4(_invMatrix).getLength();
				_ellipsoid.radius.setScalar(pivotRadius);

				// find the hit point and use the closest point on the horizon if we miss
				if (camera.isPerspectiveCamera) {
					if (!_ellipsoid.intersectRay(raycaster.ray, _vec$1)) {
						closestRaySpherePointFromRotation(raycaster.ray, pivotRadius, _vec$1);
					}
				} else {
					closestRayEllipsoidSurfacePointEstimate(raycaster.ray, _ellipsoid, _vec$1);
				}
				_vec$1.applyMatrix4(tilesGroup.worldMatrix);

				// get the point directions
				_center.setFromMatrixPosition(tilesGroup.worldMatrix);
				pivotDir.subVectors(pivotPoint, _center).normalize();
				newPivotDir.subVectors(_vec$1, _center).normalize();

				// construct the rotation
				_quaternion$1.setFromUnitVectors(newPivotDir, pivotDir);
				makeRotateAroundPoint(_center, _quaternion$1, _rotMatrix);

				// apply the rotation
				camera.worldMatrix.premultiply(_rotMatrix);
				camera.worldMatrix.decompose(camera.position, camera.quaternion, _vec$1);
				if (pointerTracker.getMoveDistance() / deltaTime < 2 * window.devicePixelRatio) {
					this.inertiaStableFrames++;
				} else {
					this.globeInertia.copy(_quaternion$1);
					this.globeInertiaFactor = 1 / deltaTime;
					this.inertiaStableFrames = 0;
				}
			}
			this._alignCameraUp(this.up);
		}

		// disable rotation once we're outside the control transition
		_updateRotation(...args) {
			if (this._rotationMode === 1 || this._isNearControls()) {
				this._rotationMode = 1;
				super._updateRotation(...args);
			} else {
				this.pivotMesh.visible = false;
				this._rotationMode = -1;
			}
			this._alignCameraUp(this.up);
		}
		_updateZoom() {
			const {
				zoomDelta,
				ellipsoid,
				zoomSpeed,
				zoomPoint,
				camera,
				maxZoom,
				state
			} = this;
			if (state !== ZOOM && zoomDelta === 0) {
				return;
			}

			// reset momentum
			this.rotationInertia.set(0, 0);
			this.dragInertia.set(0, 0, 0);
			this.globeInertia.identity();
			this.globeInertiaFactor = 0;

			// used to scale the tilt transitions based on zoom intensity
			const deltaAlpha = t3d.MathUtils.clamp(t3d.MathUtils.mapLinear(Math.abs(zoomDelta), 0, 20, 0, 1), 0, 1);
			if (this._isNearControls() || zoomDelta > 0) {
				this._updateZoomDirection();

				// When zooming try to tilt the camera towards the center of the planet to avoid the globe
				// spinning as you zoom out from the horizon
				if (zoomDelta < 0 && (this.zoomPointSet || this._updateZoomPoint())) {
					// get the forward vector and vector toward the center of the ellipsoid
					_forward.set(0, 0, -1).transformDirection(camera.worldMatrix).normalize();
					_toCenter.copy(this.up).multiplyScalar(-1);

					// Calculate alpha values to use to scale the amount of tilt that occurs as the camera moves.
					// Scales based on mouse position near the horizon and current tilt.
					this.getUpDirection(zoomPoint, _zoomPointUp);
					const upAlpha = t3d.MathUtils.clamp(t3d.MathUtils.mapLinear(-_zoomPointUp.dot(_toCenter), 1, 0.95, 0, 1), 0, 1);
					const forwardAlpha = 1 - _forward.dot(_toCenter);
					const cameraAlpha = camera.isOrthographicCamera ? 0.05 : 1;
					const adjustedDeltaAlpha = t3d.MathUtils.clamp(deltaAlpha * 3, 0, 1);

					// apply scale
					const alpha = Math.min(upAlpha * forwardAlpha * cameraAlpha * adjustedDeltaAlpha, 0.1);
					_toCenter.lerpVectors(_forward, _toCenter, alpha).normalize();

					// perform rotation
					_quaternion$1.setFromUnitVectors(_forward, _toCenter);
					makeRotateAroundPoint(zoomPoint, _quaternion$1, _rotMatrix);
					camera.worldMatrix.premultiply(_rotMatrix);
					camera.worldMatrix.decompose(camera.position, camera.quaternion, _toCenter);

					// update zoom direction
					this.zoomDirection.subVectors(zoomPoint, camera.position).normalize();
				}
				super._updateZoom();
			} else if (camera.isPerspectiveCamera) {
				// orient the camera to focus on the earth during the zoom
				const transitionDistance = this._getPerspectiveTransitionDistance();
				const maxDistance = this._getMaxPerspectiveDistance();
				const distanceAlpha = t3d.MathUtils.mapLinear(this.getDistanceToCenter(), transitionDistance, maxDistance, 0, 1);
				this._tiltTowardsCenter(t3d.MathUtils.lerp(0, 0.4, distanceAlpha * deltaAlpha));
				this._alignCameraUpToNorth(t3d.MathUtils.lerp(0, 0.2, distanceAlpha * deltaAlpha));

				// calculate zoom in a similar way to environment controls so
				// the zoom speeds are comparable
				const dist = this.getDistanceToCenter() - ellipsoid.radius.x;
				const scale = zoomDelta * dist * zoomSpeed * 0.0025;
				const clampedScale = Math.max(scale, Math.min(this.getDistanceToCenter() - maxDistance, 0));

				// zoom out directly from the globe center
				this.getVectorToCenter(_vec$1).normalize();
				this.camera.position.addScaledVector(_vec$1, clampedScale);
				this.camera.updateMatrix();
				this.zoomDelta = 0;
			} else {
				const transitionZoom = this._getOrthographicTransitionZoom();
				const minZoom = this._getMinOrthographicZoom();
				const distanceAlpha = t3d.MathUtils.mapLinear(camera.zoom, transitionZoom, minZoom, 0, 1);
				this._tiltTowardsCenter(t3d.MathUtils.lerp(0, 0.4, distanceAlpha * deltaAlpha));
				this._alignCameraUpToNorth(t3d.MathUtils.lerp(0, 0.2, distanceAlpha * deltaAlpha));
				const scale = this.zoomDelta;
				const normalizedDelta = Math.pow(0.95, Math.abs(scale * 0.05));
				const scaleFactor = scale > 0 ? 1 / Math.abs(normalizedDelta) : normalizedDelta;
				const maxScaleFactor = minZoom / camera.zoom;
				const clampedScaleFactor = Math.max(scaleFactor * zoomSpeed, Math.min(maxScaleFactor, 1));
				camera.zoom = Math.min(maxZoom, camera.zoom * clampedScaleFactor);
				camera.updateProjectionMatrix();
				this.zoomDelta = 0;
				this.zoomDirectionSet = false;
			}
		}

		// tilt the camera to align with north
		_alignCameraUpToNorth(alpha) {
			const {
				tilesGroup
			} = this;
			_globalUp.set(0, 0, 1).transformDirection(tilesGroup.worldMatrix);
			this._alignCameraUp(_globalUp, alpha);
		}

		// tilt the camera to align with the provided "up" value
		_alignCameraUp(up, alpha = null) {
			const {
				camera
			} = this;
			_forward.set(0, 0, -1).transformDirection(camera.worldMatrix);
			_right.set(-1, 0, 0).transformDirection(camera.worldMatrix);
			_targetRight.crossVectors(up, _forward);

			// compute the alpha based on how far away from boresight the up vector is
			// so we can ease into the correct orientation
			if (alpha === null) {
				alpha = 1 - Math.abs(_forward.dot(up));
				alpha = t3d.MathUtils.mapLinear(alpha, 0, 1, -0.01, 1);
				alpha = t3d.MathUtils.clamp(alpha, 0, 1) ** 2;
			}
			_targetRight.lerp(_right, 1 - alpha).normalize();
			_quaternion$1.setFromUnitVectors(_right, _targetRight);
			camera.quaternion.premultiply(_quaternion$1);
			camera.updateMatrix();
		}

		// tilt the camera to look at the center of the globe
		_tiltTowardsCenter(alpha) {
			const {
				camera,
				tilesGroup
			} = this;
			_forward.set(0, 0, -1).transformDirection(camera.worldMatrix).normalize();
			_vec$1.setFromMatrixPosition(tilesGroup.worldMatrix).sub(camera.position).normalize();
			_vec$1.lerp(_forward, 1 - alpha).normalize();
			_quaternion$1.setFromUnitVectors(_forward, _vec$1);
			camera.quaternion.premultiply(_quaternion$1);
			camera.updateMatrix();
		}

		// returns the perspective camera transition distance can move to based on globe size and fov
		_getPerspectiveTransitionDistance() {
			const {
				camera,
				ellipsoid
			} = this;
			if (!camera.isPerspectiveCamera) {
				throw new Error();
			}

			// When the smallest fov spans 65% of the ellipsoid then we use the near controls
			const ellipsoidRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
			const fovHoriz = 2 * Math.atan(Math.tan(t3d.MathUtils.DEG2RAD * camera.fov * 0.5) * camera.aspect);
			const distVert = ellipsoidRadius / Math.tan(t3d.MathUtils.DEG2RAD * camera.fov * 0.5);
			const distHoriz = ellipsoidRadius / Math.tan(fovHoriz * 0.5);
			const dist = Math.max(distVert, distHoriz);
			return dist;
		}

		// returns the max distance the perspective camera can move to based on globe size and fov
		_getMaxPerspectiveDistance() {
			const {
				camera,
				ellipsoid
			} = this;
			if (!camera.isPerspectiveCamera) {
				throw new Error();
			}

			// allow for zooming out such that the ellipsoid is half the size of the largest fov
			const ellipsoidRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
			const fovHoriz = 2 * Math.atan(Math.tan(t3d.MathUtils.DEG2RAD * camera.fov * 0.5) * camera.aspect);
			const distVert = ellipsoidRadius / Math.tan(t3d.MathUtils.DEG2RAD * camera.fov * 0.5);
			const distHoriz = ellipsoidRadius / Math.tan(fovHoriz * 0.5);
			const dist = 2 * Math.max(distVert, distHoriz);
			return dist;
		}

		// returns the transition threshold for orthographic zoom based on the globe size and camera settings
		_getOrthographicTransitionZoom() {
			const {
				camera,
				ellipsoid
			} = this;
			if (!camera.isOrthographicCamera) {
				throw new Error();
			}
			const orthoHeight = camera.top - camera.bottom;
			const orthoWidth = camera.right - camera.left;
			const orthoSize = Math.max(orthoHeight, orthoWidth);
			const ellipsoidRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
			const ellipsoidDiameter = 2 * ellipsoidRadius;
			return 2 * orthoSize / ellipsoidDiameter;
		}

		// returns the minimum allowed orthographic zoom based on the globe size and camera settings
		_getMinOrthographicZoom() {
			const {
				camera,
				ellipsoid
			} = this;
			if (!camera.isOrthographicCamera) {
				throw new Error();
			}
			const orthoHeight = camera.top - camera.bottom;
			const orthoWidth = camera.right - camera.left;
			const orthoSize = Math.min(orthoHeight, orthoWidth);
			const ellipsoidRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
			const ellipsoidDiameter = 2 * ellipsoidRadius;
			return 0.7 * orthoSize / ellipsoidDiameter;
		}

		// returns the "virtual position" of the orthographic based on where it is and
		// where it's looking primarily so we can reasonably position the camera object
		// in space and derive a reasonable "up" value.
		_getVirtualOrthoCameraPosition(target, camera = this.camera) {
			const {
				tilesGroup,
				ellipsoid
			} = this;
			if (!camera.isOrthographicCamera) {
				throw new Error();
			}

			// get ray in globe coordinate frame
			_ray.origin.copy(camera.position);
			_ray.direction.set(0, 0, -1).transformDirection(camera.worldMatrix);
			_invMatrix.copy(tilesGroup.worldMatrix).invert();
			_ray.applyMatrix4(_invMatrix);

			// get the closest point to the ray on the globe in the global coordinate frame
			closestRayEllipsoidSurfacePointEstimate(_ray, ellipsoid, _pos$3);
			_pos$3.applyMatrix4(tilesGroup.worldMatrix);

			// get ortho camera info
			const orthoHeight = camera.top - camera.bottom;
			const orthoWidth = camera.right - camera.left;
			const orthoSize = Math.max(orthoHeight, orthoWidth) / camera.zoom;
			_forward.set(0, 0, -1).transformDirection(camera.worldMatrix);

			// ensure we move the camera exactly along the forward vector to avoid shifting
			// the camera in other directions due to floating point error
			const dist = _pos$3.sub(camera.position).dot(_forward);
			target.copy(camera.position).addScaledVector(_forward, dist - orthoSize * 4);
		}
		_isNearControls() {
			const {
				camera
			} = this;
			if (camera.isPerspectiveCamera) {
				return this.getDistanceToCenter() < this._getPerspectiveTransitionDistance();
			} else {
				return camera.zoom > this._getOrthographicTransitionZoom();
			}
		}
		_raycast(raycaster) {
			const result = super._raycast(raycaster);
			if (result === null) {
				// if there was no hit then fallback to intersecting the ellipsoid.
				const {
					ellipsoid,
					tilesGroup
				} = this;
				_invMatrix.copy(tilesGroup.worldMatrix).invert();
				_ray.copy(raycaster.ray).applyMatrix4(_invMatrix);
				const point = ellipsoid.intersectRay(_ray, _vec$1);
				if (point !== null) {
					return {
						point: point.clone().applyMatrix4(tilesGroup.worldMatrix)
					};
				} else {
					return null;
				}
			} else {
				return result;
			}
		}
	}

	class LoaderBase {
		constructor() {
			this.fetchOptions = {};
			this.workingPath = '';
		}
		load(...args) {
			console.warn('Loader: "load" function has been deprecated in favor of "loadAsync".');
			return this.loadAsync(...args);
		}
		loadAsync(url) {
			return fetch(url, this.fetchOptions).then(res => {
				if (!res.ok) {
					throw new Error(`Failed to load file "${url}" with status ${res.status} : ${res.statusText}`);
				}
				return res.arrayBuffer();
			}).then(buffer => {
				if (this.workingPath === '') {
					this.workingPath = this.workingPathForURL(url);
				}
				return this.parse(buffer);
			});
		}
		resolveExternalURL(url) {
			if (/^[^\\/]/.test(url) && !/^http/.test(url)) {
				return this.workingPath + '/' + url;
			} else {
				return url;
			}
		}
		workingPathForURL(url) {
			const splits = url.split(/[\\/]/g);
			splits.pop();
			const workingPath = splits.join('/');
			return workingPath + '/';
		}
		parse(buffer) {
			throw new Error('LoaderBase: Parse not implemented.');
		}
	}

	const utf8decoder = new TextDecoder();
	function arrayToString(array) {
		return utf8decoder.decode(array);
	}

	/**
	 * Structure almost identical to Cesium, also the comments and the names are kept
	 * https://github.com/CesiumGS/cesium/blob/0a69f67b393ba194eefb7254600811c4b712ddc0/packages/engine/Source/Scene/Implicit3DTileContent.js
	 */
	function isOctreeSubdivision(tile) {
		return tile.__implicitRoot.implicitTiling.subdivisionScheme === 'OCTREE';
	}
	function getBoundsDivider(tile) {
		return isOctreeSubdivision(tile) ? 8 : 4;
	}
	function getSubtreeCoordinates(tile, parentTile) {
		if (!parentTile) {
			return [0, 0, 0];
		}
		const x = 2 * parentTile.__x + tile.__subtreeIdx % 2;
		const y = 2 * parentTile.__y + Math.floor(tile.__subtreeIdx / 2) % 2;
		const z = isOctreeSubdivision(tile) ? 2 * parentTile.__z + Math.floor(tile.__subtreeIdx / 4) % 2 : 0;
		return [x, y, z];
	}
	class SubtreeTile {
		constructor(parentTile, childMortonIndex) {
			this.parent = parentTile;
			this.children = [];
			this.__level = parentTile.__level + 1;
			this.__implicitRoot = parentTile.__implicitRoot;
			// Index inside the tree
			this.__subtreeIdx = childMortonIndex;
			[this.__x, this.__y, this.__z] = getSubtreeCoordinates(this, parentTile);
		}
		static copy(tile) {
			const copyTile = {};
			copyTile.children = [];
			copyTile.__level = tile.__level;
			copyTile.__implicitRoot = tile.__implicitRoot;
			// Index inside the tree
			copyTile.__subtreeIdx = tile.__subtreeIdx;
			[copyTile.__x, copyTile.__y, copyTile.__z] = [tile.__x, tile.__y, tile.__z];
			copyTile.boundingVolume = tile.boundingVolume;
			copyTile.geometricError = tile.geometricError;
			return copyTile;
		}
	}
	class SUBTREELoader extends LoaderBase {
		constructor(tile) {
			super();
			this.tile = tile;
			this.rootTile = tile.__implicitRoot; // The implicit root tile
			this.workingPath = null;
		}

		/**
		 * A helper object for storing the two parts of the subtree binary
		 *
		 * @typedef {object} Subtree
		 * @property {number} version
		 * @property {JSON} subtreeJson
		 * @property {ArrayBuffer} subtreeByte
		 * @private
		 */

		/**
		 *
		 * @param buffer
		 * @return {Subtree}
		 */
		parseBuffer(buffer) {
			const dataView = new DataView(buffer);
			let offset = 0;
			// 16-byte header
			// 4 bytes
			const magic = readMagicBytes(dataView);
			console.assert(magic === 'subt', 'SUBTREELoader: The magic bytes equal "subt".');
			offset += 4;
			// 4 bytes
			const version = dataView.getUint32(offset, true);
			console.assert(version === 1, 'SUBTREELoader: The version listed in the header is "1".');
			offset += 4;
			// From Cesium
			// Read the bottom 32 bits of the 64-bit byte length.
			// This is ok for now because:
			// 1) not all browsers have native 64-bit operations
			// 2) the data is well under 4GB
			// 8 bytes
			const jsonLength = dataView.getUint32(offset, true);
			offset += 8;
			// 8 bytes
			const byteLength = dataView.getUint32(offset, true);
			offset += 8;
			const subtreeJson = JSON.parse(arrayToString(new Uint8Array(buffer, offset, jsonLength)));
			offset += jsonLength;
			const subtreeByte = buffer.slice(offset, offset + byteLength);
			return {
				version,
				subtreeJson,
				subtreeByte
			};
		}
		async parse(buffer) {
			// todo here : handle json
			const subtree = this.parseBuffer(buffer);
			const subtreeJson = subtree.subtreeJson;

			// TODO Handle metadata
			/*
			 const subtreeMetadata = subtreeJson.subtreeMetadata;
			 subtree._metadata = subtreeMetadata;
			*/

			/*
				Tile availability indicates which tiles exist within the subtree
				Content availability indicates which tiles have associated content resources
				Child subtree availability indicates what subtrees are reachable from this subtree
			*/

			// After identifying how availability is stored, put the results in this new array for consistent processing later
			subtreeJson.contentAvailabilityHeaders = [].concat(subtreeJson.contentAvailability);
			const bufferHeaders = this.preprocessBuffers(subtreeJson.buffers);
			const bufferViewHeaders = this.preprocessBufferViews(subtreeJson.bufferViews, bufferHeaders);

			// Buffers and buffer views are inactive until explicitly marked active.
			// This way we can avoid fetching buffers that will not be used.
			this.markActiveBufferViews(subtreeJson, bufferViewHeaders);

			// Await the active buffers. If a buffer is external (isExternal === true),
			// fetch it from its URI.
			const buffersU8 = await this.requestActiveBuffers(bufferHeaders, subtree.subtreeByte);
			const bufferViewsU8 = this.parseActiveBufferViews(bufferViewHeaders, buffersU8);
			this.parseAvailability(subtree, subtreeJson, bufferViewsU8);
			this.expandSubtree(this.tile, subtree);
		}

		/**
		 * Determine which buffer views need to be loaded into memory. This includes:
		 *
		 * <ul>
		 * <li>The tile availability bitstream (if a bitstream is defined)</li>
		 * <li>The content availability bitstream(s) (if a bitstream is defined)</li>
		 * <li>The child subtree availability bitstream (if a bitstream is defined)</li>
		 * </ul>
		 *
		 * <p>
		 * This function modifies the buffer view headers' isActive flags in place.
		 * </p>
		 *
		 * @param {JSON} subtreeJson The JSON chunk from the subtree
		 * @param {BufferViewHeader[]} bufferViewHeaders The preprocessed buffer view headers
		 * @private
		 */
		markActiveBufferViews(subtreeJson, bufferViewHeaders) {
			let header;
			const tileAvailabilityHeader = subtreeJson.tileAvailability;
			// Check for bitstream first, which is part of the current schema.
			// bufferView is the name of the bitstream from an older schema.
			if (!isNaN(tileAvailabilityHeader.bitstream)) {
				header = bufferViewHeaders[tileAvailabilityHeader.bitstream];
			} else if (!isNaN(tileAvailabilityHeader.bufferView)) {
				header = bufferViewHeaders[tileAvailabilityHeader.bufferView];
			}
			if (header) {
				header.isActive = true;
				header.bufferHeader.isActive = true;
			}
			const contentAvailabilityHeaders = subtreeJson.contentAvailabilityHeaders;
			for (let i = 0; i < contentAvailabilityHeaders.length; i++) {
				header = undefined;
				if (!isNaN(contentAvailabilityHeaders[i].bitstream)) {
					header = bufferViewHeaders[contentAvailabilityHeaders[i].bitstream];
				} else if (!isNaN(contentAvailabilityHeaders[i].bufferView)) {
					header = bufferViewHeaders[contentAvailabilityHeaders[i].bufferView];
				}
				if (header) {
					header.isActive = true;
					header.bufferHeader.isActive = true;
				}
			}
			header = undefined;
			const childSubtreeAvailabilityHeader = subtreeJson.childSubtreeAvailability;
			if (!isNaN(childSubtreeAvailabilityHeader.bitstream)) {
				header = bufferViewHeaders[childSubtreeAvailabilityHeader.bitstream];
			} else if (!isNaN(childSubtreeAvailabilityHeader.bufferView)) {
				header = bufferViewHeaders[childSubtreeAvailabilityHeader.bufferView];
			}
			if (header) {
				header.isActive = true;
				header.bufferHeader.isActive = true;
			}
		}

		/**
		 * Go through the list of buffers and gather all the active ones into
		 * a dictionary.
		 * <p>
		 * The results are put into a dictionary object. The keys are indices of
		 * buffers, and the values are Uint8Arrays of the contents. Only buffers
		 * marked with the isActive flag are fetched.
		 * </p>
		 * <p>
		 * The internal buffer (the subtree's binary chunk) is also stored in this
		 * dictionary if it is marked active.
		 * </p>
		 * @param {BufferHeader[]} bufferHeaders The preprocessed buffer headers
		 * @param {ArrayBuffer} internalBuffer The binary chunk of the subtree file
		 * @returns {object} buffersU8 A dictionary of buffer index to a Uint8Array of its contents.
		 * @private
		 */
		async requestActiveBuffers(bufferHeaders, internalBuffer) {
			const promises = [];
			for (let i = 0; i < bufferHeaders.length; i++) {
				const bufferHeader = bufferHeaders[i];
				// If the buffer is not active, resolve with undefined.
				if (!bufferHeader.isActive) {
					promises.push(Promise.resolve());
				} else if (bufferHeader.isExternal) {
					// Get the absolute URI of the external buffer.
					const url = this.parseImplicitURIBuffer(this.tile, this.rootTile.implicitTiling.subtrees.uri, bufferHeader.uri);
					const fetchPromise = fetch(url, this.fetchOptions).then(response => {
						if (!response.ok) {
							throw new Error(`SUBTREELoader: Failed to load external buffer from ${bufferHeader.uri} with error code ${response.status}.`);
						}
						return response.arrayBuffer();
					}).then(arrayBuffer => new Uint8Array(arrayBuffer));
					promises.push(fetchPromise);
				} else {
					promises.push(Promise.resolve(new Uint8Array(internalBuffer)));
				}
			}
			const bufferResults = await Promise.all(promises);
			const buffersU8 = {};
			for (let i = 0; i < bufferResults.length; i++) {
				const result = bufferResults[i];
				if (result) {
					buffersU8[i] = result;
				}
			}
			return buffersU8;
		}

		/**
		 * Go through the list of buffer views, and if they are marked as active,
		 * extract a subarray from one of the active buffers.
		 *
		 * @param {BufferViewHeader[]} bufferViewHeaders
		 * @param {object} buffersU8 A dictionary of buffer index to a Uint8Array of its contents.
		 * @returns {object} A dictionary of buffer view index to a Uint8Array of its contents.
		 * @private
		 */
		parseActiveBufferViews(bufferViewHeaders, buffersU8) {
			const bufferViewsU8 = {};
			for (let i = 0; i < bufferViewHeaders.length; i++) {
				const bufferViewHeader = bufferViewHeaders[i];
				if (!bufferViewHeader.isActive) {
					continue;
				}
				const start = bufferViewHeader.byteOffset;
				const end = start + bufferViewHeader.byteLength;
				const buffer = buffersU8[bufferViewHeader.buffer];
				bufferViewsU8[i] = buffer.slice(start, end);
			}
			return bufferViewsU8;
		}

		/**
		 * A buffer header is the JSON header from the subtree JSON chunk plus
		 * a couple extra boolean flags for easy reference.
		 *
		 * Buffers are assumed inactive until explicitly marked active. This is used
		 * to avoid fetching unneeded buffers.
		 *
		 * @typedef {object} BufferHeader
		 * @property {boolean} isActive Whether this buffer is currently used.
		 * @property {string} [uri] The URI of the buffer (external buffers only)
		 * @property {number} byteLength The byte length of the buffer, including any padding contained within.
		 * @private
		 */

		/**
		 * Iterate over the list of buffers from the subtree JSON and add the isActive field for easier parsing later.
		 * This modifies the objects in place.
		 * @param {Object[]} [bufferHeaders=[]] The JSON from subtreeJson.buffers.
		 * @returns {BufferHeader[]} The same array of headers with additional fields.
		 * @private
		 */
		preprocessBuffers(bufferHeaders = []) {
			for (let i = 0; i < bufferHeaders.length; i++) {
				const bufferHeader = bufferHeaders[i];
				bufferHeader.isActive = false;
				bufferHeader.isExternal = !!bufferHeader.uri;
			}
			return bufferHeaders;
		}

		/**
		 * A buffer view header is the JSON header from the subtree JSON chunk plus
		 * the isActive flag and a reference to the header for the underlying buffer.
		 *
		 * @typedef {object} BufferViewHeader
		 * @property {BufferHeader} bufferHeader A reference to the header for the underlying buffer
		 * @property {boolean} isActive Whether this bufferView is currently used.
		 * @property {number} buffer The index of the underlying buffer.
		 * @property {number} byteOffset The start byte of the bufferView within the buffer.
		 * @property {number} byteLength The length of the bufferView. No padding is included in this length.
		 * @private
		 */

		/**
		 * Iterate the list of buffer views from the subtree JSON and add the
		 * isActive flag. Also save a reference to the bufferHeader.
		 *
		 * @param {Object[]} [bufferViewHeaders=[]] The JSON from subtree.bufferViews.
		 * @param {BufferHeader[]} bufferHeaders The preprocessed buffer headers.
		 * @returns {BufferViewHeader[]} The same array of bufferView headers with additional fields.
		 * @private
		 */
		preprocessBufferViews(bufferViewHeaders = [], bufferHeaders) {
			for (let i = 0; i < bufferViewHeaders.length; i++) {
				const bufferViewHeader = bufferViewHeaders[i];
				bufferViewHeader.bufferHeader = bufferHeaders[bufferViewHeader.buffer];
				bufferViewHeader.isActive = false;
				// Keep the external flag for potential use in requestActiveBuffers
				bufferViewHeader.isExternal = bufferViewHeader.bufferHeader.isExternal;
			}
			return bufferViewHeaders;
		}

		/**
		 * Parse the three availability bitstreams and store them in the subtree.
		 *
		 * @param {Subtree} subtree The subtree to modify.
		 * @param {Object} subtreeJson The subtree JSON.
		 * @param {Object} bufferViewsU8 A dictionary of buffer view index to a Uint8Array of its contents.
		 * @private
		 */
		parseAvailability(subtree, subtreeJson, bufferViewsU8) {
			const branchingFactor = getBoundsDivider(this.rootTile);
			const subtreeLevels = this.rootTile.implicitTiling.subtreeLevels;
			const tileAvailabilityBits = (Math.pow(branchingFactor, subtreeLevels) - 1) / (branchingFactor - 1);
			const childSubtreeBits = Math.pow(branchingFactor, subtreeLevels);
			subtree._tileAvailability = this.parseAvailabilityBitstream(subtreeJson.tileAvailability, bufferViewsU8, tileAvailabilityBits);
			subtree._contentAvailabilityBitstreams = [];
			for (let i = 0; i < subtreeJson.contentAvailabilityHeaders.length; i++) {
				const bitstream = this.parseAvailabilityBitstream(subtreeJson.contentAvailabilityHeaders[i], bufferViewsU8,
				// content availability has the same length as tile availability.
				tileAvailabilityBits);
				subtree._contentAvailabilityBitstreams.push(bitstream);
			}
			subtree._childSubtreeAvailability = this.parseAvailabilityBitstream(subtreeJson.childSubtreeAvailability, bufferViewsU8, childSubtreeBits);
		}

		/**
		 * Given the JSON describing an availability bitstream, turn it into an
		 * in-memory representation using an object. This handles bitstreams from a bufferView.
		 *
		 * @param {Object} availabilityJson A JSON object representing the availability.
		 * @param {Object} bufferViewsU8 A dictionary of buffer view index to its Uint8Array contents.
		 * @param {number} lengthBits The length of the availability bitstream in bits.
		 * @returns {object}
		 * @private
		 */
		parseAvailabilityBitstream(availabilityJson, bufferViewsU8, lengthBits) {
			if (!isNaN(availabilityJson.constant)) {
				return {
					constant: Boolean(availabilityJson.constant),
					lengthBits: lengthBits
				};
			}
			let bufferView;
			// Check for bitstream first, which is part of the current schema.
			// bufferView is the name of the bitstream from an older schema.
			if (!isNaN(availabilityJson.bitstream)) {
				bufferView = bufferViewsU8[availabilityJson.bitstream];
			} else if (!isNaN(availabilityJson.bufferView)) {
				bufferView = bufferViewsU8[availabilityJson.bufferView];
			}
			return {
				bitstream: bufferView,
				lengthBits: lengthBits
			};
		}

		/**
		 * Expand a single subtree tile. This transcodes the subtree into
		 * a tree of {@link SubtreeTile}. The root of this tree is stored in
		 * the placeholder tile's children array. This method also creates
		 * tiles for the child subtrees to be lazily expanded as needed.
		 *
		 * @param {Object | SubtreeTile} subtreeRoot The first node of the subtree.
		 * @param {Subtree} subtree The parsed subtree.
		 * @private
		 */
		expandSubtree(subtreeRoot, subtree) {
			// TODO If multiple contents were supported then this tile could contain both renderable and un renderable content.
			const contentTile = SubtreeTile.copy(subtreeRoot);
			// If the subtree root tile has content, then create a placeholder child with cloned parameters
			// Todo Multiple contents not handled, keep the first content found
			for (let i = 0; subtree && i < subtree._contentAvailabilityBitstreams.length; i++) {
				if (subtree && this.getBit(subtree._contentAvailabilityBitstreams[i], 0)) {
					// Create a child holding the content uri, this child is similar to its parent and doesn't have any children.
					contentTile.content = {
						uri: this.parseImplicitURI(subtreeRoot, this.rootTile.content.uri)
					};
					break;
				}
			}
			subtreeRoot.children.push(contentTile);
			// Creating each leaf inside the current subtree.
			const bottomRow = this.transcodeSubtreeTiles(contentTile, subtree);
			// For each child subtree, create a tile containing the uri of the next subtree to fetch.
			const childSubtrees = this.listChildSubtrees(subtree, bottomRow);
			for (let i = 0; i < childSubtrees.length; i++) {
				const subtreeLocator = childSubtrees[i];
				const leafTile = subtreeLocator.tile;
				const subtreeTile = this.deriveChildTile(null, leafTile, null, subtreeLocator.childMortonIndex);
				// Assign subtree uri as content.
				subtreeTile.content = {
					uri: this.parseImplicitURI(subtreeTile, this.rootTile.implicitTiling.subtrees.uri)
				};
				leafTile.children.push(subtreeTile);
			}
		}

		/**
		 * Transcode the implicitly defined tiles within this subtree and generate
		 * explicit {@link SubtreeTile} objects. This function only transcodes tiles,
		 * child subtrees are handled separately.
		 *
		 * @param {Object | SubtreeTile} subtreeRoot The root of the current subtree.
		 * @param {Subtree} subtree The subtree to get availability information.
		 * @returns {Array} The bottom row of transcoded tiles. This is helpful for processing child subtrees.
		 * @private
		 */
		transcodeSubtreeTiles(subtreeRoot, subtree) {
			// Sliding window over the levels of the tree.
			// Each row is branchingFactor * length of previous row.
			// Tiles within a row are ordered by Morton index.
			let parentRow = [subtreeRoot];
			let currentRow = [];
			for (let level = 1; level < this.rootTile.implicitTiling.subtreeLevels; level++) {
				const branchingFactor = getBoundsDivider(this.rootTile);
				const levelOffset = (Math.pow(branchingFactor, level) - 1) / (branchingFactor - 1);
				const numberOfChildren = branchingFactor * parentRow.length;
				for (let childMortonIndex = 0; childMortonIndex < numberOfChildren; childMortonIndex++) {
					const childBitIndex = levelOffset + childMortonIndex;
					const parentMortonIndex = childMortonIndex >> Math.log2(branchingFactor);
					const parentTile = parentRow[parentMortonIndex];
					// Check if tile is available.
					if (!this.getBit(subtree._tileAvailability, childBitIndex)) {
						currentRow.push(undefined);
						continue;
					}
					// Create a tile and add it as a child.
					const childTile = this.deriveChildTile(subtree, parentTile, childBitIndex, childMortonIndex);
					parentTile.children.push(childTile);
					currentRow.push(childTile);
				}
				parentRow = currentRow;
				currentRow = [];
			}
			return parentRow;
		}

		/**
		 * Given a parent tile and information about which child to create, derive
		 * the properties of the child tile implicitly.
		 * <p>
		 * This creates a real tile for rendering.
		 * </p>
		 *
		 * @param {Subtree} subtree The subtree the child tile belongs to.
		 * @param {Object | SubtreeTile} parentTile The parent of the new child tile.
		 * @param {number} childBitIndex The index of the child tile within the tile's availability information.
		 * @param {number} childMortonIndex The morton index of the child tile relative to its parent.
		 * @returns {SubtreeTile} The new child tile.
		 * @private
		 */
		deriveChildTile(subtree, parentTile, childBitIndex, childMortonIndex) {
			const subtreeTile = new SubtreeTile(parentTile, childMortonIndex);
			subtreeTile.boundingVolume = this.getTileBoundingVolume(subtreeTile);
			subtreeTile.geometricError = this.getGeometricError(subtreeTile);
			// Todo Multiple contents not handled, keep the first found content.
			for (let i = 0; subtree && i < subtree._contentAvailabilityBitstreams.length; i++) {
				if (subtree && this.getBit(subtree._contentAvailabilityBitstreams[i], childBitIndex)) {
					subtreeTile.content = {
						uri: this.parseImplicitURI(subtreeTile, this.rootTile.content.uri)
					};
					break;
				}
			}
			return subtreeTile;
		}

		/**
		 * Get a bit from the bitstream as a Boolean. If the bitstream
		 * is a constant, the constant value is returned instead.
		 *
		 * @param {ParsedBitstream} object
		 * @param {number} index The integer index of the bit.
		 * @returns {boolean} The value of the bit.
		 * @private
		 */
		getBit(object, index) {
			if (index < 0 || index >= object.lengthBits) {
				throw new Error('Bit index out of bounds.');
			}
			if (object.constant !== undefined) {
				return object.constant;
			}
			// byteIndex is floor(index / 8)
			const byteIndex = index >> 3;
			const bitIndex = index % 8;
			return (new Uint8Array(object.bitstream)[byteIndex] >> bitIndex & 1) === 1;
		}

		/**
		 * //TODO Adapt for Sphere
		 * To maintain numerical stability during this subdivision process,
		 * the actual bounding volumes should not be computed progressively by subdividing a non-root tile volume.
		 * Instead, the exact bounding volumes are computed directly for a given level.
		 * @param {Object | SubtreeTile} tile
		 * @return {Object} object containing the bounding volume.
		 */
		getTileBoundingVolume(tile) {
			const boundingVolume = {};
			if (this.rootTile.boundingVolume.region) {
				const region = [...this.rootTile.boundingVolume.region];
				const minX = region[0];
				const maxX = region[2];
				const minY = region[1];
				const maxY = region[3];
				const sizeX = (maxX - minX) / Math.pow(2, tile.__level);
				const sizeY = (maxY - minY) / Math.pow(2, tile.__level);
				region[0] = minX + sizeX * tile.__x; // west
				region[2] = minX + sizeX * (tile.__x + 1); // east
				region[1] = minY + sizeY * tile.__y; // south
				region[3] = minY + sizeY * (tile.__y + 1); // north
				for (let k = 0; k < 4; k++) {
					const coord = region[k];
					if (coord < -Math.PI) {
						region[k] += 2 * Math.PI;
					} else if (coord > Math.PI) {
						region[k] -= 2 * Math.PI;
					}
				}
				// Also divide the height in the case of octree.
				if (isOctreeSubdivision(tile)) {
					const minZ = region[4];
					const maxZ = region[5];
					const sizeZ = (maxZ - minZ) / Math.pow(2, tile.__level);
					region[4] = minZ + sizeZ * tile.__z; // minimum height
					region[5] = minZ + sizeZ * (tile.__z + 1); // maximum height
				}
				boundingVolume.region = region;
			}
			if (this.rootTile.boundingVolume.box) {
				// 0-2: center of the box
				// 3-5: x axis direction and half length
				// 6-8: y axis direction and half length
				// 9-11: z axis direction and half length
				const box = [...this.rootTile.boundingVolume.box];
				const cellSteps = 2 ** tile.__level - 1;
				const scale = Math.pow(2, -tile.__level);
				const axisNumber = isOctreeSubdivision(tile) ? 3 : 2;
				for (let i = 0; i < axisNumber; i++) {
					// scale the bounds axes
					box[3 + i * 3 + 0] *= scale;
					box[3 + i * 3 + 1] *= scale;
					box[3 + i * 3 + 2] *= scale;
					// axis vector
					const x = box[3 + i * 3 + 0];
					const y = box[3 + i * 3 + 1];
					const z = box[3 + i * 3 + 2];
					// adjust the center by the x, y and z axes
					const axisOffset = i === 0 ? tile.__x : i === 1 ? tile.__y : tile.__z;
					box[0] += 2 * x * (-0.5 * cellSteps + axisOffset);
					box[1] += 2 * y * (-0.5 * cellSteps + axisOffset);
					box[2] += 2 * z * (-0.5 * cellSteps + axisOffset);
				}
				boundingVolume.box = box;
			}
			return boundingVolume;
		}

		/**
		 * Each child’s geometricError is half of its parent’s geometricError.
		 * @param {Object | SubtreeTile} tile
		 * @return {number}
		 */
		getGeometricError(tile) {
			return this.rootTile.geometricError / Math.pow(2, tile.__level);
		}

		/**
		 * Determine what child subtrees exist and return a list of information.
		 *
		 * @param {Object} subtree The subtree for looking up availability.
		 * @param {Array} bottomRow The bottom row of tiles in a transcoded subtree.
		 * @returns {[]} A list of identifiers for the child subtrees.
		 * @private
		 */
		listChildSubtrees(subtree, bottomRow) {
			const results = [];
			const branchingFactor = getBoundsDivider(this.rootTile);
			for (let i = 0; i < bottomRow.length; i++) {
				const leafTile = bottomRow[i];
				if (leafTile === undefined) {
					continue;
				}
				for (let j = 0; j < branchingFactor; j++) {
					const index = i * branchingFactor + j;
					if (this.getBit(subtree._childSubtreeAvailability, index)) {
						results.push({
							tile: leafTile,
							childMortonIndex: index
						});
					}
				}
			}
			return results;
		}
		/**
		 * Replaces placeholder tokens in a URI template with the corresponding tile properties.
		 *
		 * The URI template should contain the tokens:
		 * - `{level}` for the tile's subdivision level.
		 * - `{x}` for the tile's x-coordinate.
		 * - `{y}` for the tile's y-coordinate.
		 * - `{z}` for the tile's z-coordinate.
		 *
		 * @param {Object} tile - The tile object containing properties __level, __x, __y, and __z.
		 * @param {string} uri - The URI template string with placeholders.
		 * @returns {string} The URI with placeholders replaced by the tile's properties.
		 */
		parseImplicitURI(tile, uri) {
			uri = uri.replace('{level}', tile.__level);
			uri = uri.replace('{x}', tile.__x);
			uri = uri.replace('{y}', tile.__y);
			uri = uri.replace('{z}', tile.__z);
			return uri;
		}

		/**
		 * Generates the full external buffer URI for a tile by combining an implicit URI with a buffer URI.
		 *
		 * First, it parses the implicit URI using the tile properties and the provided template. Then, it creates a new URL
		 * relative to the tile's base path, removes the last path segment, and appends the buffer URI.
		 *
		 * @param {Object} tile - The tile object that contains properties:
		 *	 - __level: the subdivision level,
		 *	 - __x, __y, __z: the tile coordinates,
		 * @param {string} uri - The URI template string with placeholders for the tile (e.g., `{level}`, `{x}`, `{y}`, `{z}`).
		 * @param {string} bufUri - The buffer file name to append (e.g., "0_1.bin").
		 * @returns {string} The full external buffer URI.
		 */
		parseImplicitURIBuffer(tile, uri, bufUri) {
			// Generate the base tile URI by replacing placeholders
			const subUri = this.parseImplicitURI(tile, uri);

			// Create a URL object relative to the tile's base path
			const url = new URL(subUri, this.workingPath + '/');

			// Remove the last path segment
			url.pathname = url.pathname.substring(0, url.pathname.lastIndexOf('/'));

			// Construct the final URL with the buffer URI appended
			return new URL(url.pathname + '/' + bufUri, this.workingPath + '/').toString();
		}
	}

	class ImplicitTilingPlugin {
		constructor() {
			this.name = 'IMPLICIT_TILING_PLUGIN';
		}
		init(tiles) {
			this.tiles = tiles;
		}
		preprocessNode(tile, tileSetDir, parentTile) {
			if (tile.implicitTiling) {
				tile.__hasUnrenderableContent = true;
				tile.__hasRenderableContent = false;

				// Declare some properties
				tile.__subtreeIdx = 0; // Idx of the tile in its subtree
				tile.__implicitRoot = tile; // Keep this tile as an Implicit Root Tile

				// Coords of the tile
				tile.__x = 0;
				tile.__y = 0;
				tile.__z = 0;
				tile.__level = 0;
			} else if (/.subtree$/i.test(tile.content?.uri)) {
				// Handling content uri pointing to a subtree file
				tile.__hasUnrenderableContent = true;
				tile.__hasRenderableContent = false;
			}
		}
		parseTile(buffer, tile, extension) {
			if (/^subtree$/i.test(extension)) {
				const loader = new SUBTREELoader(tile);
				loader.workingPath = tile.__basePath;
				loader.fetchOptions = this.tiles.fetchOptions;
				return loader.parse(buffer);
			}
		}
		preprocessURL(url, tile) {
			if (tile && tile.implicitTiling) {
				const implicitUri = tile.implicitTiling.subtrees.uri.replace('{level}', tile.__level).replace('{x}', tile.__x).replace('{y}', tile.__y).replace('{z}', tile.__z);
				return new URL(implicitUri, tile.__basePath + '/').toString();
			}
			return url;
		}
		disposeTile(tile) {
			if (/.subtree$/i.test(tile.content?.uri)) {
				// TODO: ideally the plugin doesn't need to know about children being processed
				tile.children.length = 0;
				tile.__childrenProcessed = 0;
			}
		}
	}

	class FadeManager {
		constructor() {
			this.duration = 250;
			this.fadeCount = 0;
			this._lastTick = -1;
			this._fadeState = new Map();
			this.onFadeComplete = null;
			this.onFadeStart = null;
			this.onFadeSetComplete = null;
			this.onFadeSetStart = null;
		}

		// delete the object from the fade, reset the material data
		deleteObject(object) {
			if (!object) {
				return;
			}
			this.completeFade(object);
		}

		// Ensure we're storing a fade timer for the provided object
		// Returns whether a new state had to be added
		guaranteeState(object) {
			const fadeState = this._fadeState;
			if (fadeState.has(object)) {
				return false;
			}
			const state = {
				fadeInTarget: 0,
				fadeOutTarget: 0,
				fadeIn: 0,
				fadeOut: 0
			};
			fadeState.set(object, state);
			return true;
		}

		// Force the fade to complete in the direction it is already trending
		completeFade(object) {
			const fadeState = this._fadeState;
			if (!fadeState.has(object)) {
				return;
			}
			const visible = fadeState.get(object).fadeOutTarget === 0;
			fadeState.delete(object);

			// fire events
			this.fadeCount--;
			if (this.onFadeComplete) {
				this.onFadeComplete(object, visible);
			}
			if (this.fadeCount === 0 && this.onFadeSetComplete) {
				this.onFadeSetComplete();
			}
		}
		completeAllFades() {
			this._fadeState.forEach((value, key) => {
				this.completeFade(key);
			});
		}
		forEachObject(cb) {
			this._fadeState.forEach((info, object) => {
				cb(object, info);
			});
		}

		// Fade the object in
		fadeIn(object) {
			const noState = this.guaranteeState(object);
			const state = this._fadeState.get(object);
			state.fadeInTarget = 1;
			state.fadeOutTarget = 0;
			state.fadeOut = 0;

			// Fire events
			if (noState) {
				this.fadeCount++;
				if (this.fadeCount === 1 && this.onFadeSetStart) {
					this.onFadeSetStart();
				}
				if (this.onFadeStart) {
					this.onFadeStart(object);
				}
			}
		}

		// Fade the object out
		fadeOut(object) {
			const noState = this.guaranteeState(object);
			const state = this._fadeState.get(object);
			state.fadeOutTarget = 1;

			// Fire events and initialize state
			if (noState) {
				state.fadeInTarget = 1;
				state.fadeIn = 1;
				this.fadeCount++;
				if (this.fadeCount === 1 && this.onFadeSetStart) {
					this.onFadeSetStart();
				}
				if (this.onFadeStart) {
					this.onFadeStart(object);
				}
			}
		}
		isFading(object) {
			return this._fadeState.has(object);
		}
		isFadingOut(object) {
			const state = this._fadeState.get(object);
			return state && state.fadeOutTarget === 1;
		}

		// Tick the fade timer for each actively fading object
		update() {
			// clamp delta in case duration is really small or 0
			const time = window.performance.now();
			if (this._lastTick === -1) {
				this._lastTick = time;
			}
			const delta = t3d.MathUtils.clamp((time - this._lastTick) / this.duration, 0, 1);
			this._lastTick = time;
			const fadeState = this._fadeState;
			fadeState.forEach((state, object) => {
				// tick the fade values
				const {
					fadeOutTarget,
					fadeInTarget
				} = state;
				let {
					fadeOut,
					fadeIn
				} = state;
				const fadeInSign = Math.sign(fadeInTarget - fadeIn);
				fadeIn = t3d.MathUtils.clamp(fadeIn + fadeInSign * delta, 0, 1);
				const fadeOutSign = Math.sign(fadeOutTarget - fadeOut);
				fadeOut = t3d.MathUtils.clamp(fadeOut + fadeOutSign * delta, 0, 1);
				state.fadeIn = fadeIn;
				state.fadeOut = fadeOut;

				// Check if the fade in and fade out animations are complete
				const fadeOutComplete = fadeOut === 1 || fadeOut === 0;
				const fadeInComplete = fadeIn === 1 || fadeIn === 0;

				// If they are or the fade out animation is further along than the
				// fade in animation then mark the fade as completed for this tile
				if (fadeOutComplete && fadeInComplete || fadeOut >= fadeIn) {
					this.completeFade(object);
				}
			});
		}
	}

	// Adjusts the provided material to support fading in and out using a bayer pattern.
	function wrapFadeMaterial(material) {
		material.shaderName = `${material.shaderName || material.type}_fade`;
		material.defines.FEATURE_FADE = 0;
		material.uniforms.fadeIn = 0;
		material.uniforms.fadeOut = 0;
		const fragmentShader = material.fragmentShader || (material.type === t3d.MATERIAL_TYPE.BASIC ? t3d.ShaderLib.basic_frag : t3d.ShaderLib.pbr_frag);
		material.type = t3d.MATERIAL_TYPE.SHADER;
		material.vertexShader = material.type === t3d.MATERIAL_TYPE.BASIC ? t3d.ShaderLib.basic_vert : t3d.ShaderLib.pbr_vert;
		material.fragmentShader = fragmentShader.replace(/void main\(/, value => /* glsl */`
			#if FEATURE_FADE

			// adapted from https://www.shadertoy.com/view/Mlt3z8
			float bayerDither2x2(vec2 v) {
				return mod(3.0 * v.y + 2.0 * v.x, 4.0);
			}

			float bayerDither4x4(vec2 v) {
				vec2 P1 = mod(v, 2.0);
				vec2 P2 = floor(0.5 * mod(v, 4.0));
				return 4.0 * bayerDither2x2(P1) + bayerDither2x2(P2);
			}

			uniform float fadeIn;
			uniform float fadeOut;

			#endif

			${value}
		`).replace(/#include <end_frag>/, value => /* glsl */`
			${value}

			#if FEATURE_FADE

			float bayerValue = bayerDither4x4(floor(mod(gl_FragCoord.xy, 4.0)));
			float bayerBins = 16.0;
			float dither = (0.5 + bayerValue) / bayerBins;
			if (dither >= fadeIn) {
				discard;
			}

			if (dither < fadeOut) {
				discard;
			}

			#endif
		`);
		return material.uniforms;
	}

	// Class for managing and updating extended fade parameters
	class FadeMaterialManager {
		constructor() {
			this._fadeParams = new WeakMap();
			this.fading = 0;
		}

		// Set the fade parameters for the given scene
		setFade(scene, fadeIn, fadeOut) {
			if (!scene) {
				return;
			}

			// traverse the scene and update the fade parameters of all materials
			const fadeParams = this._fadeParams;
			scene.traverse(child => {
				const material = child.material;
				if (material) {
					const params = fadeParams.get(material);
					params.fadeIn = fadeIn;
					params.fadeOut = fadeOut;
					const fadeInComplete = fadeIn === 0 || fadeIn === 1;
					const fadeOutComplete = fadeOut === 0 || fadeOut === 1;
					const value = Number(!fadeInComplete || !fadeOutComplete);
					if (material.defines.FEATURE_FADE !== value) {
						this.fading += value === 1 ? 1 : -1;
						material.defines.FEATURE_FADE = value;
						material.needsUpdate = true;
					}
				}
			});
		}

		// initialize materials in the object
		prepareScene(scene) {
			scene.traverse(child => {
				if (child.material) {
					this.prepareMaterial(child.material);
				}
			});
		}

		// delete the object from the fade, reset the material data
		deleteScene(scene) {
			if (!scene) {
				return;
			}

			// revert the materials
			const fadeParams = this._fadeParams;
			scene.traverse(child => {
				const material = child.material;
				if (material) {
					fadeParams.delete(material);
					material.defines.FEATURE_FADE = false;
					material.needsUpdate = true;
				}
			});
		}

		// initialize the material
		prepareMaterial(material) {
			const fadeParams = this._fadeParams;
			if (fadeParams.has(material)) {
				return;
			}
			fadeParams.set(material, wrapFadeMaterial(material));
		}
	}

	class TilesFadePlugin {
		get fadeDuration() {
			return this._fadeManager.duration;
		}
		set fadeDuration(value) {
			this._fadeManager.duration = Number(value);
		}
		get fadingTiles() {
			return this._fadeManager.fadeCount;
		}
		constructor(options) {
			options = {
				maximumFadeOutTiles: 50,
				fadeRootTiles: false,
				fadeDuration: 250,
				...options
			};
			this.name = 'FADE_TILES_PLUGIN';
			this.priority = -2;
			this.tiles = null;
			this.batchedMesh = null;
			this._fadeManager = new FadeManager();
			this._fadeMaterialManager = new FadeMaterialManager();
			this._prevCameraTransforms = null;
			this._fadingOutCount = 0;
			this.maximumFadeOutTiles = options.maximumFadeOutTiles;
			this.fadeRootTiles = options.fadeRootTiles;
			this.fadeDuration = options.fadeDuration;
		}
		init(tiles) {
			// event callback initialization
			this._onLoadModel = ({
				scene
			}) => {
				// initialize all the scene materials to fade
				this._fadeMaterialManager.prepareScene(scene);
			};
			this._onDisposeModel = ({
				tile,
				scene
			}) => {
				// delete the fade info from the managers on disposal of model
				this._fadeManager.deleteObject(tile);
				this._fadeMaterialManager.deleteScene(scene);
			};
			this._onAddCamera = ({
				camera
			}) => {
				// track the camera transform
				this._prevCameraTransforms.set(camera, new t3d.Matrix4());
			};
			this._onDeleteCamera = ({
				camera
			}) => {
				// remove the camera transform
				this._prevCameraTransforms.delete(camera);
			};
			this._onTileVisibilityChange = ({
				tile,
				visible
			}) => {
				// this function gets fired _after_ all set visible callbacks including the batched meshes

				// revert the scene and fade to the initial state when toggling
				const scene = tile.cached.scene;
				if (scene) {
					scene.visible = true;
				}
			};
			this._onUpdateBefore = () => {
				onUpdateBefore.call(this);
			};
			this._onUpdateAfter = () => {
				onUpdateAfter.call(this);
			};
			tiles.addEventListener('load-model', this._onLoadModel);
			tiles.addEventListener('dispose-model', this._onDisposeModel);
			tiles.addEventListener('add-camera', this._onAddCamera);
			tiles.addEventListener('delete-camera', this._onDeleteCamera);
			tiles.addEventListener('update-before', this._onUpdateBefore);
			tiles.addEventListener('update-after', this._onUpdateAfter);
			tiles.addEventListener('tile-visibility-change', this._onTileVisibilityChange);

			// initialize fade manager
			const fadeManager = this._fadeManager;
			fadeManager.onFadeSetStart = () => {
				tiles.dispatchEvent({
					type: 'fade-start'
				});
				tiles.dispatchEvent({
					type: 'needs-render'
				});
			};
			fadeManager.onFadeSetComplete = () => {
				tiles.dispatchEvent({
					type: 'fade-end'
				});
				tiles.dispatchEvent({
					type: 'needs-render'
				});
			};
			fadeManager.onFadeComplete = (tile, visible) => {
				// mark the fade as finished and reset the fade parameters
				this._fadeMaterialManager.setFade(tile.cached.scene, 0, 0);
				if (!visible) {
					// now that the tile is hidden we can run the built-in setTileVisible function for the tile
					tiles.invokeOnePlugin(plugin => plugin !== this && plugin.setTileVisible && plugin.setTileVisible(tile, false));
					this._fadingOutCount--;
				}
			};

			// initialize the state based on what's already present
			const prevCameraTransforms = new Map();
			tiles.$cameras._cameras.forEach(camera => {
				prevCameraTransforms.set(camera, new t3d.Matrix4());
			});
			tiles.forEachLoadedModel((scene, tile) => {
				this._onLoadModel({
					scene
				});
			});
			this.tiles = tiles;
			this._fadeManager = fadeManager;
			this._prevCameraTransforms = prevCameraTransforms;
		}

		// callback for fading to prevent tiles from being removed until the fade effect has completed
		setTileVisible(tile, visible) {
			const fadeManager = this._fadeManager;

			// track the fade state
			const wasFading = fadeManager.isFading(tile);
			if (fadeManager.isFadingOut(tile)) {
				this._fadingOutCount--;
			}

			// trigger any necessary fades
			if (!visible) {
				this._fadingOutCount++;
				fadeManager.fadeOut(tile);
			} else {
				// if this is a root renderable tile and this is the first time rendering in
				// then pop it in
				const isRootRenderableTile = tile.__depthFromRenderedParent === 1;
				if (isRootRenderableTile) {
					if (tile[HAS_POPPED_IN] || this.fadeRootTiles) {
						this._fadeManager.fadeIn(tile);
					}
					tile[HAS_POPPED_IN] = true;
				} else {
					this._fadeManager.fadeIn(tile);
				}
			}

			// if a tile was already fading then it's already marked as visible and in the scene
			if (wasFading) {
				return true;
			}

			// cancel the visibility change trigger because we're fading and will call this after
			// fade completes.
			const isFading = this._fadeManager.isFading(tile);
			if (!visible && isFading) {
				return true;
			}
			return false;
		}
		dispose() {
			const tiles = this.tiles;
			this._fadeManager.completeAllFades();
			if (this.batchedMesh !== null) {
				this._onBatchedMeshDispose();
			}
			tiles.removeEventListener('load-model', this._onLoadModel);
			tiles.removeEventListener('dispose-model', this._onDisposeModel);
			tiles.removeEventListener('add-camera', this._onAddCamera);
			tiles.removeEventListener('delete-camera', this._onDeleteCamera);
			tiles.removeEventListener('update-before', this._onUpdateBefore);
			tiles.removeEventListener('update-after', this._onUpdateAfter);
			tiles.removeEventListener('tile-visibility-change', this._onTileVisibilityChange);
			tiles.forEachLoadedModel((scene, tile) => {
				this._fadeManager.deleteObject(tile);
				if (scene) {
					scene.visible = true; // TODO
				}
			});
		}
	}
	const HAS_POPPED_IN = Symbol('HAS_POPPED_IN');
	const _fromPos = new t3d.Vector3();
	const _toPos = new t3d.Vector3();
	const _fromQuat = new t3d.Quaternion();
	const _toQuat = new t3d.Quaternion();
	const _scale = new t3d.Vector3();
	function onUpdateBefore() {
		const fadeManager = this._fadeManager;
		const tiles = this.tiles;

		// store the tiles renderer state before the tiles update so we can check
		// whether fading started or stopped completely
		this._fadingBefore = fadeManager.fadeCount;
		this._displayActiveTiles = tiles.displayActiveTiles;

		// we need to display all active tiles in this case so we don't fade tiles in
		// when moving from off screen
		tiles.displayActiveTiles = true;
	}
	function onUpdateAfter() {
		const fadeManager = this._fadeManager;
		const fadeMaterialManager = this._fadeMaterialManager;
		const displayActiveTiles = this._displayActiveTiles;
		const fadingBefore = this._fadingBefore;
		const prevCameraTransforms = this._prevCameraTransforms;
		const {
			tiles,
			maximumFadeOutTiles
		} = this;
		const cameras = tiles.$cameras._cameras;

		// reset the active tiles flag
		tiles.displayActiveTiles = displayActiveTiles;

		// update fade step
		fadeManager.update();

		// fire an event
		const fadingAfter = fadeManager.fadeCount;
		if (fadingBefore !== 0 && fadingAfter !== 0) {
			tiles.dispatchEvent({
				type: 'fade-change'
			});
			tiles.dispatchEvent({
				type: 'needs-render'
			});
		}

		// update the visibility of tiles based on visibility since we must use
		// the active tiles for rendering fade
		if (!displayActiveTiles) {
			tiles.visibleTiles.forEach(t => {
				// if a tile is fading out then it may not be traversed and thus will not have
				// the frustum flag set correctly.
				const scene = t.cached.scene;
				if (scene) {
					scene.visible = t.__inFrustum;
				}
			});
		}
		if (maximumFadeOutTiles < this._fadingOutCount) {
			// determine whether all the rendering cameras are moving
			// quickly so we can adjust how tiles fade accordingly
			let isMovingFast = true;
			cameras.forEach(camera => {
				if (!prevCameraTransforms.has(camera)) {
					return;
				}
				const currMatrix = camera.worldMatrix;
				const prevMatrix = prevCameraTransforms.get(camera);
				currMatrix.decompose(_toPos, _toQuat, _scale);
				prevMatrix.decompose(_fromPos, _fromQuat, _scale);
				const angleTo = _toQuat.angleTo(_fromQuat);
				const positionTo = _toPos.distanceTo(_fromPos);

				// if rotation is moving > 0.25 radians per frame or position is moving > 0.1 units
				// then we are considering the camera to be moving too fast to notice a faster / abrupt fade
				isMovingFast = isMovingFast && (angleTo > 0.25 || positionTo > 0.1);
			});
			if (isMovingFast) {
				fadeManager.completeAllFades();
			}
		}

		// track the camera movement so we can use it for next frame
		cameras.forEach(camera => {
			prevCameraTransforms.get(camera).copy(camera.worldMatrix);
		});

		// update the fade state for each tile
		fadeManager.forEachObject((tile, {
			fadeIn,
			fadeOut
		}) => {
			// prevent faded tiles from being unloaded
			const scene = tile.cached.scene;
			const isFadingOut = fadeManager.isFadingOut(tile);
			tiles.markTileUsed(tile);
			if (scene) {
				fadeMaterialManager.setFade(scene, fadeIn, fadeOut);
				if (isFadingOut) {
					scene.visible = true;
				}
			}
		});
	}

	class GoogleAttributionsManager {
		constructor() {
			this.creditsCount = {};
		}
		_adjustAttributions(line, add) {
			const creditsCount = this.creditsCount;
			const tokens = line.split(/;/g);
			for (let i = 0, l = tokens.length; i < l; i++) {
				const t = tokens[i];
				if (!(t in creditsCount)) {
					creditsCount[t] = 0;
				}
				creditsCount[t] += add ? 1 : -1;
				if (creditsCount[t] <= 0) {
					delete creditsCount[t];
				}
			}
		}
		addAttributions(line) {
			this._adjustAttributions(line, true);
		}
		removeAttributions(line) {
			this._adjustAttributions(line, false);
		}
		toString() {
			// attribution guidelines: https://developers.google.com/maps/documentation/tile/create-renderer#display-attributions

			const sortedByCount = Object.entries(this.creditsCount).sort((a, b) => {
				const countA = a[1];
				const countB = b[1];
				return countB - countA; // Descending order
			});
			return sortedByCount.map(pair => pair[0]).join('; ');
		}
	}

	class GoogleCloudAuthPlugin {
		constructor({
			apiToken,
			autoRefreshToken = false,
			logoUrl = null,
			useRecommendedSettings = true
		}) {
			this.name = 'GOOGLE_CLOUD_AUTH_PLUGIN';
			this.priority = -Infinity;
			this.apiToken = apiToken;
			this.autoRefreshToken = autoRefreshToken;
			this.useRecommendedSettings = useRecommendedSettings;
			this.logoUrl = logoUrl;
			this.sessionToken = null;
			this.tiles = null;
			this._onLoadCallback = null;
			this._visibilityChangeCallback = null;
			this._tokenRefreshPromise = null;
			this._attributionsManager = new GoogleAttributionsManager();
			this._logoAttribution = {
				value: '',
				type: 'image',
				collapsible: false
			};
			this._attribution = {
				value: '',
				type: 'string',
				collapsible: true
			};
		}
		init(tiles) {
			if (tiles == null) {
				return;
			}

			// reset the tiles in case this plugin was removed and re-added
			tiles.resetFailedTiles();
			if (tiles.rootURL == null) {
				tiles.rootURL = 'https://tile.googleapis.com/v1/3dtiles/root.json';
			}
			if (this.useRecommendedSettings) {
				// This plugin changes below values to be more efficient for the photorealistic tiles
				tiles.parseQueue.maxJobs = 10;
				tiles.downloadQueue.maxJobs = 30;
				tiles.errorTarget = 20;
			}
			this.tiles = tiles;
			this._onLoadCallback = ({
				tileSet
			}) => {
				// the first tile set loaded will be the root
				this.sessionToken = getSessionToken(tileSet.root);

				// clear the callback once the root is loaded
				tiles.removeEventListener('load-tile-set', this._onLoadCallback);
			};
			this._visibilityChangeCallback = ({
				tile,
				visible
			}) => {
				// TODO
				// const copyright = tile.cached.metadata.asset.copyright || '';
				// if (visible) {
				// 	this._attributionsManager.addAttributions(copyright);
				// } else {
				// 	this._attributionsManager.removeAttributions(copyright);
				// }
			};
			tiles.addEventListener('load-tile-set', this._onLoadCallback);
			tiles.addEventListener('tile-visibility-change', this._visibilityChangeCallback);
		}
		getAttributions(target) {
			if (this.tiles.visibleTiles.size > 0) {
				if (this.logoUrl) {
					this._logoAttribution.value = this.logoUrl;
					target.push(this._logoAttribution);
				}
				this._attribution.value = this._attributionsManager.toString();
				target.push(this._attribution);
			}
		}
		preprocessURL(uri) {
			uri = new URL(uri);
			if (/^http/.test(uri.protocol)) {
				uri.searchParams.append('key', this.apiToken);
				if (this.sessionToken !== null) {
					uri.searchParams.append('session', this.sessionToken);
				}
			}
			return uri.toString();
		}
		dispose() {
			const {
				tiles
			} = this;
			tiles.removeEventListener('load-tile-set', this._onLoadCallback);
			tiles.removeEventListener('tile-visibility-change', this._visibilityChangeCallback);
		}
		async fetchData(uri, options) {
			// wait for the token to refresh if loading
			if (this._tokenRefreshPromise !== null) {
				await this._tokenRefreshPromise;
				uri = this.preprocessURL(uri);
			}
			const res = await fetch(uri, options);
			if (res.status >= 400 && res.status <= 499 && this.autoRefreshToken) {
				await this._refreshToken(options);
				return fetch(this.preprocessURL(uri), options);
			} else {
				return res;
			}
		}
		_refreshToken(options) {
			if (this._tokenRefreshPromise === null) {
				// refetch the root if the token has expired
				const rootURL = new URL(this.tiles.rootURL);
				rootURL.searchParams.append('key', this.apiToken);
				this._tokenRefreshPromise = fetch(rootURL, options).then(res => res.json()).then(res => {
					this.sessionToken = getSessionToken(res.root);
					this._tokenRefreshPromise = null;
				});

				// dispatch an error if we fail to refresh the token
				this._tokenRefreshPromise.catch(error => {
					this.tiles.dispatchEvent({
						type: 'load-error',
						tile: null,
						error,
						url: rootURL
					});
				});
			}
			return this._tokenRefreshPromise;
		}
	}
	function getSessionToken(root) {
		let sessionToken = null;
		traverseSet(root, tile => {
			if (tile.content && tile.content.uri) {
				const [, params] = tile.content.uri.split('?');
				sessionToken = new URLSearchParams(params).get('session');
				return true;
			}
			return false;
		});
		return sessionToken;
	}

	class GeometryUtils {
		/**
		 * @param {Geometry} geometry
		 */
		static computeNormals(geometry) {
			const index = geometry.index;
			const attributes = geometry.attributes;
			const positionAttribute = attributes.a_Position;
			if (positionAttribute === undefined) {
				return;
			}
			let normalAttribute = attributes.a_Normal;
			if (normalAttribute === undefined) {
				normalAttribute = new t3d.Attribute(new t3d.Buffer(new Float32Array(positionAttribute.buffer.count * 3), 3));
				geometry.addAttribute('a_Normal', normalAttribute);
			} else {
				for (let i = 0; i < normalAttribute.buffer.array.length; i++) {
					normalAttribute.buffer.array[i] = 0; // reset existing normals to zero
				}
				normalAttribute.buffer.version++;
			}
			const pA = new t3d.Vector3(),
				pB = new t3d.Vector3(),
				pC = new t3d.Vector3();
			const nA = new t3d.Vector3(),
				nB = new t3d.Vector3(),
				nC = new t3d.Vector3();
			const cb = new t3d.Vector3(),
				ab = new t3d.Vector3();
			if (index) {
				// indexed elements
				for (let i = 0, il = index.buffer.count; i < il; i += 3) {
					const vA = index.buffer.array[i + 0];
					const vB = index.buffer.array[i + 1];
					const vC = index.buffer.array[i + 2];
					pA.fromArray(positionAttribute.buffer.array, vA * 3);
					pB.fromArray(positionAttribute.buffer.array, vB * 3);
					pC.fromArray(positionAttribute.buffer.array, vC * 3);
					cb.subVectors(pC, pB);
					ab.subVectors(pA, pB);
					cb.cross(ab);
					nA.fromArray(normalAttribute.buffer.array, vA * 3);
					nB.fromArray(normalAttribute.buffer.array, vB * 3);
					nC.fromArray(normalAttribute.buffer.array, vC * 3);
					nA.add(cb);
					nB.add(cb);
					nC.add(cb);
					nA.toArray(normalAttribute.buffer.array, vA * 3);
					nB.toArray(normalAttribute.buffer.array, vB * 3);
					nC.toArray(normalAttribute.buffer.array, vC * 3);
				}
			} else {
				// non-indexed elements (unconnected triangle soup)
				for (let i = 0, il = positionAttribute.buffer.count * 3; i < il; i += 9) {
					pA.fromArray(positionAttribute.buffer.array, i + 0);
					pB.fromArray(positionAttribute.buffer.array, i + 3);
					pC.fromArray(positionAttribute.buffer.array, i + 6);
					cb.subVectors(pC, pB);
					ab.subVectors(pA, pB);
					cb.cross(ab);
					cb.toArray(normalAttribute.buffer.array, i + 0);
					cb.toArray(normalAttribute.buffer.array, i + 3);
					cb.toArray(normalAttribute.buffer.array, i + 6);
				}
			}
			this.normalizeNormals(geometry);
		}

		/**
		 * @param {Geometry} geometry
		 */
		static normalizeNormals(geometry) {
			const normals = geometry.attributes.a_Normal.buffer;
			for (let i = 0; i < normals.array.length; i += 3) {
				_vec3_1.fromArray(normals.array, i);
				_vec3_1.normalize();
				_vec3_1.toArray(normals.array, i);
			}
		}

		/**
		 * @param {Geometry} geometry
		 */
		static computeTangents(geometry) {
			const index = geometry.index;
			const attributes = geometry.attributes;

			// based on http://www.terathon.com/code/tangent.html
			// (per vertex tangents)

			if (index === null || attributes.a_Position === undefined || attributes.a_Normal === undefined || attributes.a_Uv === undefined) {
				console.warn('GeometryUtils: .computeTangents() failed. Missing required attributes (index, a_Position, a_Normal or a_Uv)');
				return;
			}
			const indices = index.buffer.array;
			const positions = attributes.a_Position.buffer.array;
			const normals = attributes.a_Normal.buffer.array;
			const uvs = attributes.a_Uv.buffer.array;
			const nVertices = positions.length / 3;
			if (!attributes.a_Tangent) {
				geometry.addAttribute('a_Tangent', new t3d.Attribute(new t3d.Buffer(new Float32Array(4 * nVertices), 4)));
			}
			const tangents = attributes.a_Tangent.buffer.array;
			const tan1 = [],
				tan2 = [];
			for (let i = 0; i < nVertices; i++) {
				tan1[i] = new t3d.Vector3();
				tan2[i] = new t3d.Vector3();
			}
			const vA = new t3d.Vector3(),
				vB = new t3d.Vector3(),
				vC = new t3d.Vector3(),
				uvA = new t3d.Vector2(),
				uvB = new t3d.Vector2(),
				uvC = new t3d.Vector2(),
				sdir = new t3d.Vector3(),
				tdir = new t3d.Vector3();
			function handleTriangle(a, b, c) {
				vA.fromArray(positions, a * 3);
				vB.fromArray(positions, b * 3);
				vC.fromArray(positions, c * 3);
				uvA.fromArray(uvs, a * 2);
				uvB.fromArray(uvs, b * 2);
				uvC.fromArray(uvs, c * 2);
				vB.sub(vA);
				vC.sub(vA);
				uvB.sub(uvA);
				uvC.sub(uvA);
				const r = 1.0 / (uvB.x * uvC.y - uvC.x * uvB.y);

				// silently ignore degenerate uv triangles having coincident or colinear vertices

				if (!isFinite(r)) return;
				sdir.copy(vB).multiplyScalar(uvC.y).addScaledVector(vC, -uvB.y).multiplyScalar(r);
				tdir.copy(vC).multiplyScalar(uvB.x).addScaledVector(vB, -uvC.x).multiplyScalar(r);
				tan1[a].add(sdir);
				tan1[b].add(sdir);
				tan1[c].add(sdir);
				tan2[a].add(tdir);
				tan2[b].add(tdir);
				tan2[c].add(tdir);
			}
			let groups = geometry.groups;
			if (groups.length === 0) {
				groups = [{
					start: 0,
					count: indices.length
				}];
			}
			for (let i = 0, il = groups.length; i < il; i++) {
				const group = groups[i];
				const start = group.start;
				const count = group.count;
				for (let j = start, jl = start + count; j < jl; j += 3) {
					handleTriangle(indices[j + 0], indices[j + 1], indices[j + 2]);
				}
			}
			const tmp = new t3d.Vector3(),
				tmp2 = new t3d.Vector3();
			const n = new t3d.Vector3(),
				n2 = new t3d.Vector3();
			function handleVertex(v) {
				n.fromArray(normals, v * 3);
				n2.copy(n);
				const t = tan1[v];

				// Gram-Schmidt orthogonalize

				tmp.copy(t);
				tmp.sub(n.multiplyScalar(n.dot(t))).normalize();

				// Calculate handedness

				tmp2.crossVectors(n2, t);
				const test = tmp2.dot(tan2[v]);
				const w = test < 0.0 ? -1 : 1.0;
				tangents[v * 4] = tmp.x;
				tangents[v * 4 + 1] = tmp.y;
				tangents[v * 4 + 2] = tmp.z;
				tangents[v * 4 + 3] = -w; // why negative?
			}
			for (let i = 0, il = groups.length; i < il; i++) {
				const group = groups[i];
				const start = group.start;
				const count = group.count;
				for (let j = start, jl = start + count; j < jl; j += 3) {
					handleVertex(indices[j + 0]);
					handleVertex(indices[j + 1]);
					handleVertex(indices[j + 2]);
				}
			}
		}

		/**
		 * @param {Array<Geometry>} geometries
		 * @param {boolean} useGroups
		 * @returns {Geometry}
		 */
		static mergeGeometries(geometries, useGroups = false) {
			const isIndexed = geometries[0].index !== null;
			const attributesUsed = new Set(Object.keys(geometries[0].attributes));
			const morphAttributesUsed = new Set(Object.keys(geometries[0].morphAttributes));
			const attributes = {};
			const morphAttributes = {};
			const mergedGeometry = new t3d.Geometry();
			let offset = 0;
			for (let i = 0; i < geometries.length; i++) {
				const geometry = geometries[i];

				// ensure that all geometries are indexed, or none

				if (isIndexed !== (geometry.index !== null)) {
					console.error('GeometryUtils: .mergeGeometries() failed with geometry at index ' + i + '. All geometries must have compatible attributes; make sure index attribute exists among all geometries, or in none of them.');
					return null;
				}

				// gather attributes, exit early if they're different

				for (const name in geometry.attributes) {
					if (!attributesUsed.has(name)) {
						console.error('GeometryUtils: .mergeGeometries() failed with geometry at index ' + i + '. All geometries must have compatible attributes; make sure "' + name + '" attribute exists among all geometries, or in none of them.');
						return null;
					}
					if (attributes[name] === undefined) attributes[name] = [];
					attributes[name].push(geometry.attributes[name]);
				}

				// gather morph attributes, exit early they're different

				for (const name in geometry.morphAttributes) {
					if (!morphAttributesUsed.has(name)) {
						console.error('GeometryUtils: .mergeGeometries() failed with geometry at index ' + i + '. .morphAttributes must be consistent throughout all geometries.');
						return null;
					}
					if (morphAttributes[name] === undefined) morphAttributes[name] = [];
					morphAttributes[name].push(geometry.morphAttributes[name]);
				}
				if (useGroups) {
					let count;
					if (isIndexed) {
						count = geometry.index.buffer.count;
					} else if (geometry.attributes.a_Position !== undefined) {
						count = geometry.attributes.a_Position.buffer.count;
					} else {
						console.error('GeometryUtils: .mergeGeometries() failed with geometry at index ' + i + '. The geometry must have either an index or an a_Position attribute');
						return null;
					}
					mergedGeometry.addGroup(offset, count, i);
					offset += count;
				}
			}

			// merge indices

			if (isIndexed) {
				let indexOffset = 0;
				const mergedIndex = [];
				for (let i = 0; i < geometries.length; i++) {
					const index = geometries[i].index;
					for (let j = 0; j < index.buffer.count; j++) {
						mergedIndex.push(index.buffer.array[j] + indexOffset);
					}
					indexOffset += geometries[i].attributes.a_Position.buffer.count;
				}
				mergedGeometry.setIndex(mergedIndex);
			}

			// merge attributes

			for (const name in attributes) {
				const mergedAttribute = this.mergeAttributes(attributes[name]);
				if (!mergedAttribute) {
					console.error('GeometryUtils: .mergeGeometries() failed while trying to merge the ' + name + ' attribute.');
					return null;
				}
				mergedGeometry.addAttribute(name, mergedAttribute);
			}

			// merge morph attributes

			for (const name in morphAttributes) {
				const numMorphTargets = morphAttributes[name][0].length;
				if (numMorphTargets === 0) break;
				mergedGeometry.morphAttributes = mergedGeometry.morphAttributes || {};
				mergedGeometry.morphAttributes[name] = [];
				for (let i = 0; i < numMorphTargets; i++) {
					const morphAttributesToMerge = [];
					for (let j = 0; j < morphAttributes[name].length; j++) {
						morphAttributesToMerge.push(morphAttributes[name][j][i]);
					}
					const mergedMorphAttribute = this.mergeAttributes(morphAttributesToMerge);
					if (!mergedMorphAttribute) {
						console.error('GeometryUtils: .mergeGeometries() failed while trying to merge the ' + name + ' morphAttribute.');
						return null;
					}
					mergedGeometry.morphAttributes[name].push(mergedMorphAttribute);
				}
			}
			return mergedGeometry;
		}

		/**
		 * @param {Array<Attribute>} attributes
		 * @returns {Attribute}
		 */
		static mergeAttributes(attributes) {
			let TypedArray;
			let size;
			let normalized;
			let arrayLength = 0;
			for (let i = 0; i < attributes.length; i++) {
				const attribute = attributes[i];
				if (attribute.buffer.stride !== attribute.size) {
					console.error('GeometryUtils: .mergeAttributes() failed. Interleaved buffer attributes are not supported.');
					return null;
				}
				if (TypedArray === undefined) TypedArray = attribute.buffer.array.constructor;
				if (TypedArray !== attribute.buffer.array.constructor) {
					console.error('GeometryUtils: .mergeAttributes() failed. Buffer.array must be of consistent array types across matching attributes.');
					return null;
				}
				if (size === undefined) size = attribute.size;
				if (size !== attribute.size) {
					console.error('GeometryUtils: .mergeAttributes() failed. Attribute.size must be consistent across matching attributes.');
					return null;
				}
				if (normalized === undefined) normalized = attribute.normalized;
				if (normalized !== attribute.normalized) {
					console.error('GeometryUtils: .mergeAttributes() failed. Attribute.normalized must be consistent across matching attributes.');
					return null;
				}
				arrayLength += attribute.buffer.array.length;
			}
			const array = new TypedArray(arrayLength);
			let offset = 0;
			for (let i = 0; i < attributes.length; i++) {
				array.set(attributes[i].buffer.array, offset);
				offset += attributes[i].buffer.array.length;
			}
			return new t3d.Attribute(new t3d.Buffer(array, size), size, 0, normalized);
		}

		/**
		 * @param {Geometry} geometry
		 * @param {Matrix4} matrix
		 * @param {boolean} updateBoundings
		 * @returns {Geometry}
		 */
		static applyMatrix4(geometry, matrix, updateBoundings) {
			let array, count, offset;
			const position = geometry.attributes['a_Position'];
			if (position !== undefined) {
				array = position.buffer.array;
				count = position.buffer.count;
				offset = position.offset;
				for (let i = 0; i < count; i++) {
					_vec3_1.fromArray(array, i * 3 + offset);
					_vec3_1.applyMatrix4(matrix);
					_vec3_1.toArray(array, i * 3 + offset);
				}
				position.buffer.version++;
			}
			const normal = geometry.attributes['a_Normal'];
			if (normal !== undefined) {
				array = normal.buffer.array;
				count = normal.buffer.count;
				offset = normal.offset;
				const normalMatrix = _mat3_1.setFromMatrix4(matrix).invert().transpose();
				for (let i = 0; i < count; i++) {
					_vec3_1.fromArray(array, i * 3 + offset);
					_vec3_1.applyMatrix3(normalMatrix).normalize();
					_vec3_1.toArray(array, i * 3 + offset);
				}
				normal.buffer.version++;
			}
			const tangent = geometry.attributes['a_Tangent'];
			if (tangent !== undefined) {
				array = tangent.buffer.array;
				count = tangent.buffer.count;
				offset = tangent.offset;
				for (let i = 0; i < count; i++) {
					_vec3_1.fromArray(array, i * 3 + offset);
					_vec3_1.transformDirection(matrix);
					_vec3_1.toArray(array, i * 3 + offset);
				}
				tangent.buffer.version++;
			}
			if (geometry.boundingBox !== null && updateBoundings) {
				geometry.computeBoundingBox();
			}
			if (geometry.boundingSphere !== null && updateBoundings) {
				geometry.computeBoundingSphere();
			}
			return geometry;
		}

		/**
		 * @param {Geometry} geometry
		 * @returns {Attribute}
		 */
		static getWireframeAttribute(geometry) {
			const indices = [];
			const geometryIndex = geometry.index;
			const geometryPosition = geometry.attributes.a_Position;
			if (!geometryPosition) {
				console.error('GeometryUtils: .getWireframeAttribute() failed. The geometry must have an a_Position attribute');
				return null;
			}
			if (geometryIndex !== null) {
				const array = geometryIndex.buffer.array;
				for (let i = 0, l = array.length; i < l; i += 3) {
					const a = array[i + 0];
					const b = array[i + 1];
					const c = array[i + 2];
					indices.push(a, b, b, c, c, a);
				}
			} else {
				const array = geometryPosition.buffer.array;
				for (let i = 0, l = array.length / 3 - 1; i < l; i += 3) {
					const a = i + 0;
					const b = i + 1;
					const c = i + 2;
					indices.push(a, b, b, c, c, a);
				}
			}
			return new t3d.Attribute(new t3d.Buffer(geometryPosition.buffer.array.length / 3 > 65536 ? new Uint32Array(indices) : new Uint16Array(indices), 1));
		}

		// deprecated since v0.2.0, add warning since v0.3.0
		static mergeBufferAttributes(attributes) {
			console.warn('GeometryUtils: mergeBufferAttributes() has been renamed to mergeAttributes().');
			return this.mergeAttributes(attributes);
		}
	}
	const _vec3_1 = new t3d.Vector3();
	const _mat3_1 = new t3d.Matrix3();

	const TILE_X$1 = Symbol('TILE_X');
	const TILE_Y$1 = Symbol('TILE_Y');
	const TILE_LEVEL$1 = Symbol('TILE_LEVEL');

	// Base class for supporting tiled images with a consistent size / resolution per tile
	class ImageFormatPlugin {
		get tiling() {
			return this.imageSource.tiling;
		}
		constructor(options = {}) {
			const {
				pixelSize = 0.01,
				center = false,
				useRecommendedSettings = true,
				imageSource = null
			} = options;
			this.priority = -10;
			this.tiles = null;

			// tiling scheme
			this.imageSource = imageSource;

			// options
			this.pixelSize = pixelSize;
			this.center = center;
			this.useRecommendedSettings = useRecommendedSettings;
		}

		// Plugin functions
		init(tiles) {
			if (this.useRecommendedSettings) {
				tiles.errorTarget = 1;
				// TODO: apply skip traversal settings here once supported, as well, for faster loading
			}
			this.tiles = tiles;
			this.imageSource.fetchOptions = tiles.fetchOptions;
			this.imageSource.fetchData = (url, options) => {
				tiles.invokeAllPlugins(plugin => url = plugin.preprocessURL ? plugin.preprocessURL(url, null) : url);
				return tiles.invokeOnePlugin(plugin => plugin !== this && plugin.fetchData && plugin.fetchData(url, options));
			};
		}
		async loadRootTileSet() {
			const {
				tiles,
				imageSource
			} = this;
			let url = tiles.rootURL;
			tiles.invokeAllPlugins(plugin => url = plugin.preprocessURL ? plugin.preprocessURL(url, null) : url);
			await imageSource.init(url);
			return this.getTileset(url);
		}
		async parseToMesh(buffer, tile, extension, uri, abortSignal) {
			// Construct texture
			const tx = tile[TILE_X$1];
			const ty = tile[TILE_Y$1];
			const level = tile[TILE_LEVEL$1];
			const texture = await this.imageSource.processBufferToTexture(buffer);
			this.imageSource.setData(tx, ty, level, texture);

			// Construct mesh
			let sx = 1,
				sy = 1;
			let x = 0,
				y = 0,
				z = 0;
			const boundingBox = tile.boundingVolume.box;
			if (boundingBox) {
				[x, y, z] = boundingBox;
				sx = boundingBox[3];
				sy = boundingBox[7];
			}

			// adjust the geometry transform itself rather than the mesh because it reduces the artifact errors
			// when using batched mesh rendering.
			const mesh = new t3d.Mesh(new t3d.PlaneGeometry(2 * sx, 2 * sy), new t3d.BasicMaterial());
			const rotation = new t3d.Euler(Math.PI / 2, 0, 0, 'XYZ');
			GeometryUtils.applyMatrix4(mesh.geometry, new t3d.Matrix4().makeRotationFromEuler(rotation), true);
			mesh.material.diffuseMap = texture;
			mesh.material.transparent = true;
			mesh.position.set(x, y, z);
			return mesh;
		}
		preprocessNode(tile) {
			// generate children
			const {
				tiling
			} = this;
			const maxLevel = tiling.maxLevel;
			const level = tile[TILE_LEVEL$1];
			if (level < maxLevel && tile.parent !== null) {
				this.expandChildren(tile);
			}
		}
		disposeTile(tile) {
			const tx = tile[TILE_X$1];
			const ty = tile[TILE_Y$1];
			const level = tile[TILE_LEVEL$1];
			this.imageSource.release(tx, ty, level);
		}

		// Local functions
		getTileset(baseUrl) {
			const {
				tiling,
				tiles
			} = this;
			const minLevel = tiling.minLevel;
			const {
				tileCountX,
				tileCountY
			} = tiling.getLevel(minLevel);

			// generate all children for the root
			const children = [];
			for (let x = 0; x < tileCountX; x++) {
				for (let y = 0; y < tileCountY; y++) {
					const child = this.createChild(x, y, minLevel);
					if (child !== null) {
						children.push(child);
					}
				}
			}

			// generate tile set
			const tileset = {
				asset: {
					version: '1.1'
				},
				geometricError: 1e5,
				root: {
					refine: 'REPLACE',
					geometricError: 1e5,
					boundingVolume: this.createBoundingVolume(0, 0, -1),
					children,
					[TILE_LEVEL$1]: -1,
					[TILE_X$1]: 0,
					[TILE_Y$1]: 0
				}
			};
			tiles.preprocessTileSet(tileset, baseUrl);
			return tileset;
		}
		getUrl(x, y, level) {
			return this.imageSource.getUrl(x, y, level);
		}
		createBoundingVolume(x, y, level) {
			const {
				center,
				pixelSize,
				tiling
			} = this;
			const {
				pixelWidth,
				pixelHeight
			} = tiling.getLevel(tiling.maxLevel);

			// calculate the world space bounds position from the range
			const [minX, minY, maxX, maxY] = level === -1 ? tiling.getFullBounds(true) : tiling.getTileBounds(x, y, level, true);
			let extentsX = (maxX - minX) / 2;
			let extentsY = (maxY - minY) / 2;
			let centerX = minX + extentsX;
			let centerY = minY + extentsY;
			if (center) {
				centerX -= 0.5;
				centerY -= 0.5;
			}

			// scale the fields
			centerX *= pixelWidth * pixelSize;
			extentsX *= pixelWidth * pixelSize;
			centerY *= pixelHeight * pixelSize;
			extentsY *= pixelHeight * pixelSize;

			// return bounding box
			return {
				box: [
				// center
				centerX, centerY, 0,
				// x, y, z half vectors
				extentsX, 0.0, 0.0, 0.0, extentsY, 0.0, 0.0, 0.0, 0.0]
			};
		}
		createChild(x, y, level) {
			const {
				pixelSize,
				tiling
			} = this;
			if (!tiling.getTileExists(x, y, level)) {
				return null;
			}

			// the scale ration of the image at this level
			const {
				pixelWidth,
				pixelHeight
			} = tiling.getLevel(tiling.maxLevel);
			const {
				pixelWidth: levelWidth,
				pixelHeight: levelHeight
			} = tiling.getLevel(level);
			const geometricError = pixelSize * (Math.max(pixelWidth / levelWidth, pixelHeight / levelHeight) - 1);

			// Generate the node
			return {
				refine: 'REPLACE',
				geometricError: geometricError,
				boundingVolume: this.createBoundingVolume(x, y, level),
				content: {
					uri: this.getUrl(x, y, level)
				},
				children: [],
				// save the tile params so we can expand later
				[TILE_X$1]: x,
				[TILE_Y$1]: y,
				[TILE_LEVEL$1]: level
			};
		}
		expandChildren(tile) {
			const level = tile[TILE_LEVEL$1];
			const x = tile[TILE_X$1];
			const y = tile[TILE_Y$1];
			for (let cx = 0; cx < 2; cx++) {
				for (let cy = 0; cy < 2; cy++) {
					const child = this.createChild(2 * x + cx, 2 * y + cy, level + 1);
					if (child) {
						tile.children.push(child);
					}
				}
			}
		}
	}

	const _v0 = /* @__PURE__ */new t3d.Vector3();
	const _v1 = /* @__PURE__ */new t3d.Vector3();
	function getCartographicToMeterDerivative(ellipsoid, lat, lon) {
		const EPS = 1e-5;
		const lonp = lon + EPS;
		let latp = lat + EPS;
		if (Math.abs(latp) > Math.PI / 2) {
			latp = latp - EPS;
		}
		ellipsoid.getCartographicToPosition(lat, lon, 0, _v0);
		ellipsoid.getCartographicToPosition(latp, lon, 0, _v1);
		const dy = _v0.distanceTo(_v1) / EPS;
		ellipsoid.getCartographicToPosition(lat, lonp, 0, _v1);
		const dx = _v0.distanceTo(_v1) / EPS;
		return [dx, dy];
	}

	const MIN_LON_VERTS = 30;
	const MIN_LAT_VERTS = 15;
	const _pos$2 = /* @__PURE__ */new t3d.Vector3();
	const _norm$1 = /* @__PURE__ */new t3d.Vector3();
	const _uv = /* @__PURE__ */new t3d.Vector2();
	const _sphere$1 = /* @__PURE__ */new t3d.Sphere();
	class EllipsoidProjectionTilesPlugin extends ImageFormatPlugin {
		get projection() {
			return this.tiling.projection;
		}
		constructor(options = {}) {
			const {
				shape = 'planar',
				endCaps = true,
				...rest
			} = options;
			super(rest);

			// options
			this.shape = shape;
			this.endCaps = endCaps;
		}

		// override the parse to mesh logic to support a region mesh
		async parseToMesh(buffer, tile, ...args) {
			const mesh = await super.parseToMesh(buffer, tile, ...args);

			// if displaying the tiles as an ellipsoid
			const {
				shape,
				projection,
				tiles,
				tiling
			} = this;
			if (shape === 'ellipsoid') {
				const ellipsoid = tiles.ellipsoid;
				const level = tile[TILE_LEVEL$1];
				const x = tile[TILE_X$1];
				const y = tile[TILE_Y$1];
				const [minU, minV, maxU, maxV] = tiling.getTileBounds(x, y, level, true);
				const [west, south, east, north] = tile.boundingVolume.region;

				// new geometry
				// default to a minimum number of vertices per degree on each axis
				const latVerts = Math.ceil((north - south) * t3d.MathUtils.RAD2DEG * 0.25);
				const lonVerts = Math.ceil((east - west) * t3d.MathUtils.RAD2DEG * 0.25);
				const yVerts = Math.max(MIN_LAT_VERTS, latVerts);
				const xVerts = Math.max(MIN_LON_VERTS, lonVerts);
				const geometry = new t3d.PlaneGeometry(1, 1, xVerts, yVerts);

				// adjust the geometry to position it at the region
				const {
					a_Position: position,
					a_Normal: normal,
					a_Uv: uv
				} = geometry.attributes;
				const vertCount = position.buffer.count;
				tile.cached.boundingVolume.getBoundingSphere(_sphere$1);
				for (let i = 0; i < vertCount; i++) {
					_pos$2.fromArray(position.buffer.array, i * 3);
					_norm$1.fromArray(normal.buffer.array, i * 3);
					_uv.fromArray(uv.buffer.array, i * 2);
					const lon = t3d.MathUtils.mapLinear(_uv.x, 0, 1, west, east);
					let lat = t3d.MathUtils.mapLinear(_uv.y, 0, 1, south, north);
					if (projection.isMercator && _uv.y !== 0 && _uv.y !== 1) {
						// ensure we have an edge loop positioned at the mercator limit
						// to avoid UV distortion as much as possible at low LoDs
						const latLimit = projection.convertProjectionToLatitude(1);
						const vStep = 1 / yVerts;
						const prevLat = t3d.MathUtils.mapLinear(_uv.y - vStep, 0, 1, south, north);
						const nextLat = t3d.MathUtils.mapLinear(_uv.y + vStep, 0, 1, south, north);
						if (lat > latLimit && prevLat < latLimit) {
							lat = latLimit;
						}
						if (lat < -latLimit && nextLat > -latLimit) {
							lat = -latLimit;
						}
					}
					ellipsoid.getCartographicToPosition(lat, lon, 0, _pos$2).sub(_sphere$1.center);
					ellipsoid.getCartographicToNormal(lat, lon, _norm$1);

					// update the geometry
					const u = t3d.MathUtils.mapLinear(projection.convertLongitudeToProjection(lon), minU, maxU, 0, 1);
					const v = t3d.MathUtils.mapLinear(projection.convertLatitudeToProjection(lat), minV, maxV, 0, 1);
					uv.buffer.array[i * 2] = u;
					uv.buffer.array[i * 2 + 1] = v;
					_pos$2.toArray(position.buffer.array, i * 3);
					_norm$1.toArray(normal.buffer.array, i * 3);
				}
				mesh.geometry = geometry;
				mesh.position.copy(_sphere$1.center);
			}
			return mesh;
		}
		createBoundingVolume(x, y, level) {
			if (this.shape === 'ellipsoid') {
				const {
					tiling,
					endCaps
				} = this;
				const isRoot = level === -1;
				const normalizedBounds = isRoot ? tiling.getFullBounds(true) : tiling.getTileBounds(x, y, level, true);
				const cartBounds = isRoot ? tiling.getFullBounds() : tiling.getTileBounds(x, y, level);
				if (endCaps) {
					// if the north side is at the edge
					if (normalizedBounds[3] === 1) {
						cartBounds[3] = Math.PI / 2;
					}

					// if the south side is at the edge
					if (normalizedBounds[1] === 0) {
						cartBounds[1] = -Math.PI / 2;
					}
				}
				return {
					region: [...cartBounds, -1, 1]
				};
			} else {
				return super.createBoundingVolume(x, y, level);
			}
		}
		preprocessNode(tile, ...rest) {
			super.preprocessNode(tile, rest);
			const {
				shape,
				projection,
				tiling
			} = this;
			if (shape === 'ellipsoid') {
				const level = tile[TILE_LEVEL$1];
				const x = tile[TILE_X$1];
				const y = tile[TILE_Y$1];

				// if this is the root node then skip calculating the geometric error
				if (level === -1) {
					tile.geometricError = 1e50;
					return parent;
				}
				const [minU, minV, maxU, maxV] = tiling.getTileBounds(x, y, level, true);
				const {
					tilePixelWidth,
					tilePixelHeight
				} = tiling.getLevel(level);
				const {
					pixelWidth,
					pixelHeight
				} = tiling.getLevel(tiling.maxLevel);

				// one pixel width in uv space
				const tileUWidth = (maxU - minU) / tilePixelWidth;
				const tileVWidth = (maxV - minV) / tilePixelHeight;
				const rootUWidth = 1 / pixelWidth;
				const rootVWidth = 1 / pixelHeight;

				// calculate the region ranges
				const [, south, east, north] = tiling.getTileBounds(x, y, level);

				// calculate the changes in lat / lon at the given point
				// find the most bowed point of the latitude range since the amount that latitude changes is
				// dependent on the Y value of the image
				const midLat = south > 0 !== north > 0 ? 0 : Math.min(Math.abs(south), Math.abs(north));
				const midV = projection.convertLatitudeToProjection(midLat);
				const lonFactor = projection.getLongitudeDerivativeAtValue(minU);
				const latFactor = projection.getLatitudeDerivativeAtValue(midV);

				// TODO: is this correct?

				// calculate the size of a pixel on the surface
				const [xDeriv, yDeriv] = getCartographicToMeterDerivative(this.tiles.ellipsoid, midLat, east);
				const tilePixelWidth2 = Math.max(tileUWidth * lonFactor * xDeriv, tileVWidth * latFactor * yDeriv);
				const rootPixelWidth = Math.max(rootUWidth * lonFactor * xDeriv, rootVWidth * latFactor * yDeriv);
				tile.geometricError = tilePixelWidth2 - rootPixelWidth;

				// if this is the root then keep the geometric error high
				if (tile.parent === null) {
					tile.geometricError = 1e50;
				}
			}
			return tile;
		}
	}

	// Class for storing and querying a certain projection scheme for an image and converting
	// between the [0, 1] image range to cartographic longitude / latitude values.
	class ProjectionScheme {
		get isMercator() {
			return this.scheme === 'EPSG:3857';
		}
		constructor(scheme = 'EPSG:4326') {
			this.scheme = scheme;
			this.tileCountX = 1;
			this.tileCountY = 1;
			this.setScheme(scheme);
		}
		setScheme(scheme) {
			this.scheme = scheme;
			switch (scheme) {
				// equirect
				case 'EPSG:4326':
					this.tileCountX = 2;
					this.tileCountY = 1;
					break;

				// mercator
				case 'EPSG:3857':
					this.tileCountX = 1;
					this.tileCountY = 1;
					break;
				default:
					throw new Error();
			}
		}
		convertProjectionToLatitude(v) {
			if (this.isMercator) {
				// https://gis.stackexchange.com/questions/447421/convert-a-point-on-a-flat-2d-web-mercator-map-image-to-a-coordinate
				const ratio = t3d.MathUtils.mapLinear(v, 0, 1, -1, 1);
				return 2 * Math.atan(Math.exp(ratio * Math.PI)) - Math.PI / 2;
			} else {
				return t3d.MathUtils.mapLinear(v, 0, 1, -Math.PI / 2, Math.PI / 2);
			}
		}
		convertProjectionToLongitude(v) {
			return t3d.MathUtils.mapLinear(v, 0, 1, -Math.PI, Math.PI);
		}
		convertLatitudeToProjection(lat) {
			if (this.isMercator) {
				// https://stackoverflow.com/questions/14329691/convert-latitude-longitude-point-to-a-pixels-x-y-on-mercator-projection
				const mercatorN = Math.log(Math.tan(Math.PI / 4 + lat / 2));
				return 1 / 2 + 1 * mercatorN / (2 * Math.PI);
			} else {
				return t3d.MathUtils.mapLinear(lat, -Math.PI / 2, Math.PI / 2, 0, 1);
			}
		}
		convertLongitudeToProjection(lon) {
			return (lon + Math.PI) / (2 * Math.PI);
		}
		getLongitudeDerivativeAtValue(value) {
			return 2 * Math.PI;
		}
		getLatitudeDerivativeAtValue(value) {
			const EPS = 1e-5;
			let yp = value - EPS;
			if (yp < 0) {
				yp = value + EPS;
			}
			if (this.isMercator) {
				// TODO: why is this 2 * Math.PI rather than Math.PI?
				return Math.abs(this.convertProjectionToLatitude(value) - this.convertProjectionToLatitude(yp)) / EPS;
			} else {
				return Math.PI;
			}
		}
		getBounds() {
			return [this.convertProjectionToLongitude(0), this.convertProjectionToLatitude(0), this.convertProjectionToLongitude(1), this.convertProjectionToLatitude(1)];
		}
	}

	function hash(...args) {
		return args.join('_');
	}

	// class for retrieving and locking data being requested
	// "fetchItem" and "disposeItem" should be implemented
	class DataCache {
		constructor() {
			this.cache = {};
			this.count = 0;
			this.cachedBytes = 0;
		}

		// overridable
		fetchItem() {}
		disposeItem() {}
		getMemoryUsage(item) {
			return 0;
		}

		// sets the data in the cache explicitly without need to load
		setData(...args) {
			const {
				cache
			} = this;
			const data = args.pop();
			const key = hash(...args);
			if (key in cache) {
				throw new Error(`DataCache: "${key}" is already present.`);
			} else {
				this.cache[key] = {
					abortController: new AbortController(),
					result: data,
					count: 1,
					bytes: this.getMemoryUsage(data)
				};
				this.count++;
				this.cachedBytes += this.cache[key].bytes;
			}
			return data;
		}

		// fetches the associated data if it doesn't exist and increments the lock counter
		lock(...args) {
			const {
				cache
			} = this;
			const key = hash(...args);
			if (key in cache) {
				cache[key].count++;
			} else {
				const abortController = new AbortController();
				const info = {
					abortController,
					result: null,
					count: 1,
					bytes: 0
				};
				info.result = this.fetchItem(...args, abortController.signal).then(res => {
					info.result = res;
					info.bytes = this.getMemoryUsage(res);
					this.cachedBytes += info.bytes;
					return res;
				});
				this.cache[key] = info;
				this.count++;
			}
			return cache[key].result;
		}

		// decrements the lock counter for the item and deletes the item if it has reached zero
		release(...args) {
			const key = hash(...args);
			this.releaseViaFullKey(key);
		}

		// get the loaded item
		get(...args) {
			const {
				cache
			} = this;
			const key = hash(...args);
			if (key in cache) {
				return cache[key].result;
			} else {
				return null;
			}
		}

		// dispose all items
		dispose() {
			const {
				cache
			} = this;
			for (const key in cache) {
				const {
					abortController
				} = cache[key];
				abortController.abort();
				this.releaseViaFullKey(key, true);
			}
			this.cache = {};
		}

		// releases an item with an optional force flag
		releaseViaFullKey(key, force = false) {
			const {
				cache
			} = this;
			if (key in cache) {
				// decrement the lock
				const info = cache[key];
				info.count--;

				// if the item is no longer being used
				if (info.count === 0 || force) {
					const disposeCallback = () => {
						// if the object isn't in the cache anymore then exit early
						if (cache[key] !== info) {
							return;
						}

						// abort any loads
						const {
							result,
							abortController
						} = info;
						abortController.abort();

						// dispose of the object even if it still is in progress
						if (result instanceof Promise) {
							// "disposeItem" will throw potentially if fetch, etc are cancelled using the abort signal
							result.then(item => this.disposeItem(item)).catch(() => {});
						} else {
							this.disposeItem(result);
						}
						delete cache[key];
						this.count--;
						this.cachedBytes -= info.bytes;
					};
					if (force) {
						// if we're forcing disposal then dispose immediately
						disposeCallback();
					} else {
						// queue for disposal in a frame here - we need to make sure we're not disposing of something twice
						// this can get called multiple times in a row to increment then decrement again.
						queueMicrotask(() => {
							if (info.count === 0) {
								disposeCallback();
							}
						});
					}
				}
				return true;
			} else {
				throw new Error('DataCache: Attempting to release key that does not exist');
			}
		}
	}

	// Class for storing and querying a tiling scheme including a bounds, origin, and negative tile indices.
	// Assumes that tiles are split into four child tiles at each level.
	function clamp(x, min, max) {
		return Math.min(Math.max(x, min), max);
	}
	class TilingScheme {
		get levelCount() {
			return this._levels.length;
		}
		get maxLevel() {
			return this.levelCount - 1;
		}
		get minLevel() {
			const levels = this._levels;
			for (let i = 0; i < levels.length; i++) {
				if (levels[i] !== null) {
					return i;
				}
			}
			return -1;
		}

		// prioritize user-set bounds over projection bounds if present
		get rootBounds() {
			return this._rootBounds ?? this.projection?.getBounds() ?? [0, 0, 1, 1];
		}
		get rootOrigin() {
			const bounds = this.rootBounds;
			return this._rootOrigin ?? [bounds[0], bounds[1]];
		}
		constructor() {
			this.flipY = false;
			this.pixelOverlap = 0;

			// The origin and bounds
			this._rootBounds = null;
			this._rootOrigin = null;
			this.projection = null;
			this._levels = [];
		}

		// build the zoom levels
		setLevel(level, options = {}) {
			const levels = this._levels;
			while (levels.length < level) {
				levels.push(null);
			}
			const {
				tilePixelWidth = 256,
				tilePixelHeight = 256,
				tileCountX = 2 ** level,
				tileCountY = 2 ** level
			} = options;
			const {
				pixelWidth = tilePixelWidth * tileCountX,
				pixelHeight = tilePixelHeight * tileCountY
			} = options;
			levels[level] = {
				tilePixelWidth,
				tilePixelHeight,
				pixelWidth,
				pixelHeight,
				tileCountX,
				tileCountY
			};
		}
		generateLevels(levels, rootTileX, rootTileY, options = {}) {
			const {
				minLevel = 0,
				tilePixelWidth = 256,
				tilePixelHeight = 256
			} = options;
			const maxLevel = levels - 1;
			const {
				pixelWidth = tilePixelWidth * rootTileX * 2 ** maxLevel,
				pixelHeight = tilePixelHeight * rootTileY * 2 ** maxLevel
			} = options;
			for (let level = minLevel; level < levels; level++) {
				const invLevel = levels - level - 1;
				const levelPixelWidth = Math.ceil(pixelWidth * 2 ** -invLevel);
				const levelPixelHeight = Math.ceil(pixelHeight * 2 ** -invLevel);
				const tileCountX = Math.ceil(levelPixelWidth / tilePixelWidth);
				const tileCountY = Math.ceil(levelPixelHeight / tilePixelHeight);
				this.setLevel(level, {
					tilePixelWidth,
					tilePixelHeight,
					pixelWidth: levelPixelWidth,
					pixelHeight: levelPixelHeight,
					tileCountX,
					tileCountY
				});
			}
		}
		getLevel(level) {
			return this._levels[level];
		}

		// bounds setters
		setOrigin(x, y) {
			this._rootOrigin = [x, y];
		}
		setBounds(minX, minY, maxX, maxY) {
			this._rootBounds = [minX, minY, maxX, maxY];
		}
		setProjection(projection) {
			this.projection = projection;
		}

		// query functions
		getTileAtPoint(bx, by, level, normalized = false, clampTiles = true) {
			const {
				projection,
				flipY
			} = this;
			const {
				tileCountX,
				tileCountY
			} = this.getLevel(level);
			const xStride = 1 / tileCountX;
			const yStride = 1 / tileCountY;
			if (projection && !normalized) {
				bx = projection.convertLongitudeToProjection(bx);
				by = projection.convertLatitudeToProjection(by);
			}
			if (clampTiles) {
				bx = clamp(bx, 0, 1);
				by = clamp(by, 0, 1);
			}
			let tx = Math.floor(bx / xStride);
			let ty = Math.floor(by / yStride);
			if (flipY) {
				ty = tileCountY - 1 - ty;
			}
			if (clampTiles) {
				tx = clamp(tx, 0, tileCountX - 1);
				ty = clamp(ty, 0, tileCountY - 1);
			}
			return [tx, ty];
		}
		getTilesInRange(minX, minY, maxX, maxY, level, normalized = false, clampTiles = true) {
			const minTile = this.getTileAtPoint(minX, minY, level, normalized, clampTiles);
			const maxTile = this.getTileAtPoint(maxX, maxY, level, normalized, clampTiles);
			if (this.flipY) {
				[minTile[1], maxTile[1]] = [maxTile[1], minTile[1]];
			}
			return [...minTile, ...maxTile];
		}
		getTileExists(x, y, level, LOG) {
			const [rminx, rminy, rmaxx, rmaxy] = this.rootBounds;
			const [tminx, tminy, tmaxx, tmaxy] = this.getTileBounds(x, y, level, LOG);
			const isDegenerate = tminx >= tmaxx || tminy >= tmaxy;
			return !isDegenerate && tminx <= rmaxx && tminy <= rmaxy && tmaxx >= rminx && tmaxy >= rminy;
		}
		getFullBounds(normalized = false) {
			const {
				projection
			} = this;
			const bounds = [...this.rootBounds];
			if (projection && normalized) {
				bounds[0] = projection.convertLongitudeToProjection(bounds[0]);
				bounds[1] = projection.convertLatitudeToProjection(bounds[1]);
				bounds[2] = projection.convertLongitudeToProjection(bounds[2]);
				bounds[3] = projection.convertLatitudeToProjection(bounds[3]);
			}
			return bounds;
		}
		getTileBounds(x, y, level, normalized = false) {
			const {
				flipY,
				pixelOverlap,
				projection
			} = this;
			const {
				tilePixelWidth,
				tilePixelHeight,
				pixelWidth,
				pixelHeight
			} = this.getLevel(level);
			let tileLeft = tilePixelWidth * x - pixelOverlap;
			let tileTop = tilePixelHeight * y - pixelOverlap;
			let tileRight = tileLeft + tilePixelWidth + pixelOverlap * 2;
			let tileBottom = tileTop + tilePixelHeight + pixelOverlap * 2;

			// clamp
			tileLeft = Math.max(tileLeft, 0);
			tileTop = Math.max(tileTop, 0);
			tileRight = Math.min(tileRight, pixelWidth);
			tileBottom = Math.min(tileBottom, pixelHeight);

			// normalized
			tileLeft = tileLeft / pixelWidth;
			tileRight = tileRight / pixelWidth;
			tileTop = tileTop / pixelHeight;
			tileBottom = tileBottom / pixelHeight;

			// invert y
			if (flipY) {
				const extents = (tileBottom - tileTop) / 2;
				const centerY = (tileTop + tileBottom) / 2;
				const invCenterY = 1.0 - centerY;
				tileTop = invCenterY - extents;
				tileBottom = invCenterY + extents;
			}
			const bounds = [tileLeft, tileTop, tileRight, tileBottom];
			if (projection && !normalized) {
				bounds[0] = projection.convertProjectionToLongitude(bounds[0]);
				bounds[1] = projection.convertProjectionToLatitude(bounds[1]);
				bounds[2] = projection.convertProjectionToLongitude(bounds[2]);
				bounds[3] = projection.convertProjectionToLatitude(bounds[3]);
			}
			return bounds;
		}
	}

	// TODO: support queries for detail at level - ie projected pixel size for geometric error mapping
	// Goes here or in "TilingScheme"?
	class TiledImageSource extends DataCache {
		constructor() {
			super();
			this.tiling = new TilingScheme();
			this.fetchOptions = {};
			this.fetchData = (...args) => fetch(...args);
		}

		// async function for initializing the tiled image set
		init(url) {}

		// helper for processing the buffer into a texture
		async processBufferToTexture(buffer) {
			const blob = new Blob([buffer]);
			const imageBitmap = await createImageBitmap(blob, {
				premultiplyAlpha: 'none',
				colorSpaceConversion: 'none',
				imageOrientation: 'flipY'
			});
			const texture = new t3d.Texture2D();
			texture.image = imageBitmap;
			texture.generateMipmaps = false;
			texture.minFilter = t3d.TEXTURE_FILTER.LINEAR;
			texture.encoding = t3d.TEXEL_ENCODING_TYPE.SRGB;
			texture.version++;
			return texture;
		}
		getMemoryUsage(tex) {
			return 0;
		}

		// fetch the item with the given key fields
		fetchItem(...args) {
			const url = this.getUrl(...args);
			return this.fetchData(url, this.fetchOptions).then(res => res.arrayBuffer()).then(buffer => this.processBufferToTexture(buffer));
		}

		// dispose of the item that was fetched
		disposeItem(texture) {
			texture.dispose();
			if (texture.image instanceof ImageBitmap) {
				texture.image.close();
			}
		}
		getUrl(...args) {}
	}

	class XYZImageSource extends TiledImageSource {
		constructor(options = {}) {
			super();
			const {
				levels = 20,
				tileDimension = 256
			} = options;
			this.tileDimension = tileDimension;
			this.levels = levels;
			this.url = null;
		}
		getUrl(x, y, level) {
			return this.url.replace('{z}', level).replace('{x}', x).replace('{y}', y);
		}
		init(url) {
			// transform the url
			const {
				tiling,
				tileDimension,
				levels
			} = this;
			tiling.flipY = true;
			tiling.setProjection(new ProjectionScheme('EPSG:3857'));
			tiling.generateLevels(levels, 1, 1, {
				tilePixelWidth: tileDimension,
				tilePixelHeight: tileDimension
			});
			this.url = url;
			return Promise.resolve();
		}
	}

	class TMSImageSource extends TiledImageSource {
		constructor() {
			super();
			this.tileSets = null;
			this.extension = null;
			this.url = null;
		}
		getUrl(x, y, level) {
			const {
				url,
				extension,
				tileSets,
				tiling
			} = this;
			return new URL(`${parseInt(tileSets[level - tiling.minLevel].href)}/${x}/${y}.${extension}`, url).toString();
		}
		init(url) {
			return this.fetchData(new URL('tilemapresource.xml', url), this.fetchOptions).then(res => res.text()).then(text => {
				const {
					tiling
				} = this;

				// elements
				const xml = new DOMParser().parseFromString(text, 'text/xml');
				const boundingBox = xml.querySelector('BoundingBox');
				const origin = xml.querySelector('Origin');
				const tileFormat = xml.querySelector('TileFormat');
				const tileSets = xml.querySelector('TileSets').querySelectorAll('TileSet');

				// tile set definitions
				const tileSetList = [...tileSets].map(ts => ({
					href: parseInt(ts.getAttribute('href')),
					unitsPerPixel: parseFloat(ts.getAttribute('units-per-pixel')),
					order: parseInt(ts.getAttribute('order'))
				})).sort((a, b) => {
					return a.order - b.order;
				});

				// bounding box
				const minX = parseFloat(boundingBox.getAttribute('minx')) * t3d.MathUtils.DEG2RAD;
				const maxX = parseFloat(boundingBox.getAttribute('maxx')) * t3d.MathUtils.DEG2RAD;
				const minY = parseFloat(boundingBox.getAttribute('miny')) * t3d.MathUtils.DEG2RAD;
				const maxY = parseFloat(boundingBox.getAttribute('maxy')) * t3d.MathUtils.DEG2RAD;

				// origin in lat / lon
				const originX = parseFloat(origin.getAttribute('x')) * t3d.MathUtils.DEG2RAD;
				const originY = parseFloat(origin.getAttribute('y')) * t3d.MathUtils.DEG2RAD;

				// image dimensions in pixels
				const tileWidth = parseInt(tileFormat.getAttribute('width'));
				const tileHeight = parseInt(tileFormat.getAttribute('height'));
				const extension = tileFormat.getAttribute('extension');
				const srs = xml.querySelector('SRS').textContent;

				// assign settings
				this.extension = extension;
				this.url = url;
				this.tileSets = tileSetList;

				// initialize tiling and projection schemes
				tiling.setProjection(new ProjectionScheme(srs));
				tiling.setOrigin(originX, originY);
				tiling.setBounds(minX, minY, maxX, maxY);
				tileSetList.forEach(({
					order
				}) => {
					tiling.setLevel(order, {
						tileCountX: tiling.projection.tileCountX * 2 ** order,
						tilePixelWidth: tileWidth,
						tilePixelHeight: tileHeight
					});
				});
			});
		}
	}

	// Support for XYZ / Slippy tile systems


	// https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
	class XYZTilesPlugin extends EllipsoidProjectionTilesPlugin {
		constructor(options = {}) {
			const {
				levels = 20,
				tileDimension = 256,
				pixelSize = 1e-5,
				...rest
			} = options;
			super({
				pixelSize,
				...rest
			});
			this.name = 'XYZ_TILES_PLUGIN';
			this.imageSource = new XYZImageSource({
				levels,
				tileDimension
			});
		}
	}

	// Support for TMS tiles
	// https://wiki.osgeo.org/wiki/Tile_Map_Service_Specification
	// NOTE: Most, if not all, TMS generation implementations do not correctly support the Origin tag
	// and tile index offsets, including CesiumJS and Ion.
	class TMSTilesPlugin extends EllipsoidProjectionTilesPlugin {
		constructor(...args) {
			super(...args);
			this.name = 'TMS_TILES_PLUGIN';
			this.imageSource = new TMSImageSource();
		}
	}

	function zigZagDecode(value) {
		return value >> 1 ^ -(value & 1);
	}
	class QuantizedMeshLoaderBase extends LoaderBase {
		constructor(...args) {
			super(...args);
			this.fetchOptions.header = {
				Accept: 'application/vnd.quantized-mesh,application/octet-stream;q=0.9'
			};
		}
		loadAsync(...args) {
			const {
				fetchOptions
			} = this;
			fetchOptions.header = fetchOptions.header || {};
			fetchOptions.header['Accept'] = 'application/vnd.quantized-mesh,application/octet-stream;q=0.9';
			fetchOptions.header['Accept'] += ';extensions=octvertexnormals-watermask-metadata';
			return super.loadAsync(...args);
		}
		parse(buffer) {
			let pointer = 0;
			const view = new DataView(buffer);
			const readFloat64 = () => {
				const result = view.getFloat64(pointer, true);
				pointer += 8;
				return result;
			};
			const readFloat32 = () => {
				const result = view.getFloat32(pointer, true);
				pointer += 4;
				return result;
			};
			const readInt = () => {
				const result = view.getUint32(pointer, true);
				pointer += 4;
				return result;
			};
			const readByte = () => {
				const result = view.getUint8(pointer);
				pointer += 1;
				return result;
			};
			const readBuffer = (count, type) => {
				const result = new type(buffer, pointer, count); // eslint-disable-line new-cap
				pointer += count * type.BYTES_PER_ELEMENT;
				return result;
			};

			// extract header
			const header = {
				center: [readFloat64(), readFloat64(), readFloat64()],
				minHeight: readFloat32(),
				maxHeight: readFloat32(),
				sphereCenter: [readFloat64(), readFloat64(), readFloat64()],
				sphereRadius: readFloat64(),
				horizonOcclusionPoint: [readFloat64(), readFloat64(), readFloat64()]
			};

			// extract vertex data
			const vertexCount = readInt();
			const uBuffer = readBuffer(vertexCount, Uint16Array);
			const vBuffer = readBuffer(vertexCount, Uint16Array);
			const hBuffer = readBuffer(vertexCount, Uint16Array);
			const uResult = new Float32Array(vertexCount);
			const vResult = new Float32Array(vertexCount);
			const hResult = new Float32Array(vertexCount);

			// decode vertex data
			let u = 0;
			let v = 0;
			let h = 0;
			const MAX_VALUE = 32767;
			for (let i = 0; i < vertexCount; ++i) {
				u += zigZagDecode(uBuffer[i]);
				v += zigZagDecode(vBuffer[i]);
				h += zigZagDecode(hBuffer[i]);
				uResult[i] = u / MAX_VALUE;
				vResult[i] = v / MAX_VALUE;
				hResult[i] = h / MAX_VALUE;
			}

			// align pointer for index data
			const is32 = vertexCount > 65536;
			const bufferType = is32 ? Uint32Array : Uint16Array;
			if (is32) {
				pointer = Math.ceil(pointer / 4) * 4;
			} else {
				pointer = Math.ceil(pointer / 2) * 2;
			}

			// extract index data
			const triangleCount = readInt();
			const indices = readBuffer(triangleCount * 3, bufferType);

			// decode the index data
			let highest = 0;
			for (let i = 0; i < indices.length; ++i) {
				const code = indices[i];
				indices[i] = highest - code;
				if (code === 0) {
					++highest;
				}
			}

			// sort functions for the edges since they are not pre-sorted
			const vSort = (a, b) => vResult[b] - vResult[a];
			const vSortReverse = (a, b) => -vSort(a, b);
			const uSort = (a, b) => uResult[a] - uResult[b];
			const uSortReverse = (a, b) => -uSort(a, b);

			// get edge indices
			const westVertexCount = readInt();
			const westIndices = readBuffer(westVertexCount, bufferType);
			westIndices.sort(vSort);
			const southVertexCount = readInt();
			const southIndices = readBuffer(southVertexCount, bufferType);
			southIndices.sort(uSort);
			const eastVertexCount = readInt();
			const eastIndices = readBuffer(eastVertexCount, bufferType);
			eastIndices.sort(vSortReverse);
			const northVertexCount = readInt();
			const northIndices = readBuffer(northVertexCount, bufferType);
			northIndices.sort(uSortReverse);
			const edgeIndices = {
				westIndices,
				southIndices,
				eastIndices,
				northIndices
			};

			// parse extensions
			const extensions = {};
			while (pointer < view.byteLength) {
				const extensionId = readByte();
				const extensionLength = readInt();
				if (extensionId === 1) {
					// oct encoded normals
					const xy = readBuffer(vertexCount * 2, Uint8Array);
					const normals = new Float32Array(vertexCount * 3);

					// https://github.com/CesiumGS/cesium/blob/baaabaa49058067c855ad050be73a9cdfe9b6ac7/packages/engine/Source/Core/AttributeCompression.js#L119-L140
					for (let i = 0; i < vertexCount; i++) {
						let x = xy[2 * i + 0] / 255 * 2 - 1;
						let y = xy[2 * i + 1] / 255 * 2 - 1;
						const z = 1.0 - (Math.abs(x) + Math.abs(y));
						if (z < 0.0) {
							const oldVX = x;
							x = (1.0 - Math.abs(y)) * signNotZero(oldVX);
							y = (1.0 - Math.abs(oldVX)) * signNotZero(y);
						}
						const len = Math.sqrt(x * x + y * y + z * z);
						normals[3 * i + 0] = x / len;
						normals[3 * i + 1] = y / len;
						normals[3 * i + 2] = z / len;
					}
					extensions['octvertexnormals'] = {
						extensionId,
						normals
					};
				} else if (extensionId === 2) {
					// water mask
					const size = extensionLength === 1 ? 1 : 256;
					const mask = readBuffer(size * size, Uint8Array);
					extensions['watermask'] = {
						extensionId,
						mask,
						size
					};
				} else if (extensionId === 4) {
					// metadata
					const jsonLength = readInt();
					const jsonBuffer = readBuffer(jsonLength, Uint8Array);
					const json = new TextDecoder().decode(jsonBuffer);
					extensions['metadata'] = {
						extensionId,
						json: JSON.parse(json)
					};
				}
			}
			return {
				header,
				indices,
				vertexData: {
					u: uResult,
					v: vResult,
					height: hResult
				},
				edgeIndices,
				extensions
			};
		}
	}
	function signNotZero(v) {
		return v < 0.0 ? -1 : 1.0;
	}

	const _norm = /* @__PURE__ */new t3d.Vector3();
	const _tri = /* @__PURE__ */new t3d.Triangle();
	const _uvh = /* @__PURE__ */new t3d.Vector3();
	const _pos$1 = /* @__PURE__ */new t3d.Vector3();
	class QuantizedMeshLoader extends QuantizedMeshLoaderBase {
		constructor(manager = t3d.DefaultLoadingManager) {
			super();
			this.manager = manager;
			this.ellipsoid = new Ellipsoid();
			this.skirtLength = 1000;
			this.smoothSkirtNormals = true;
			this.solid = false;

			// set the range of the tile
			this.minLat = -Math.PI / 2;
			this.maxLat = Math.PI / 2;
			this.minLon = -Math.PI;
			this.maxLon = Math.PI;
		}
		parse(buffer) {
			const {
				ellipsoid,
				solid,
				skirtLength,
				smoothSkirtNormals,
				minLat,
				maxLat,
				minLon,
				maxLon
			} = this;
			const {
				header,
				indices,
				vertexData,
				edgeIndices,
				extensions
			} = super.parse(buffer);
			const geometry = new t3d.Geometry();
			const material = new t3d.PBRMaterial();
			material.roughness = 1.0;
			material.metalness = 0.0;
			const mesh = new t3d.Mesh(geometry, material);
			mesh.position.set(...header.center);
			const includeNormals = 'octvertexnormals' in extensions;
			const vertexCount = vertexData.u.length;
			const positions = [];
			const uvs = [];
			const indexArr = [];
			const normals = [];
			// const groupOffset = 0;
			// const materialIndex = 0;

			// construct terrain
			for (let i = 0; i < vertexCount; i++) {
				readUVHeight(i, _uvh);
				readPosition(_uvh.x, _uvh.y, _uvh.z, _pos$1);
				uvs.push(_uvh.x, _uvh.y);
				positions.push(..._pos$1);
			}
			for (let i = 0, l = indices.length; i < l; i++) {
				indexArr.push(indices[i]);
			}
			if (includeNormals) {
				const extNormals = extensions['octvertexnormals'].normals;
				for (let i = 0, l = extNormals.length; i < l; i++) {
					normals.push(extNormals[i]);
				}
			}

			// add material group
			// geometry.addGroup(groupOffset, indices.length, materialIndex);
			// groupOffset += indices.length;
			// materialIndex++;

			// create a lower cap
			if (solid) {
				const indexOffset = positions.length / 3;
				for (let i = 0; i < vertexCount; i++) {
					readUVHeight(i, _uvh);
					readPosition(_uvh.x, _uvh.y, _uvh.z, _pos$1, -skirtLength);
					uvs.push(_uvh.x, _uvh.y);
					positions.push(..._pos$1);
				}
				for (let i = indices.length - 1; i >= 0; i--) {
					indexArr.push(indices[i] + indexOffset);
				}
				if (includeNormals) {
					const extNormals = extensions['octvertexnormals'].normals;
					for (let i = 0, l = extNormals.length; i < l; i++) {
						normals.push(-extNormals[i]);
					}
				}

				// add material group
				// geometry.addGroup(groupOffset, indices.length, materialIndex);
				// groupOffset += indices.length;
				// materialIndex++;
			}

			// construct skirts
			if (skirtLength > 0) {
				const {
					westIndices,
					eastIndices,
					southIndices,
					northIndices
				} = edgeIndices;

				// construct the indices
				let offset;

				// west
				const westStrip = constructEdgeStrip(westIndices);
				offset = positions.length / 3;
				uvs.push(...westStrip.uv);
				positions.push(...westStrip.positions);
				for (let i = 0, l = westStrip.indices.length; i < l; i++) {
					indexArr.push(westStrip.indices[i] + offset);
				}

				// east
				const eastStrip = constructEdgeStrip(eastIndices);
				offset = positions.length / 3;
				uvs.push(...eastStrip.uv);
				positions.push(...eastStrip.positions);
				for (let i = 0, l = eastStrip.indices.length; i < l; i++) {
					indexArr.push(eastStrip.indices[i] + offset);
				}

				// south
				const southStrip = constructEdgeStrip(southIndices);
				offset = positions.length / 3;
				uvs.push(...southStrip.uv);
				positions.push(...southStrip.positions);
				for (let i = 0, l = southStrip.indices.length; i < l; i++) {
					indexArr.push(southStrip.indices[i] + offset);
				}

				// north
				const northStrip = constructEdgeStrip(northIndices);
				offset = positions.length / 3;
				uvs.push(...northStrip.uv);
				positions.push(...northStrip.positions);
				for (let i = 0, l = northStrip.indices.length; i < l; i++) {
					indexArr.push(northStrip.indices[i] + offset);
				}

				// add the normals
				if (includeNormals) {
					normals.push(...westStrip.normals);
					normals.push(...eastStrip.normals);
					normals.push(...southStrip.normals);
					normals.push(...northStrip.normals);
				}

				// add material group
				// geometry.addGroup(groupOffset, indices.length, materialIndex);
				// groupOffset += indices.length;
				// materialIndex++;
			}

			// shift the positions by the center of the tile
			for (let i = 0, l = positions.length; i < l; i += 3) {
				positions[i + 0] -= header.center[0];
				positions[i + 1] -= header.center[1];
				positions[i + 2] -= header.center[2];
			}

			// generate geometry and mesh
			const indexBuffer = positions.length / 3 > 65535 ? new Uint32Array(indexArr) : new Uint16Array(indexArr);
			geometry.setIndex(new t3d.Attribute(new t3d.Buffer(indexBuffer, 1)));
			geometry.addAttribute('a_Position', new t3d.Attribute(new t3d.Buffer(new Float32Array(positions), 3)));
			geometry.addAttribute('a_Uv', new t3d.Attribute(new t3d.Buffer(new Float32Array(uvs), 2)));
			if (includeNormals) {
				geometry.addAttribute('a_Normal', new t3d.Attribute(new t3d.Buffer(new Float32Array(normals), 3)));
			}

			// generate the water texture
			if ('watermask' in extensions) {
				// invert the mask data
				// TODO: this inversion step can be a bit slow
				const {
					mask,
					size
				} = extensions['watermask'];
				const maskBuffer = new Uint8Array(2 * size * size);
				for (let i = 0, l = mask.length; i < l; i++) {
					const v = mask[i] === 255 ? 0 : 255;
					maskBuffer[2 * i + 0] = v;
					maskBuffer[2 * i + 1] = v;
				}
				const map = new t3d.Texture2D();
				map.image = {
					data: maskBuffer,
					width: size,
					height: size
				};
				map.flipY = true;
				map.format = t3d.PIXEL_FORMAT.RG;
				map.type = t3d.PIXEL_TYPE.UNSIGNED_BYTE;
				map.minFilter = t3d.TEXTURE_FILTER.LINEAR_MIPMAP_LINEAR;
				map.magFilter = t3d.TEXTURE_FILTER.LINEAR;
				map.version++;
				material.roughnessMap = map;
			}

			// set metadata
			mesh.userData.minHeight = header.minHeight;
			mesh.userData.maxHeight = header.maxHeight;
			if ('metadata' in extensions) {
				mesh.userData.metadata = extensions['metadata'].json;
			}
			return mesh;
			function readUVHeight(index, target) {
				target.x = vertexData.u[index];
				target.y = vertexData.v[index];
				target.z = vertexData.height[index];
				return target;
			}
			function readPosition(u, v, h, target, heightOffset = 0) {
				const height = t3d.MathUtils.lerp(header.minHeight, header.maxHeight, h);
				const lon = t3d.MathUtils.lerp(minLon, maxLon, u);
				const lat = t3d.MathUtils.lerp(minLat, maxLat, v);
				ellipsoid.getCartographicToPosition(lat, lon, height + heightOffset, target);
				return target;
			}
			function constructEdgeStrip(indices) {
				const topUvs = [];
				const topPos = [];
				const botUvs = [];
				const botPos = [];
				const sideIndices = [];
				for (let i = 0, l = indices.length; i < l; i++) {
					readUVHeight(indices[i], _uvh);
					topUvs.push(_uvh.x, _uvh.y);
					botUvs.push(_uvh.x, _uvh.y);
					readPosition(_uvh.x, _uvh.y, _uvh.z, _pos$1);
					topPos.push(..._pos$1);
					readPosition(_uvh.x, _uvh.y, _uvh.z, _pos$1, -skirtLength);
					botPos.push(..._pos$1);
				}
				const triCount = indices.length - 1;
				for (let i = 0; i < triCount; i++) {
					const t0 = i;
					const t1 = i + 1;
					const b0 = i + indices.length;
					const b1 = i + indices.length + 1;
					sideIndices.push(t0, b0, t1);
					sideIndices.push(t1, b0, b1);
				}
				let normals = null;
				if (includeNormals) {
					const total = (topPos.length + botPos.length) / 3;
					if (smoothSkirtNormals) {
						normals = new Array(total * 3);
						const extNormals = extensions['octvertexnormals'].normals;
						const botOffset = normals.length / 2;
						for (let i = 0, l = total / 2; i < l; i++) {
							const index = indices[i];
							const i3 = 3 * i;
							const nx = extNormals[3 * index + 0];
							const ny = extNormals[3 * index + 1];
							const nz = extNormals[3 * index + 2];
							normals[i3 + 0] = nx;
							normals[i3 + 1] = ny;
							normals[i3 + 2] = nz;
							normals[botOffset + i3 + 0] = nx;
							normals[botOffset + i3 + 1] = ny;
							normals[botOffset + i3 + 2] = nz;
						}
					} else {
						normals = [];
						_tri.a.fromArray(topPos, 0);
						_tri.b.fromArray(botPos, 0);
						_tri.c.fromArray(topPos, 3);
						_tri.getNormal(_norm);
						for (let i = 0; i < total; i++) {
							normals.push(..._norm);
						}
					}
				}
				return {
					uv: [...topUvs, ...botUvs],
					positions: [...topPos, ...botPos],
					indices: sideIndices,
					normals
				};
			}
		}

		// generates a child mesh in the given quadrant using the same settings as the loader
		clipToQuadrant(mesh, left, bottom) {
			// scratch vectors
			const _uv0 = new t3d.Vector3();
			const _uv1 = new t3d.Vector3();
			const _pos0 = new t3d.Vector3();
			const _pos1 = new t3d.Vector3();
			const _pos2 = new t3d.Vector3();
			const _pos3 = new t3d.Vector3();
			const _temp = new t3d.Vector3();
			const _temp2 = new t3d.Vector3();
			const _cart = {};

			// helper variables
			const SPLIT_VALUE = 0.5;
			const triPool = new TrianglePool();
			const vertNames = ['a', 'b', 'c'];
			const {
				ellipsoid,
				skirtLength,
				solid,
				smoothSkirtNormals
			} = this;

			// source geometry
			const sourceGeometry = mesh.geometry;
			const normal = sourceGeometry.attributes.a_Normal;
			const index = sourceGeometry.index;

			// geometry data
			let nextIndex = 0;
			const vertToNewIndexMap = {};
			const newPosition = [];
			const newNormal = normal ? [] : null;
			const newUv = [];
			const newIndex = [];

			// uv offsets
			const xUvOffset = left ? 0 : -0.5;
			const yUvOffset = bottom ? 0 : -0.5;

			// iterate over each group separately to retain the group information
			const geometry = new t3d.Geometry();
			const capGroup = sourceGeometry.groups[0];

			// construct the cap geometry
			// const newStart = newIndex.length;
			// const materialIndex = 0;
			for (let i = capGroup.start / 3; i < (capGroup.start + capGroup.count) / 3; i++) {
				const i0 = index.getX(i * 3 + 0);
				const i1 = index.getX(i * 3 + 1);
				const i2 = index.getX(i * 3 + 2);
				const tri = triPool.get();
				tri.setFromAttributeAndIndices(sourceGeometry, i0, i1, i2);

				// split the triangle by the first axis
				const xResult = [];
				splitTriangle(tri, 'x', left, xResult);

				// split the triangles by the second axis
				const yResult = [];
				for (let t = 0, l = xResult.length; t < l; t++) {
					splitTriangle(xResult[t], 'y', bottom, yResult);
				}

				// save the geometry
				const {
					minLat,
					maxLat,
					minLon,
					maxLon,
					ellipsoid
				} = this;
				for (let t = 0, l = yResult.length; t < l; t++) {
					const tri = yResult[t];
					vertNames.forEach(n => {
						const uv = tri.uv[n];
						if (uv.x !== SPLIT_VALUE && uv.y !== SPLIT_VALUE) {
							return;
						}
						const point = tri.position[n];
						const lat = t3d.MathUtils.lerp(minLat, maxLat, uv.y);
						const lon = t3d.MathUtils.lerp(minLon, maxLon, uv.x);
						point.add(mesh.position);
						ellipsoid.getPositionToCartographic(point, _cart);
						ellipsoid.getCartographicToPosition(lat, lon, _cart.height, point);
						point.sub(mesh.position);
					});
					pushVertex(tri.position.a, tri.uv.a, tri.normal.a);
					pushVertex(tri.position.b, tri.uv.b, tri.normal.b);
					pushVertex(tri.position.c, tri.uv.c, tri.normal.c);
				}
				triPool.reset();
			}

			// geometry.addGroup(newStart, newIndex.length - newStart, materialIndex);
			// materialIndex++;

			// construct bottom cap
			const capTriangles = newIndex.length / 3;
			if (solid) {
				// newStart = newIndex.length;
				for (let i = capTriangles * 3 - 1; i >= 0; i--) {
					const index = newIndex[i];
					_temp.fromArray(newPosition, index * 3).add(mesh.position);
					ellipsoid.getPositionToNormal(_temp, _temp);
					_pos0.fromArray(newPosition, index * 3).addScaledVector(_temp, -skirtLength);
					_uv0.fromArray(newUv, index * 2);
					_temp.fromArray(newNormal, index * 3);
					pushVertex(_pos0, _uv0, _temp);
				}

				// geometry.addGroup(newStart, newIndex.length - newStart, materialIndex);
				// materialIndex++;
			}

			// construct the skirt
			if (skirtLength > 0) {
				// TODO: this seems to have some problematic cases at the root tiles near the poles
				// newStart = newIndex.length;
				for (let i = 0; i < capTriangles; i++) {
					const triOffset = 3 * i;
					for (let e = 0; e < 3; e++) {
						const ne = (e + 1) % 3;
						const i0 = newIndex[triOffset + e];
						const i1 = newIndex[triOffset + ne];
						_uv0.fromArray(newUv, i0 * 2);
						_uv1.fromArray(newUv, i1 * 2);

						// find the vertices that lie on the edge
						if (_uv0.x === _uv1.x && (_uv0.x === 0 || _uv0.x === SPLIT_VALUE || _uv0.x === 1.0) || _uv0.y === _uv1.y && (_uv0.y === 0 || _uv0.y === SPLIT_VALUE || _uv0.y === 1.0)) {
							_pos0.fromArray(newPosition, i0 * 3);
							_pos1.fromArray(newPosition, i1 * 3);
							const u0 = _pos0;
							const u1 = _pos1;
							const b0 = _pos2.copy(_pos0);
							const b1 = _pos3.copy(_pos1);
							_temp.copy(b0).add(mesh.position);
							ellipsoid.getPositionToNormal(_temp, _temp);
							b0.addScaledVector(_temp, -skirtLength);
							_temp.copy(b1).add(mesh.position);
							ellipsoid.getPositionToNormal(_temp, _temp);
							b1.addScaledVector(_temp, -skirtLength);
							if (smoothSkirtNormals && newNormal) {
								_temp.fromArray(newNormal, i0 * 3);
								_temp2.fromArray(newNormal, i1 * 3);
							} else {
								_temp.subVectors(u0, u1);
								_temp2.subVectors(u0, b0).cross(_temp).normalize();
								_temp.copy(_temp2);
							}
							pushVertex(u1, _uv1, _temp2);
							pushVertex(u0, _uv0, _temp);
							pushVertex(b0, _uv0, _temp);
							pushVertex(u1, _uv1, _temp2);
							pushVertex(b0, _uv0, _temp);
							pushVertex(b1, _uv1, _temp2);
						}
					}
				}

				// geometry.addGroup(newStart, newIndex.length - newStart, materialIndex);
				// materialIndex++;
			}

			// offset the uvs
			for (let i = 0, l = newUv.length; i < l; i += 2) {
				newUv[i] = (newUv[i] + xUvOffset) * 2.0;
				newUv[i + 1] = (newUv[i + 1] + yUvOffset) * 2.0;
			}

			// new geometry
			const indexBuffer = newPosition.length / 3 > 65535 ? new Uint32Array(newIndex) : new Uint16Array(newIndex);
			geometry.setIndex(new t3d.Attribute(new t3d.Buffer(indexBuffer, 1)));
			geometry.setAttribute('a_Position', new t3d.Attribute(new t3d.Buffer(new Float32Array(newPosition), 3)));
			geometry.setAttribute('a_Uv', new t3d.Attribute(new t3d.Buffer(new Float32Array(newUv), 2)));
			if (normal) {
				geometry.setAttribute('a_Normal', new t3d.Attribute(new t3d.Buffer(new Float32Array(newNormal), 3)));
			}

			// new mesh
			const result = new t3d.Mesh(geometry, mesh.material.clone());
			result.position.copy(mesh.position);
			result.quaternion.copy(mesh.quaternion);
			result.scale.copy(mesh.scale);
			result.userData.minHeight = mesh.userData.minHeight;
			result.userData.maxHeight = mesh.userData.maxHeight;
			return result;
			function splitTriangle(tri, axis, negativeSide, target) {
				// TODO: clean up, add scratch variables, optimize
				const edgeIndices = [];
				const edges = [];
				const lerpValues = [];
				for (let i = 0; i < 3; i++) {
					const v = vertNames[i];
					const nv = vertNames[(i + 1) % 3];
					const p = tri.uv[v];
					const np = tri.uv[nv];
					const pValue = p[axis];
					const npValue = np[axis];

					// if the uv values span across the halfway divide
					if (pValue < SPLIT_VALUE !== npValue < SPLIT_VALUE || pValue === SPLIT_VALUE) {
						edgeIndices.push(i);
						edges.push([v, nv]);
						lerpValues.push(t3d.MathUtils.mapLinear(SPLIT_VALUE, pValue, npValue, 0, 1));
					}
				}
				if (edgeIndices.length !== 2) {
					const minBound = Math.min(tri.uv.a[axis], tri.uv.b[axis], tri.uv.c[axis]);
					if (minBound < SPLIT_VALUE === negativeSide) {
						target.push(tri);
					}
				} else if (edgeIndices.length === 2) {
					// TODO: how can we determine which triangles actually need to be added here ahead of time
					const tri0 = triPool.get();
					const tri1 = triPool.get();
					const tri2 = triPool.get();
					const sequential = (edgeIndices[0] + 1) % 3 === edgeIndices[1];
					if (sequential) {
						tri0.lerpVertex(tri, edges[0][0], edges[0][1], lerpValues[0], 'a');
						tri0.copyVertex(tri, edges[0][1], 'b');
						tri0.lerpVertex(tri, edges[1][0], edges[1][1], lerpValues[1], 'c');
						tri0.uv.a[axis] = SPLIT_VALUE;
						tri0.uv.c[axis] = SPLIT_VALUE;
						tri1.lerpVertex(tri, edges[0][0], edges[0][1], lerpValues[0], 'a');
						tri1.copyVertex(tri, edges[1][1], 'b');
						tri1.copyVertex(tri, edges[0][0], 'c');
						tri1.uv.a[axis] = SPLIT_VALUE;
						tri2.lerpVertex(tri, edges[0][0], edges[0][1], lerpValues[0], 'a');
						tri2.lerpVertex(tri, edges[1][0], edges[1][1], lerpValues[1], 'b');
						tri2.copyVertex(tri, edges[1][1], 'c');
						tri2.uv.a[axis] = SPLIT_VALUE;
						tri2.uv.b[axis] = SPLIT_VALUE;
					} else {
						tri0.lerpVertex(tri, edges[0][0], edges[0][1], lerpValues[0], 'a');
						tri0.lerpVertex(tri, edges[1][0], edges[1][1], lerpValues[1], 'b');
						tri0.copyVertex(tri, edges[0][0], 'c');
						tri0.uv.a[axis] = SPLIT_VALUE;
						tri0.uv.b[axis] = SPLIT_VALUE;
						tri1.lerpVertex(tri, edges[0][0], edges[0][1], lerpValues[0], 'a');
						tri1.copyVertex(tri, edges[0][1], 'b');
						tri1.lerpVertex(tri, edges[1][0], edges[1][1], lerpValues[1], 'c');
						tri1.uv.a[axis] = SPLIT_VALUE;
						tri1.uv.c[axis] = SPLIT_VALUE;
						tri2.copyVertex(tri, edges[0][1], 'a');
						tri2.copyVertex(tri, edges[1][0], 'b');
						tri2.lerpVertex(tri, edges[1][0], edges[1][1], lerpValues[1], 'c');
						tri2.uv.c[axis] = SPLIT_VALUE;
					}
					let minBound;
					minBound = Math.min(tri0.uv.a[axis], tri0.uv.b[axis], tri0.uv.c[axis]);
					if (minBound < SPLIT_VALUE === negativeSide) {
						target.push(tri0);
					}
					minBound = Math.min(tri1.uv.a[axis], tri1.uv.b[axis], tri1.uv.c[axis]);
					if (minBound < SPLIT_VALUE === negativeSide) {
						target.push(tri1);
					}
					minBound = Math.min(tri2.uv.a[axis], tri2.uv.b[axis], tri2.uv.c[axis]);
					if (minBound < SPLIT_VALUE === negativeSide) {
						target.push(tri2);
					}
				}
			}

			// hash the vertex for index generation
			function hashVertex(x, y, z) {
				const scalar = 1e5;
				const additive = 0.5;
				const hx = ~~(x * scalar + additive);
				const hy = ~~(y * scalar + additive);
				const hz = ~~(z * scalar + additive);
				return `${hx}_${hy}_${hz}`;
			}

			// add the vertex to the geometry
			function pushVertex(pos, uv, norm) {
				let hash = hashVertex(pos.x, pos.y, pos.z);
				if (newNormal) {
					hash += `_${hashVertex(norm.x, norm.y, norm.z)}`;
				}
				if (!(hash in vertToNewIndexMap)) {
					vertToNewIndexMap[hash] = nextIndex;
					nextIndex++;
					newPosition.push(pos.x, pos.y, pos.z);
					newUv.push(uv.x, uv.y);
					if (newNormal) {
						newNormal.push(norm.x, norm.y, norm.z);
					}
				}
				const index = vertToNewIndexMap[hash];
				newIndex.push(index);
				return index;
			}
		}
	}

	// Pool of reusable triangles
	class TrianglePool {
		constructor() {
			this.pool = [];
			this.index = 0;
		}
		get() {
			if (this.index >= this.pool.length) {
				const tri = new AttributeTriangle();
				this.pool.push(tri);
			}
			const res = this.pool[this.index];
			this.index++;
			return res;
		}
		reset() {
			this.index = 0;
		}
	}

	// Set of triangle definitions for quantized mesh attributes
	class AttributeTriangle {
		constructor() {
			this.position = new t3d.Triangle();
			this.uv = new t3d.Triangle();
			this.normal = new t3d.Triangle();
		}
		setFromAttributeAndIndices(geometry, i0, i1, i2) {
			this.position.setFromAttributeAndIndices(geometry.attributes.a_Position, i0, i1, i2);
			this.uv.setFromAttributeAndIndices(geometry.attributes.a_Uv, i0, i1, i2);
			if (geometry.attributes.a_Normal) {
				this.normal.setFromAttributeAndIndices(geometry.attributes.a_Normal, i0, i1, i2);
			}
		}
		lerpVertex(other, e0, e1, alpha, targetVertex) {
			this.position[targetVertex].lerpVectors(other.position[e0], other.position[e1], alpha);
			this.uv[targetVertex].lerpVectors(other.uv[e0], other.uv[e1], alpha);
			this.normal[targetVertex].lerpVectors(other.normal[e0], other.normal[e1], alpha);
		}
		copyVertex(other, fromVertex, targetVertex) {
			this.position[targetVertex].copy(other.position[fromVertex]);
			this.uv[targetVertex].copy(other.uv[fromVertex]);
			this.normal[targetVertex].copy(other.normal[fromVertex]);
		}
	}

	const TILE_X = Symbol('TILE_X');
	const TILE_Y = Symbol('TILE_Y');
	const TILE_LEVEL = Symbol('TILE_LEVEL');
	const TILE_AVAILABLE = Symbol('TILE_AVAILABLE');

	// We don't know the height ranges for the tile set on load so assume a large range and
	// adjust it once the tiles have actually loaded based on the min and max height
	const INITIAL_HEIGHT_RANGE = 1e4;
	const _vec = /* @__PURE__ */new t3d.Vector3();

	// Checks if the given tile is available
	function isTileAvailable(available, level, x, y) {
		if (level < available.length) {
			// TODO: consider a binary search
			const availableSet = available[level];
			for (let i = 0, l = availableSet.length; i < l; i++) {
				const {
					startX,
					startY,
					endX,
					endY
				} = availableSet[i];
				if (x >= startX && x <= endX && y >= startY && y <= endY) {
					return true;
				}
			}
		}
		return false;
	}

	// Calculates the max level that can be loaded.
	function getMaxLevel(layer) {
		const {
			available = null,
			maxzoom = null
		} = layer;
		return maxzoom === null ? available.length - 1 : maxzoom;
	}

	// Calculates whether metadata availability is present - returns -1 if not.
	function getMetadataAvailability(layer) {
		const {
			metadataAvailability = -1
		} = layer;
		return metadataAvailability;
	}

	// Calculates whether the given tile should have metadata availability
	function getTileHasMetadata(tile, layer) {
		const level = tile[TILE_LEVEL];
		const metadataAvailability = getMetadataAvailability(layer);
		const maxLevel = getMaxLevel(layer);
		return level < maxLevel && metadataAvailability !== -1 && level % metadataAvailability === 0;
	}

	// Constructs the url for the given tile content
	function getContentUrl(x, y, level, version, layer) {
		return layer.tiles[0].replace(/{\s*z\s*}/g, level).replace(/{\s*x\s*}/g, x).replace(/{\s*y\s*}/g, y).replace(/{\s*version\s*}/g, version);
	}
	class QuantizedMeshPlugin {
		constructor(options = {}) {
			const {
				useRecommendedSettings = true,
				skirtLength = null,
				smoothSkirtNormals = true,
				solid = false
			} = options;
			this.name = 'QUANTIZED_MESH_PLUGIN';
			this.tiles = null;
			this.layer = null;
			this.useRecommendedSettings = useRecommendedSettings;
			this.skirtLength = skirtLength;
			this.smoothSkirtNormals = smoothSkirtNormals;
			this.solid = solid;
			this.attribution = null;
			this.tiling = new TilingScheme();
			this.projection = new ProjectionScheme();
		}

		// Plugin function
		init(tiles) {
			// TODO: should we avoid setting this globally?
			tiles.fetchOptions.headers = tiles.fetchOptions.headers || {};
			tiles.fetchOptions.headers.Accept = 'application/vnd.quantized-mesh,application/octet-stream;q=0.9';
			if (this.useRecommendedSettings) {
				tiles.errorTarget = 2;
			}
			this.tiles = tiles;
		}
		loadRootTileSet() {
			const {
				tiles
			} = this;

			// initialize href to resolve the root in case it's specified as a relative url
			let url = new URL('layer.json', new URL(tiles.rootURL, location.href));
			tiles.invokeAllPlugins(plugin => url = plugin.preprocessURL ? plugin.preprocessURL(url, null) : url);
			return tiles.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(url, this.tiles.fetchOptions)).then(res => res.json()).then(json => {
				this.layer = json;
				const {
					projection: layerProjection = 'EPSG:4326',
					extensions = [],
					attribution = '',
					available = null
				} = json;
				const {
					tiling,
					tiles,
					projection
				} = this;

				// attribution
				if (attribution) {
					this.attribution = {
						value: attribution,
						type: 'string',
						collapsible: true
					};
				}

				// extensions
				if (extensions.length > 0) {
					tiles.fetchOptions.headers['Accept'] += `;extensions=${extensions.join('-')}`;
				}

				// initialize tiling, projection
				projection.setScheme(layerProjection);
				const {
					tileCountX,
					tileCountY
				} = projection;
				tiling.setProjection(projection);
				tiling.generateLevels(getMaxLevel(json) + 1, tileCountX, tileCountY);

				// initialize children
				const children = [];
				for (let x = 0; x < tileCountX; x++) {
					const child = this.createChild(0, x, 0, available);
					if (child) {
						children.push(child);
					}
				}

				// produce the tile set root
				const tileset = {
					asset: {
						version: '1.1'
					},
					geometricError: Infinity,
					root: {
						refine: 'REPLACE',
						geometricError: Infinity,
						boundingVolume: {
							region: [...this.tiling.getFullBounds(), -INITIAL_HEIGHT_RANGE, INITIAL_HEIGHT_RANGE]
						},
						children: children,
						[TILE_AVAILABLE]: available,
						[TILE_LEVEL]: -1
					}
				};
				let baseUrl = tiles.rootURL;
				tiles.invokeAllPlugins(plugin => baseUrl = plugin.preprocessURL ? plugin.preprocessURL(baseUrl, null) : baseUrl);
				tiles.preprocessTileSet(tileset, baseUrl);
				return tileset;
			});
		}
		async parseToMesh(buffer, tile, extension, uri) {
			const {
				skirtLength,
				solid,
				smoothSkirtNormals,
				tiles
			} = this;

			// set up loader
			const ellipsoid = tiles.ellipsoid;
			const loader = new QuantizedMeshLoader(tiles.manager);
			loader.ellipsoid.copy(ellipsoid);
			loader.solid = solid;
			loader.smoothSkirtNormals = smoothSkirtNormals;
			loader.skirtLength = skirtLength === null ? tile.geometricError : skirtLength;

			// split the parent tile if needed
			let result;
			if (extension === 'tile_split') {
				// split the parent tile
				const searchParams = new URL(uri).searchParams;
				const left = searchParams.get('left') === 'true';
				const bottom = searchParams.get('bottom') === 'true';
				const [west, south, east, north] = tile.parent.boundingVolume.region;
				loader.minLat = south;
				loader.maxLat = north;
				loader.minLon = west;
				loader.maxLon = east;
				result = loader.clipToQuadrant(tile.parent.cached.scene, left, bottom);
			} else {
				const [west, south, east, north] = tile.boundingVolume.region;
				loader.minLat = south;
				loader.maxLat = north;
				loader.minLon = west;
				loader.maxLon = east;

				// parse the tile data
				result = loader.parse(buffer);
			}

			// adjust the bounding region to be more accurate based on the contents of the terrain file
			// NOTE: The debug region bounds are only created after the tile is first shown so the debug
			// region bounding volume will have the correct dimensions.
			const {
				minHeight,
				maxHeight,
				metadata
			} = result.userData;
			tile.boundingVolume.region[4] = minHeight;
			tile.boundingVolume.region[5] = maxHeight;
			tile.cached.boundingVolume.setRegionData(ellipsoid, ...tile.boundingVolume.region);

			// use the geometric error value if it's present
			if (metadata) {
				if ('geometricerror' in metadata) {
					tile.geometricError = metadata.geometricerror;
				}

				// if the tile hasn't been expanded yet and isn't in the queue to do so then
				// mark it for expansion again
				const hasMetadata = getTileHasMetadata(tile, this.layer);
				if (hasMetadata && 'available' in metadata && tile.children.length === 0) {
					// add an offset to account for the current and previous layers
					tile[TILE_AVAILABLE] = [...new Array(tile[TILE_LEVEL] + 1).fill(null), ...metadata.available];
				}
			}

			// NOTE: we expand children only once the parent mesh data is loaded to ensure the mesh
			// data is ready for clipping. It's possible that this child data gets to the parse stage
			// first, otherwise, while the parent is still downloading.
			// Ideally we would be able to guarantee parents are loaded first but this is an odd case.
			this.expandChildren(tile);
			return result;
		}
		getAttributions(target) {
			if (this.attribution) {
				target.push(this.attribution);
			}
		}

		// Local functions
		createChild(level, x, y, available) {
			const {
				tiles,
				layer,
				tiling,
				projection
			} = this;
			const ellipsoid = tiles.ellipsoid;
			const isAvailable = available === null || isTileAvailable(available, level, x, y);
			const url = getContentUrl(x, y, level, 1, layer);
			const region = [...tiling.getTileBounds(x, y, level), -INITIAL_HEIGHT_RANGE, INITIAL_HEIGHT_RANGE];
			const [/* east */ /* minHeight */, south,, north,, maxHeight] = region;
			const midLat = south > 0 !== north > 0 ? 0 : Math.min(Math.abs(south), Math.abs(north));

			// get the projected perimeter
			ellipsoid.getCartographicToPosition(midLat, 0, maxHeight, _vec);
			_vec.z = 0;

			// https://github.com/CesiumGS/cesium/blob/53889cbed2a91d38e0fae4b6f2dcf6783632fc92/packages/engine/Source/Scene/QuadtreeTileProvider.js#L24-L31
			// Implicit quantized mesh tile error halves with every layer
			const tileCountX = projection.tileCountX;
			const maxRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
			const rootGeometricError = maxRadius * 2 * Math.PI * 0.25 / (65 * tileCountX);
			const geometricError = rootGeometricError / 2 ** level;

			// Create the child
			const tile = {
				[TILE_AVAILABLE]: null,
				[TILE_LEVEL]: level,
				[TILE_X]: x,
				[TILE_Y]: y,
				refine: 'REPLACE',
				geometricError: geometricError,
				boundingVolume: {
					region
				},
				content: isAvailable ? {
					uri: url
				} : null,
				children: []
			};

			// if we're relying on tile metadata availability then skip storing the tile metadata
			if (!getTileHasMetadata(tile, layer)) {
				tile[TILE_AVAILABLE] = available;
			}
			return tile;
		}
		expandChildren(tile) {
			const level = tile[TILE_LEVEL];
			const x = tile[TILE_X];
			const y = tile[TILE_Y];
			const available = tile[TILE_AVAILABLE];
			let hasChildren = false;
			for (let cx = 0; cx < 2; cx++) {
				for (let cy = 0; cy < 2; cy++) {
					const child = this.createChild(level + 1, 2 * x + cx, 2 * y + cy, available);
					if (child.content !== null) {
						tile.children.push(child);
						hasChildren = true;
					} else {
						tile.children.push(child);
						child.content = {
							uri: `tile.tile_split?bottom=${cy === 0}&left=${cx === 0}`
						};
					}
				}
			}
			if (!hasChildren) {
				tile.children.length = 0;
			}
		}
		fetchData(uri, options) {
			// if this is our custom url indicating a tile split then return fake response
			if (/tile_split/.test(uri)) {
				return new ArrayBuffer();
			}
		}
		disposeTile(tile) {
			// dispose of the generated children past the metadata layer to avoid accumulating too much
			if (getTileHasMetadata(tile, this.layer)) {
				tile.children.length = 0;
				tile.__childrenProcessed = 0;
				tile[TILE_AVAILABLE] = null;
			}
			tile.children.length = 0;
			tile.__childrenProcessed = 0;
		}
	}

	class CesiumIonAuthPlugin {
		constructor({
			apiToken,
			assetId = null,
			autoRefreshToken = false,
			useRecommendedSettings = true
		}) {
			this.name = 'CESIUM_ION_AUTH_PLUGIN';
			this.priority = -Infinity;
			this.apiToken = apiToken;
			this.assetId = assetId;
			this.autoRefreshToken = autoRefreshToken;
			this.useRecommendedSettings = useRecommendedSettings;
			this.tiles = null;
			this.endpointURL = null;
			this._bearerToken = null;
			this._tileSetVersion = -1;
			this._tokenRefreshPromise = null;
			this._attributions = [];
			this._disposed = false;
		}
		init(tiles) {
			if (this.assetId !== null) {
				tiles.rootURL = `https://api.cesium.com/v1/assets/${this.assetId}/endpoint`;
			}
			this.tiles = tiles;
			this.endpointURL = tiles.rootURL;

			// reset the tiles in case this plugin was removed and re-added
			tiles.resetFailedTiles();
		}
		loadRootTileSet() {
			// ensure we have an up-to-date token and root url, then trigger the internal
			// root tile set load function
			return this._refreshToken().then(() => {
				return this.tiles.invokeOnePlugin(plugin => plugin !== this && plugin.loadRootTileSet && plugin.loadRootTileSet());
			});
		}
		preprocessURL(uri) {
			uri = new URL(uri);
			if (/^http/.test(uri.protocol) && this._tileSetVersion != -1) {
				uri.searchParams.append('v', this._tileSetVersion);
			}
			return uri.toString();
		}
		fetchData(uri, options) {
			const tiles = this.tiles;
			if (tiles.getPluginByName('GOOGLE_CLOUD_AUTH_PLUGIN') !== null) {
				return null;
			} else {
				return Promise.resolve().then(async () => {
					// wait for the token to refresh if loading
					if (this._tokenRefreshPromise !== null) {
						await this._tokenRefreshPromise;
						uri = this.preprocessURL(uri);
					}
					const res = await fetch(uri, options);
					if (res.status >= 400 && res.status <= 499 && this.autoRefreshToken) {
						await this._refreshToken(options);
						return fetch(this.preprocessURL(uri), options);
					} else {
						return res;
					}
				});
			}
		}
		getAttributions(target) {
			if (this.tiles.visibleTiles.size > 0) {
				target.push(...this._attributions);
			}
		}
		_refreshToken(options) {
			if (this._tokenRefreshPromise === null) {
				// construct the url to fetch the endpoint
				const url = new URL(this.endpointURL);
				url.searchParams.append('access_token', this.apiToken);
				this._tokenRefreshPromise = fetch(url, options).then(res => {
					if (this._disposed) {
						return null;
					}
					if (!res.ok) {
						throw new Error(`CesiumIonAuthPlugin: Failed to load data with error code ${res.status}`);
					}
					return res.json();
				}).then(json => {
					if (this._disposed) {
						return null;
					}
					const tiles = this.tiles;
					if ('externalType' in json) {
						const url = new URL(json.options.url);
						tiles.rootURL = json.options.url;

						// if the tile set is "external" then assume it's a google API tile set
						tiles.registerPlugin(new GoogleCloudAuthPlugin({
							apiToken: url.searchParams.get('key'),
							autoRefreshToken: this.autoRefreshToken,
							useRecommendedSettings: this.useRecommendedSettings
						}));
					} else {
						// GLTF
						// CZML
						// KML
						// GEOJSON
						if (json.type === 'TERRAIN' && tiles.getPluginByName('QUANTIZED_MESH_PLUGIN') === null) {
							tiles.registerPlugin(new QuantizedMeshPlugin({
								useRecommendedSettings: this.useRecommendedSettings
							}));
						} else if (json.type === 'IMAGERY' && tiles.getPluginByName('TMS_TILES_PLUGIN') === null) {
							tiles.registerPlugin(new TMSTilesPlugin({
								useRecommendedSettings: this.useRecommendedSettings,
								shape: 'ellipsoid'
							}));
						}
						tiles.rootURL = json.url;
						tiles.fetchOptions.headers = tiles.fetchOptions.headers || {};
						tiles.fetchOptions.headers.Authorization = `Bearer ${json.accessToken}`;

						// save the version key if present
						if (url.searchParams.has('v') && this._tileSetVersion === -1) {
							const url = new URL(json.url);
							this._tileSetVersion = url.searchParams.get('v');
						}
						this._bearerToken = json.accessToken;
						if (json.attributions) {
							this._attributions = json.attributions.map(att => ({
								value: att.html,
								type: 'html',
								collapsible: att.collapsible
							}));
						}
					}
					this._tokenRefreshPromise = null;
					return json;
				});

				// dispatch an error if we fail to refresh the token
				this._tokenRefreshPromise.catch(error => {
					this.tiles.dispatchEvent({
						type: 'load-error',
						tile: null,
						error,
						url
					});
				});
			}
			return this._tokenRefreshPromise;
		}
		dispose() {
			this._disposed = true;
		}
	}

	class Box3Helper extends t3d.Mesh {
		constructor(box, color = 0xffff00) {
			const indices = new Uint16Array([0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7]);
			const positions = [1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1];
			const geometry = new t3d.Geometry();
			geometry.setIndex(new t3d.Attribute(new t3d.Buffer(indices, 1)));
			geometry.addAttribute('a_Position', new t3d.Attribute(new t3d.Buffer(new Float32Array(positions), 3)));

			// Skip update bounding box
			// Because we may not want to consider Box3Helper's bounding box
			geometry.computeBoundingSphere();
			const material = new t3d.LineMaterial();
			material.diffuse.setHex(color);
			super(geometry, material);
			this.box = box;
		}
		updateMatrix(force) {
			const box = this.box;
			box.getCenter(this.position);
			if (box.isEmpty()) {
				this.scale.multiplyScalar(0);
			} else {
				box.getSize(this.scale);
				this.scale.multiplyScalar(0.5);
			}
			super.updateMatrix(force);
		}
	}
	Box3Helper.prototype.isBox3Helper = true;

	class SphereHelper extends t3d.Mesh {
		constructor(sphere, color = 0xffff00, angleSteps = 40) {
			const positions = [];
			for (let i = 0; i < 3; i++) {
				const axis1 = axes[i];
				const axis2 = axes[(i + 1) % 3];
				_vector.set(0, 0, 0);
				for (let a = 0; a < angleSteps; a++) {
					let angle;
					angle = 2 * Math.PI * a / (angleSteps - 1);
					_vector[axis1] = Math.sin(angle);
					_vector[axis2] = Math.cos(angle);
					positions.push(_vector.x, _vector.y, _vector.z);
					angle = 2 * Math.PI * (a + 1) / (angleSteps - 1);
					_vector[axis1] = Math.sin(angle);
					_vector[axis2] = Math.cos(angle);
					positions.push(_vector.x, _vector.y, _vector.z);
				}
			}
			const geometry = new t3d.Geometry();
			geometry.addAttribute('a_Position', new t3d.Attribute(new t3d.Buffer(new Float32Array(positions), 3)));
			geometry.computeBoundingSphere();
			const lineMaterial = new t3d.LineMaterial();
			lineMaterial.diffuse.setHex(color);
			super(geometry, lineMaterial);
			this.sphere = sphere;
		}
		updateMatrix(force) {
			const sphere = this.sphere;
			this.position.copy(sphere.center);
			if (sphere.isEmpty()) {
				this.scale.setScalar(0);
			} else {
				this.scale.setScalar(sphere.radius);
			}
			super.updateMatrix(force);
		}
	}
	SphereHelper.prototype.isSphereHelper = true;
	const _vector = new t3d.Vector3();
	const axes = ['x', 'y', 'z'];

	/**
	 * EdgesBuilder is a helper for building edges geometry from given triangles data.
	 */
	const EdgesBuilder = {
		/**
		 * @param {Array} bufferArray - Flat buffer array containing vertex positions.
		 * @param {Array} [indices] - Flat buffer array of indices, must be multiple of 3.
		 * @param {object} [options={}] - The options object.
		 * @param {number} [options.thresholdAngle=1] - An edge is only rendered if the angle (in degrees) between the face normals of the adjoining faces exceeds this value.
		 * @param {number} [options.stride=3] - The number of values of the array that should be associated with a particular vertex.
		 * @param {number} [options.offset=0] - The offset in the buffer array where the position starts.
		 * @returns {object} - The edges geometry data.
		 */
		getGeometryData: function (bufferArray, indices, options = {}) {
			const thresholdAngle = options.thresholdAngle !== undefined ? options.thresholdAngle : 1;
			const stride = options.stride !== undefined ? options.stride : 3;
			const offset = options.offset !== undefined ? options.offset : 0;
			let i, j, l, key, face;
			const result = [];

			/** merge vertices */

			const verticesMap = {};
			const unique = [],
				changes = [];
			const precisionPoints = 4; // number of decimal points, e.g. 4 for epsilon of 0.0001
			const precision = Math.pow(10, precisionPoints);
			let offsetIndex, x, y, z;
			l = bufferArray.length / stride;
			for (i = 0; i < l; i++) {
				offsetIndex = i * stride + offset;
				x = bufferArray[offsetIndex + 0];
				y = bufferArray[offsetIndex + 1];
				z = bufferArray[offsetIndex + 2];
				key = Math.round(x * precision) + '_' + Math.round(y * precision) + '_' + Math.round(z * precision);
				if (verticesMap[key] === undefined) {
					verticesMap[key] = i;
					unique.push(x, y, z);
					changes[i] = unique.length / 3 - 1;
				} else {
					changes[i] = changes[verticesMap[key]];
				}
			}

			/** get faces	(vertices and normal) */

			const faces = [];
			if (indices) {
				l = indices.length / 3;
				for (i = 0; i < l; i++) {
					face = {
						i: [0, 0, 0],
						n: [1, 1, 1]
					};
					face.i[0] = changes[indices[i * 3 + 0]];
					face.i[1] = changes[indices[i * 3 + 1]];
					face.i[2] = changes[indices[i * 3 + 2]];
					computeFaceNormal(face, unique);
					faces.push(face);
				}
			} else {
				for (i = 0; i < l; i++) {
					face = {
						i: [0, 0, 0],
						n: [1, 1, 1]
					};
					face.i[0] = changes[i * 3 + 0];
					face.i[1] = changes[i * 3 + 1];
					face.i[2] = changes[i * 3 + 2];
					computeFaceNormal(face, unique);
					faces.push(face);
				}
			}

			/**
			 * get edges { index1: edge[ 0 ], index2: edge[ 1 ], face1: i, face2: undefined }
			 */
			let edge1, edge2;
			const edge = [0, 0],
				edges = {};
			for (i = 0, l = faces.length; i < l; i++) {
				face = faces[i];
				for (j = 0; j < 3; j++) {
					edge1 = face.i[j];
					edge2 = face.i[(j + 1) % 3];
					edge[0] = Math.min(edge1, edge2);
					edge[1] = Math.max(edge1, edge2);
					key = edge[0] + ',' + edge[1];
					if (edges[key] === undefined) {
						edges[key] = {
							index1: edge[0],
							index2: edge[1],
							face1: i,
							face2: undefined
						};
					} else {
						edges[key].face2 = i;
					}
				}
			}

			/** edges filter */
			const thresholdDot = Math.cos(DEG2RAD * thresholdAngle);
			for (key in edges) {
				const e = edges[key];
				// an edge is only rendered if the angle (in degrees) between the face normals of the adjoining faces exceeds this value. default = 1 degree.
				if (e.face2 === undefined || dot(faces[e.face1].n, faces[e.face2].n) <= thresholdDot) {
					result.push(unique[e.index1 * 3 + 0], unique[e.index1 * 3 + 1], unique[e.index1 * 3 + 2], unique[e.index2 * 3 + 0], unique[e.index2 * 3 + 1], unique[e.index2 * 3 + 2]);
				}
			}

			/** return */
			return {
				positions: result
			};
		}
	};
	function dot(v1, v2) {
		return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	}
	function computeFaceNormal(face, buffer) {
		const vAX = buffer[face.i[0] * 3 + 0];
		const vAY = buffer[face.i[0] * 3 + 1];
		const vAZ = buffer[face.i[0] * 3 + 2];
		const vBX = buffer[face.i[1] * 3 + 0];
		const vBY = buffer[face.i[1] * 3 + 1];
		const vBZ = buffer[face.i[1] * 3 + 2];
		const vCX = buffer[face.i[2] * 3 + 0];
		const vCY = buffer[face.i[2] * 3 + 1];
		const vCZ = buffer[face.i[2] * 3 + 2];
		const cbX = vCX - vBX; // ax
		const cbY = vCY - vBY; // ay
		const cbZ = vCZ - vBZ; // az

		const abX = vAX - vBX; // bx
		const abY = vAY - vBY; // by
		const abZ = vAZ - vBZ; // bz

		let nX = cbY * abZ - cbZ * abY;
		let nY = cbZ * abX - cbX * abZ;
		let nZ = cbX * abY - cbY * abX;
		const nLen = Math.sqrt(nX * nX + nY * nY + nZ * nZ);
		nX /= nLen;
		nY /= nLen;
		nZ /= nLen;
		face.n[0] = nX;
		face.n[1] = nY;
		face.n[2] = nZ;
	}
	const DEG2RAD = Math.PI / 180;

	const _pos = new t3d.Vector3();
	function getRegionGeometry(ellipsoidRegion) {
		// retrieve the relevant fields
		const {
			latRange,
			lonRange,
			heightRange
		} = ellipsoidRegion;
		const {
			x: latStart,
			y: latEnd
		} = latRange;
		const {
			x: lonStart,
			y: lonEnd
		} = lonRange;
		const {
			x: heightStart,
			y: heightEnd
		} = heightRange;

		// get the attributes
		const geometry = new t3d.BoxGeometry(1, 1, 1, 32, 32);
		const {
			a_Position: position
		} = geometry.attributes;

		// perturb the position buffer into an ellipsoid region
		for (let i = 0, l = position.buffer.count; i < l; i++) {
			_pos.fromArray(position.buffer.array, i * 3);
			const lat = t3d.MathUtils.mapLinear(_pos.x, -0.5, 0.5, latStart, latEnd);
			const lon = t3d.MathUtils.mapLinear(_pos.y, -0.5, 0.5, lonStart, lonEnd);
			let height = heightStart;
			if (_pos.z < 0) {
				height = heightEnd;
			}
			ellipsoidRegion.getCartographicToPosition(lat, lon, height, _pos);
			_pos.toArray(position.buffer.array, i * 3);
		}
		return geometry;
	}
	class EllipsoidRegionHelper extends t3d.Mesh {
		constructor(ellipsoidRegion = new EllipsoidRegion(), color = 0xffff00) {
			super(new t3d.Geometry(), new t3d.LineMaterial());
			this.ellipsoidRegion = ellipsoidRegion;
			this.material.diffuse.setHex(color);
			this.update();
			this.raycast = () => {}; // disable raycasting
		}
		update() {
			this.geometry.dispose();
			const regionGeometry = getRegionGeometry(this.ellipsoidRegion);
			const {
				positions
			} = EdgesBuilder.getGeometryData(regionGeometry.attributes.a_Position.buffer.array, regionGeometry.index.buffer.array, {
				thresholdAngle: 80
			});
			const geometry = new t3d.Geometry();
			geometry.addAttribute('a_Position', new t3d.Attribute(new t3d.Buffer(new Float32Array(positions), 3)));
			// geometry.computeBoundingBox();
			// geometry.computeBoundingSphere();

			this.geometry = geometry;
			this.geometry.computeBoundingBox();
			this.geometry.computeBoundingSphere();
		}
		dispose() {
			this.geometry.dispose();
			this.material.dispose();
		}
	}

	const ORIGINAL_MATERIAL = Symbol('ORIGINAL_MATERIAL');
	const HAS_RANDOM_COLOR = Symbol('HAS_RANDOM_COLOR');
	const HAS_RANDOM_NODE_COLOR = Symbol('HAS_RANDOM_NODE_COLOR');
	const LOAD_TIME = Symbol('LOAD_TIME');
	const PARENT_BOUND_REF_COUNT = Symbol('PARENT_BOUND_REF_COUNT');
	const _sphere = /* @__PURE__ */new t3d.Sphere();
	const emptyRaycast = () => {};
	const colors = {};

	// Return a consistant random color for an index
	function getIndexedRandomColor(index) {
		if (!colors[index]) {
			const h = Math.random();
			const s = 0.5 + Math.random() * 0.5;
			const l = 0.375 + Math.random() * 0.25;
			colors[index] = new t3d.Color3().setHSL(h, s, l);
		}
		return colors[index];
	}

	// color modes
	const NONE = 0;
	const SCREEN_ERROR = 1;
	const GEOMETRIC_ERROR = 2;
	const DISTANCE = 3;
	const DEPTH = 4;
	const RELATIVE_DEPTH = 5;
	const IS_LEAF = 6;
	const RANDOM_COLOR = 7;
	const RANDOM_NODE_COLOR = 8;
	const CUSTOM_COLOR = 9;
	const LOAD_ORDER = 10;
	const ColorModes = Object.freeze({
		NONE,
		SCREEN_ERROR,
		GEOMETRIC_ERROR,
		DISTANCE,
		DEPTH,
		RELATIVE_DEPTH,
		IS_LEAF,
		RANDOM_COLOR,
		RANDOM_NODE_COLOR,
		CUSTOM_COLOR,
		LOAD_ORDER
	});
	class DebugTilesPlugin {
		static get ColorModes() {
			return ColorModes;
		}
		get unlit() {
			return this._unlit;
		}
		set unlit(v) {
			if (v !== this._unlit) {
				this._unlit = v;
				this.materialsNeedUpdate = true;
			}
		}
		get colorMode() {
			return this._colorMode;
		}
		set colorMode(v) {
			if (v !== this._colorMode) {
				this._colorMode = v;
				this.materialsNeedUpdate = true;
			}
		}
		constructor(options) {
			options = {
				displayParentBounds: false,
				displayBoxBounds: false,
				displaySphereBounds: false,
				displayRegionBounds: false,
				colorMode: NONE,
				maxDebugDepth: -1,
				maxDebugDistance: -1,
				maxDebugError: -1,
				customColorCallback: null,
				unlit: false,
				enabled: true,
				...options
			};
			this.name = 'DEBUG_TILES_PLUGIN';
			this.tiles = null;
			this._colorMode = null;
			this._unlit = null;
			this.materialsNeedUpdate = false;
			this.extremeDebugDepth = -1;
			this.extremeDebugError = -1;
			this.boxGroup = null;
			this.sphereGroup = null;
			this.regionGroup = null;

			// options
			this._enabled = options.enabled;
			this._displayParentBounds = options.displayParentBounds;
			this.displayBoxBounds = options.displayBoxBounds;
			this.displaySphereBounds = options.displaySphereBounds;
			this.displayRegionBounds = options.displayRegionBounds;
			this.colorMode = options.colorMode;
			this.maxDebugDepth = options.maxDebugDepth;
			this.maxDebugDistance = options.maxDebugDistance;
			this.maxDebugError = options.maxDebugError;
			this.customColorCallback = options.customColorCallback;
			this.unlit = options.unlit;
			this.getDebugColor = (value, target) => {
				target.setRGB(value, value, value);
			};
		}
		get enabled() {
			return this._enabled;
		}
		set enabled(v) {
			if (v !== this._enabled && this.tiles !== null) {
				if (v) {
					this.init(this.tiles);
				} else {
					this.dispose();
				}
			}
			this._enabled = v;
		}
		get displayParentBounds() {
			return this._displayParentBounds;
		}
		set displayParentBounds(v) {
			if (this._displayParentBounds !== v) {
				this._displayParentBounds = v;
				if (!v) {
					// Reset all ref counts
					traverseSet(this.tiles.root, null, tile => {
						tile[PARENT_BOUND_REF_COUNT] = null;
						this._onTileVisibilityChange(tile, tile.__visible);
					});
				} else {
					// Initialize ref count for existing tiles
					this.tiles.traverse(tile => {
						if (tile.__visible) {
							this._onTileVisibilityChange(tile, true);
						}
					});
				}
			}
		}

		// initialize the groups for displaying helpers, register events, and initialize existing tiles
		init(tiles) {
			this.tiles = tiles;
			if (!this.enabled) {
				return;
			}

			// initialize groups
			this.boxGroup = new t3d.Object3D();
			this.boxGroup.name = 'DebugTilesPlugin.boxGroup';
			tiles.add(this.boxGroup);
			this.boxGroup.updateMatrix();
			this.sphereGroup = new t3d.Object3D();
			this.sphereGroup.name = 'DebugTilesPlugin.sphereGroup';
			tiles.add(this.sphereGroup);
			this.sphereGroup.updateMatrix();
			this.regionGroup = new t3d.Object3D();
			this.regionGroup.name = 'DebugTilesPlugin.regionGroup';
			tiles.add(this.regionGroup);
			this.regionGroup.updateMatrix();

			// register events
			this._onLoadTileSetCB = () => {
				this._initExtremes();
			};
			this._onLoadModelCB = ({
				scene,
				tile
			}) => {
				this._onLoadModel(scene, tile);
			};
			this._onDisposeModelCB = ({
				tile
			}) => {
				this._onDisposeModel(tile);
			};
			this._onUpdateAfterCB = () => {
				this._onUpdateAfter();
			};
			this._onTileVisibilityChangeCB = ({
				scene,
				tile,
				visible
			}) => {
				this._onTileVisibilityChange(tile, visible);
			};
			tiles.addEventListener('load-tile-set', this._onLoadTileSetCB);
			tiles.addEventListener('load-model', this._onLoadModelCB);
			tiles.addEventListener('dispose-model', this._onDisposeModelCB);
			tiles.addEventListener('update-after', this._onUpdateAfterCB);
			tiles.addEventListener('tile-visibility-change', this._onTileVisibilityChangeCB);
			this._initExtremes();

			// initialize an already-loaded tiles
			tiles.traverse(tile => {
				if (tile.cached.scene) {
					this._onLoadModel(tile.cached.scene, tile);
				}
			});
			tiles.visibleTiles.forEach(tile => {
				this._onTileVisibilityChange(tile, true);
			});
		}
		getTileInformationFromActiveObject(object) {
			// Find which tile this scene is associated with. This is slow and
			// intended for debug purposes only.
			let targetTile = null;
			const activeTiles = this.tiles.activeTiles;
			activeTiles.forEach(tile => {
				if (targetTile) {
					return true;
				}
				const scene = tile.cached.scene;
				if (scene) {
					scene.traverse(c => {
						if (c === object) {
							targetTile = tile;
						}
					});
				}
			});
			if (targetTile) {
				return {
					distanceToCamera: targetTile.__distanceFromCamera,
					geometricError: targetTile.geometricError,
					screenSpaceError: targetTile.__error,
					depth: targetTile.__depth,
					isLeaf: targetTile.__isLeaf
				};
			} else {
				return null;
			}
		}
		_initExtremes() {
			if (!(this.tiles && this.tiles.root)) {
				return;
			}

			// initialize the extreme values of the hierarchy
			let maxDepth = -1;
			let maxError = -1;

			// Note that we are not using this.tiles.traverse()
			// as we don't want to pay the cost of preprocessing tiles.
			traverseSet(this.tiles.root, null, (tile, _, depth) => {
				maxDepth = Math.max(maxDepth, depth);
				maxError = Math.max(maxError, tile.geometricError);
			});
			this.extremeDebugDepth = maxDepth;
			this.extremeDebugError = maxError;
		}
		_onUpdateAfter() {
			const {
				tiles,
				colorMode
			} = this;
			if (!tiles.root) {
				return;
			}
			if (this.materialsNeedUpdate) {
				tiles.forEachLoadedModel(scene => {
					this._updateMaterial(scene);
				});
				this.materialsNeedUpdate = false;
			}

			// set box or sphere visibility
			this.boxGroup.visible = this.displayBoxBounds;
			this.sphereGroup.visible = this.displaySphereBounds;
			this.regionGroup.visible = this.displayRegionBounds;

			// get max values to use for materials
			let maxDepth = -1;
			if (this.maxDebugDepth === -1) {
				maxDepth = this.extremeDebugDepth;
			} else {
				maxDepth = this.maxDebugDepth;
			}
			let maxError = -1;
			if (this.maxDebugError === -1) {
				maxError = this.extremeDebugError;
			} else {
				maxError = this.maxDebugError;
			}
			let maxDistance = -1;
			if (this.maxDebugDistance === -1) {
				tiles.getBoundingSphere(_sphere);
				maxDistance = _sphere.radius;
			} else {
				maxDistance = this.maxDebugDistance;
			}
			const {
				errorTarget,
				visibleTiles
			} = tiles;
			let sortedTiles;
			if (colorMode === LOAD_ORDER) {
				sortedTiles = Array.from(visibleTiles).sort((a, b) => {
					return a[LOAD_TIME] - b[LOAD_TIME];
				});
			}

			// update plugins
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
					if (colorMode === RANDOM_NODE_COLOR) {
						h = Math.random();
						s = 0.5 + Math.random() * 0.5;
						l = 0.375 + Math.random() * 0.25;
					}
					if (c.material) {
						if (colorMode !== RANDOM_COLOR) {
							delete c.material[HAS_RANDOM_COLOR];
						}
						if (colorMode !== RANDOM_NODE_COLOR) {
							delete c.material[HAS_RANDOM_NODE_COLOR];
						}
						switch (colorMode) {
							case DEPTH:
								{
									const val = tile.__depth / maxDepth;
									this.getDebugColor(val, c.material.diffuse);
									break;
								}
							case RELATIVE_DEPTH:
								{
									const val = tile.__depthFromRenderedParent / maxDepth;
									this.getDebugColor(val, c.material.diffuse);
									break;
								}
							case SCREEN_ERROR:
								{
									const val = tile.__error / errorTarget;
									if (val > 1.0) {
										c.material.diffuse.setRGB(1.0, 0.0, 0.0);
									} else {
										this.getDebugColor(val, c.material.diffuse);
									}
									break;
								}
							case GEOMETRIC_ERROR:
								{
									const val = Math.min(tile.geometricError / maxError, 1);
									this.getDebugColor(val, c.material.diffuse);
									break;
								}
							case DISTANCE:
								{
									// We don't update the distance if the geometric error is 0.0 so
									// it will always be black.
									const val = Math.min(tile.__distanceFromCamera / maxDistance, 1);
									this.getDebugColor(val, c.material.diffuse);
									break;
								}
							case IS_LEAF:
								{
									if (!tile.children || tile.children.length === 0) {
										this.getDebugColor(1.0, c.material.diffuse);
									} else {
										this.getDebugColor(0.0, c.material.diffuse);
									}
									break;
								}
							case RANDOM_NODE_COLOR:
								{
									if (!c.material[HAS_RANDOM_NODE_COLOR]) {
										c.material.diffuse.setHSL(h, s, l);
										c.material[HAS_RANDOM_NODE_COLOR] = true;
									}
									break;
								}
							case RANDOM_COLOR:
								{
									if (!c.material[HAS_RANDOM_COLOR]) {
										c.material.diffuse.setHSL(h, s, l);
										c.material[HAS_RANDOM_COLOR] = true;
									}
									break;
								}
							case CUSTOM_COLOR:
								{
									if (this.customColorCallback) {
										this.customColorCallback(tile, c);
									} else {
										console.warn('DebugTilesPlugin: customColorCallback not defined');
									}
									break;
								}
							case LOAD_ORDER:
								{
									const value = sortedTiles.indexOf(tile);
									this.getDebugColor(value / (sortedTiles.length - 1), c.material.diffuse);
									break;
								}
						}
					}
				});
			});
		}
		_onTileVisibilityChange(tile, visible) {
			if (this.displayParentBounds) {
				traverseAncestors(tile, current => {
					if (current[PARENT_BOUND_REF_COUNT] == null) {
						current[PARENT_BOUND_REF_COUNT] = 0;
					}
					if (visible) {
						current[PARENT_BOUND_REF_COUNT]++;
					} else if (current[PARENT_BOUND_REF_COUNT] > 0) {
						current[PARENT_BOUND_REF_COUNT]--;
					}
					const tileVisible = current === tile && visible || this.displayParentBounds && current[PARENT_BOUND_REF_COUNT] > 0;
					this._updateBoundHelper(current, tileVisible);
				});
			} else {
				this._updateBoundHelper(tile, visible);
			}
		}
		_createBoundHelper(tile) {
			const tiles = this.tiles;
			const cached = tile.cached;
			const {
				sphere,
				obb,
				region
			} = cached.boundingVolume;
			if (obb) {
				// Create debug bounding box
				// In some cases the bounding box may have a scale of 0 in one dimension resulting
				// in the NaNs in an extracted rotation so we disable matrix updates instead.
				const boxHelperGroup = new t3d.Object3D();
				boxHelperGroup.name = 'DebugTilesPlugin.boxHelperGroup';
				boxHelperGroup.matrix.copy(obb._originBoxTransform);
				boxHelperGroup.matrixAutoUpdate = false;
				boxHelperGroup.matrixNeedsUpdate = false;
				const boxHelper = new Box3Helper(obb._originBox, getIndexedRandomColor(tile.__depth).getHex());
				boxHelper.raycast = emptyRaycast;
				boxHelperGroup.add(boxHelper);
				cached.boxHelperGroup = boxHelperGroup;
				if (tiles.visibleTiles.has(tile) && this.displayBoxBounds) {
					this.boxGroup.add(boxHelperGroup);
					boxHelperGroup.updateMatrix(true);
				}
			}
			if (sphere) {
				// Create debug bounding sphere
				const sphereHelper = new SphereHelper(sphere, getIndexedRandomColor(tile.__depth).getHex());
				sphereHelper.raycast = emptyRaycast;
				cached.sphereHelper = sphereHelper;
				if (tiles.visibleTiles.has(tile) && this.displaySphereBounds) {
					this.sphereGroup.add(sphereHelper);
					sphereHelper.updateMatrix(true);
				}
			}
			if (region) {
				// Create debug bounding region
				const regionHelper = new EllipsoidRegionHelper(region, getIndexedRandomColor(tile.__depth).getHex());
				regionHelper.raycast = emptyRaycast;
				cached.regionHelper = regionHelper;
				if (tiles.visibleTiles.has(tile) && this.displayRegionBounds) {
					this.regionGroup.add(regionHelper);
					regionHelper.updateMatrix(true);
				}
			}
		}
		_updateHelperMaterial(tile, material) {
			if (tile.__visible || !this.displayParentBounds) {
				material.opacity = 1;
			} else {
				material.opacity = 0.2;
			}
			material.transparent = material.opacity < 1;
		}
		_updateBoundHelper(tile, visible) {
			const cached = tile.cached;
			if (!cached) {
				return;
			}
			const sphereGroup = this.sphereGroup;
			const boxGroup = this.boxGroup;
			const regionGroup = this.regionGroup;
			if (visible && cached.boxHelperGroup == null && cached.sphereHelper == null && cached.regionHelper == null) {
				this._createBoundHelper(tile);
			}
			const boxHelperGroup = cached.boxHelperGroup;
			const sphereHelper = cached.sphereHelper;
			const regionHelper = cached.regionHelper;
			if (!visible) {
				if (boxHelperGroup) {
					boxGroup.remove(boxHelperGroup);
				}
				if (sphereHelper) {
					sphereGroup.remove(sphereHelper);
				}
				if (regionHelper) {
					regionGroup.remove(regionHelper);
				}
			} else {
				// TODO: consider updating the volumes based on the bounding regions here in case they've been changed
				if (boxHelperGroup) {
					boxGroup.add(boxHelperGroup);
					boxHelperGroup.updateMatrix(true);
					this._updateHelperMaterial(tile, boxHelperGroup.children[0].material);
				}
				if (sphereHelper) {
					sphereGroup.add(sphereHelper);
					sphereHelper.updateMatrix(true);
					this._updateHelperMaterial(tile, sphereHelper.material);
				}
				if (regionHelper) {
					regionGroup.add(regionHelper);
					regionHelper.updateMatrix(true);
					this._updateHelperMaterial(tile, regionHelper.material);
				}
			}
		}
		_updateMaterial(scene) {
			// update the materials for debug rendering
			const {
				colorMode,
				unlit
			} = this;
			scene.traverse(c => {
				if (!c.material) {
					return;
				}
				const currMaterial = c.material;
				const originalMaterial = c[ORIGINAL_MATERIAL];

				// dispose the previous material
				if (currMaterial !== originalMaterial) {
					currMaterial.dispose();
				}

				// assign the new material
				if (colorMode !== NONE || unlit) {
					if (c.material.drawMode === t3d.DRAW_MODE.POINTS) {
						const pointsMaterial = new t3d.PointsMaterial();
						pointsMaterial.size = originalMaterial.size;
						pointsMaterial.sizeAttenuation = originalMaterial.sizeAttenuation;
						c.material = pointsMaterial;
					} else if (c.geometry.instanceCount >= 0) {
						c.material = unlit ? new InstancedBasicMaterial() : new InstancedPBRMaterial();
						c.material.metalness = 0.0;
						c.material.roughness = 1.0;
						c.material.shading = unlit ? t3d.SHADING_TYPE.SMOOTH_SHADING : t3d.SHADING_TYPE.FLAT_SHADING;
					} else {
						c.material = unlit ? new t3d.BasicMaterial() : new t3d.PBRMaterial();
						c.material.metalness = 0.0;
						c.material.roughness = 1.0;
						c.material.shading = unlit ? t3d.SHADING_TYPE.SMOOTH_SHADING : t3d.SHADING_TYPE.FLAT_SHADING;
					}
					c.material.envMap = undefined;

					// if no debug rendering is happening then assign the material properties
					if (colorMode === NONE) {
						c.material.diffuseMap = originalMaterial.diffuseMap;
						c.material.diffuse.setRGB(originalMaterial.diffuse);
					}
				} else {
					c.material = originalMaterial;
				}
			});
		}
		_onLoadModel(scene, tile) {
			tile[LOAD_TIME] = performance.now();

			// Cache the original materials
			scene.traverse(c => {
				const material = c.material;
				if (material) {
					c[ORIGINAL_MATERIAL] = material;
				}
			});

			// Update the materials to align with the settings
			this._updateMaterial(scene);
		}
		_onDisposeModel(tile) {
			const cached = tile.cached;
			if (cached.boxHelperGroup) {
				cached.boxHelperGroup.children[0].geometry.dispose();
				delete cached.boxHelperGroup;
			}
			if (cached.sphereHelper) {
				cached.sphereHelper.geometry.dispose();
				delete cached.sphereHelper;
			}
			if (cached.regionHelper) {
				cached.regionHelper.geometry.dispose();
				delete cached.regionHelper;
			}
		}
		dispose() {
			if (!this.enabled) {
				return;
			}
			const tiles = this.tiles;
			tiles.removeEventListener('load-tile-set', this._onLoadTileSetCB);
			tiles.removeEventListener('load-model', this._onLoadModelCB);
			tiles.removeEventListener('dispose-model', this._onDisposeModelCB);
			tiles.removeEventListener('update-after', this._onUpdateAfterCB);
			tiles.removeEventListener('tile-visibility-change', this._onTileVisibilityChangeCB);

			// reset all materials
			this.colorMode = NONE;
			this.unlit = false;
			tiles.forEachLoadedModel(scene => {
				this._updateMaterial(scene);
			});

			// dispose of all helper objects
			tiles.traverse(tile => {
				this._onDisposeModel(tile);
			});
			this.boxGroup?.removeFromParent();
			this.sphereGroup?.removeFromParent();
			this.regionGroup?.removeFromParent();
		}
	}

	class ReorientationPlugin {
		constructor(options) {
			options = {
				up: '+z',
				recenter: true,
				lat: null,
				lon: null,
				height: 0,
				...options
			};
			this.tiles = null;
			this.up = options.up.toLowerCase().replace(/\s+/, '');
			this.lat = options.lat;
			this.lon = options.lon;
			this.height = options.height;
			this.recenter = options.recenter;
			this._callback = null;
		}
		init(tiles) {
			this.tiles = tiles;
			this._callback = () => {
				const {
					up,
					lat,
					lon,
					height,
					recenter
				} = this;
				if (lat !== null && lon !== null) {
					// if the latitude and longitude are provided then remove the position offset
					this.transformLatLonHeightToOrigin(lat, lon, height);
				} else {
					const {
						ellipsoid
					} = tiles;
					const minRadii = Math.min(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
					tiles.getBoundingSphere(sphere);
					if (sphere.center.getLength() > minRadii * 0.5) {
						// otherwise see if this is possibly a tile set on the surface of the globe based on the positioning
						const cart = {};
						ellipsoid.getPositionToCartographic(sphere.center, cart);
						this.transformLatLonHeightToOrigin(cart.lat, cart.lon, cart.height);
					} else {
						// lastly fall back to orienting the up direction to +Y
						tiles.euler.set(0, 0, 0);
						switch (up) {
							case 'x':
							case '+x':
								tiles.euler.z = Math.PI / 2;
								break;
							case '-x':
								tiles.euler.z = -Math.PI / 2;
								break;
							case 'y':
							case '+y':
								break;
							case '-y':
								tiles.euler.z = Math.PI;
								break;
							case 'z':
							case '+z':
								tiles.euler.x = -Math.PI / 2;
								break;
							case '-z':
								tiles.euler.x = Math.PI / 2;
								break;
						}
						tiles.position.copy(sphere.center).applyEuler(tiles.euler).multiplyScalar(-1);
					}
				}
				if (!recenter) {
					tiles.position.setScalar(0);
				}
				tiles.removeEventListener('load-tile-set', this._callback);
			};
			tiles.addEventListener('load-tile-set', this._callback);
		}
		transformLatLonHeightToOrigin(lat, lon, height = 0) {
			const tiles = this.tiles;
			const {
				ellipsoid
			} = tiles;

			// get ENU orientation (Z facing north and X facing west) and position
			ellipsoid.getRotationMatrixFromAzElRoll(lat, lon, 0, 0, 0, tiles.matrix, OBJECT_FRAME);
			ellipsoid.getCartographicToPosition(lat, lon, height, vec);

			// adjust the tiles matrix
			tiles.matrix.setPosition(vec).invert().decompose(tiles.position, tiles.quaternion, tiles.scale);
			tiles.updateMatrix();
		}
		dispose() {
			const tiles = this.tiles;
			tiles.position.setScalar(0);
			tiles.quaternion.identity();
			tiles.scale.set(1, 1, 1);
			tiles.addEventListener('load-tile-set', this._callback);
		}
	}
	const sphere = new t3d.Sphere();
	const vec = new t3d.Vector3();

	class LoadParser {
		static parse(context, loader) {
			const extension = getUrlExtension(context.url);
			let pr = null;
			if (extension === 'gltf') {
				pr = loader.loadFile(context.url).then(buffer => {
					context.options.buffer = buffer;
				});
			} else {
				pr = loader.loadFile(context.url, 'arraybuffer').then(buffer => {
					context.options.buffer = buffer;
				});
			}
			return pr;
		}
	}

	const _quaternion = new t3d.Quaternion();
	t3d.Vector3.prototype.applyEuler = function (euler) {
		return this.applyQuaternion(_quaternion.setFromEuler(euler));
	};
	t3d.Vector3.prototype.applyAxisAngle = function (axis, angle) {
		return this.applyQuaternion(_quaternion.setFromAxisAngle(axis, angle));
	};
	t3d.Triangle.prototype.setFromAttributeAndIndices = function (attribute, i0, i1, i2) {
		const array = attribute.buffer.array;
		const itemSize = attribute.size;
		const offset = attribute.offset;
		this.a.fromArray(array, i0 * itemSize + offset);
		this.b.fromArray(array, i1 * itemSize + offset);
		this.c.fromArray(array, i2 * itemSize + offset);
	};
	t3d.Object3D.prototype.removeFromParent = function () {
		const parent = this.parent;
		if (parent !== null) {
			parent.remove(this);
		}
		return this;
	};
	t3d.MathUtils.DEG2RAD = Math.PI / 180;
	t3d.MathUtils.RAD2DEG = 180 / Math.PI;
	let oldMethod;
	oldMethod = t3d.Camera.prototype.setOrtho;
	t3d.Camera.prototype.setOrtho = function (left, right, bottom, top, near, far) {
		this.left = left;
		this.right = right;
		this.bottom = bottom;
		this.top = top;
		this.near = near;
		this.far = far;
		this.zoom = 1;
		this.isPerspectiveCamera = false;
		this.isOrthographicCamera = true;
		oldMethod.call(this, left, right, bottom, top, near, far);
	};
	oldMethod = t3d.Camera.prototype.setPerspective;
	t3d.Camera.prototype.setPerspective = function (fov, aspect, near, far) {
		this.fov = fov;
		this.aspect = aspect;
		this.near = near;
		this.far = far;
		this.isPerspectiveCamera = true;
		this.isOrthographicCamera = false;
		oldMethod.call(this, fov, aspect, near, far);
	};
	t3d.Camera.prototype.updateProjectionMatrix = function () {
		if (this.isOrthographicCamera) {
			this.setOrtho(this.left, this.right, this.bottom, this.top, this.near, this.far);
		} else if (this.isPerspectiveCamera) {
			this.setPerspective(this.fov, this.aspect, this.near, this.far);
		}
		return this;
	};

	exports.B3DMLoader = B3DMLoader;
	exports.CMPTLoader = CMPTLoader;
	exports.CesiumIonAuthPlugin = CesiumIonAuthPlugin;
	exports.DebugLoadParser = LoadParser;
	exports.DebugTilesPlugin = DebugTilesPlugin;
	exports.EnvironmentControls = EnvironmentControls;
	exports.GlobeControls = GlobeControls;
	exports.I3DMLoader = I3DMLoader;
	exports.ImplicitTilingPlugin = ImplicitTilingPlugin;
	exports.OBB = OBB;
	exports.PNTSLoader = PNTSLoader;
	exports.QuantizedMeshPlugin = QuantizedMeshPlugin;
	exports.ReorientationPlugin = ReorientationPlugin;
	exports.TMSTilesPlugin = TMSTilesPlugin;
	exports.TileGLTFLoader = TileGLTFLoader;
	exports.Tiles3D = Tiles3D;
	exports.TilesFadePlugin = TilesFadePlugin;
	exports.XYZTilesPlugin = XYZTilesPlugin;

}));
