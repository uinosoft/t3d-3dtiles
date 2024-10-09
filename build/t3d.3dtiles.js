// t3d-3dtiles
(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports, require('t3d')) :
	typeof define === 'function' && define.amd ? define(['exports', 't3d'], factory) :
	(global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.t3d = global.t3d || {}, global.t3d));
})(this, (function (exports, t3d) { 'use strict';

	/**
	 * Returns the file extension of the path component of a URL
	 * @param {string} url
	 * @returns {string} null if no extension found
	 */
	const getUrlExtension = url => {
		let parsedUrl;
		try {
			parsedUrl = new URL(url, 'http://fakehost.com/');
		} catch (_) {
			// Ignore invalid URLs
			return null;
		}
		const filename = parsedUrl.pathname.split('/').pop();
		const dotIndex = filename.lastIndexOf('.');
		if (dotIndex === -1 || dotIndex === filename.length - 1) {
			// Has no extension or has trailing . character
			return null;
		}
		const extension = filename.substring(dotIndex + 1);
		return extension;
	};
	const traverseSet = (tile, beforeCb = null, afterCb = null, parent = null, depth = 0) => {
		if (beforeCb && beforeCb(tile, parent, depth)) {
			if (afterCb) {
				afterCb(tile, parent, depth);
			}
			return;
		}
		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			traverseSet(children[i], beforeCb, afterCb, tile, depth + 1);
		}
		if (afterCb) {
			afterCb(tile, parent, depth);
		}
	};

	const raycastTraverse = (tile, tiles3D, ray, intersects, localRay = null) => {
		const {
			activeTiles
		} = tiles3D;
		const boundingVolume = tile.cached.boundingVolume;

		// reuse the ray when traversing the tree
		if (localRay === null) {
			localRay = _ray_1$1;
			_mat4_1$1.copy(tiles3D.worldMatrix).inverse();
			localRay.copy(ray).applyMatrix4(_mat4_1$1);
		}
		if (!tile.__used || !boundingVolume.intersectsRay(localRay)) {
			return;
		}
		if (activeTiles.has(tile)) {
			_intersectTileScene(tile, ray, intersects);
		}
		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			raycastTraverse(children[i], tiles3D, ray, intersects, localRay);
		}
	};

	// Returns the closest hit when traversing the tree
	const raycastTraverseFirstHit = (tile, tiles3D, ray, localRay = null) => {
		const {
			activeTiles
		} = tiles3D;

		// reuse the ray when traversing the tree
		if (localRay === null) {
			localRay = _ray_1$1;
			_mat4_1$1.copy(tiles3D.worldMatrix).inverse();
			localRay.copy(ray).applyMatrix4(_mat4_1$1);
		}

		// get a set of intersections so we intersect the nearest one first
		const array = [];
		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			const child = children[i];
			if (!child.__used) {
				continue;
			}
			const boundingVolume = child.cached.boundingVolume;
			if (boundingVolume.intersectRay(localRay, _vec3_1$6)) {
				_vec3_1$6.applyMatrix4(tiles3D.worldMatrix);
				array.push({
					distance: _vec3_1$6.distanceToSquared(ray.origin),
					tile: child
				});
			}
		}

		// sort them by ascending distance
		array.sort(distanceSort);
		let bestHit = null;
		let bestHitDistSq = Infinity;

		// If the root is active make sure we've checked it
		if (activeTiles.has(tile)) {
			_intersectTileScene(tile, ray, _hitArray);
			if (_hitArray.length > 0) {
				if (_hitArray.length > 1) {
					_hitArray.sort(distanceSort);
				}
				const hit = _hitArray[0];
				_hitArray.length = 0;
				bestHit = hit;
				bestHitDistSq = hit.distance * hit.distance;
			}
		}

		// traverse until we find the best hit and early out if a tile bounds
		// couldn't possible include a best hit
		for (let i = 0, l = array.length; i < l; i++) {
			const data = array[i];
			const distanceSquared = data.distance;
			const tile = data.tile;
			if (distanceSquared > bestHitDistSq) {
				break;
			}
			const hit = raycastTraverseFirstHit(tile, tiles3D, ray, localRay);
			if (hit) {
				const hitDistSq = hit.distance * hit.distance;
				if (hitDistSq < bestHitDistSq) {
					bestHit = hit;
					bestHitDistSq = hitDistSq;
				}
			}
		}
		return bestHit;
	};
	const distanceSort = (a, b) => {
		return a.distance - b.distance;
	};
	const _intersectTileScene = (tile, ray, intersects) => {
		const scene = tile.cached.scene;
		const lengthBefore = intersects.length;
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
	};
	const _mat4_1$1 = new t3d.Matrix4();
	const _ray_1$1 = new t3d.Ray();
	const _vec3_1$6 = new t3d.Vector3();
	const _hitArray = [];

	/**
	 * State of the request.
	 *
	 * @enum {Number}
	 */
	const RequestState = {
		UNLOADED: 0,
		LOADING: 1,
		PARSING: 2,
		LOADED: 3,
		FAILED: 4
	};
	var RequestState$1 = Object.freeze(RequestState);

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
			const invScaleX = _vec3_1$5.setFromMatrixColumn(_mat4_1, 0).getLength();
			const invScaleY = _vec3_1$5.setFromMatrixColumn(_mat4_1, 1).getLength();
			const invScaleZ = _vec3_1$5.setFromMatrixColumn(_mat4_1, 2).getLength();
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
	const _vec3_1$5 = new t3d.Vector3();

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
			_mat3_1.setFromMatrix4(matrix);
			const invSX = 1 / sx;
			const invSY = 1 / sy;
			const invSZ = 1 / sz;
			_mat3_1.elements[0] *= invSX;
			_mat3_1.elements[1] *= invSX;
			_mat3_1.elements[2] *= invSX;
			_mat3_1.elements[3] *= invSY;
			_mat3_1.elements[4] *= invSY;
			_mat3_1.elements[5] *= invSY;
			_mat3_1.elements[6] *= invSZ;
			_mat3_1.elements[7] *= invSZ;
			_mat3_1.elements[8] *= invSZ;
			this.rotation.multiply(_mat3_1);
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
	const _mat3_1 = new t3d.Matrix3();
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
			this.radius = radius;
		}
		getCartographicToPosition(lat, lon, height, target) {
			// From Cesium function Ellipsoid.cartographicToCartesian
			// https://github.com/CesiumGS/cesium/blob/665ec32e813d5d6fe906ec3e87187f6c38ed5e49/packages/engine/Source/Core/Ellipsoid.js#L396
			this.getCartographicToNormal(lat, lon, _norm);
			const radius = this.radius;
			_vec3_1$2.copy(_norm);
			_vec3_1$2.x *= radius.x ** 2;
			_vec3_1$2.y *= radius.y ** 2;
			_vec3_1$2.z *= radius.z ** 2;
			const gamma = Math.sqrt(_norm.dot(_vec3_1$2));
			_vec3_1$2.multiplyScalar(1 / gamma);
			return target.copy(_vec3_1$2).addScaledVector(_norm, height);
		}
		getCartographicToNormal(lat, lon, target) {
			_spherical.set(1, latitudeToSphericalPhi(lat), lon);
			target.setFromSpherical(_spherical).normalize();

			// swap frame from the t3d.js frame to the geo coord frame
			swapToGeoFrame(target);
			return target;
		}
	}
	const _norm = new t3d.Vector3();
	const _spherical = new t3d.Spherical();
	const _vec3_1$2 = new t3d.Vector3();

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
				const midLat = mapLinear(0.5, 0, 1, latRange.x, latRange.y);
				const midLon = mapLinear(0.5, 0, 1, lonRange.x, lonRange.y);

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
			_invMatrix.setFromMatrix3(target.rotation).inverse();
			const points = this._getPoints(true);

			// get the center of the region
			_center.set(0, 0, 0);
			for (let i = 0, l = points.length; i < l; i++) {
				_center.add(points[i]);
			}
			_center.multiplyScalar(1 / points.length);
			for (let i = 0, l = points.length; i < l; i++) {
				points[i].sub(_center).applyMatrix4(_invMatrix).add(_center);
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
			const midLat = mapLinear(0.5, 0, 1, latRange.x, latRange.y);
			const midLon = mapLinear(0.5, 0, 1, lonRange.x, lonRange.y);
			const lonOffset = Math.floor(lonRange.x / HALF_PI) * HALF_PI;
			const latlon = [[-PI / 2, 0], [PI / 2, 0], [0, lonOffset], [0, lonOffset + PI / 2], [0, lonOffset + PI], [0, lonOffset + 3 * PI / 2], [latRange.x, lonRange.y], [latRange.y, lonRange.y], [latRange.x, lonRange.x], [latRange.y, lonRange.x], [0, lonRange.x], [0, lonRange.y], [midLat, midLon], [latRange.x, midLon], [latRange.y, midLon], [midLat, lonRange.x], [midLat, lonRange.y]];
			const target = [];
			const total = latlon.length;
			for (let z = 0; z <= 1; z++) {
				const height = mapLinear(z, 0, 1, heightRange.x, heightRange.y);
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
	const _center = new t3d.Vector3();
	const _invMatrix = new t3d.Matrix4();
	const PI = Math.PI;
	const HALF_PI = PI / 2;

	// Linear mapping from range <a1, a2> to range <b1, b2>
	function mapLinear(x, a1, a2, b1, b2) {
		return b1 + (x - a1) * (b2 - b1) / (a2 - a1);
	}
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
			obb.setFromCenterAndAxes(_vec3_4.set(data[0], data[1], data[2]), _vec3_1$1.set(data[3], data[4], data[5]), _vec3_2.set(data[6], data[7], data[8]), _vec3_3.set(data[9], data[10], data[11])).applyMatrix4(transform);
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
		setRegionData(west, south, east, north, minHeight, maxHeight) {
			const region = new EllipsoidRegion(WGS84_RADIUS.clone(), new t3d.Vector2(south, north), new t3d.Vector2(west, east), new t3d.Vector2(minHeight, maxHeight));
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
				if (ray.intersectSphere(sphere, _vec3_1$1)) {
					sphereDistSq = sphere.containsPoint(ray.origin) ? 0 : ray.origin.distanceToSquared(_vec3_1$1);
				}
			}
			if (obb) {
				if (obb.intersectRay(ray, _vec3_1$1)) {
					obbDistSq = obb.containsPoint(ray.origin) ? 0 : ray.origin.distanceToSquared(_vec3_1$1);
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
	const WGS84_RADIUS = new t3d.Vector3(6378137, 6378137, 6356752.3142451793);
	const _vec3_1$1 = new t3d.Vector3();
	const _vec3_2 = new t3d.Vector3();
	const _vec3_3 = new t3d.Vector3();
	const _vec3_4 = new t3d.Vector3();

	class LRUCache {
		constructor({
			maxSize = 800,
			minSize = 600,
			unloadPercent = 0.05,
			unloadPriorityCallback = defaultPriorityCallback$1
		}) {
			// options
			this.maxSize = maxSize;
			this.minSize = minSize;
			this.unloadPercent = unloadPercent;

			// "itemSet" doubles as both the list of the full set of items currently
			// stored in the cache (keys) as well as a map to the time the item was last
			// used so it can be sorted appropriately.
			this.itemSet = new Map();
			this.itemList = [];
			this.usedSet = new Set();
			this.callbacks = new Map();
			this.scheduled = false;
			this.unloadPriorityCallback = unloadPriorityCallback;
		}
		isFull() {
			return this.itemSet.size >= this.maxSize;
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
			itemList.push(item);
			usedSet.add(item);
			itemSet.set(item, Date.now());
			callbacks.set(item, removeCb);
			return true;
		}
		remove(item) {
			const usedSet = this.usedSet;
			const itemSet = this.itemSet;
			const itemList = this.itemList;
			const callbacks = this.callbacks;
			if (itemSet.has(item)) {
				callbacks.get(item)(item);
				const index = itemList.indexOf(item);
				itemList.splice(index, 1);
				usedSet.delete(item);
				itemSet.delete(item);
				callbacks.delete(item);
				return true;
			}
			return false;
		}
		markUsed(item) {
			const itemSet = this.itemSet;
			const usedSet = this.usedSet;
			if (!itemSet.has(item) || usedSet.has(item)) {
				return false;
			}
			itemSet.set(item, Date.now());
			usedSet.add(item);
			return true;
		}
		markAllUnused() {
			this.usedSet.clear();
		}
		unloadToMinSize() {
			const unloadPercent = this.unloadPercent;
			const targetSize = this.minSize;
			const itemList = this.itemList;
			const itemSet = this.itemSet;
			const usedSet = this.usedSet;
			const callbacks = this.callbacks;
			const unused = itemList.length - usedSet.size;
			const excess = itemList.length - targetSize;
			const unloadPriorityCallback = this.unloadPriorityCallback;
			if (excess <= 0 || unused <= 0) {
				return false;
			}

			// used items should be at the end of the array
			itemList.sort((a, b) => {
				const usedA = usedSet.has(a);
				const usedB = usedSet.has(b);
				if (usedA && usedB) {
					// If they're both used then don't bother moving them
					return 0;
				} else if (!usedA && !usedB) {
					// Use the sort function otherwise
					// higher priority should be further to the left
					return unloadPriorityCallback(itemSet, b) - unloadPriorityCallback(itemSet, a);
				} else {
					// If one is used and the other is not move the used one towards the end of the array
					return usedA ? 1 : -1;
				}
			});

			// address corner cases where the minSize might be zero or smaller than maxSize - minSize,
			// which would result in a very small or no items being unloaded.
			const unusedExcess = Math.min(excess, unused);
			const maxUnload = Math.max(targetSize * unloadPercent, unusedExcess * unloadPercent);
			let nodesToUnload = Math.min(maxUnload, unused);
			nodesToUnload = Math.ceil(nodesToUnload);
			const removedItems = itemList.splice(0, nodesToUnload);
			for (let i = 0, l = removedItems.length; i < l; i++) {
				const item = removedItems[i];
				callbacks.get(item)(item);
				itemSet.delete(item);
				callbacks.delete(item);
			}
			return true;
		}
		scheduleUnload(markAllUnused = true) {
			if (this.scheduled) {
				return false;
			}
			this.scheduled = true;
			enqueueMicrotask(() => {
				this.scheduled = false;
				this.unloadToMinSize();
				if (markAllUnused) {
					this.markAllUnused();
				}
			});
		}
	}
	const defaultPriorityCallback$1 = (map, key) => {
		return map.get(key);
	};
	const enqueueMicrotask = callback => {
		Promise.resolve().then(callback);
	};

	class PriorityQueue {
		constructor({
			maxJobs = 6,
			autoUpdate = true,
			priorityCallback = defaultPriorityCallback
		}) {
			// options
			this.maxJobs = maxJobs;
			this.autoUpdate = autoUpdate;
			this.items = [];
			this.callbacks = new Map();
			this.currJobs = 0;
			this.scheduled = false;
			this.priorityCallback = priorityCallback;
			this._runjobs = () => {
				this.tryRunJobs();
				this.scheduled = false;
			};
		}

		// Customizable scheduling callback. Default using requestAnimationFrame()
		schedulingCallback(func) {
			requestAnimationFrame(func);
		}
		sort() {
			const priorityCallback = this.priorityCallback;
			const items = this.items;
			items.sort(priorityCallback);
		}
		add(item, callback) {
			return new Promise((resolve, reject) => {
				const prCallback = (...args) => callback(...args).then(resolve).catch(reject);
				const items = this.items;
				const callbacks = this.callbacks;
				items.push(item);
				callbacks.set(item, prCallback);
				if (this.autoUpdate) {
					this.scheduleJobRun();
				}
			});
		}
		remove(item) {
			const items = this.items;
			const callbacks = this.callbacks;
			const index = items.indexOf(item);
			if (index !== -1) {
				items.splice(index, 1);
				callbacks.delete(item);
			}
		}
		tryRunJobs() {
			this.sort();
			const items = this.items;
			const callbacks = this.callbacks;
			const maxJobs = this.maxJobs;
			let currJobs = this.currJobs;
			while (maxJobs > currJobs && items.length > 0) {
				currJobs++;
				const item = items.pop();
				const callback = callbacks.get(item);
				callbacks.delete(item);
				callback(item).then(() => {
					this.currJobs--;
					if (this.autoUpdate) {
						this.scheduleJobRun();
					}
				}).catch(() => {
					this.currJobs--;
					if (this.autoUpdate) {
						this.scheduleJobRun();
					}
				});
			}
			this.currJobs = currJobs;
		}
		scheduleJobRun() {
			if (!this.scheduled) {
				this.schedulingCallback(this._runjobs.bind(this));
				this.scheduled = true;
			}
		}
	}
	const defaultPriorityCallback = () => {
		throw new Error('PriorityQueue: PriorityCallback function not defined.');
	};

	class TilesLoader {
		constructor() {
			this.lruCache = new LRUCache({
				unloadPriorityCallback: lruPriorityCallback
			});
			this.downloadQueue = new PriorityQueue({
				maxJobs: 4,
				priorityCallback
			});
			this.parseQueue = new PriorityQueue({
				maxJobs: 1,
				priorityCallback
			});
		}
		fetchTileSet(url, parent = null, fetchOptions = {}) {
			return fetch(url, fetchOptions).then(res => {
				if (res.ok) {
					return res.json();
				} else {
					throw new Error(`TilesLoader: Failed to load tileset "${url}" with status ${res.status} : ${res.statusText}`);
				}
			}).then(json => {
				const version = json.asset.version;
				const [major, minor] = version.split('.').map(v => parseInt(v));
				console.assert(major <= 1, 'TilesLoader: asset.version is expected to be a 1.x or a compatible version.');
				if (major === 1 && minor > 0) {
					console.warn('TilesLoader: tiles versions at 1.1 or higher have limited support. Some new extensions and features may not be supported.');
				}

				// remove trailing slash and last path-segment from the URL
				const basePath = url.replace(/\/[^/]*\/?$/, '');
				traverseSet(json.root, (node, parent) => preprocessTile(node, parent, basePath), null, parent, parent ? parent.__depth : 0);
				return json;
			});
		}
		requestTileContents(tile, tiles3D) {
			// If the tile is already being loaded then don't
			// start it again.
			if (tile.__loadingState !== RequestState$1.UNLOADED) {
				return;
			}
			const stats = tiles3D.stats;
			const lruCache = this.lruCache;
			const downloadQueue = this.downloadQueue;
			const parseQueue = this.parseQueue;
			const isExternalTileSet = tile.__externalTileSet;
			lruCache.add(tile, t => {
				if (t.__loadingState === RequestState$1.LOADING) {
					// Stop the load if it's started
					t.__loadAbort.abort();
					t.__loadAbort = null;
				} else if (isExternalTileSet) {
					t.children.length = 0;
				} else {
					tiles3D.$disposeTile(t);
				}

				// Decrement stats
				if (t.__loadingState === RequestState$1.LOADING) {
					stats.downloading--;
				} else if (t.__loadingState === RequestState$1.PARSING) {
					stats.parsing--;
				}
				t.__loadingState = RequestState$1.UNLOADED;
				t.__loadIndex++;
				downloadQueue.remove(t);
				parseQueue.remove(t);
			});

			// Track a new load index so we avoid the condition where this load is stopped and
			// another begins soon after so we don't parse twice.
			tile.__loadIndex++;
			const loadIndex = tile.__loadIndex;
			const controller = new AbortController();
			const signal = controller.signal;
			stats.downloading++;
			tile.__loadAbort = controller;
			tile.__loadingState = RequestState$1.LOADING;
			const errorCallback = e => {
				// if it has been unloaded then the tile has been disposed
				if (tile.__loadIndex !== loadIndex) {
					return;
				}
				if (e.name !== 'AbortError') {
					downloadQueue.remove(tile);
					parseQueue.remove(tile);
					if (tile.__loadingState === RequestState$1.PARSING) {
						stats.parsing--;
					} else if (tile.__loadingState === RequestState$1.LOADING) {
						stats.downloading--;
					}
					stats.failed++;
					tile.__loadingState = RequestState$1.FAILED;

					// Handle fetch bug for switching examples in index.html.
					// https://stackoverflow.com/questions/12009423/what-does-status-canceled-for-a-resource-mean-in-chrome-developer-tools
					if (e.message !== 'Failed to fetch') {
						console.error(`TilesLoader: Failed to load tile at url "${tile.content.uri}".`);
						console.error(e);
					}
				} else {
					lruCache.remove(tile);
				}
			};
			if (isExternalTileSet) {
				downloadQueue.add(tile, tileCb => {
					// if it has been unloaded then the tile has been disposed
					if (tileCb.__loadIndex !== loadIndex) {
						return Promise.resolve();
					}
					const preprocessURL = tiles3D.preprocessURL;
					const fetchOptions = tiles3D.fetchOptions;
					const uri = preprocessURL ? preprocessURL(tileCb.content.uri) : tileCb.content.uri;
					return this.fetchTileSet(uri, tileCb, Object.assign({
						signal
					}, fetchOptions));
				}).then(json => {
					// if it has been unloaded then the tile has been disposed
					if (tile.__loadIndex !== loadIndex) {
						return;
					}
					stats.downloading--;
					tile.__loadAbort = null;
					tile.__loadingState = RequestState$1.LOADED;
					tile.children.push(json.root);
				}).catch(errorCallback);
			} else {
				downloadQueue.add(tile, downloadTile => {
					if (downloadTile.__loadIndex !== loadIndex) {
						return Promise.resolve();
					}
					const preprocessURL = tiles3D.preprocessURL;
					const fetchOptions = tiles3D.fetchOptions;
					const uri = preprocessURL ? preprocessURL(downloadTile.content.uri) : downloadTile.content.uri;
					return fetch(uri, Object.assign({
						signal
					}, fetchOptions));
				}).then(res => {
					if (tile.__loadIndex !== loadIndex) {
						return;
					}
					if (res.ok) {
						const extension = getUrlExtension(res.url);
						if (extension === 'gltf') {
							return res.json();
						} else {
							return res.arrayBuffer();
						}
					} else {
						throw new Error(`Failed to load model with error code ${res.status}`);
					}
				}).then(buffer => {
					// if it has been unloaded then the tile has been disposed
					if (tile.__loadIndex !== loadIndex) {
						return;
					}
					stats.downloading--;
					stats.parsing++;
					tile.__loadAbort = null;
					tile.__loadingState = RequestState$1.PARSING;
					return parseQueue.add(tile, parseTile => {
						// if it has been unloaded then the tile has been disposed
						if (parseTile.__loadIndex !== loadIndex) {
							return Promise.resolve();
						}
						const uri = parseTile.content.uri;
						const extension = getUrlExtension(uri);
						return tiles3D.$parseTile(buffer, parseTile, extension);
					});
				}).then(() => {
					// if it has been unloaded then the tile has been disposed
					if (tile.__loadIndex !== loadIndex) {
						return;
					}
					stats.parsing--;
					tile.__loadingState = RequestState$1.LOADED;
					if (tile.__wasSetVisible) {
						tiles3D.$setTileVisible(tile, true);
					}
					if (tile.__wasSetActive) {
						tiles3D.$setTileActive(tile, true);
					}
				}).catch(errorCallback);
			}
		}
	}

	/**
	 * Function for sorting the evicted LRU items. We should evict the shallowest depth first.
	 * @param {Tile} tile
	 * @returns number
	 */
	const lruPriorityCallback = (map, tile) => 1 / (tile.__depthFromRenderedParent + 1);

	/**
	 * Function for provided to sort all tiles for prioritizing loading.
	 *
	 * @param {Tile} a
	 * @param {Tile} b
	 * @returns number
	 */
	const priorityCallback = (a, b) => {
		if (a.__depth !== b.__depth) {
			// load shallower tiles first
			return a.__depth > b.__depth ? -1 : 1;
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
	const preprocessTile = (tile, parentTile, tileSetDir) => {
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
			if (tile.content.uri) {
				// tile content uri has to be interpreted relative to the tileset.json
				// tile.content.uri = new URL( tile.content.uri, tileSetDir + '/' ).toString();
				tile.content.uri = `${tileSetDir}/${tile.content.uri}`;
			}

			// NOTE: fix for some cases where tile provide the bounding volume
			// but volumes are not present.
			if (tile.content.boundingVolume && !('box' in tile.content.boundingVolume || 'sphere' in tile.content.boundingVolume || 'region' in tile.content.boundingVolume)) {
				delete tile.content.boundingVolume;
			}
		}
		tile.parent = parentTile;
		tile.children = tile.children || [];
		const uri = tile.content && tile.content.uri;
		if (uri) {
			// "content" should only indicate loadable meshes, not external tile sets
			const extension = getUrlExtension(tile.content.uri);
			const isExternalTileSet = Boolean(extension && extension.toLowerCase() === 'json');
			tile.__externalTileSet = isExternalTileSet;
			tile.__contentEmpty = isExternalTileSet;
		} else {
			tile.__externalTileSet = false;
			tile.__contentEmpty = true;
		}

		// Expected to be set during calculateError()
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
		tile.__loadingState = RequestState$1.UNLOADED;
		tile.__loadIndex = 0;
		tile.__loadAbort = null;
		tile.__depthFromRenderedParent = -1;
		if (parentTile === null) {
			tile.__depth = 0;
			tile.refine = tile.refine || 'REPLACE';
		} else {
			tile.__depth = parentTile.__depth + 1;
			tile.refine = tile.refine || parentTile.refine;
		}

		//

		const transform = new t3d.Matrix4();
		if (tile.transform) {
			transform.fromArray(tile.transform);
		}
		if (parentTile) {
			transform.premultiply(parentTile.cached.transform);
		}
		const transformInverse = new t3d.Matrix4().copy(transform).inverse();
		const transformScale = _vec3_1.setFromMatrixScale(transform);
		const uniformScale = Math.max(transformScale.x, transformScale.y, transformScale.z);
		let geometricError = tile.geometricError * uniformScale;
		const boundingVolume = new TileBoundingVolume();
		if ('box' in tile.boundingVolume) {
			boundingVolume.setOBBData(tile.boundingVolume.box, transform);
		}
		if ('sphere' in tile.boundingVolume) {
			boundingVolume.setSphereData(tile.boundingVolume.sphere, transform);
		}
		if ('region' in tile.boundingVolume) {
			boundingVolume.setRegionData(...tile.boundingVolume.region);
			geometricError = tile.geometricError;
		}
		tile.cached = {
			loadIndex: 0,
			transform,
			transformInverse,
			geometricError,
			// geometric error applied tile transform scale

			boundingVolume,
			active: false,
			inFrustum: [],
			scene: null,
			geometry: null,
			material: null,
			featureTable: null,
			batchTable: null
		};
	};
	const _vec3_1 = new t3d.Vector3();

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
		OPAQUE: 'OPAQUE',
		MASK: 'MASK',
		BLEND: 'BLEND'
	};
	const PATH_PROPERTIES = {
		scale: 'scale',
		translation: 'position',
		rotation: 'quaternion',
		weights: 'morphTargetInfluences'
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
		TRIANGLES: 4,
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

	let MaterialParser$2 = class MaterialParser {
		static parse(context, loader) {
			const {
				gltf,
				textures
			} = context;
			if (!gltf.materials) return;
			const transformExt = loader.extensions.get('KHR_texture_transform');
			const unlitExt = loader.extensions.get('KHR_materials_unlit');
			const pbrSpecularGlossinessExt = loader.extensions.get('KHR_materials_pbrSpecularGlossiness');
			const clearcoatExt = loader.extensions.get('KHR_materials_clearcoat');
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
				const {
					KHR_materials_unlit,
					KHR_materials_pbrSpecularGlossiness,
					KHR_materials_clearcoat
				} = extensions;
				let material = null;
				if (KHR_materials_unlit && unlitExt) {
					material = unlitExt.getMaterial();
				} else if (KHR_materials_pbrSpecularGlossiness && pbrSpecularGlossinessExt) {
					material = pbrSpecularGlossinessExt.getMaterial();
					pbrSpecularGlossinessExt.parseParams(material, KHR_materials_pbrSpecularGlossiness, textures, transformExt);
				} else if (KHR_materials_clearcoat && clearcoatExt) {
					material = clearcoatExt.getMaterial();
					clearcoatExt.parseParams(material, KHR_materials_clearcoat, textures);
				} else {
					material = new t3d.PBRMaterial();
				}
				material.name = name;
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
		static parse(context) {
			const {
				gltf,
				nodes,
				accessors
			} = context;
			const {
				animations
			} = gltf;
			if (!animations) return;
			const animationClips = animations.map((gltfAnimation, index) => {
				const {
					channels,
					samplers,
					name = `animation_${index}`
				} = gltfAnimation;
				const tracks = [];
				let duration = 0;
				for (let i = 0; i < channels.length; i++) {
					const gltfChannel = channels[i];
					const gltfSampler = samplers[gltfChannel.sampler];
					if (gltfSampler) {
						const target = gltfChannel.target;
						const name = target.node !== undefined ? target.node : target.id; // Note: target.id is deprecated.
						const inputAccessor = accessors[gltfSampler.input];
						const outputAccessor = accessors[gltfSampler.output];
						const node = nodes[name];
						if (!node) continue;
						node.updateMatrix();
						node.matrixAutoUpdate = true;
						let TypedKeyframeTrack;
						switch (PATH_PROPERTIES[target.path]) {
							case PATH_PROPERTIES.rotation:
								TypedKeyframeTrack = t3d.QuaternionKeyframeTrack;
								break;
							case PATH_PROPERTIES.weights:
								TypedKeyframeTrack = t3d.NumberKeyframeTrack;
								break;
							case PATH_PROPERTIES.position:
							case PATH_PROPERTIES.scale:
							default:
								TypedKeyframeTrack = t3d.VectorKeyframeTrack;
								break;
						}
						if (!TypedKeyframeTrack) {
							continue;
						}
						const input = new inputAccessor.buffer.array.constructor(inputAccessor.buffer.array);
						const output = new Float32Array(outputAccessor.buffer.array);
						if (outputAccessor.normalized) {
							const scale = GLTFUtils.getNormalizedComponentScale(outputAccessor.buffer.array.constructor);
							for (let j = 0, jl = output.length; j < jl; j++) {
								output[j] *= scale;
							}
						}
						const targetNodes = [];
						if (PATH_PROPERTIES[target.path] === PATH_PROPERTIES.weights) {
							// Node may be a Object3D (glTF mesh with several primitives) or a Mesh.
							node.traverse(function (object) {
								if (object.isMesh && object.morphTargetInfluences) {
									targetNodes.push(object);
								}
							});
						} else {
							targetNodes.push(node);
						}
						for (let j = 0, jl = targetNodes.length; j < jl; j++) {
							const interpolant = getInterpolant(gltfSampler.interpolation, TypedKeyframeTrack === t3d.QuaternionKeyframeTrack);
							const track = new TypedKeyframeTrack(targetNodes[j], PATH_PROPERTIES[target.path], input, output, interpolant);
							tracks.push(track);
						}
						const maxTime = input[input.length - 1];
						if (duration < maxTime) {
							duration = maxTime;
						}
					}
				}
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
	 * Clearcoat Materials Extension
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_materials_clearcoat
	 */
	class KHR_materials_clearcoat {
		static getMaterial() {
			return new t3d.PBRMaterial();
		}
		static parseParams(material, extension, textures) {
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
			let offsetX = 0,
				offsetY = 0,
				repeatX = 1,
				repeatY = 1,
				rotation = 0;
			if (extDef.offset !== undefined) {
				offsetX = extDef.offset[0];
				offsetY = extDef.offset[1];
			}
			if (extDef.rotation !== undefined) {
				rotation = extDef.rotation;
			}
			if (extDef.scale !== undefined) {
				repeatX = extDef.scale[0];
				repeatY = extDef.scale[1];
			}
			const matrix = material[mapType + 'Transform'];
			if (matrix) {
				matrix.setUvTransform(offsetX, offsetY, repeatX, repeatY, rotation, 0, 0);
			}

			// If texCoord is present, it overrides the texture's texCoord
			if (extDef.texCoord !== undefined) {
				material[mapType + 'Coord'] = extDef.texCoord;
			}
		}
	}

	const DefaultParsePipeline = [IndexParser$1, ReferenceParser, Validator, BufferParser, BufferViewParser, ImageParser, TextureParser, MaterialParser$2, AccessorParser, PrimitiveParser$1, NodeParser, SkinParser, SceneParser, AnimationParser];
	const DefaultExtensions = new Map([['EXT_meshopt_compression', EXT_meshopt_compression], ['KHR_draco_mesh_compression', KHR_draco_mesh_compression], ['KHR_lights_punctual', KHR_lights_punctual], ['KHR_materials_clearcoat', KHR_materials_clearcoat], ['KHR_materials_pbrSpecularGlossiness', KHR_materials_pbrSpecularGlossiness], ['KHR_materials_unlit', KHR_materials_unlit], ['KHR_mesh_quantization', {}],
	// This is supported by default
	['KHR_texture_basisu', KHR_texture_basisu], ['KHR_texture_transform', KHR_texture_transform]]);
	class GLTFLoader {
		constructor(manager = t3d.DefaultLoadingManager, parsers = DefaultParsePipeline, extensions = DefaultExtensions) {
			this.manager = manager;

			// If ture, loading manager will dispatch progress for every buffer and image.
			// otherwise, loading manager will only dispatch progress for the whole gltf resource.
			this.detailLoadProgress = true;

			// If set false, need add Promise.catch to catch errors.
			this.autoLogError = true;
			this.extensions = new Map(extensions);
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

	let MaterialParser$1 = class MaterialParser {
		static parse(context, loader) {
			const {
				gltf,
				textures
			} = context;
			if (!gltf.materials) return;
			const transformExt = loader.extensions.get('KHR_texture_transform');
			const unlitExt = loader.extensions.get('KHR_materials_unlit');
			const pbrSpecularGlossinessExt = loader.extensions.get('KHR_materials_pbrSpecularGlossiness');
			const clearcoatExt = loader.extensions.get('KHR_materials_clearcoat');
			const techniquesExt = loader.extensions.get('KHR_techniques_webgl');
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
				const {
					KHR_materials_unlit,
					KHR_materials_pbrSpecularGlossiness,
					KHR_materials_clearcoat,
					KHR_techniques_webgl
				} = extensions;
				let material = null;
				if (KHR_materials_unlit && unlitExt) {
					material = unlitExt.getMaterial();
				} else if (KHR_materials_pbrSpecularGlossiness && pbrSpecularGlossinessExt) {
					material = pbrSpecularGlossinessExt.getMaterial();
					pbrSpecularGlossinessExt.parseParams(material, KHR_materials_pbrSpecularGlossiness, textures, transformExt);
				} else if (KHR_materials_clearcoat && clearcoatExt) {
					material = clearcoatExt.getMaterial();
					clearcoatExt.parseParams(material, KHR_materials_clearcoat, textures);
				} else if (KHR_techniques_webgl && techniquesExt) {
					// @parser-modification - add KHR_techniques_webgl
					material = techniquesExt.getMaterial();
					techniquesExt.parseParams(material, KHR_techniques_webgl, textures);
				} else {
					material = new t3d.PBRMaterial();
				}
				material.name = name;
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
			ReferenceParser, Validator, BufferParser, BufferViewParser, ImageParser, TextureParser, MaterialParser$1,
			// replace MaterialParser
			AccessorParser, PrimitiveParser$1, NodeParser, SkinParser, SceneParser, AnimationParser, B3DMRootParser // insert B3DMRootParser
			]);
			this.extensions.set('KHR_techniques_webgl', KHR_techniques_webgl);
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

	const instancingParsVert = `
		#ifdef USE_INSTANCING
				attribute mat4 instanceMatrix;
		#endif
`;
	const instancingPositionVert = `
		#ifdef USE_INSTANCING
				transformed = (instanceMatrix * vec4(transformed, 1.0)).xyz;
		#endif
`;
	const instancingNormalVert = `
		#ifdef USE_INSTANCING
				#ifdef USE_INSTANCING
						objectNormal = (transposeMat4(inverseMat4(instanceMatrix)) * vec4(objectNormal, 0.0)).xyz;
				#endif

				#ifdef USE_TANGENT
						objectTangent = (transposeMat4(inverseMat4(instanceMatrix)) * vec4(objectTangent, 0.0)).xyz;
				#endif
		#endif
`;

	class InstancedPBRMaterial extends t3d.PBRMaterial {
		constructor() {
			super();
			this.type = t3d.MATERIAL_TYPE.SHADER;
			this.shaderName = 'TILE_I_PBR';
			this.vertexShader = vertexShader$1;
			this.fragmentShader = t3d.ShaderLib.pbr_frag;
			this.defines.USE_INSTANCING = true;
		}
	}
	InstancedPBRMaterial.prototype.isInstancedPBRMaterial = true;
	let vertexShader$1 = t3d.ShaderLib.pbr_vert;
	vertexShader$1 = vertexShader$1.replace('#include <logdepthbuf_pars_vert>', `
		#include <logdepthbuf_pars_vert>
		${instancingParsVert}
`);
	vertexShader$1 = vertexShader$1.replace('#include <pvm_vert>', `
		${instancingPositionVert}
		#include <pvm_vert>
`);
	vertexShader$1 = vertexShader$1.replace('#include <normal_vert>', `
		${instancingNormalVert}
		#include <normal_vert>
`);

	class MaterialParser {
		static parse(context, loader) {
			const {
				gltf,
				textures
			} = context;
			if (!gltf.materials) return;
			const transformExt = loader.extensions.get('KHR_texture_transform');
			const unlitExt = loader.extensions.get('KHR_materials_unlit');
			const pbrSpecularGlossinessExt = loader.extensions.get('KHR_materials_pbrSpecularGlossiness');
			const clearcoatExt = loader.extensions.get('KHR_materials_clearcoat');
			const techniquesExt = loader.extensions.get('KHR_techniques_webgl');
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
				const {
					KHR_materials_unlit,
					KHR_materials_pbrSpecularGlossiness,
					KHR_materials_clearcoat,
					KHR_techniques_webgl
				} = extensions;
				let material = null;
				if (KHR_materials_unlit && unlitExt) {
					material = unlitExt.getMaterial();
				} else if (KHR_materials_pbrSpecularGlossiness && pbrSpecularGlossinessExt) {
					material = pbrSpecularGlossinessExt.getMaterial();
					pbrSpecularGlossinessExt.parseParams(material, KHR_materials_pbrSpecularGlossiness, textures, transformExt);
				} else if (KHR_materials_clearcoat && clearcoatExt) {
					material = clearcoatExt.getMaterial();
					clearcoatExt.parseParams(material, KHR_materials_clearcoat, textures);
				} else if (KHR_techniques_webgl && techniquesExt) {
					// @parser-modification - add KHR_techniques_webgl
					material = techniquesExt.getMaterial();
					techniquesExt.parseParams(material, KHR_techniques_webgl, textures);
				} else {
					material = new InstancedPBRMaterial(); // @parser-modification - instanced materials
				}
				material.name = name;
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

	class InstancedBasicMaterial extends t3d.BasicMaterial {
		constructor() {
			super();
			this.type = t3d.MATERIAL_TYPE.SHADER;
			this.shaderName = 'TILE_I_BASIC';
			this.vertexShader = vertexShader;
			this.fragmentShader = t3d.ShaderLib.basic_frag;
			this.defines.USE_INSTANCING = true;
		}
	}
	InstancedBasicMaterial.prototype.isInstancedBasicMaterial = true;
	let vertexShader = t3d.ShaderLib.basic_vert;
	vertexShader = vertexShader.replace('#include <logdepthbuf_pars_vert>', `
		#include <logdepthbuf_pars_vert>
		${instancingParsVert}
`);
	vertexShader = vertexShader.replace('#include <pvm_vert>', `
		${instancingPositionVert}
		#include <pvm_vert>
`);
	vertexShader = vertexShader.replace('#include <normal_vert>', `
		${instancingNormalVert}
		#include <normal_vert>
`);

	class KHR_materials_unlit_i {
		static getMaterial() {
			return new InstancedBasicMaterial();
		}
	}

	class KHR_materials_pbrSpecularGlossiness_i extends KHR_materials_pbrSpecularGlossiness {
		static getMaterial() {
			const material = new InstancedPBRMaterial();
			material.specular = new t3d.Color3(0x111111);
			return material;
		}
	}

	/**
	 * Clearcoat Materials Extension
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_materials_clearcoat
	 */
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
		return url.substring(url.lastIndexOf('.') + 1) === 'glb';
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

	class ModelLoader {
		constructor(manager) {
			const b3dmLoader = new B3DMLoader(manager);
			const i3dmLoader = new I3DMLoader(manager);
			const pntsLoader = new PNTSLoader(manager);
			const cmptLoader = new CMPTLoader(manager);
			const gltfLoader = new TileGLTFLoader(manager);
			this._loaders = new Map([['b3dm', b3dmLoader], ['i3dm', i3dmLoader], ['pnts', pntsLoader], ['cmpt', cmptLoader], ['gltf', gltfLoader], ['glb', gltfLoader]]);
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
	const mat4_1 = new t3d.Matrix4();
	const X_AXIS = new t3d.Vector3(1, 0, 0);
	const Y_AXIS = new t3d.Vector3(0, 1, 0);

	const schedulingTiles = (tile, tiles3D) => {
		determineFrustumSet(tile, tiles3D);
		markUsedSetLeaves(tile, tiles3D);
		skipTraversal(tile, tiles3D);
		toggleTiles(tile, tiles3D);
		tiles3D.$tilesLoader.lruCache.scheduleUnload();
	};
	const determineFrustumSet = (tile, tiles3D) => {
		const stats = tiles3D.stats;
		const frameCount = tiles3D.frameCount;
		const errorTarget = tiles3D.errorTarget;
		const maxDepth = tiles3D.maxDepth;
		const loadSiblings = tiles3D.loadSiblings;
		const lruCache = tiles3D.$tilesLoader.lruCache;
		const stopAtEmptyTiles = tiles3D.stopAtEmptyTiles;
		_resetFrameState(tile, frameCount);

		// Early out if this tile is not within view.
		const inFrustum = _tileInView(tile, tiles3D);
		if (inFrustum === false) {
			return false;
		}
		tile.__used = true;
		lruCache.markUsed(tile);
		tile.__inFrustum = true;
		stats.inFrustum++;

		// Early out if this tile has less error than we're targeting but don't stop
		// at an external tile set.
		if ((stopAtEmptyTiles || !tile.__contentEmpty) && !tile.__externalTileSet) {
			// compute the _error and __distanceFromCamera fields
			_calculateError(tile, tiles3D);
			const error = tile.__error;
			if (error <= errorTarget) {
				return true;
			}

			// Early out if we've reached the maximum allowed depth.
			if (maxDepth > 0 && tile.__depth + 1 >= maxDepth) {
				return true;
			}
		}

		// Traverse children and see if any children are in view.
		let anyChildrenUsed = false;
		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			const c = children[i];
			const r = determineFrustumSet(c, tiles3D);
			anyChildrenUsed = anyChildrenUsed || r;
		}

		// If there are children within view and we are loading siblings then mark
		// all sibling tiles as used, as well.
		if (anyChildrenUsed && loadSiblings) {
			for (let i = 0, l = children.length; i < l; i++) {
				const c = children[i];
				_recursivelyMarkUsed(c, frameCount, lruCache);
			}
		}
		return true;
	};

	// Traverse and mark the tiles that are at the leaf nodes of the "used" tree.
	const markUsedSetLeaves = (tile, tiles3D) => {
		const stats = tiles3D.stats;
		const frameCount = tiles3D.frameCount;
		if (!_isUsedThisFrame(tile, frameCount)) {
			return;
		}
		stats.used++;

		// This tile is a leaf if none of the children had been used.
		const children = tile.children;
		let anyChildrenUsed = false;
		for (let i = 0, l = children.length; i < l; i++) {
			const c = children[i];
			anyChildrenUsed = anyChildrenUsed || _isUsedThisFrame(c, frameCount);
		}
		if (!anyChildrenUsed) {
			// TODO: This isn't necessarily right because it's possible that a parent tile is considered in the
			// frustum while the child tiles are not, making them unused. If all children have loaded and were properly
			// considered to be in the used set then we shouldn't set ourselves to a leaf here.
			tile.__isLeaf = true;
		} else {
			let childrenWereVisible = false;
			let allChildrenLoaded = true;
			for (let i = 0, l = children.length; i < l; i++) {
				const c = children[i];
				markUsedSetLeaves(c, tiles3D);
				childrenWereVisible = childrenWereVisible || c.__wasSetVisible || c.__childrenWereVisible;
				if (_isUsedThisFrame(c, frameCount)) {
					const childLoaded = c.__allChildrenLoaded || !c.__contentEmpty && _isDownloadFinished(c.__loadingState) || c.__externalTileSet && c.__loadingState === RequestState$1.FAILED;
					allChildrenLoaded = allChildrenLoaded && childLoaded;
				}
			}
			tile.__childrenWereVisible = childrenWereVisible;
			tile.__allChildrenLoaded = allChildrenLoaded;
		}
	};

	// Skip past tiles we consider unrenderable because they are outside the error threshold.
	const skipTraversal = (tile, tiles3D) => {
		const stats = tiles3D.stats;
		const frameCount = tiles3D.frameCount;
		if (!_isUsedThisFrame(tile, frameCount)) {
			return;
		}
		const parent = tile.parent;
		const parentDepthToParent = parent ? parent.__depthFromRenderedParent : -1;
		tile.__depthFromRenderedParent = parentDepthToParent;

		// Request the tile contents or mark it as visible if we've found a leaf.
		const lruCache = tiles3D.$tilesLoader.lruCache;
		if (tile.__isLeaf) {
			tile.__depthFromRenderedParent++;
			if (tile.__loadingState === RequestState$1.LOADED) {
				if (tile.__inFrustum) {
					tile.__visible = true;
					stats.visible++;
				}
				tile.__active = true;
				stats.active++;
			} else if (!lruCache.isFull() && (!tile.__contentEmpty || tile.__externalTileSet)) {
				tiles3D.$tilesLoader.requestTileContents(tile, tiles3D);
			}
			return;
		}
		const errorRequirement = (tiles3D.errorTarget + 1) * tiles3D.errorThreshold;
		const meetsSSE = tile.__error <= errorRequirement;
		const includeTile = meetsSSE || tile.refine === 'ADD';
		const hasModel = !tile.__contentEmpty;
		const hasContent = hasModel || tile.__externalTileSet;
		const loadedContent = _isDownloadFinished(tile.__loadingState) && hasContent;
		const childrenWereVisible = tile.__childrenWereVisible;
		const children = tile.children;
		const allChildrenHaveContent = tile.__allChildrenLoaded;

		// Increment the relative depth of the node to the nearest rendered parent if it has content
		// and is being rendered.
		if (includeTile && hasModel) {
			tile.__depthFromRenderedParent++;
		}

		// If we've met the SSE requirements and we can load content then fire a fetch.
		if (includeTile && !loadedContent && !lruCache.isFull() && hasContent) {
			tiles3D.$tilesLoader.requestTileContents(tile, tiles3D);
		}

		// Only mark this tile as visible if it meets the screen space error requirements, has loaded content, not
		// all children have loaded yet, and if no children were visible last frame. We want to keep children visible
		// that _were_ visible to avoid a pop in level of detail as the camera moves around and parent / sibling tiles
		// load in.

		// Skip the tile entirely if there's no content to load
		if (meetsSSE && !allChildrenHaveContent && !childrenWereVisible && loadedContent || tile.refine === 'ADD' && loadedContent) {
			if (tile.__inFrustum) {
				tile.__visible = true;
				stats.visible++;
			}
			tile.__active = true;
			stats.active++;
		}

		// If we're additive then don't stop the traversal here because it doesn't matter whether the children load in
		// at the same rate.
		if (tile.refine !== 'ADD' && meetsSSE && !allChildrenHaveContent && loadedContent) {
			// load the child content if we've found that we've been loaded so we can move down to the next tile
			// layer when the data has loaded.
			for (let i = 0, l = children.length; i < l; i++) {
				const c = children[i];
				if (_isUsedThisFrame(c, frameCount) && !lruCache.isFull()) {
					c.__depthFromRenderedParent = tile.__depthFromRenderedParent + 1;
					_recursivelyLoadTiles(c, c.__depthFromRenderedParent, tiles3D);
				}
			}
		} else {
			for (let i = 0, l = children.length; i < l; i++) {
				const c = children[i];
				if (_isUsedThisFrame(c, frameCount)) {
					skipTraversal(c, tiles3D);
				}
			}
		}
	};

	// Final traverse to toggle tile visibility.
	const toggleTiles = (tile, tiles3D) => {
		const frameCount = tiles3D.frameCount;
		const isUsed = _isUsedThisFrame(tile, frameCount);
		if (isUsed || tile.__usedLastFrame) {
			let setActive = false;
			let setVisible = false;
			if (isUsed) {
				// enable visibility if active due to shadows
				setActive = tile.__active;
				if (tiles3D.displayActiveTiles) {
					setVisible = tile.__active || tile.__visible;
				} else {
					setVisible = tile.__visible;
				}
			}

			// If the active or visible state changed then call the functions.
			if (!tile.__contentEmpty && tile.__loadingState === RequestState$1.LOADED) {
				if (tile.__wasSetActive !== setActive) {
					tiles3D.$setTileActive(tile, setActive);
				}
				if (tile.__wasSetVisible !== setVisible) {
					tiles3D.$setTileVisible(tile, setVisible);
				}
			}
			tile.__wasSetActive = setActive;
			tile.__wasSetVisible = setVisible;
			tile.__usedLastFrame = isUsed;
			const children = tile.children;
			for (let i = 0, l = children.length; i < l; i++) {
				const c = children[i];
				toggleTiles(c, tiles3D);
			}
		}
	};

	// Resets the frame frame information for the given tile
	const _resetFrameState = (tile, frameCount) => {
		if (tile.__lastFrameVisited !== frameCount) {
			tile.__lastFrameVisited = frameCount;
			tile.__used = false;
			tile.__inFrustum = false;
			tile.__isLeaf = false;
			tile.__visible = false;
			tile.__active = false;
			tile.__error = Infinity;
			tile.__distanceFromCamera = Infinity;
			tile.__childrenWereVisible = false;
			tile.__allChildrenLoaded = false;
		}
	};
	const _tileInView = (tile, tiles3D) => {
		const cameraInfos = tiles3D.$cameras.getInfos();
		const cached = tile.cached;
		const boundingVolume = cached.boundingVolume;
		const inFrustum = cached.inFrustum;
		let inView = false;
		for (let i = 0, l = cameraInfos.length; i < l; i++) {
			// Track which camera frustums this tile is in so we can use it
			// to ignore the error calculations for cameras that can't see it
			const frustum = cameraInfos[i].frustum;
			if (boundingVolume.intersectsFrustum(frustum)) {
				inView = true;
				inFrustum[i] = true;
			} else {
				inFrustum[i] = false;
			}
		}
		return inView;
	};
	const _calculateError = (tile, tiles3D) => {
		const cameraInfos = tiles3D.$cameras.getInfos();
		const cached = tile.cached;
		const inFrustum = cached.inFrustum;
		const boundingVolume = cached.boundingVolume;
		let maxError = -Infinity;
		let minDistance = Infinity;
		for (let i = 0, l = cameraInfos.length; i < l; i++) {
			if (!inFrustum[i]) {
				continue;
			}

			// transform camera position into local frame of the tile bounding box
			const info = cameraInfos[i];
			const invScale = info.invScale;
			let error;
			if (info.isOrthographic) {
				const pixelSize = info.pixelSize;
				error = cached.geometricError / (pixelSize * invScale);
			} else {
				const distance = boundingVolume.distanceToPoint(info.position);
				const scaledDistance = distance * invScale;
				const sseDenominator = info.sseDenominator;
				error = cached.geometricError / (scaledDistance * sseDenominator);
				minDistance = Math.min(minDistance, scaledDistance);
			}
			maxError = Math.max(maxError, error);
		}
		tile.__distanceFromCamera = minDistance;
		tile.__error = maxError;
	};

	// Recursively mark tiles used down to the next tile with content
	const _recursivelyMarkUsed = (tile, frameCount, lruCache) => {
		_resetFrameState(tile, frameCount);
		tile.__used = true;
		lruCache.markUsed(tile);
		if (tile.__contentEmpty) {
			const children = tile.children;
			for (let i = 0, l = children.length; i < l; i++) {
				_recursivelyMarkUsed(children[i], frameCount, lruCache);
			}
		}
	};

	// Checks whether this tile was last used on the given frame.
	const _isUsedThisFrame = (tile, frameCount) => {
		return tile.__lastFrameVisited === frameCount && tile.__used;
	};
	function _isDownloadFinished(value) {
		return value === RequestState$1.LOADED || value === RequestState$1.FAILED;
	}
	const _recursivelyLoadTiles = (tile, depthFromRenderedParent, tiles3D) => {
		// Try to load any external tile set children if the external tile set has loaded.
		const doTraverse = tile.__contentEmpty && (!tile.__externalTileSet || _isDownloadFinished(tile.__loadingState));
		if (doTraverse) {
			const children = tile.children;
			for (let i = 0, l = children.length; i < l; i++) {
				// don't increment depth to rendered parent here because we should treat
				// the next layer of rendered children as just a single depth away for the
				// sake of sorting.
				const child = children[i];
				child.__depthFromRenderedParent = depthFromRenderedParent;
				_recursivelyLoadTiles(child, depthFromRenderedParent, tiles3D);
			}
		} else {
			tiles3D.$tilesLoader.requestTileContents(tile, tiles3D);
		}
	};

	class Tiles3D extends t3d.Object3D {
		constructor(url, manager = new t3d.LoadingManager()) {
			super();

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

			this._rootURL = url;
			this._rootTileSet = null;
			this._autoDisableRendererCulling = true;
			this.$cameras = new CameraList();
			this.$tilesLoader = new TilesLoader();
			this.$modelLoader = new ModelLoader(manager);
			this.$events = new t3d.EventDispatcher();
		}
		get rootURL() {
			return this._rootURL;
		}
		get rootTileSet() {
			const rootTileSet = this._rootTileSet;
			if (!rootTileSet) {
				const url = this._rootURL;
				this._rootTileSet = this.$tilesLoader.fetchTileSet(this.preprocessURL ? this.preprocessURL(url) : url, null, this.fetchOptions).then(json => {
					this._rootTileSet = json;
				}).then(json => {
					// Push this onto the end of the event stack to ensure this runs
					// after the base renderer has placed the provided json where it
					// needs to be placed and is ready for an update.
					Promise.resolve().then(() => {
						// TODO dispatch event only if this is the root tileset for now, we can
						// dispatch this event for all tilesets in the future
						_TileSetLoadedEvent.json = json;
						_TileSetLoadedEvent.url = url;
						this.$events.dispatchEvent(_TileSetLoadedEvent);
					});
				}).catch(err => {
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
		update() {
			const rootTile = this.root;
			if (!rootTile) return;
			const {
				stats
			} = this;
			stats.inFrustum = 0;
			stats.used = 0;
			stats.active = 0;
			stats.visible = 0;
			this.frameCount++;
			this.$cameras.updateInfos(this.worldMatrix);
			schedulingTiles(rootTile, this);
		}
		addCamera(camera) {
			return this.$cameras.add(camera);
		}
		removeCamera(camera) {
			return this.$cameras.remove(camera);
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
				if (tile.__loadingState === RequestState$1.FAILED) {
					tile.__loadingState = RequestState$1.UNLOADED;
				}
			});
			stats.failed = 0;
		}
		dispose() {
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
		$parseTile(buffer, tile, extension) {
			return this.$modelLoader.loadTileContent(buffer, tile, extension, this).then(scene => {
				scene.traverse(c => {
					c[INITIAL_FRUSTUM_CULLED] = c.frustumCulled; // store initial value
					c.frustumCulled = c.frustumCulled && !this._autoDisableRendererCulling;
				});
				_TileLoadedEvent.scene = scene;
				_TileLoadedEvent.tile = tile;
				this.$events.dispatchEvent(_TileLoadedEvent);
			});
		}
		$setTileVisible(tile, visible) {
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
			_TileVisibilityChangedEvent.scene = scene;
			_TileVisibilityChangedEvent.tile = tile;
			_TileVisibilityChangedEvent.visible = visible;
			this.$events.dispatchEvent(_TileVisibilityChangedEvent);
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
				_TileDisposedEvent.scene = cached.scene;
				_TileDisposedEvent.tile = tile;
				this.$events.dispatchEvent(_TileDisposedEvent);
				cached.scene = null;
				cached.materials = null;
				cached.textures = null;
				cached.geometry = null;
			}
			tile._loadIndex++;
		}
	}
	const INITIAL_FRUSTUM_CULLED = Symbol('INITIAL_FRUSTUM_CULLED');
	const _TileSetLoadedEvent = {
		type: 'TileSetLoaded',
		json: null,
		url: null
	};
	const _TileLoadedEvent = {
		type: 'TileLoaded',
		scene: null,
		tile: null
	};
	const _TileDisposedEvent = {
		type: 'TileDisposed',
		scene: null,
		tile: null
	};
	const _TileVisibilityChangedEvent = {
		type: 'TileVisibilityChanged',
		scene: null,
		tile: null,
		visible: false
	};
	const tempMat = new t3d.Matrix4();
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

	exports.B3DMLoader = B3DMLoader;
	exports.CMPTLoader = CMPTLoader;
	exports.DebugLoadParser = LoadParser;
	exports.I3DMLoader = I3DMLoader;
	exports.InstancedBasicMaterial = InstancedBasicMaterial;
	exports.InstancedPBRMaterial = InstancedPBRMaterial;
	exports.OBB = OBB;
	exports.PNTSLoader = PNTSLoader;
	exports.TileGLTFLoader = TileGLTFLoader;
	exports.Tiles3D = Tiles3D;

}));
