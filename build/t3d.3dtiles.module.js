// t3d-3dtiles
import { Vector3, Matrix3, Box3, Ray, Plane, Matrix4, Spherical, Sphere, Euler, MathUtils, Vector2, Frustum, PBRMaterial, ShaderLib, MATERIAL_TYPE, TEXEL_ENCODING_TYPE, DRAW_SIDE, Geometry, PointsMaterial, Material, BasicMaterial, VERTEX_COLOR, SHADING_TYPE, Quaternion, Attribute, Buffer, Color3, Mesh, Object3D, LoadingManager, EventDispatcher, PlaneGeometry, ShaderMaterial, DRAW_MODE, Camera } from 't3d';
import { GLTFLoader } from 't3d/addons/loaders/glTF/GLTFLoader.js';
import { ReferenceParser } from 't3d/addons/loaders/glTF/parsers/ReferenceParser.js';
import { Validator } from 't3d/addons/loaders/glTF/parsers/Validator.js';
import { BufferParser } from 't3d/addons/loaders/glTF/parsers/BufferParser.js';
import { BufferViewParser } from 't3d/addons/loaders/glTF/parsers/BufferViewParser.js';
import { ImageParser } from 't3d/addons/loaders/glTF/parsers/ImageParser.js';
import { TextureParser } from 't3d/addons/loaders/glTF/parsers/TextureParser.js';
import { MaterialParser as MaterialParser$1 } from 't3d/addons/loaders/glTF/parsers/MaterialParser.js';
import { AccessorParser } from 't3d/addons/loaders/glTF/parsers/AccessorParser.js';
import { PrimitiveParser as PrimitiveParser$1 } from 't3d/addons/loaders/glTF/parsers/PrimitiveParser.js';
import { NodeParser } from 't3d/addons/loaders/glTF/parsers/NodeParser.js';
import { SkinParser } from 't3d/addons/loaders/glTF/parsers/SkinParser.js';
import { SceneParser } from 't3d/addons/loaders/glTF/parsers/SceneParser.js';
import { AnimationParser } from 't3d/addons/loaders/glTF/parsers/AnimationParser.js';
import { GLTFUtils } from 't3d/addons/loaders/glTF/GLTFUtils.js';
import { ALPHA_MODES, ATTRIBUTES, ACCESSOR_COMPONENT_TYPES, WEBGL_DRAW_MODES } from 't3d/addons/loaders/glTF/Constants.js';
import { KHR_materials_pbrSpecularGlossiness } from 't3d/addons/loaders/glTF/extensions/KHR_materials_pbrSpecularGlossiness.js';
import { KHR_materials_clearcoat } from 't3d/addons/loaders/glTF/extensions/KHR_materials_clearcoat.js';
import { Raycaster } from 't3d/addons/Raycaster.js';
import { Box3Helper } from 't3d/addons/objects/Box3Helper.js';
import { SphereHelper } from 't3d/addons/objects/SphereHelper.js';

/**
 * An oriented bounding box.
 */
class OBB {

	/**
     * Create a new OBB.
     * @param {Box3} [box] - The axis-aligned bounding box.
     * @param {Matrix3} [rotation] - The rotation matrix.
     */
	constructor(box = new Box3(), rotation = new Matrix3()) {
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

		this.rotation.set(
			_vec3_1$4.x, _vec3_2$1.x, _vec3_3$1.x,
			_vec3_1$4.y, _vec3_2$1.y, _vec3_3$1.y,
			_vec3_1$4.z, _vec3_2$1.z, _vec3_3$1.z
		);

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

		_mat3_1$1.setFromMatrix4(matrix);

		const invSX = 1 / sx;
		const invSY = 1 / sy;
		const invSZ = 1 / sz;

		_mat3_1$1.elements[0] *= invSX;
		_mat3_1$1.elements[1] *= invSX;
		_mat3_1$1.elements[2] *= invSX;

		_mat3_1$1.elements[3] *= invSY;
		_mat3_1$1.elements[4] *= invSY;
		_mat3_1$1.elements[5] *= invSY;

		_mat3_1$1.elements[6] *= invSZ;
		_mat3_1$1.elements[7] *= invSZ;
		_mat3_1$1.elements[8] *= invSZ;

		this.rotation.multiply(_mat3_1$1);

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
					points[index]
						.set(
							x < 0 ? min.x : max.x,
							y < 0 ? min.y : max.y,
							z < 0 ? min.z : max.z
						)
						.applyMatrix3(this.rotation)
						.add(center);
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

		return Math.abs(v.dot(obbStruct.u[0])) <= obbStruct.e[0] &&
				Math.abs(v.dot(obbStruct.u[1])) <= obbStruct.e[1] &&
				Math.abs(v.dot(obbStruct.u[2])) <= obbStruct.e[2];
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
		return closestPoint.distanceToSquared(sphere.center) <= (sphere.radius * sphere.radius);
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

		transform.set(
			e[0], e[3], e[6], center.x,
			e[1], e[4], e[7], center.y,
			e[2], e[5], e[8], center.z,
			0, 0, 0, 1
		);
	}

}

const closestPoint = new Vector3();

const _vec3_1$4 = new Vector3();
const _vec3_2$1 = new Vector3();
const _vec3_3$1 = new Vector3();
const _mat3_1$1 = new Matrix3();

const R = [[], [], []];
const AbsR = [[], [], []];
const t = [];

const _halfSize = new Vector3();

const a = {
	c: new Vector3(), // center
	u: [new Vector3(), new Vector3(), new Vector3()], // basis vectors
	e: [] // half width
};

const b = {
	c: new Vector3(), // center
	u: [new Vector3(), new Vector3(), new Vector3()], // basis vectors
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
		_vec3_1$3.copy(point).applyMatrix4(this._originBoxTransformInverse);
		return this.box.containsPoint(_vec3_1$3);
	}

	intersectsRay(ray) {
		_ray_1$1.copy(ray).applyMatrix4(this._originBoxTransformInverse);
		return _ray_1$1.intersectsBox(this._originBox);
	}

	intersectRay(ray, target) {
		_ray_1$1.copy(ray).applyMatrix4(this._originBoxTransformInverse);

		if (_ray_1$1.intersectBox(this._originBox, target)) {
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

const _ray_1$1 = new Ray();
const _vec3_1$3 = new Vector3();

// Cesium / 3D tiles Spheroid:
// - Up is Z at 90 degrees latitude
// - 0, 0 latitude, longitude is X axis
//      Z
//      |
//      |
//      .----- Y
//     /
//   X


// t3d.js Spherical Coordinates
// - Up is Y at 90 degrees latitude
// - 0, 0 latitude, longitude is Z
//      Y
//      |
//      |
//      .----- X
//     /
//   Z

function swapToGeoFrame(target) {
	const { x, y, z } = target;
	target.x = z;
	target.y = x;
	target.z = y;
}

function latitudeToSphericalPhi(latitude) {
	return -latitude + Math.PI / 2;
}

class Ellipsoid {

	constructor(radius = new Vector3(1, 1, 1)) {
		this.name = '';
		this.radius = radius;
	}

	intersectRay(ray, target) {
		_matrix$1.makeScale(...this.radius.toArray([])).invert();
		_sphere$1.center.set(0, 0, 0);
		_sphere$1.radius = 1;

		_ray$3.copy(ray).applyMatrix4(_matrix$1);
		if (_ray$3.intersectSphere(_sphere$1, target)) {
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
		this.getEastNorthUpAxes(lat, lon, _vecX, _vecY, _vecZ, _pos$1);
		return target.makeBasis(_vecX, _vecY, _vecZ).setPosition(_pos$1);
	}

	getEastNorthUpAxes(lat, lon, vecEast, vecNorth, vecUp, point = _pos$1) {
		this.getCartographicToPosition(lat, lon, 0, point);
		this.getCartographicToNormal(lat, lon, vecUp);		// up
		vecEast.set(-point.y, point.x, 0).normalize();		// east
		vecNorth.crossVectors(vecUp, vecEast).normalize();	// north
	}

	getRotationMatrixFromAzElRoll(lat, lon, az, el, roll, target, frame = ENU_FRAME) {
		this.getEastNorthUpFrame(lat, lon, _matrix$1);
		_euler.set(el, roll, -az, 'ZXY');

		target
			.makeRotationFromEuler(_euler)
			.premultiply(_matrix$1)
			.setPosition(0, 0, 0);

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
		this.getCartographicToNormal(lat, lon, _norm);

		const radius = this.radius;
		_vec$4.copy(_norm);
		_vec$4.x *= radius.x ** 2;
		_vec$4.y *= radius.y ** 2;
		_vec$4.z *= radius.z ** 2;

		const gamma = Math.sqrt(_norm.dot(_vec$4));
		_vec$4.multiplyScalar(1 / gamma);

		return target.copy(_vec$4).addScaledVector(_norm, height);
	}

	getPositionToCartographic(pos, target) {
		// From Cesium function Ellipsoid.cartesianToCartographic
		// https://github.com/CesiumGS/cesium/blob/665ec32e813d5d6fe906ec3e87187f6c38ed5e49/packages/engine/Source/Core/Ellipsoid.js#L463
		this.getPositionToSurfacePoint(pos, _vec$4);
		this.getPositionToNormal(pos, _norm);

		const heightDelta = _vec2$1.subVectors(pos, _vec$4);

		target.lon = Math.atan2(_norm.y, _norm.x);
		target.lat = Math.asin(_norm.z);
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
		const invRadiusSqX = 1 / (radius.x ** 2);
		const invRadiusSqY = 1 / (radius.y ** 2);
		const invRadiusSqZ = 1 / (radius.z ** 2);

		const x2 = pos.x * pos.x * invRadiusSqX;
		const y2 = pos.y * pos.y * invRadiusSqY;
		const z2 = pos.z * pos.z * invRadiusSqZ;

		// Compute the squared ellipsoid norm.
		const squaredNorm = x2 + y2 + z2;
		const ratio = Math.sqrt(1.0 / squaredNorm);

		// As an initial approximation, assume that the radial intersection is the projection point.
		const intersection = _vec$4.copy(pos).multiplyScalar(ratio);
		if (squaredNorm < CENTER_EPS) {
			return !isFinite(ratio) ? null : target.copy(intersection);
		}

		// Use the gradient at the intersection point in place of the true unit normal.
		// The difference in magnitude will be absorbed in the multiplier.
		const gradient = _vec2$1.set(
			intersection.x * invRadiusSqX * 2.0,
			intersection.y * invRadiusSqY * 2.0,
			intersection.z * invRadiusSqZ * 2.0
		);

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
			denominator =
				x2 * xMultiplier3 * invRadiusSqX +
				y2 * yMultiplier3 * invRadiusSqY +
				z2 * zMultiplier3 * invRadiusSqZ;

			const derivative = -2 * denominator;
			correction = func / derivative;
		} while (Math.abs(func) > EPSILON12);

		return target.set(
			pos.x * xMultiplier,
			pos.y * yMultiplier,
			pos.z * zMultiplier
		);
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
		const eSquared = 1 - (semiMinorAxis ** 2 / semiMajorAxis ** 2);
		const phi = latitude * MathUtils.DEG2RAD;

		const sinPhiSquared = Math.sin(phi) ** 2;
		const N = semiMajorAxis / Math.sqrt(1 - eSquared * sinPhiSquared);
		return N;
	}

	getPositionElevation(pos) {
		// logic from "getPositionToCartographic"
		this.getPositionToSurfacePoint(pos, _vec$4);

		const heightDelta = _vec2$1.subVectors(pos, _vec$4);
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

const _spherical = new Spherical();
const _norm = new Vector3();
const _vec$4 = new Vector3();
const _vec2$1 = new Vector3();
const _matrix$1 = new Matrix4();
const _matrix2 = new Matrix4();
const _sphere$1 = new Sphere();
const _euler = new Euler();

const _vecX = new Vector3();
const _vecY = new Vector3();
const _vecZ = new Vector3();
const _pos$1 = new Vector3();

const _ray$3 = new Ray();

const EPSILON12 = 1e-12;
const CENTER_EPS = 0.1;

const ENU_FRAME = 0;
const CAMERA_FRAME = 1;
const OBJECT_FRAME = 2;

class EllipsoidRegion extends Ellipsoid {

	constructor(
		radius = new Vector3(1, 1, 1),
		latRange = new Vector2(-HALF_PI, HALF_PI),
		lonRange = new Vector2(0, 2 * PI),
		heightRange = new Vector2(0, 1)
	) {
		super(radius);

		this.latRange = latRange;
		this.lonRange = lonRange;
		this.heightRange = heightRange;
	}

	// refer to https://github.com/CesiumGS/cesium/blob/1.119/packages/engine/Source/Core/OrientedBoundingBox.js#L343
	// refer to https://github.com/NASA-AMMOS/3DTilesRendererJS/blob/c46d59c674e9ac1e652b9e9b65849bf02a645a6a/src/three/math/EllipsoidRegion.js#L120
	getOrientedBoundingBox(target) {
		resetPool();

		const { latRange, lonRange } = this;

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

		target.rotation.set(
			_orthoX.x, _orthoY.x, _orthoZ.x,
			_orthoX.y, _orthoY.y, _orthoZ.y,
			_orthoX.z, _orthoY.z, _orthoZ.z
		);

		// transform the points into the local frame
		_invMatrix$1.setFromMatrix3(target.rotation).inverse();

		const points = this._getPoints(true);

		// get the center of the region
		_center$1.set(0, 0, 0);
		for (let i = 0, l = points.length; i < l; i++) {
			_center$1.add(points[i]);
		}
		_center$1.multiplyScalar(1 / points.length);

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
		const { latRange, lonRange, heightRange } = this;

		const midLat = mapLinear(0.5, 0, 1, latRange.x, latRange.y);
		const midLon = mapLinear(0.5, 0, 1, lonRange.x, lonRange.y);

		const lonOffset = Math.floor(lonRange.x / HALF_PI) * HALF_PI;

		const latlon = [
			[-PI / 2, 0],
			[PI / 2, 0],
			[0, lonOffset],
			[0, lonOffset + PI / 2],
			[0, lonOffset + PI],
			[0, lonOffset + 3 * PI / 2],

			[latRange.x, lonRange.y],
			[latRange.y, lonRange.y],
			[latRange.x, lonRange.x],
			[latRange.y, lonRange.x],

			[0, lonRange.x],
			[0, lonRange.y],

			[midLat, midLon],
			[latRange.x, midLon],
			[latRange.y, midLon],
			[midLat, lonRange.x],
			[midLat, lonRange.y]
		];

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

const _orthoX = new Vector3();
const _orthoY = new Vector3();
const _orthoZ = new Vector3();
const _center$1 = new Vector3();
const _invMatrix$1 = new Matrix4();

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
		return new Vector3();
	}

	if (!_pointsPool[_poolIndex]) {
		_pointsPool[_poolIndex] = new Vector3();
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

		obb.setFromCenterAndAxes(
			_vec3_4.set(data[0], data[1], data[2]),
			_vec3_1$2.set(data[3], data[4], data[5]),
			_vec3_2.set(data[6], data[7], data[8]),
			_vec3_3.set(data[9], data[10], data[11])
		).applyMatrix4(transform);
		obb.updateCache();

		this.obb = obb;
	}

	setSphereData(data, transform) {
		const sphere = new Sphere();

		sphere.center.set(data[0], data[1], data[2]);
		sphere.radius = data[3];
		sphere.applyMatrix4(transform);

		this.sphere = sphere;
	}

	setRegionData(ellipsoid, west, south, east, north, minHeight, maxHeight) {
		const region = new EllipsoidRegion(
			ellipsoid.radius.clone(),
			new Vector2(south, north),
			new Vector2(west, east),
			new Vector2(minHeight, maxHeight)
		);
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

const _vec3_1$2 = new Vector3();
const _vec3_2 = new Vector3();
const _vec3_3 = new Vector3();
const _vec3_4 = new Vector3();

/**
 * Returns the file extension of the path component of a URL
 * @param {string} url
 * @returns {string} null if no extension found
 */
const getUrlExtension = url => {
	if (!url) {
		return null;
	}

	const filename = url
		.replace(/[a-z]+:\/\/[^/]+/i, '') 	// remove origin
		.replace(/\?.*$/i, '') 				// remove query
		.replace(/.*\//g, ''); 				// remove path

	const lastPeriod = filename.lastIndexOf('.');
	if (lastPeriod === -1) {
		return null;
	}

	return filename.substring(lastPeriod + 1) || null;
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

const raycastTraverse = (tile, tiles3D, ray, intersects, localRay = null) => {
	const { activeTiles } = tiles3D;
	const boundingVolume = tile.cached.boundingVolume;

	// reuse the ray when traversing the tree
	if (localRay === null) {
		localRay = _ray_1;
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
	const { activeTiles } = tiles3D;

	// reuse the ray when traversing the tree
	if (localRay === null) {
		localRay = _ray_1;
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

		if (boundingVolume.intersectRay(localRay, _vec3_1$1)) {
			_vec3_1$1.applyMatrix4(tiles3D.worldMatrix);
			array.push({
				distance: _vec3_1$1.distanceToSquared(ray.origin),
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

const _mat4_1$1 = new Matrix4();
const _ray_1 = new Ray();
const _vec3_1$1 = new Vector3();

const _hitArray = [];

class FastFrustum extends Frustum {

	constructor() {
		super();

		this.points = new Array(8).fill().map(() => new Vector3());
	}

	updateCache() {
		const { planes, points } = this;
		const planeIntersections = [
			[planes[0], planes[3], planes[4]], // Near top left
			[planes[1], planes[3], planes[4]], // Near top right
			[planes[0], planes[2], planes[4]], // Near bottom left
			[planes[1], planes[2], planes[4]], // Near bottom right
			[planes[0], planes[3], planes[5]], // Far top left
			[planes[1], planes[3], planes[5]], // Far top right
			[planes[0], planes[2], planes[5]], // Far bottom left
			[planes[1], planes[2], planes[5]] // Far bottom right
		];

		planeIntersections.forEach((planes, index) => {
			findIntersectionPoint(planes[0], planes[1], planes[2], points[index]);
		});
	}

}

const _mat3_1 = new Matrix3();

// Solve a system of equations to find the point where the three planes intersect
function findIntersectionPoint(plane1, plane2, plane3, target) {
	// Create the matrix A using the normals of the planes as rows
	const A = _mat3_1.set(
		plane1.normal.x, plane1.normal.y, plane1.normal.z,
		plane2.normal.x, plane2.normal.y, plane2.normal.z,
		plane3.normal.x, plane3.normal.y, plane3.normal.z
	);

	// Create the vector B using the constants of the planes
	target.set(-plane1.constant, -plane2.constant, -plane3.constant);

	// Solve for X by applying the inverse matrix to B
	target.applyMatrix3(A.inverse());

	return target;
}

class CameraList {

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
		const { itemSet, loadedSet } = this;
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
			while (
				this.cachedBytes - removedBytes > maxBytesSize ||
				itemList.length - removedNodes > maxSize
			) {
				const item = itemList[removedNodes];
				const bytes = bytesMap.get(item) || 0;
				if (
					usedSet.has(item) && loadedSet.has(item) ||
					this.cachedBytes - removedBytes - bytes < maxBytesSize &&
					itemList.length - removedNodes <= maxSize
				) {
					break;
				}

				removedBytes += bytes;
				removedNodes++;
			}

			// evict up to the min node or bytes size, keeping one more item over the min bytes limit
			// so we're meeting it
			while (
				removedBytes < bytesToUnload ||
				removedNodes < nodesToUnload
			) {
				const item = itemList[removedNodes];
				const bytes = bytesMap.get(item) || 0;
				if (
					usedSet.has(item) ||
					this.cachedBytes - removedBytes - bytes < minBytesSize &&
					removedNodes >= nodesToUnload
				) {
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
			const { callback, resolve, reject } = callbacks.get(item);
			callbacks.delete(item);

			let result;
			try {
				result = callback(item);
			} catch (err) {
				reject(err);
				completedCallback();
			}

			if (result instanceof Promise) {
				result
					.then(resolve)
					.catch(reject)
					.finally(completedCallback);
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

		const magic =
            String.fromCharCode(dataView.getUint8(0)) +
            String.fromCharCode(dataView.getUint8(1)) +
            String.fromCharCode(dataView.getUint8(2)) +
            String.fromCharCode(dataView.getUint8(3));

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
			const { buffer, binOffset, binLength } = this;
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
		const { header, options } = context;
		const buffer = options.buffer;

		const featureTableStart = header.magic === 'i3dm' ? 32 : 28;
		const featureTableEnd = featureTableStart + header.featureTableJSONByteLength + header.featureTableBinaryByteLength;
		const batchTableStart = featureTableEnd;
		const batchTableEnd = batchTableStart + header.batchTableJSONByteLength + header.batchTableBinaryByteLength;

		// parse the feature table

		const featureTableBuffer = buffer.slice(
			featureTableStart,
			featureTableEnd
		);

		const featureTable = new FeatureTable(
			featureTableBuffer,
			0,
			header.featureTableJSONByteLength,
			header.featureTableBinaryByteLength
		);

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

		const batchTableBuffer = buffer.slice(
			batchTableStart,
			batchTableEnd
		);

		const batchTable = new BatchTable(
			batchTableBuffer,
			batchSize,
			0,
			header.batchTableJSONByteLength,
			header.batchTableBinaryByteLength
		);

		// output the tables to the context

		context.featureTable = featureTable;
		context.batchTable = batchTable;
		context.batchTableEnd = batchTableEnd;
	}

}

class B3DMParser {

	static parse(context, loader) {
		const glbBytes = new Uint8Array(
			context.options.buffer,
			context.batchTableEnd,
			context.header.byteLength - context.batchTableEnd
		);

		const glbData = GLTFUtils.parseGLB(glbBytes.slice().buffer);

		context.gltf = glbData.gltf;
		context.buffers = glbData.buffers;
	}

}

class B3DMRootParser {

	static parse(context, loader) {
		const { root, featureTable, options } = context;

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
		return new PBRMaterial();
	}

	static parseParams(material, extension, textures) {
		const { values } = extension;
		const { u_diffuse } = values;

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
		super(manager, [
			HeaderParser, // insert HeaderParser
			TableParser, // insert TableParser
			B3DMParser, // insert B3DMParser
			ReferenceParser,
			Validator,
			BufferParser,
			BufferViewParser,
			ImageParser,
			TextureParser,
			MaterialParser$1,
			AccessorParser,
			PrimitiveParser$1,
			NodeParser,
			SkinParser,
			SceneParser,
			AnimationParser,
			B3DMRootParser // insert B3DMRootParser
		]);

		this.extensions.set('KHR_techniques_webgl', KHR_techniques_webgl);

		this.autoParseConfig.materials.push('KHR_techniques_webgl');
	}

}

class I3DMParser {

	static parse(context, loader) {
		const bodyBytes = new Uint8Array(
			context.options.buffer,
			context.batchTableEnd,
			context.header.byteLength - context.batchTableEnd
		);

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

class InstancedPBRMaterial extends PBRMaterial {

	constructor() {
		super();
		this.type = MATERIAL_TYPE.SHADER;
		this.shaderName = 'TILE_I_PBR';
		this.vertexShader = vertexShader$1;
		this.fragmentShader = ShaderLib.pbr_frag;
		this.defines.USE_INSTANCING = true;
	}

}

InstancedPBRMaterial.prototype.isInstancedPBRMaterial = true;

let vertexShader$1 = ShaderLib.pbr_vert;

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
		const { gltf, textures } = context;

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

			const { KHR_materials_unlit, KHR_materials_pbrSpecularGlossiness } = extensions;

			if (pbrMetallicRoughness) {
				const { baseColorFactor, baseColorTexture, metallicFactor, roughnessFactor, metallicRoughnessTexture } = pbrMetallicRoughness;

				if (Array.isArray(baseColorFactor)) {
					material.diffuse.fromArray(baseColorFactor);
					material.opacity = (baseColorFactor[3] !== undefined) ? baseColorFactor[3] : 1;
				}

				if (baseColorTexture) {
					material.diffuseMap = textures[baseColorTexture.index];
					material.diffuseMapCoord = baseColorTexture.texCoord || 0;
					if (material.diffuseMap) {
						material.diffuseMap.encoding = TEXEL_ENCODING_TYPE.SRGB;
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
					material.emissiveMap.encoding = TEXEL_ENCODING_TYPE.SRGB;
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

			material.side = doubleSided === true ? DRAW_SIDE.DOUBLE : DRAW_SIDE.FRONT;

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
		const { gltf, accessors, materials, bufferViews } = context;

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
				const { KHR_draco_mesh_compression } = extensions;

				let geometryPromise;

				const geometryKey = createGeometryKey(gltfPrimitive);
				if (geometryPromiseCache.has(geometryKey)) {
					geometryPromise = geometryPromiseCache.get(geometryKey);
				} else {
					if (KHR_draco_mesh_compression && dracoExt) {
						geometryPromise = dracoExt.getGeometry(KHR_draco_mesh_compression, bufferViews, gltfPrimitive.attributes, gltf.accessors, loader.getDRACOLoader());
					} else {
						geometryPromise = Promise.resolve(new Geometry());
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
						material: material === undefined ? new InstancedPBRMaterial() : materials[material], // @parser-modification - instanced materials
						weights: (Object.keys(geometry.morphAttributes).length > 0 && gltfMesh.weights) ? gltfMesh.weights.slice(0) : undefined,
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
	const { attributes, indices, targets } = gltfPrimitive;

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

	const { boundingBox, boundingSphere } = geometry;
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
	let { geometry, material, skinned, mode } = primitive;

	// If the material will be modified later on, clone it now.
	const useVertexTangents = geometry.attributes[ATTRIBUTES.TANGENT] !== undefined;
	const useVertexColors = geometry.attributes[ATTRIBUTES.COLOR_0] !== undefined;
	const useFlatShading = geometry.attributes[ATTRIBUTES.NORMAL] === undefined;
	const useSkinning = skinned;

	if (mode === WEBGL_DRAW_MODES.POINTS) {
		const cacheKey = 'PointsMaterial:' + material.id;
		let pointsMaterial = materialCache.get(cacheKey);
		if (!pointsMaterial) {
			pointsMaterial = new PointsMaterial();
			Material.prototype.copy.call(pointsMaterial, material);
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
			basicMaterial = new BasicMaterial();
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
					cachedMaterial.vertexColors = VERTEX_COLOR.RGB;
				} else if (geometry.attributes[ATTRIBUTES.COLOR_0].size === 4) {
					cachedMaterial.vertexColors = VERTEX_COLOR.RGBA;
				} else {
					console.warn('Illegal vertex color size: ' + geometry.attributes[ATTRIBUTES.COLOR_0].size);
				}
			}

			if (useFlatShading) {
				cachedMaterial.shading = SHADING_TYPE.FLAT_SHADING;
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
		geometryKey = 'draco:' + dracoExtension.bufferView
				+ ':' + dracoExtension.indices
				+ ':' + createAttributesKey(dracoExtension.attributes);
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
		const { featureTable, root, options } = context;

		const INSTANCES_LENGTH = featureTable.getData('INSTANCES_LENGTH');
		const POSITION = featureTable.getData('POSITION', INSTANCES_LENGTH, 'FLOAT', 'VEC3');
		const NORMAL_UP = featureTable.getData('NORMAL_UP', INSTANCES_LENGTH, 'FLOAT', 'VEC3');
		const NORMAL_RIGHT = featureTable.getData('NORMAL_RIGHT', INSTANCES_LENGTH, 'FLOAT', 'VEC3');
		const SCALE = featureTable.getData('SCALE', INSTANCES_LENGTH, 'FLOAT', 'SCALAR');
		const SCALE_NON_UNIFORM = featureTable.getData('SCALE_NON_UNIFORM', INSTANCES_LENGTH, 'FLOAT', 'VEC3');

		// check unsupported features

		[
			// Global Properties
			'QUANTIZED_VOLUME_OFFSET',
			'QUANTIZED_VOLUME_SCALE',
			'EAST_NORTH_UP',

			// Per-instance Properties
			'POSITION_QUANTIZED',
			'NORMAL_UP_OCT32P',
			'NORMAL_RIGHT_OCT32P'
		].forEach(feature => {
			if (feature in featureTable.header) {
				console.warn(`I3DMLoader: Unsupported FeatureTable feature "${feature}" detected.`);
			}
		});

		// set instance matrix for all geometries

		const averageVector = new Vector3();
		for (let i = 0; i < INSTANCES_LENGTH; i++) {
			averageVector.x += POSITION[i * 3 + 0] / INSTANCES_LENGTH;
			averageVector.y += POSITION[i * 3 + 1] / INSTANCES_LENGTH;
			averageVector.z += POSITION[i * 3 + 2] / INSTANCES_LENGTH;
		}

		const instances = [];

		root.traverse(child => {
			if (child.isMesh) {
				const { geometry } = child;
				geometry.instanceCount = INSTANCES_LENGTH;

				const instanceMatrix = new Attribute(new Buffer(new Float32Array(INSTANCES_LENGTH * 16), 16), 16);
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

				tempMat$1.set(
					tempRight.x, tempUp.x, tempFwd.x, 0,
					tempRight.y, tempUp.y, tempFwd.y, 0,
					tempRight.z, tempUp.z, tempFwd.z, 0,
					0, 0, 0, 1
				);

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
				const { geometry } = instances[j];
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

const tempFwd = new Vector3();
const tempUp = new Vector3();
const tempRight = new Vector3();
const tempPos = new Vector3();
const tempQuat = new Quaternion();
const tempSca = new Vector3();
const tempMat$1 = new Matrix4();

class KHR_techniques_webgl_i extends KHR_techniques_webgl {

	static getMaterial() {
		return new InstancedPBRMaterial();
	}

}

class InstancedBasicMaterial extends BasicMaterial {

	constructor() {
		super();
		this.type = MATERIAL_TYPE.SHADER;
		this.shaderName = 'TILE_I_BASIC';
		this.vertexShader = vertexShader;
		this.fragmentShader = ShaderLib.basic_frag;
		this.defines.USE_INSTANCING = true;
	}

}

InstancedBasicMaterial.prototype.isInstancedBasicMaterial = true;

let vertexShader = ShaderLib.basic_vert;

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
		const material = new InstancedPBRMaterial(); // fallback to InstancedPBRMaterial
		material.specular = new Color3(0x111111);
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
		super(manager, [
			HeaderParser, // insert HeaderParser
			TableParser, // insert TableParser
			I3DMParser, // insert I3DMParser
			ReferenceParser,
			Validator,
			BufferParser,
			BufferViewParser,
			ImageParser,
			TextureParser,
			MaterialParser, // replace MaterialParser
			AccessorParser,
			PrimitiveParser, // replace PrimitiveParser
			NodeParser,
			SkinParser,
			SceneParser,
			AnimationParser,
			I3DMRootParser // insert I3DMSetupParser
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
		const { featureTable } = context;

		const POINTS_LENGTH = featureTable.getData('POINTS_LENGTH');
		const POSITION = featureTable.getData('POSITION', POINTS_LENGTH, 'FLOAT', 'VEC3');
		const RGB = featureTable.getData('RGB', POINTS_LENGTH, 'UNSIGNED_BYTE', 'VEC3');
		const RGBA = featureTable.getData('RGBA', POINTS_LENGTH, 'UNSIGNED_BYTE', 'VEC4');

		// check unsupported features

		[
			// Global Properties
			'QUANTIZED_VOLUME_OFFSET',
			'QUANTIZED_VOLUME_SCALE',
			'CONSTANT_RGBA',
			'BATCH_LENGTH',

			// Per-point Properties
			'POSITION_QUANTIZED',
			'RGB565',
			'NORMAL',
			'NORMAL_OCT16P',
			'BATCH_ID'
		].forEach(feature => {
			if (feature in featureTable.header) {
				console.warn(`PNTSLoader: Unsupported FeatureTable feature "${feature}" detected.`);
			}
		});

		// generate root

		const geometry = new Geometry();
		geometry.addAttribute('a_Position', new Attribute(new Buffer(POSITION, 3), 3, 0, true));
		geometry.computeBoundingBox();
		geometry.computeBoundingSphere();

		const material = new PointsMaterial();
		material.size = 2;
		material.sizeAttenuation = false;

		if (RGB !== null) {
			geometry.addAttribute('a_Color', new Attribute(new Buffer(RGB, 3), 3, 0, true));
			material.vertexColors = VERTEX_COLOR.RGB;
		} else if (RGBA !== null) {
			geometry.addAttribute('a_Color', new Attribute(new Buffer(RGBA, 4), 4, 0, true));
			material.vertexColors = VERTEX_COLOR.RGBA;
		}

		const root = new Mesh(geometry, material);

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
		super(manager, [
			HeaderParser,
			TableParser,
			PNTSRootParser
		]);
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

			const tileMagic =
				String.fromCharCode(tileView.getUint8(0)) +
				String.fromCharCode(tileView.getUint8(1)) +
				String.fromCharCode(tileView.getUint8(2)) +
				String.fromCharCode(tileView.getUint8(3));

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
		const { tiles, options, path } = context;

		const adjustmentTransform = options.adjustmentTransform;

		const promises = [];

		for (const i in tiles) {
			const { type, buffer } = tiles[i];

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
			const group = new Object3D();

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
		super(manager, [
			HeaderParser,
			CMPTParser,
			CMPTRootParser
		]);

		const b3dmLoader = new B3DMLoader(manager);
		const i3dmLoader = new I3DMLoader(manager);
		const pntsLoader = new PNTSLoader(manager);

		this._loaders = new Map([
			['b3dm', b3dmLoader],
			['i3dm', i3dmLoader],
			['pnts', pntsLoader]
		]);
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
		const { url, options } = context;
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

class ModelLoader {

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
	// return tile.__childrenProcessed === tile.children.length;
	return true; // TODO: implement this
}

// Resets the frame frame information for the given tile
const resetFrameState = (tile, renderer) => {
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
};

// Recursively mark tiles used down to the next tile with content
const recursivelyMarkUsed = (tile, renderer) => {
	renderer.ensureChildrenArePreprocessed(tile);

	resetFrameState(tile, renderer);
	markUsed(tile, renderer);

	// don't traverse if the children have not been processed, yet
	if (!tile.__hasRenderableContent && areChildrenProcessed()) {
		const children = tile.children;
		for (let i = 0, l = children.length; i < l; i++) {
			recursivelyMarkUsed(children[i], renderer);
		}
	}
};

// Recursively traverses to the next tiles with unloaded renderable content to load them
function recursivelyLoadNextRenderableTiles(tile, renderer) {
	renderer.ensureChildrenArePreprocessed(tile);

	// exit the recursion if the tile hasn't been used this frame
	if (isUsedThisFrame(tile, renderer.frameCount)) {
		// queue this tile to download content
		if (tile.__hasContent && tile.__loadingState === UNLOADED && !renderer.lruCache.isFull()) {
			renderer.requestTileContents(tile);
		}

		{
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

	return true;
}

const markUsedTiles = (tile, renderer) => {
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
};

// Traverse and mark the tiles that are at the leaf nodes of the "used" tree.
const markUsedSetLeaves = (tile, renderer) => {
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
				const childLoaded =
                    c.__allChildrenLoaded ||
                    (c.__hasRenderableContent && isDownloadFinished(c.__loadingState)) ||
                    (c.__hasUnrenderableContent && c.__loadingState === FAILED);
				allChildrenLoaded = allChildrenLoaded && childLoaded;
			}
		}
		tile.__childrenWereVisible = childrenWereVisible;
		tile.__allChildrenLoaded = allChildrenLoaded;
	}
};

// Skip past tiles we consider unrenderable because they are outside the error threshold.
const markVisibleTiles = (tile, renderer) => {
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
			renderer.requestTileContents(tile);
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
		renderer.requestTileContents(tile);
	}

	// Only mark this tile as visible if it meets the screen space error requirements, has loaded content, not
	// all children have loaded yet, and if no children were visible last frame. We want to keep children visible
	// that _were_ visible to avoid a pop in level of detail as the camera moves around and parent / sibling tiles
	// load in.

	// Skip the tile entirely if there's no content to load
	if (
		(meetsSSE && !allChildrenLoaded && !childrenWereVisible && loadedContent)
        || (tile.refine === 'ADD' && loadedContent)
	) {
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
};

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

const WGS84_ELLIPSOID = new Ellipsoid(new Vector3(WGS84_RADIUS, WGS84_RADIUS, WGS84_HEIGHT));
WGS84_ELLIPSOID.name = 'WGS84 Earth';

const _updateBeforeEvent = { type: 'update-before' };
const _updateAfterEvent = { type: 'update-after' };

const PLUGIN_REGISTERED = Symbol('PLUGIN_REGISTERED');

const INITIAL_FRUSTUM_CULLED = Symbol('INITIAL_FRUSTUM_CULLED');

const tempMat = new Matrix4();
const viewErrorTarget = {
	inView: false,
	error: Infinity
};

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

class Tiles3D extends Object3D {

	get root() {
		const rootTileSet = this.rootTileSet;
		return rootTileSet ? rootTileSet.root : null;
	}

	constructor(url, manager = new LoadingManager()) {
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
		this.usedSet = new Set();

		this.rootLoadingState = UNLOADED;

		// internals

		this.rootURL = url;
		this.rootTileSet = null;

		this._autoDisableRendererCulling = true;

		this.plugins = [];

		const lruCache = new LRUCache();
		lruCache.unloadPriorityCallback = lruPriorityCallback;

		const downloadQueue = new PriorityQueue();
		downloadQueue.maxJobs = 10;
		downloadQueue.priorityCallback = priorityCallback;

		const parseQueue = new PriorityQueue();
		parseQueue.maxJobs = 1;
		parseQueue.priorityCallback = priorityCallback;

		this.lruCache = lruCache;
		this.downloadQueue = downloadQueue;
		this.parseQueue = parseQueue;

		this.$cameras = new CameraList();
		this.$modelLoader = new ModelLoader(manager);
		this.$events = new EventDispatcher();

		this.lruCache.computeMemoryUsageCallback = tile => tile.cached.bytesUsed ?? null;
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

	markTileUsed(tile) {
		// save the tile in a separate "used set" so we can mark it as unused
		// before the next tile set traversal
		this.usedSet.add(tile);
		this.lruCache.markUsed(tile);
	}

	// Public API
	update() {
		const { lruCache, usedSet, stats, root } = this;

		if (this.rootLoadingState === UNLOADED) {
			this.rootLoadingState = LOADING;
			this.invokeOnePlugin(plugin => plugin.loadRootTileSet && plugin.loadRootTileSet())
				.then(root => {
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
				})
				.catch(error => {
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

		this.lruCache.scheduleUnload();

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

	dispatchEvent(event) {
		this.$events.dispatchEvent(event);
	}

	fetchData(url, options) {
		return fetch(url, options);
	}

	parseTile(buffer, tile, extension) {
		return this.$modelLoader.loadTileContent(buffer, tile, extension, this)
			.then(scene => {
				scene.traverse(c => {
					c[INITIAL_FRUSTUM_CULLED] = c.frustumCulled; // store initial value
					c.frustumCulled = c.frustumCulled && !this._autoDisableRendererCulling;
				});

				this.dispatchEvent({
					type: 'load-model',
					scene,
					tile
				});
			});
	}

	disposeTile(tile) {
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

		tile._loadIndex++;
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
			if (
				tile.content.boundingVolume &&
			!(
				'box' in tile.content.boundingVolume ||
				'sphere' in tile.content.boundingVolume ||
				'region' in tile.content.boundingVolume
			)
			) {
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
			tile.__depthFromRenderedParent = (tile.__hasRenderableContent ? 1 : 0);
			tile.refine = tile.refine || 'REPLACE';
		} else {
			// increment the "depth from parent" when we encounter a new tile with content
			tile.__depth = parentTile.__depth + 1;
			tile.__depthFromRenderedParent = parentTile.__depthFromRenderedParent + (tile.__hasRenderableContent ? 1 : 0);
			tile.refine = tile.refine || parentTile.refine;
		}

		tile.__basePath = tileSetDir;

		tile.__lastFrameVisited = -1;

		tile.__loadIndex = 0; // TODO remove this
		tile.__loadAbort = null; // TODO remove this

		this.invokeAllPlugins(plugin => {
			plugin !== this && plugin.preprocessNode && plugin.preprocessNode(tile, tileSetDir, parentTile);
		});

		// cached

		const transform = new Matrix4();
		if (tile.transform) {
			transform.fromArray(tile.transform);
		}

		if (parentTile) {
			transform.premultiply(parentTile.cached.transform);
		}

		const transformInverse = new Matrix4().copy(transform).inverse();
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

			// TODO remove this

			inFrustum: [],

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
		// TODO
	}

	// Private Functions
	preprocessTileSet(json, url, parent = null) {
		const version = json.asset.version;
		const [major, minor] = version.split('.').map(v => parseInt(v));
		console.assert(
			major <= 1,
			'Tiles3D: asset.version is expected to be a 1.x or a compatible version.'
		);

		if (major === 1 && minor > 0) {
			console.warn('Tiles3D: tiles versions at 1.1 or higher have limited support. Some new extensions and features may not be supported.');
		}

		// remove the last file path path-segment from the URL including the trailing slash
		let basePath = url.replace(/\/[^/]*$/, '');
		basePath = new URL(basePath, window.location.href).toString();

		traverseSet(
			json.root,
			(node, parent) => this.preprocessNode(node, basePath, parent),
			null,
			parent,
			parent ? parent.__depth : 0
		);
	}

	loadRootTileSet() {
		// transform the url
		let processedUrl = this.rootURL;
		this.invokeAllPlugins(plugin => processedUrl = plugin.preprocessURL ? plugin.preprocessURL(processedUrl, null) : processedUrl);

		// load the tile set root
		const pr = this
			.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(processedUrl, this.fetchOptions))
			.then(res => {
				if (res.ok) {
					return res.json();
				} else {
					throw new Error(`Tiles3D: Failed to load tileset "${processedUrl}" with status ${res.status} : ${res.statusText}`);
				}
			})
			.then(root => {
				const { extensions = {} } = root;

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

		const stats = this.stats;
		const lruCache = this.lruCache;
		const downloadQueue = this.downloadQueue;
		const parseQueue = this.parseQueue;

		const isExternalTileSet = tile.__hasUnrenderableContent;

		lruCache.add(tile, t => {
			if (t.__loadingState === LOADING) {
				// Stop the load if it's started
				t.__loadAbort.abort();
				t.__loadAbort = null;
			} else if (isExternalTileSet) {
				t.children.length = 0;
			} else {
				this.disposeTile(t);
			}

			// Decrement stats
			if (t.__loadingState === LOADING) {
				stats.downloading--;
			} else if (t.__loadingState === PARSING) {
				stats.parsing--;
			}

			t.__loadingState = UNLOADED;
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
		tile.__loadingState = LOADING;

		const errorCallback = e => {
			// if it has been unloaded then the tile has been disposed
			if (tile.__loadIndex !== loadIndex) {
				return;
			}

			if (e.name !== 'AbortError') {
				downloadQueue.remove(tile);
				parseQueue.remove(tile);

				if (tile.__loadingState === PARSING) {
					stats.parsing--;
				} else if (tile.__loadingState === LOADING) {
					stats.downloading--;
				}

				stats.failed++;
				tile.__loadingState = FAILED;

				// Handle fetch bug for switching examples in index.html.
				// https://stackoverflow.com/questions/12009423/what-does-status-canceled-for-a-resource-mean-in-chrome-developer-tools
				if (e.message !== 'Failed to fetch') {
					console.error(`Tiles3D: Failed to load tile at url "${tile.content.uri}".`);
					console.error(e);
				}
			} else {
				lruCache.remove(tile);
			}
		};

		let uri = new URL(tile.content.uri, tile.__basePath + '/').toString();
		this.invokeAllPlugins(plugin => uri = plugin.preprocessURL ? plugin.preprocessURL(uri, tile) : uri);

		if (isExternalTileSet) {
			downloadQueue.add(tile, tileCb => {
				// if it has been unloaded then the tile has been disposed
				if (tileCb.__loadIndex !== loadIndex) {
					return Promise.resolve();
				}

				return this.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(uri, { ...this.fetchOptions, signal }));
			}).then(res => {
				if (res.ok) {
					return res.json();
				} else {
					throw new Error(`Tiles3D: Failed to load tileset "${uri}" with status ${res.status} : ${res.statusText}`);
				}
			}).then(json => {
				this.preprocessTileSet(json, uri, tile);
				return json;
			}).then(json => {
				// if it has been unloaded then the tile has been disposed
				if (tile.__loadIndex !== loadIndex) {
					return;
				}

				stats.downloading--;
				tile.__loadAbort = null;
				tile.__loadingState = LOADED;

				tile.children.push(json.root);

				this.dispatchEvent({
					type: 'load-tile-set',
					tileSet: json,
					url: uri
				});
			}).catch(errorCallback);
		} else {
			downloadQueue.add(tile, tileCb => {
				if (tileCb.__loadIndex !== loadIndex) {
					return Promise.resolve();
				}

				return this.invokeOnePlugin(plugin => plugin.fetchData && plugin.fetchData(uri, { ...this.fetchOptions, signal }));
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
				tile.__loadingState = PARSING;

				return parseQueue.add(tile, parseTile => {
					// if it has been unloaded then the tile has been disposed
					if (parseTile.__loadIndex !== loadIndex) {
						return Promise.resolve();
					}

					const uri = parseTile.content.uri;
					const extension = getUrlExtension(uri);

					return this.parseTile(buffer, parseTile, extension);
				});
			}).then(() => {
				// if it has been unloaded then the tile has been disposed
				if (tile.__loadIndex !== loadIndex) {
					return;
				}

				stats.parsing--;
				tile.__loadingState = LOADED;

				if (tile.__wasSetVisible) {
					this.invokeOnePlugin(plugin => plugin.setTileVisible && plugin.setTileVisible(tile, true));
				}

				if (tile.__wasSetActive) {
					this.setTileActive(tile, true);
				}
			}).catch(errorCallback);
		}
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

	addCamera(camera) {
		const success = this.$cameras.add(camera);
		if (success) {
			this.dispatchEvent({ type: 'add-camera', camera });
		}
		return success;
	}

	removeCamera(camera) {
		const success = this.$cameras.remove(camera);
		if (success) {
			this.dispatchEvent({ type: 'delete-camera', camera });
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

class PivotPointMesh extends Mesh {

	constructor() {
		super(new PlaneGeometry(0, 0), new PivotMaterial());
		this.renderOrder = Infinity;
	}

}

class PivotMaterial extends ShaderMaterial {

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

const _vec$3 = new Vector2();
const _vec2 = new Vector2();

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
		this.hoverPosition = new Vector2();
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
		this.hoverPosition = new Vector2();
		this.hoverSet = false;
	}

	// The pointers can be set multiple times per frame so track whether the pointer has
	// been set this frame or not so we don't overwrite the previous position and lose information
	// about pointer movement
	updateFrame() {
		const { previousPositions, pointerPositions } = this;
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
		const position = new Vector2();
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
		this.getCenterPoint(_vec$3);
		this.getPreviousCenterPoint(_vec2);

		return _vec$3.sub(_vec2).getLength();
	}

	getTouchPointerDistance(pointerPositions = this.pointerPositions) {
		if (this.getPointerCount() <= 1 || this.getPointerType() === 'mouse') {
			return 0;
		}

		const { pointerOrder } = this;
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

const _matrix = new Matrix4();
const _ray$2 = new Ray();
const _vec$2 = new Vector3();

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
	target.x = ((clientX - element.offsetLeft) / element.clientWidth) * 2 - 1;
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

		_vec$2.set(0, 0, 0);
		_ray$2.closestPointToPoint(_vec$2, target).normalize();

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
	target
		.copy(ray.origin)
		.multiplyScalar(-1)
		.normalize();

	// get the normal of the plane the ray and origin lie in
	const rotationVec = _vec$2
		.crossVectors(target, ray.direction)
		.normalize();

	// rotate the camera direction by angle and scale it to the surface
	target
		.multiplyScalar(-1)
		.applyAxisAngle(rotationVec, -theta)
		.normalize()
		.multiplyScalar(radius);
}


// custom version of set raycaster from camera that relies on the underlying matrices
// so the ray origin is position at the camera near clip.
function setRaycasterFromCamera(raycaster, coords, camera) {
	const ray = raycaster instanceof Ray ? raycaster : raycaster.ray;
	const { origin, direction } = ray;

	// get the origin and direction of the frustum ray
	origin
		.set(coords.x, coords.y, -1)
		.unproject(camera);

	direction
		.set(coords.x, coords.y, 1)
		.unproject(camera)
		.sub(origin);

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

const _rotMatrix$1 = /* @__PURE__ */ new Matrix4();
const _delta = /* @__PURE__ */ new Vector3();
const _vec$1 = /* @__PURE__ */ new Vector3();
const _forward$1 = /* @__PURE__ */ new Vector3();
const _right$1 = /* @__PURE__ */ new Vector3();
const _rotationAxis = /* @__PURE__ */ new Vector3();
const _quaternion$2 = /* @__PURE__ */ new Quaternion();
const _plane = /* @__PURE__ */ new Plane();
const _localUp = /* @__PURE__ */ new Vector3();
const _mouseBefore = /* @__PURE__ */ new Vector3();
const _mouseAfter = /* @__PURE__ */ new Vector3();
const _identityQuat = /* @__PURE__ */ new Quaternion();
const _ray$1 = /* @__PURE__ */ new Ray();

const _zoomPointPointer = /* @__PURE__ */ new Vector2();
const _pointer$1 = /* @__PURE__ */ new Vector2();
const _prevPointer = /* @__PURE__ */ new Vector2();
const _deltaPointer = /* @__PURE__ */ new Vector2();
const _centerPoint = /* @__PURE__ */ new Vector2();
const _startCenterPoint = /* @__PURE__ */ new Vector2();

const _changeEvent = { type: 'change' };
const _startEvent = { type: 'start' };
const _endEvent = { type: 'end' };

class EnvironmentControls extends EventDispatcher {

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

		this.fallbackPlane = new Plane(new Vector3(0, 1, 0), 0);
		this.useFallbackPlane = true;

		// settings for GlobeControls
		this.reorientOnDrag = true;
		this.scaleZoomOrientationAtEdges = false;

		// internal state
		this.state = NONE$1;
		this.pointerTracker = new PointerTracker();
		this.needsUpdate = false;
		this.actionHeightOffset = 0;

		this.pivotPoint = new Vector3();

		// used for zoom
		this.zoomDirectionSet = false;
		this.zoomPointSet = false;
		this.zoomDirection = new Vector3();
		this.zoomPoint = new Vector3();
		this.zoomDelta = 0;

		// fields used for inertia
		this.rotationInertiaPivot = new Vector3();
		this.rotationInertia = new Vector2();
		this.dragInertia = new Vector3();
		this.inertiaTargetDistance = Infinity; 		// track the distance from the camera that we want to use to calculate the inertia end threshold
		this.inertiaStableFrames = 0; 				// the number of frames that the camera has not moved while the user is interacting

		// circular pivot mesh
		this.pivotMesh = new PivotPointMesh();
		this.pivotMesh.raycast = () => {};
		this.pivotMesh.scale.setScalar(0.25);

		// raycaster
		this.raycaster = new Raycaster();

		this.up = new Vector3(0, 1, 0);

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
				if (
					pointerTracker.getPointerCount() === 2 ||
					pointerTracker.isRightClicked() ||
					pointerTracker.isLeftClicked() && e.shiftKey
				) {
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

			const { pointerTracker } = this;
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

			const { pointerTracker } = this;
			pointerTracker.deletePointer(e);

			if (
				pointerTracker.getPointerType() === 'touch' &&
				pointerTracker.getPointerCount() === 0
			) {
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

			const { pointerTracker } = this;
			pointerTracker.setHoverEvent(e);
			pointerTracker.updatePointer(e);

			// TODO: do we need events here?
			this.dispatchEvent(_startEvent);

			let delta;
			switch (e.deltaMode) {
				case 2: // Pages
					delta = e.deltaY * 800;
					break;
				case 1: // Lines
					delta = e.deltaY * 40;
					break;
				case 0: // Pixels
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

			const { pointerTracker } = this;
			if (e.buttons !== pointerTracker.getPointerButtons()) {
				pointerTracker.deletePointer(e);
				this.resetState();
			}
		};

		domElement.addEventListener('contextmenu', contextMenuCallback);
		domElement.addEventListener('pointerdown', pointerdownCallback);
		domElement.addEventListener('pointermove', pointermoveCallback);
		domElement.addEventListener('pointerup', pointerupCallback);
		domElement.addEventListener('wheel', wheelCallback, { passive: false });
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
		const { camera, raycaster } = this;
		if (result !== null) {
			_vec$1.copy(result).project(camera);
			if (_vec$1.x < -1 || _vec$1.x > 1 || _vec$1.y < -1 || _vec$1.y > 1) {
				result = null;
			}
		}

		// default to the raycast hit if we have not result or the hit is closer to the camera
		// set a ray in the local ellipsoid frame
		setRaycasterFromCamera(raycaster, { x: 0, y: 0 }, camera);

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
				this.inertiaTargetDistance = _vec$1.copy(this.pivotPoint).sub(camera.position).dot(_forward$1);
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
			const { actionHeightOffset } = this;
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
		const { adjustHeight, cameraRadius } = this;
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
			setRaycasterFromCamera(_ray$1, _vec$1.set(0, 0, -1), camera);
			_ray$1.applyMatrix4(camera.viewMatrix);
			_ray$1.direction.normalize();
			_ray$1.recast(-_ray$1.direction.dot(_ray$1.origin)).at(stableDistance / _ray$1.direction.z, _vec$1);
			_vec$1.applyMatrix4(camera.worldMatrix);

			setRaycasterFromCamera(_ray$1, _delta.set(pixelThreshold, pixelThreshold, -1), camera);
			_ray$1.applyMatrix4(camera.viewMatrix);
			_ray$1.direction.normalize();
			_ray$1.recast(-_ray$1.direction.dot(_ray$1.origin)).at(stableDistance / _ray$1.direction.z, _delta);
			_delta.applyMatrix4(camera.worldMatrix);

			// get implied angle
			_vec$1.sub(pivotPoint).normalize();
			_delta.sub(pivotPoint).normalize();

			// calculate the rotation threshold
			const threshold = _vec$1.angleTo(_delta) / deltaTime;
			rotationInertia.multiplyScalar(factor);
			if (rotationInertia.getLengthSquared() < threshold ** 2 || !enableDamping) {
				rotationInertia.set(0, 0);
			}
		}

		// scale the residual translation motion
		if (dragInertia.getLengthSquared() > 0) {
			// calculate two screen points at 1 pixel apart in our notional resolution so we can stop when the delta is ~ 1 pixel
			// projected into world space
			setRaycasterFromCamera(_ray$1, _vec$1.set(0, 0, -1), camera);
			_ray$1.applyMatrix4(camera.viewMatrix);
			_ray$1.direction.normalize();
			_ray$1.recast(-_ray$1.direction.dot(_ray$1.origin)).at(stableDistance / _ray$1.direction.z, _vec$1);
			_vec$1.applyMatrix4(camera.worldMatrix);

			setRaycasterFromCamera(_ray$1, _delta.set(pixelThreshold, pixelThreshold, -1), camera);
			_ray$1.applyMatrix4(camera.viewMatrix);
			_ray$1.direction.normalize();
			_ray$1.recast(-_ray$1.direction.dot(_ray$1.origin)).at(stableDistance / _ray$1.direction.z, _delta);
			_delta.applyMatrix4(camera.worldMatrix);

			// calculate movement threshold
			const threshold = _vec$1.distanceTo(_delta) / deltaTime;
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
		const { rotationInertia, dragInertia } = this;
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
		if (!pointerTracker.getLatestPoint(_pointer$1) || (scale === 0 && state !== ZOOM)) {
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
			const finalZoomDirection = _vec$1.copy(zoomDirection);

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

		const { domElement, raycaster, camera, zoomDirection, pointerTracker } = this;
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
		const { raycaster } = this;
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

				_rotationAxis
					.crossVectors(raycaster.ray.direction, up)
					.normalize();

				raycaster.ray.direction
					.copy(up)
					.applyAxisAngle(_rotationAxis, angle)
					.multiplyScalar(-1);
			}

			// TODO: dragging causes the camera to rise because we're getting "pushed" up by lower resolution tiles and
			// don't lower back down. We should maintain a target height above tiles where possible
			// prevent the drag from inverting

			// if we drag to a point that's near the edge of the earth then we want to prevent it
			// from wrapping around and causing unexpected rotations
			this.getUpDirection(pivotPoint, _localUp);
			if (Math.abs(raycaster.ray.direction.dot(_localUp)) < DRAG_UP_THRESHOLD) {
				const angle = Math.acos(DRAG_UP_THRESHOLD);

				_rotationAxis
					.crossVectors(raycaster.ray.direction, _localUp)
					.normalize();

				raycaster.ray.direction
					.copy(_localUp)
					.applyAxisAngle(_rotationAxis, angle)
					.multiplyScalar(-1);
			}

			// find the point on the plane that we should drag to
			if (raycaster.ray.intersectPlane(_plane, _vec$1)) {
				_delta.subVectors(pivotPoint, _vec$1);
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
		_vec$1.crossVectors(_localUp, _forward$1).normalize();
		_right$1.set(1, 0, 0).transformDirection(camera.worldMatrix).normalize();
		const sign = Math.sign(_vec$1.dot(_right$1));
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
		camera.worldMatrix.decompose(camera.position, camera.quaternion, _vec$1);
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
			this.getUpDirection(zoomPoint, _vec$1);

			if (scaleZoomOrientationAtEdges) {
				let amt = Math.max(_vec$1.dot(up) - 0.6, 0) / 0.4;
				amt = MathUtils.mapLinear(amt, 0, 0.5, 0, 1);
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
			camera.worldMatrix.decompose(camera.position, camera.quaternion, _vec$1);

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
				camera.worldMatrix.decompose(camera.position, camera.quaternion, _vec$1);
			}
		}

		up.copy(newUp);
		camera.updateMatrix();
	}

	_raycast(raycaster) {
		const { scene, useFallbackPlane, fallbackPlane } = this;
		const result = raycaster.intersectObject(scene, true)[0] || null;
		if (result) {
			return result;
		} else if (useFallbackPlane) {
			// if we don't hit any geometry then try to intersect the fallback
			// plane so the camera can still be manipulated
			const plane = fallbackPlane;
			if (raycaster.ray.intersectPlane(plane, _vec$1)) {
				const planeHit = {
					point: _vec$1.clone(),
					distance: raycaster.ray.origin.distanceTo(_vec$1)
				};

				return planeHit;
			}
		}

		return null;
	}

}

const _invMatrix = /* @__PURE__ */ new Matrix4();
const _rotMatrix = /* @__PURE__ */ new Matrix4();
const _pos = /* @__PURE__ */ new Vector3();
const _vec = /* @__PURE__ */ new Vector3();
const _center = /* @__PURE__ */ new Vector3();
const _forward = /* @__PURE__ */ new Vector3();
const _right = /* @__PURE__ */ new Vector3();
const _targetRight = /* @__PURE__ */ new Vector3();
const _globalUp = /* @__PURE__ */ new Vector3();
const _quaternion$1 = /* @__PURE__ */ new Quaternion();
const _zoomPointUp = /* @__PURE__ */ new Vector3();
const _toCenter = /* @__PURE__ */ new Vector3();
const _ray = /* @__PURE__ */ new Ray();
const _ellipsoid = /* @__PURE__ */ new Ellipsoid();
const _latLon = {};

const _pointer = new Vector2();
const MIN_ELEVATION = 400;

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

		this.globeInertia = new Quaternion();
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
		const { camera, tilesGroup, ellipsoid } = this;

		// get camera values
		_forward.set(0, 0, -1).transformDirection(camera.worldMatrix);

		// set a ray in the local ellipsoid frame
		_ray.origin.copy(camera.position);
		_ray.direction.copy(_forward);
		_invMatrix.copy(tilesGroup.worldMatrix).invert();
		_ray.applyMatrix4(_invMatrix);

		// get the estimated closest point
		closestRayEllipsoidSurfacePointEstimate(_ray, ellipsoid, _vec);
		_vec.applyMatrix4(tilesGroup.worldMatrix);

		// use the closest point if no pivot was provided or it's closer
		if (
			super.getPivotPoint(target) === null ||
			target.distanceTo(_ray.origin) > _vec.distanceTo(_ray.origin)
		) {
			target.copy(_vec);
		}

		return target;
	}

	// get the vector to the center of the provided globe
	getVectorToCenter(target) {
		const { tilesGroup, camera } = this;
		return target
			.setFromMatrixPosition(tilesGroup.worldMatrix)
			.sub(camera.position);
	}

	// get the distance to the center of the globe
	getDistanceToCenter() {
		return this
			.getVectorToCenter(_vec)
			.getLength();
	}

	getUpDirection(point, target) {
		// get the "up" direction based on the wgs84 ellipsoid
		const { tilesGroup, ellipsoid } = this;
		_invMatrix.copy(tilesGroup.worldMatrix).invert();
		_vec.copy(point).applyMatrix4(_invMatrix);

		ellipsoid.getPositionToNormal(_vec, target);
		target.transformDirection(tilesGroup.worldMatrix);
	}

	getCameraUpDirection(target) {
		const { tilesGroup, ellipsoid, camera } = this;
		if (camera.isOrthographicCamera) {
			this._getVirtualOrthoCameraPosition(_vec);

			_invMatrix.copy(tilesGroup.worldMatrix).invert();
			_vec.applyMatrix4(_invMatrix);

			ellipsoid.getPositionToNormal(_vec, target);
			target.transformDirection(tilesGroup.worldMatrix);
		} else {
			this.getUpDirection(camera.position, target);
		}
	}

	update(deltaTime = 64 / 1000) {
		if (!this.enabled || !this.tilesGroup || !this.camera || deltaTime === 0) {
			return;
		}

		const { camera, pivotMesh } = this;

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

		const { tilesGroup, ellipsoid, nearMargin, farMargin } = this;
		const maxRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
		if (camera.isPerspectiveCamera) {
			// adjust the clip planes
			const distanceToCenter = _vec
				.setFromMatrixPosition(tilesGroup.worldMatrix)
				.sub(camera.position).getLength();

			// update the projection matrix
			// interpolate from the 25% radius margin around the globe down to the surface
			// so we can avoid z fighting when near value is too far at a high altitude
			const margin = nearMargin * maxRadius;
			const alpha = MathUtils.clamp((distanceToCenter - maxRadius) / margin, 0, 1);
			const minNear = MathUtils.lerp(1, 1000, alpha);
			camera.near = Math.max(minNear, distanceToCenter - maxRadius - margin);

			// update the far plane to the horizon distance
			_invMatrix.copy(tilesGroup.worldMatrix).invert();
			_pos.copy(camera.position).applyMatrix4(_invMatrix);
			ellipsoid.getPositionToCartographic(_pos, _latLon);

			// use a minimum elevation for computing the horizon distance to avoid the far clip
			// plane approaching zero as the camera goes to or below sea level.
			const elevation = Math.max(ellipsoid.getPositionElevation(_pos), MIN_ELEVATION);
			const horizonDistance = ellipsoid.calculateHorizonDistance(_latLon.lat, elevation);

			// extend the horizon distance by 2.5 to handle cases where geometry extends above the horizon
			camera.far = horizonDistance * 2.5 + 0.1 + maxRadius * farMargin;
			camera.updateProjectionMatrix();
		} else {
			this._getVirtualOrthoCameraPosition(camera.position, camera);
			camera.updateMatrix();

			_invMatrix.copy(camera.worldMatrix).invert();
			_vec.setFromMatrixPosition(tilesGroup.worldMatrix).applyMatrix4(_invMatrix);

			const distanceToCenter = -_vec.z;
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
			setRaycasterFromCamera(_ray, _vec.set(0, 0, -1), camera);
			_ray.applyMatrix4(camera.viewMatrix);
			_ray.direction.normalize();
			_ray.recast(-_ray.direction.dot(_ray.origin)).at(stableDistance / _ray.direction.z, _vec);
			_vec.applyMatrix4(camera.worldMatrix);

			setRaycasterFromCamera(_ray, _pos.set(pixelThreshold, pixelThreshold, -1), camera);
			_ray.applyMatrix4(camera.viewMatrix);
			_ray.direction.normalize();
			_ray.recast(-_ray.direction.dot(_ray.origin)).at(stableDistance / _ray.direction.z, _pos);
			_pos.applyMatrix4(camera.worldMatrix);

			// get implied angle
			_vec.sub(_center).normalize();
			_pos.sub(_center).normalize();

			this.globeInertiaFactor *= factor;
			const threshold = _vec.angleTo(_pos) / deltaTime;
			const globeAngle = 2 * Math.acos(globeInertia.w) * this.globeInertiaFactor;
			if (globeAngle < threshold || !enableDamping) {
				this.globeInertiaFactor = 0;
				globeInertia.identity();
			}
		}

		if (this.globeInertiaFactor !== 0) {
			// ensure our w component is non-one if the xyz values are
			// non zero to ensure we can animate
			if (
				globeInertia.w === 1 && (
					globeInertia.x !== 0 ||
					globeInertia.y !== 0 ||
					globeInertia.z !== 0
				)
			) {
				globeInertia.w = Math.min(globeInertia.w, 1 - 1e-9);
			}

			// construct the rotation matrix
			_center.setFromMatrixPosition(tilesGroup.worldMatrix);
			_quaternion$1.identity().slerp(globeInertia, this.globeInertiaFactor * deltaTime);
			makeRotateAroundPoint(_center, _quaternion$1, _rotMatrix);

			// apply the rotation
			camera.worldMatrix.premultiply(_rotMatrix);
			camera.worldMatrix.decompose(camera.position, camera.quaternion, _vec);
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
			const pivotDir = _pos;
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
			const pivotRadius = _vec.copy(pivotPoint).applyMatrix4(_invMatrix).getLength();
			_ellipsoid.radius.setScalar(pivotRadius);

			// find the hit point and use the closest point on the horizon if we miss
			if (camera.isPerspectiveCamera) {
				if (!_ellipsoid.intersectRay(raycaster.ray, _vec)) {
					closestRaySpherePointFromRotation(raycaster.ray, pivotRadius, _vec);
				}
			} else {
				closestRayEllipsoidSurfacePointEstimate(raycaster.ray, _ellipsoid, _vec);
			}
			_vec.applyMatrix4(tilesGroup.worldMatrix);

			// get the point directions
			_center.setFromMatrixPosition(tilesGroup.worldMatrix);
			pivotDir.subVectors(pivotPoint, _center).normalize();
			newPivotDir.subVectors(_vec, _center).normalize();

			// construct the rotation
			_quaternion$1.setFromUnitVectors(newPivotDir, pivotDir);
			makeRotateAroundPoint(_center, _quaternion$1, _rotMatrix);

			// apply the rotation
			camera.worldMatrix.premultiply(_rotMatrix);
			camera.worldMatrix.decompose(camera.position, camera.quaternion, _vec);

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
		const { zoomDelta, ellipsoid, zoomSpeed, zoomPoint, camera, maxZoom, state } = this;

		if (state !== ZOOM && zoomDelta === 0) {
			return;
		}

		// reset momentum
		this.rotationInertia.set(0, 0);
		this.dragInertia.set(0, 0, 0);
		this.globeInertia.identity();
		this.globeInertiaFactor = 0;

		// used to scale the tilt transitions based on zoom intensity
		const deltaAlpha = MathUtils.clamp(MathUtils.mapLinear(Math.abs(zoomDelta), 0, 20, 0, 1), 0, 1);
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
				const upAlpha = MathUtils.clamp(MathUtils.mapLinear(-_zoomPointUp.dot(_toCenter), 1, 0.95, 0, 1), 0, 1);
				const forwardAlpha = 1 - _forward.dot(_toCenter);
				const cameraAlpha = camera.isOrthographicCamera ? 0.05 : 1;
				const adjustedDeltaAlpha = MathUtils.clamp(deltaAlpha * 3, 0, 1);

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
			const distanceAlpha = MathUtils.mapLinear(this.getDistanceToCenter(), transitionDistance, maxDistance, 0, 1);
			this._tiltTowardsCenter(MathUtils.lerp(0, 0.4, distanceAlpha * deltaAlpha));
			this._alignCameraUpToNorth(MathUtils.lerp(0, 0.2, distanceAlpha * deltaAlpha));

			// calculate zoom in a similar way to environment controls so
			// the zoom speeds are comparable
			const dist = this.getDistanceToCenter() - ellipsoid.radius.x;
			const scale = zoomDelta * dist * zoomSpeed * 0.0025;
			const clampedScale = Math.max(scale, Math.min(this.getDistanceToCenter() - maxDistance, 0));

			// zoom out directly from the globe center
			this.getVectorToCenter(_vec).normalize();
			this.camera.position.addScaledVector(_vec, clampedScale);
			this.camera.updateMatrix();

			this.zoomDelta = 0;
		} else {
			const transitionZoom = this._getOrthographicTransitionZoom();
			const minZoom = this._getMinOrthographicZoom();
			const distanceAlpha = MathUtils.mapLinear(camera.zoom, transitionZoom, minZoom, 0, 1);
			this._tiltTowardsCenter(MathUtils.lerp(0, 0.4, distanceAlpha * deltaAlpha));
			this._alignCameraUpToNorth(MathUtils.lerp(0, 0.2, distanceAlpha * deltaAlpha));

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
		const { tilesGroup } = this;
		_globalUp.set(0, 0, 1).transformDirection(tilesGroup.worldMatrix);
		this._alignCameraUp(_globalUp, alpha);
	}

	// tilt the camera to align with the provided "up" value
	_alignCameraUp(up, alpha = null) {
		const { camera } = this;
		_forward.set(0, 0, -1).transformDirection(camera.worldMatrix);
		_right.set(-1, 0, 0).transformDirection(camera.worldMatrix);
		_targetRight.crossVectors(up, _forward);

		// compute the alpha based on how far away from boresight the up vector is
		// so we can ease into the correct orientation
		if (alpha === null) {
			alpha = 1 - Math.abs(_forward.dot(up));
			alpha = MathUtils.mapLinear(alpha, 0, 1, -0.01, 1);
			alpha = MathUtils.clamp(alpha, 0, 1) ** 2;
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
		_vec.setFromMatrixPosition(tilesGroup.worldMatrix).sub(camera.position).normalize();
		_vec.lerp(_forward, 1 - alpha).normalize();

		_quaternion$1.setFromUnitVectors(_forward, _vec);
		camera.quaternion.premultiply(_quaternion$1);
		camera.updateMatrix();
	}

	// returns the perspective camera transition distance can move to based on globe size and fov
	_getPerspectiveTransitionDistance() {
		const { camera, ellipsoid } = this;
		if (!camera.isPerspectiveCamera) {
			throw new Error();
		}

		// When the smallest fov spans 65% of the ellipsoid then we use the near controls
		const ellipsoidRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
		const fovHoriz = 2 * Math.atan(Math.tan(MathUtils.DEG2RAD * camera.fov * 0.5) * camera.aspect);
		const distVert = ellipsoidRadius / Math.tan(MathUtils.DEG2RAD * camera.fov * 0.5);
		const distHoriz = ellipsoidRadius / Math.tan(fovHoriz * 0.5);
		const dist = Math.max(distVert, distHoriz);

		return dist;
	}

	// returns the max distance the perspective camera can move to based on globe size and fov
	_getMaxPerspectiveDistance() {
		const { camera, ellipsoid } = this;
		if (!camera.isPerspectiveCamera) {
			throw new Error();
		}

		// allow for zooming out such that the ellipsoid is half the size of the largest fov
		const ellipsoidRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
		const fovHoriz = 2 * Math.atan(Math.tan(MathUtils.DEG2RAD * camera.fov * 0.5) * camera.aspect);
		const distVert = ellipsoidRadius / Math.tan(MathUtils.DEG2RAD * camera.fov * 0.5);
		const distHoriz = ellipsoidRadius / Math.tan(fovHoriz * 0.5);
		const dist = 2 * Math.max(distVert, distHoriz);

		return dist;
	}

	// returns the transition threshold for orthographic zoom based on the globe size and camera settings
	_getOrthographicTransitionZoom() {
		const { camera, ellipsoid } = this;
		if (!camera.isOrthographicCamera) {
			throw new Error();
		}

		const orthoHeight = (camera.top - camera.bottom);
		const orthoWidth = (camera.right - camera.left);
		const orthoSize = Math.max(orthoHeight, orthoWidth);
		const ellipsoidRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
		const ellipsoidDiameter = 2 * ellipsoidRadius;
		return 2 * orthoSize / ellipsoidDiameter;
	}

	// returns the minimum allowed orthographic zoom based on the globe size and camera settings
	_getMinOrthographicZoom() {
		const { camera, ellipsoid } = this;
		if (!camera.isOrthographicCamera) {
			throw new Error();
		}

		const orthoHeight = (camera.top - camera.bottom);
		const orthoWidth = (camera.right - camera.left);
		const orthoSize = Math.min(orthoHeight, orthoWidth);
		const ellipsoidRadius = Math.max(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
		const ellipsoidDiameter = 2 * ellipsoidRadius;
		return 0.7 * orthoSize / ellipsoidDiameter;
	}

	// returns the "virtual position" of the orthographic based on where it is and
	// where it's looking primarily so we can reasonably position the camera object
	// in space and derive a reasonable "up" value.
	_getVirtualOrthoCameraPosition(target, camera = this.camera) {
		const { tilesGroup, ellipsoid } = this;
		if (!camera.isOrthographicCamera) {
			throw new Error();
		}

		// get ray in globe coordinate frame
		_ray.origin.copy(camera.position);
		_ray.direction.set(0, 0, -1).transformDirection(camera.worldMatrix);
		_invMatrix.copy(tilesGroup.worldMatrix).invert();
		_ray.applyMatrix4(_invMatrix);

		// get the closest point to the ray on the globe in the global coordinate frame
		closestRayEllipsoidSurfacePointEstimate(_ray, ellipsoid, _pos);
		_pos.applyMatrix4(tilesGroup.worldMatrix);

		// get ortho camera info
		const orthoHeight = (camera.top - camera.bottom);
		const orthoWidth = (camera.right - camera.left);
		const orthoSize = Math.max(orthoHeight, orthoWidth) / camera.zoom;
		_forward.set(0, 0, -1).transformDirection(camera.worldMatrix);

		// ensure we move the camera exactly along the forward vector to avoid shifting
		// the camera in other directions due to floating point error
		const dist = _pos.sub(camera.position).dot(_forward);
		target.copy(camera.position).addScaledVector(_forward, dist - orthoSize * 4);
	}

	_isNearControls() {
		const { camera } = this;
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
			const { ellipsoid, tilesGroup } = this;
			_invMatrix.copy(tilesGroup.worldMatrix).invert();
			_ray.copy(raycaster.ray).applyMatrix4(_invMatrix);

			const point = ellipsoid.intersectRay(_ray, _vec);
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

		const delta = MathUtils.clamp((time - this._lastTick) / this.duration, 0, 1);
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
			fadeIn = MathUtils.clamp(fadeIn + fadeInSign * delta, 0, 1);

			const fadeOutSign = Math.sign(fadeOutTarget - fadeOut);
			fadeOut = MathUtils.clamp(fadeOut + fadeOutSign * delta, 0, 1);

			state.fadeIn = fadeIn;
			state.fadeOut = fadeOut;

			// Check if the fade in and fade out animations are complete
			const fadeOutComplete = fadeOut === 1 || fadeOut === 0;
			const fadeInComplete = fadeIn === 1 || fadeIn === 0;

			// If they are or the fade out animation is further along than the
			// fade in animation then mark the fade as completed for this tile
			if ((fadeOutComplete && fadeInComplete) || fadeOut >= fadeIn) {
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

	const fragmentShader = material.fragmentShader ||
		(material.type === MATERIAL_TYPE.BASIC ?
			ShaderLib.basic_frag : ShaderLib.pbr_frag);

	material.type = MATERIAL_TYPE.SHADER;

	material.vertexShader = material.type === MATERIAL_TYPE.BASIC ?
		ShaderLib.basic_vert : ShaderLib.pbr_vert;
	material.fragmentShader = fragmentShader
		.replace(/void main\(/, value => /* glsl */`
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
		`)
		.replace(/#include <end_frag>/, value => /* glsl */`
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
		this._onLoadModel = ({ scene }) => {
			// initialize all the scene materials to fade
			this._fadeMaterialManager.prepareScene(scene);
		};
		this._onDisposeModel = ({ tile, scene }) => {
			// delete the fade info from the managers on disposal of model
			this._fadeManager.deleteObject(tile);
			this._fadeMaterialManager.deleteScene(scene);
		};
		this._onAddCamera = ({ camera }) => {
			// track the camera transform
			this._prevCameraTransforms.set(camera, new Matrix4());
		};
		this._onDeleteCamera = ({ camera }) => {
			// remove the camera transform
			this._prevCameraTransforms.delete(camera);
		};
		this._onTileVisibilityChange = ({ tile, visible }) => {
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
			tiles.dispatchEvent({ type: 'fade-start' });
			tiles.dispatchEvent({ type: 'needs-render' });
		};

		fadeManager.onFadeSetComplete = () => {
			tiles.dispatchEvent({ type: 'fade-end' });
			tiles.dispatchEvent({ type: 'needs-render' });
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
			prevCameraTransforms.set(camera, new Matrix4());
		});

		tiles.forEachLoadedModel((scene, tile) => {
			this._onLoadModel({ scene });
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
const _fromPos = new Vector3();
const _toPos = new Vector3();
const _fromQuat = new Quaternion();
const _toQuat = new Quaternion();
const _scale = new Vector3();

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
	const { tiles, maximumFadeOutTiles } = this;
	const cameras = tiles.$cameras._cameras;

	// reset the active tiles flag
	tiles.displayActiveTiles = displayActiveTiles;

	// update fade step
	fadeManager.update();

	// fire an event
	const fadingAfter = fadeManager.fadeCount;
	if (fadingBefore !== 0 && fadingAfter !== 0) {
		tiles.dispatchEvent({ type: 'fade-change' });
		tiles.dispatchEvent({ type: 'needs-render' });
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
	fadeManager.forEachObject((tile, { fadeIn, fadeOut }) => {
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

	constructor({ apiToken, autoRefreshToken = false, logoUrl = null, useRecommendedSettings = true }) {
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
		this._onLoadCallback = ({ tileSet }) => {
			// the first tile set loaded will be the root
			this.sessionToken = getSessionToken(tileSet.root);

			// clear the callback once the root is loaded
			tiles.removeEventListener('load-tile-set', this._onLoadCallback);
		};

		this._visibilityChangeCallback = ({ tile, visible }) => {
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
		const { tiles } = this;
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
			this._tokenRefreshPromise = fetch(rootURL, options)
				.then(res => res.json())
				.then(res => {
					this.sessionToken = getSessionToken(res.root);
					this._tokenRefreshPromise = null;
				});

			// dispatch an error if we fail to refresh the token
			this._tokenRefreshPromise
				.catch(error => {
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

class CesiumIonAuthPlugin {

	constructor({ apiToken, assetId = null, autoRefreshToken = false, useRecommendedSettings = true }) {
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
		return this._refreshToken()
			.then(() => {
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

			this._tokenRefreshPromise = fetch(url, options)
				.then(res => {
					if (this._disposed) {
						return null;
					}

					if (!res.ok) {
						throw new Error(`CesiumIonAuthPlugin: Failed to load data with error code ${res.status}`);
					}

					return res.json();
				})
				.then(json => {
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
						if (json.type === 'TERRAIN' && tiles.getPluginByName('QUANTIZED_MESH_PLUGIN') === null) ; else if (json.type === 'IMAGERY' && tiles.getPluginByName('TMS_TILES_PLUGIN') === null) ;

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
			this._tokenRefreshPromise
				.catch(error => {
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

class DebugTilesPlugin {

	static get ColorModes() {
		return ColorModes;
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
			...options
		};

		this.name = 'DEBUG_TILES_PLUGIN';
		this.tiles = null;

		this._enabled = true;

		this.extremeDebugDepth = -1;
		this.extremeDebugError = -1;
		this.boxGroup = null;
		this.sphereGroup = null;
		this.regionGroup = null;

		// options
		this._displayParentBounds = options.displayParentBounds;
		this.displayBoxBounds = options.displayBoxBounds;
		this.displaySphereBounds = options.displaySphereBounds;
		this.displayRegionBounds = options.displayRegionBounds;
		this.colorMode = options.colorMode;
		this.maxDebugDepth = options.maxDebugDepth;
		this.maxDebugDistance = options.maxDebugDistance;
		this.maxDebugError = options.maxDebugError;
		this.customColorCallback = options.customColorCallback;

		this.getDebugColor = (value, target) => {
			target.setRGB(value, value, value);
		};
	}

	get enabled() {
		return this._enabled;
	}

	set enabled(v) {
		if (v !== this._enabled) {
			this._enabled = v;

			if (this._enabled) {
				if (this.tiles) {
					this.init(this.tiles);
				}
			} else {
				this.dispose();
			}
		}
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

		// initialize groups
		this.boxGroup = new Object3D();
		this.boxGroup.name = 'DebugTilesPlugin.boxGroup';
		tiles.add(this.boxGroup);
		this.boxGroup.updateMatrix();

		this.sphereGroup = new Object3D();
		this.sphereGroup.name = 'DebugTilesPlugin.sphereGroup';
		tiles.add(this.sphereGroup);
		this.sphereGroup.updateMatrix();

		this.regionGroup = new Object3D();
		this.regionGroup.name = 'DebugTilesPlugin.regionGroup';
		tiles.add(this.regionGroup);
		this.regionGroup.updateMatrix();

		// register events
		this._onLoadTileSetCB = () => {
			this._initExtremes();
		};

		this._onLoadModelCB = ({ scene, tile }) => {
			this._onLoadModel(scene, tile);
		};

		this._onDisposeModelCB = ({ tile }) => {
			this._onDisposeModel(tile);
		};

		this._onUpdateAfterCB = () => {
			this._onUpdateAfter();
		};

		this._onTileVisibilityChangeCB = ({ scene, tile, visible }) => {
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
		const tiles = this.tiles;

		if (!tiles.root) {
			return;
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

		const errorTarget = tiles.errorTarget;
		const colorMode = this.colorMode;
		const visibleTiles = tiles.visibleTiles;
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

				const currMaterial = c.material;
				if (currMaterial) {
					// Reset the material if needed
					const originalMaterial = c[ORIGINAL_MATERIAL];

					if (colorMode === NONE && currMaterial !== originalMaterial) {
						c.material.dispose();
						c.material = c[ORIGINAL_MATERIAL];
					} else if (colorMode !== NONE && currMaterial === originalMaterial) {
						if (c.material.drawMode === DRAW_MODE.POINTS) {
							const pointsMaterial = new PointsMaterial();
							pointsMaterial.size = originalMaterial.size;
							pointsMaterial.sizeAttenuation = originalMaterial.sizeAttenuation;
							c.material = pointsMaterial;
						} else {
							if (c.material.isInstancedPBRMaterial) {
								c.material = new InstancedPBRMaterial();
								c.material.metalness = 0.0;
								c.material.roughness = 1.0;
							} else if (c.material.isInstancedBasicMaterial) {
								c.material = new InstancedBasicMaterial();
							} else {
								c.material = new PBRMaterial();
								c.material.metalness = 0.0;
								c.material.roughness = 1.0;
							}

							c.material.shading = SHADING_TYPE.FLAT_SHADING;
							c.material.envMap = undefined;
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
							const val = Math.min(tile.geometricError / maxError, 1);
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
								console.warn('DebugTilesPlugin: customColorCallback not defined');
							}
							break;
						}
						case LOAD_ORDER: {
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

				const tileVisible = (current === tile && visible) || (this.displayParentBounds && current[PARENT_BOUND_REF_COUNT] > 0);

				this._updateBoundHelper(current, tileVisible);
			});
		} else {
			this._updateBoundHelper(tile, visible);
		}
	}

	_createBoundHelper(tile) {
		const tiles = this.tiles;
		const cached = tile.cached;
		const { sphere, obb, region } = cached.boundingVolume;
		if (obb) {
			// Create debug bounding box
			// In some cases the bounding box may have a scale of 0 in one dimension resulting
			// in the NaNs in an extracted rotation so we disable matrix updates instead.
			const boxHelperGroup = new Object3D();
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

		if (visible && (cached.boxHelperGroup == null && cached.sphereHelper == null && cached.regionHelper == null)) {
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

	_onLoadModel(scene, tile) {
		tile[LOAD_TIME] = performance.now();

		// Cache the original materials
		scene.traverse(c => {
			const material = c.material;
			if (material) {
				c[ORIGINAL_MATERIAL] = material;
			}
		});
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
		const tiles = this.tiles;

		if (tiles) {
			tiles.removeEventListener('load-tile-set', this._onLoadTileSetCB);
			tiles.removeEventListener('load-model', this._onLoadModelCB);
			tiles.removeEventListener('dispose-model', this._onDisposeModelCB);
			tiles.removeEventListener('update-after', this._onUpdateAfterCB);
			tiles.removeEventListener('tile-visibility-change', this._onTileVisibilityChangeCB);

			// reset all materials
			this.colorMode = NONE;
			this._onUpdateAfter();

			// dispose of all helper objects
			tiles.traverse(tile => {
				this._onDisposeModel(tile);
			});
		}

		this.boxGroup?.removeFromParent();
		this.sphereGroup?.removeFromParent();
		this.regionGroup?.removeFromParent();
	}

}

const ORIGINAL_MATERIAL = Symbol('ORIGINAL_MATERIAL');
const HAS_RANDOM_COLOR = Symbol('HAS_RANDOM_COLOR');
const HAS_RANDOM_NODE_COLOR = Symbol('HAS_RANDOM_NODE_COLOR');
const LOAD_TIME = Symbol('LOAD_TIME');
const PARENT_BOUND_REF_COUNT = Symbol('PARENT_BOUND_REF_COUNT');

const _sphere = new Sphere();
const emptyRaycast = () => {};
const colors = {};

// Return a consistant random color for an index
function getIndexedRandomColor(index) {
	if (!colors[index]) {
		const h = Math.random();
		const s = 0.5 + Math.random() * 0.5;
		const l = 0.375 + Math.random() * 0.25;

		colors[index] = new Color3().setHSL(h, s, l);
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
			const { up, lat, lon, height, recenter } = this;

			if (lat !== null && lon !== null) {
				// if the latitude and longitude are provided then remove the position offset
				this.transformLatLonHeightToOrigin(lat, lon, height);
			} else {
				const { ellipsoid } = tiles;
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
						case 'x': case '+x':
							tiles.euler.z = Math.PI / 2;
							break;
						case '-x':
							tiles.euler.z = -Math.PI / 2;
							break;

						case 'y': case '+y':
							break;
						case '-y':
							tiles.euler.z = Math.PI;
							break;

						case 'z': case '+z':
							tiles.euler.x = -Math.PI / 2;
							break;
						case '-z':
							tiles.euler.x = Math.PI / 2;
							break;
					}

					tiles.position
						.copy(sphere.center)
						.applyEuler(tiles.euler)
						.multiplyScalar(-1);
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
		const { ellipsoid } = tiles;

		// get ENU orientation (Z facing north and X facing west) and position
		ellipsoid.getRotationMatrixFromAzElRoll(lat, lon, 0, 0, 0, tiles.matrix, OBJECT_FRAME);
		ellipsoid.getCartographicToPosition(lat, lon, height, vec);

		// adjust the tiles matrix
		tiles.matrix
			.setPosition(vec)
			.invert()
			.decompose(tiles.position, tiles.quaternion, tiles.scale);
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

const sphere = new Sphere();
const vec = new Vector3();

class LoadParser {

	static parse(context, loader) {
		const extension = getUrlExtension(context.url);
		let pr = null;
		if (extension === 'gltf') {
			pr = loader.loadFile(context.url)
				.then(buffer => {
					context.options.buffer = buffer;
				});
		} else {
			pr = loader.loadFile(context.url, 'arraybuffer')
				.then(buffer => {
					context.options.buffer = buffer;
				});
		}
		return pr;
	}

}

const _quaternion = new Quaternion();
const _vector = new Vector3();

if (!Matrix4.prototype.makeScale) {
	Matrix4.prototype.makeScale = function(x, y, z) {
		return this.set(
			x, 0, 0, 0,
			0, y, 0, 0,
			0, 0, z, 0,
			0, 0, 0, 1
		);
	};
}

Matrix4.prototype.makeBasis = function(xAxis, yAxis, zAxis) {
	this.set(
		xAxis.x, yAxis.x, zAxis.x, 0,
		xAxis.y, yAxis.y, zAxis.y, 0,
		xAxis.z, yAxis.z, zAxis.z, 0,
		0, 0, 0, 1
	);

	return this;
};

Matrix4.prototype.setPosition = function(x, y, z) {
	const te = this.elements;

	if (x.isVector3) {
		te[12] = x.x;
		te[13] = x.y;
		te[14] = x.z;
	} else {
		te[12] = x;
		te[13] = y;
		te[14] = z;
	}

	return this;
};

Matrix4.prototype.makeRotationFromEuler = function(euler) {
	const te = this.elements;

	const x = euler.x, y = euler.y, z = euler.z;
	const a = Math.cos(x), b = Math.sin(x);
	const c = Math.cos(y), d = Math.sin(y);
	const e = Math.cos(z), f = Math.sin(z);

	if (euler.order === 'XYZ') {
		const ae = a * e, af = a * f, be = b * e, bf = b * f;

		te[0] = c * e;
		te[4] = -c * f;
		te[8] = d;

		te[1] = af + be * d;
		te[5] = ae - bf * d;
		te[9] = -b * c;

		te[2] = bf - ae * d;
		te[6] = be + af * d;
		te[10] = a * c;
	} else if (euler.order === 'YXZ') {
		const ce = c * e, cf = c * f, de = d * e, df = d * f;

		te[0] = ce + df * b;
		te[4] = de * b - cf;
		te[8] = a * d;

		te[1] = a * f;
		te[5] = a * e;
		te[9] = -b;

		te[2] = cf * b - de;
		te[6] = df + ce * b;
		te[10] = a * c;
	} else if (euler.order === 'ZXY') {
		const ce = c * e, cf = c * f, de = d * e, df = d * f;

		te[0] = ce - df * b;
		te[4] = -a * f;
		te[8] = de + cf * b;

		te[1] = cf + de * b;
		te[5] = a * e;
		te[9] = df - ce * b;

		te[2] = -a * d;
		te[6] = b;
		te[10] = a * c;
	} else if (euler.order === 'ZYX') {
		const ae = a * e, af = a * f, be = b * e, bf = b * f;

		te[0] = c * e;
		te[4] = be * d - af;
		te[8] = ae * d + bf;

		te[1] = c * f;
		te[5] = bf * d + ae;
		te[9] = af * d - be;

		te[2] = -d;
		te[6] = b * c;
		te[10] = a * c;
	} else if (euler.order === 'YZX') {
		const ac = a * c, ad = a * d, bc = b * c, bd = b * d;

		te[0] = c * e;
		te[4] = bd - ac * f;
		te[8] = bc * f + ad;

		te[1] = f;
		te[5] = a * e;
		te[9] = -b * e;

		te[2] = -d * e;
		te[6] = ad * f + bc;
		te[10] = ac - bd * f;
	} else if (euler.order === 'XZY') {
		const ac = a * c, ad = a * d, bc = b * c, bd = b * d;

		te[0] = c * e;
		te[4] = -f;
		te[8] = d * e;

		te[1] = ac * f + bd;
		te[5] = a * e;
		te[9] = ad * f - bc;

		te[2] = bc * f - ad;
		te[6] = b * e;
		te[10] = bd * f + ac;
	}

	// bottom row
	te[3] = 0;
	te[7] = 0;
	te[11] = 0;

	// last column
	te[12] = 0;
	te[13] = 0;
	te[14] = 0;
	te[15] = 1;

	return this;
};

Matrix4.prototype.invert = function() {
	return this.getInverse(this);
};

Vector3.prototype.applyEuler = function(euler) {
	return this.applyQuaternion(_quaternion.setFromEuler(euler));
};

Vector3.prototype.applyAxisAngle = function(axis, angle) {
	return this.applyQuaternion(_quaternion.setFromAxisAngle(axis, angle));
};

Vector3.prototype.angleTo = function(v) {
	const denominator = Math.sqrt(this.getLengthSquared() * v.getLengthSquared());

	if (denominator === 0) return Math.PI / 2;

	const theta = this.dot(v) / denominator;

	// clamp, to handle numerical problems

	return Math.acos(MathUtils.clamp(theta, -1, 1));
};

Vector3.prototype.isVector3 = true;

Quaternion.prototype.angleTo = function(q) {
	return 2 * Math.acos(Math.abs(MathUtils.clamp(this.dot(q), -1, 1)));
};

Quaternion.prototype.identity = function() {
	return this.set(0, 0, 0, 1);
};

if (!Quaternion.prototype.slerp) {
	Quaternion.prototype.slerp = function(q, t) {
		this.slerpQuaternions(this, q, t);
		return this;
	};
}


Object3D.prototype.removeFromParent = function() {
	const parent = this.parent;

	if (parent !== null) {
		parent.remove(this);
	}

	return this;
};

Ray.prototype.closestPointToPoint = function(point, target) {
	target.subVectors(point, this.origin);

	const directionDistance = target.dot(this.direction);

	if (directionDistance < 0) {
		return target.copy(this.origin);
	}

	return target.copy(this.origin).addScaledVector(this.direction, directionDistance);
};

Ray.prototype.recast = function(t) {
	this.origin.copy(this.at(t, _vector));

	return this;
};

MathUtils.mapLinear = function(x, a1, a2, b1, b2) {
	return b1 + (x - a1) * (b2 - b1) / (a2 - a1);
};

MathUtils.DEG2RAD = Math.PI / 180;

MathUtils.lerp = function(x, y, t) {
	return x + (y - x) * t;
};

let oldMethod;

oldMethod = Camera.prototype.setOrtho;
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

	oldMethod.call(this, left, right, bottom, top, near, far);
};

oldMethod = Camera.prototype.setPerspective;
Camera.prototype.setPerspective = function(fov, aspect, near, far) {
	this.fov = fov;
	this.aspect = aspect;
	this.near = near;
	this.far = far;

	this.isPerspectiveCamera = true;
	this.isOrthographicCamera = false;

	oldMethod.call(this, fov, aspect, near, far);
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

export { B3DMLoader, CMPTLoader, CesiumIonAuthPlugin, LoadParser as DebugLoadParser, DebugTilesPlugin, EnvironmentControls, GlobeControls, I3DMLoader, InstancedBasicMaterial, InstancedPBRMaterial, OBB, PNTSLoader, ReorientationPlugin, TileGLTFLoader, Tiles3D, TilesFadePlugin };
