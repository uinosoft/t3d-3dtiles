import { Box3, Matrix3, Vector3 } from 't3d';

/**
 * An oriented bounding box.
 */
export class OBB {

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
		_vec3_1.copy(axisX);
		_vec3_2.copy(axisY);
		_vec3_3.copy(axisZ);

		const scaleX = _vec3_1.getLength();
		const scaleY = _vec3_2.getLength();
		const scaleZ = _vec3_3.getLength();

		_vec3_1.normalize();
		_vec3_2.normalize();
		_vec3_3.normalize();

		// handle the case where the box has a dimension of 0 in one axis
		if (scaleX === 0) {
			_vec3_1.crossVectors(_vec3_2, _vec3_3);
		}

		if (scaleY === 0) {
			_vec3_2.crossVectors(_vec3_1, _vec3_3);
		}

		if (scaleZ === 0) {
			_vec3_3.crossVectors(_vec3_1, _vec3_2);
		}

		this.rotation.set(
			_vec3_1.x, _vec3_2.x, _vec3_3.x,
			_vec3_1.y, _vec3_2.y, _vec3_3.y,
			_vec3_1.z, _vec3_2.z, _vec3_3.z
		);

		const halfSize = _vec3_1.set(scaleX, scaleY, scaleZ);
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

		let sx = _vec3_1.set(e[0], e[1], e[2]).getLength();
		const sy = _vec3_1.set(e[4], e[5], e[6]).getLength();
		const sz = _vec3_1.set(e[8], e[9], e[10]).getLength();

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

		const center = this.box.getCenter(_vec3_1);
		const halfSize = this.box.getSize(_vec3_2).multiplyScalar(0.5);

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
		const center = this.box.getCenter(_vec3_1);
		const min = _vec3_2.subVectors(this.box.min, center);
		const max = _vec3_3.subVectors(this.box.max, center);

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
		const center = this.box.getCenter(_vec3_1);
		const worldMin = _vec3_2.subVectors(this.box.min, center).applyMatrix3(this.rotation).add(center);
		const worldMax = _vec3_3.subVectors(this.box.max, center).applyMatrix3(this.rotation).add(center);

		_vec3_1.set(0, 0, 1).applyMatrix3(this.rotation).normalize();
		planes[0].setFromNormalAndCoplanarPoint(_vec3_1, worldMin);
		planes[1].setFromNormalAndCoplanarPoint(_vec3_1, worldMax);
		planes[1].normal.negate();
		planes[1].constant *= -1;

		_vec3_1.set(0, 1, 0).applyMatrix3(this.rotation).normalize();
		planes[2].setFromNormalAndCoplanarPoint(_vec3_1, worldMin);
		planes[3].setFromNormalAndCoplanarPoint(_vec3_1, worldMax);
		planes[3].normal.negate();
		planes[3].constant *= -1;

		_vec3_1.set(1, 0, 0).applyMatrix3(this.rotation).normalize();
		planes[4].setFromNormalAndCoplanarPoint(_vec3_1, worldMin);
		planes[5].setFromNormalAndCoplanarPoint(_vec3_1, worldMax);
		planes[5].normal.negate();
		planes[5].constant *= -1;

		return planes;
	}

	containsPoint(point) {
		const obbStruct = getOBBStruct(this, a);

		const v = _vec3_1.subVectors(point, obbStruct.c);

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

		const v = _vec3_1.subVectors(point, obbStruct.c);

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

		const v1 = _vec3_1.subVectors(b.c, a.c);

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
		const center = this.box.getCenter(_vec3_1);

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

const _vec3_1 = new Vector3();
const _vec3_2 = new Vector3();
const _vec3_3 = new Vector3();
const _mat3_1 = new Matrix3();

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