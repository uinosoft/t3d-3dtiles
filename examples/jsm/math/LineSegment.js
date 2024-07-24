import { Vector3, Matrix4, Triangle, Ray } from 't3d';

export class LineSegment {

	constructor(startPoint = new Vector3(), endPoint = new Vector3()) {
		this.startPoint = startPoint;
		this.endPoint = endPoint;
	}

	setPoints(startPoint, endPoint) {
		this.startPoint.copy(startPoint);
		this.endPoint.copy(endPoint);
		return this;
	}

	copy(lineSegment) {
		this.startPoint.copy(lineSegment.startPoint);
		this.endPoint.copy(lineSegment.endPoint);
		return this;
	}

	applyMatrix4(m) {
		this.startPoint.applyMatrix4(m);
		this.endPoint.applyMatrix4(m);
		return this;
	}

	getLength() {
		return this.startPoint.distanceTo(this.endPoint);
	}

	getLengthSquared() {
		return this.startPoint.distanceToSquared(this.endPoint);
	}

	getCenter(target) {
		return target.addVectors(this.startPoint, this.endPoint).multiplyScalar(0.5);
	}

	getOBB(radius, target) {
		const center = this.getCenter(_vec3_1);
		const halfSize = _vec3_2.set(radius, radius, this.getLength() * 0.5 + radius); // expand box size by radius

		const rotationMatrix4 = _mat4_1.identity().lookAtRH(this.startPoint, this.endPoint, _up);

		target.box.min.subVectors(center, halfSize);
		target.box.max.addVectors(center, halfSize);
		target.rotation.setFromMatrix4(rotationMatrix4);

		return target;
	}

	closestLineToLineSegment(lineSegment, target) {
		const lineVec = _lineVec1.subVectors(this.endPoint, this.startPoint);
		const lineVec2 = _lineVec2.subVectors(lineSegment.endPoint, lineSegment.startPoint);
		const r = _lineVec3.subVectors(this.startPoint, lineSegment.startPoint);
		const a = lineVec.getLengthSquared();
		const invA = 1.0 / a;
		const e = lineVec2.getLengthSquared();
		const f = lineVec2.dot(r);
		const EPSILON = 1E-5;
		let s, t;

		if (a <= EPSILON && e <= EPSILON) {
			// Both segments degenerate into points
			target.startPoint.copy(this.startPoint);
			target.endPoint.copy(lineSegment.startPoint);
			return target;
		}

		if (a <= EPSILON) {
			// First segment degenerates into a point
			s = 0.0;
			t = f / e;
			t = Math.min(Math.max(t, 0.0), 1.0);
		} else {
			const c = lineVec.dot(r);
			if (e <= EPSILON) {
				// Second segment degenerates into a point
				t = 0.0;
				s = Math.min(Math.max(-c * invA, 0.0), 1.0);
			} else {
				// The general nondegenerate case starts here
				const b = lineVec.dot(lineVec2);
				const denom = a * e - b * b;
				// If segments not parallel, compute closest point on L1 to L2 and clamp to segment S1. Else pick arbitrary s (here 0)
				if (denom !== 0.0) { s = Math.min(Math.max((b * f - c * e) / denom, 0.0), 1.0) } else { s = 0.0 }
				// Compute point on L2 closest to S1(s) using
				// t = Dot((P1 + D1 * s) - P2, D2) / Dot(D2, D2) = (b * s + f) / e
				t = (b * s + f) / e;
				// If t in [0,1] done. Else clamp t, recompute s for the new value of t using s = Dot((P2 + D2 * t) - P1, D1) / Dot(D1, D1) = (t * b - c) / a and clamp s to [0, 1]
				if (t < 0.0) {
					t = 0.0;
					s = Math.min(Math.max(-c * invA, 0.0), 1.0);
				} else if (t > 1.0) {
					t = 1.0;
					s = Math.min(Math.max((b - c) * invA, 0.0), 1.0);
				}
			}
		}

		target.startPoint.addVectors(this.startPoint, lineVec.multiplyScalar(s));
		target.endPoint.addVectors(lineSegment.startPoint, lineVec2.multiplyScalar(t));
		return target;
	}

	intersectTriangle(triangle, target) {
		const lineVec = _lineVec0.subVectors(this.endPoint, this.startPoint).normalize();
		const ray = _ray.set(this.startPoint, lineVec);
		const backfaceCulling = false;

		const intersect = ray.intersectTriangle(triangle.a, triangle.b, triangle.c, backfaceCulling, target);
		if (intersect === null) return null;

		const distance = this.startPoint.distanceTo(target);
		const t = distance / this.getLength();
		if (t >= 0.0 && t <= 1.0) {
			return true;
		}
		return false;
	}

	closestLineToTriangle(triangle, target) {
		let min = Infinity, d;

		const intersectWithTriangle = this.intersectTriangle(triangle, _vec3_1);
		if (intersectWithTriangle) {
			target.startPoint.copy(_vec3_1);
			target.endPoint.copy(_vec3_1);
			return target;
		} else {
			const isPointAintri = Triangle.containsPoint(this.startPoint, triangle.a, triangle.b, triangle.c);
			const isPointBintri = Triangle.containsPoint(this.endPoint, triangle.a, triangle.b, triangle.c);

			// segment endpoint A and plane of triangle (when A projects inside V0V1V2)
			let computed = false;
			let a = NaN, b = NaN, c = NaN, nd = NaN;
			if (isPointAintri) {
				const lineVec1 = _lineVec1.subVectors(triangle.b, triangle.a);
				const lineVec2 = _lineVec2.subVectors(triangle.c, triangle.a);
				a = lineVec1.y * lineVec2.z - lineVec1.z * lineVec2.y;
				b = lineVec1.z * lineVec2.x - lineVec1.x * lineVec2.z;
				c = lineVec1.x * lineVec2.y - lineVec1.y * lineVec2.x;
				const lineVec3 = _lineVec3.set(a, b, c);
				computed = true;
				nd = -lineVec3.normalize().dot(triangle.a);

				d = lineVec3.dot(this.startPoint) + nd;
				const l = d;
				d *= d;
				if (d < min) {
					min = d;
					target.startPoint.copy(this.startPoint);
					target.endPoint.subVectors(this.startPoint, lineVec3.multiplyScalar(l));
				}
			}

			// segment endpoint B and plane of triangle (when B projects inside V0V1V2)
			if (isPointBintri) {
				const lineVec3 = _lineVec3;
				if (!computed) {
					const lineVec1 = _lineVec1.subVectors(triangle.b, triangle.a);
					const lineVec2 = _lineVec2.subVectors(triangle.c, triangle.a);
					a = lineVec1.y * lineVec2.z - lineVec1.z * lineVec2.y;
					b = lineVec1.z * lineVec2.x - lineVec1.x * lineVec2.z;
					c = lineVec1.x * lineVec2.y - lineVec1.y * lineVec2.x;
					lineVec3.set(a, b, c);
					computed = true;
					nd = -lineVec3.normalize().dot(triangle.a);
				}
				d = lineVec3.dot(this.endPoint) + nd;
				const l = d;
				d *= d;
				if (d < min) {
					min = d;
					target.startPoint.copy(this.endPoint);
					target.endPoint.subVectors(this.endPoint, lineVec3.multiplyScalar(l));
				}
			}

			if (isPointBintri && isPointAintri) {
				return target;
			}

			// AB -> V0V1
			_lineSegmentEdge.setPoints(triangle.a, triangle.b);
			this.closestLineToLineSegment(_lineSegmentEdge, _tampTarget);
			min = _tampTarget.getLength();
			target.startPoint.copy(_tampTarget.startPoint);
			target.endPoint.copy(_tampTarget.endPoint);

			// AB -> V1V2
			_lineSegmentEdge.setPoints(triangle.b, triangle.c);
			this.closestLineToLineSegment(_lineSegmentEdge, _tampTarget);
			d = _tampTarget.getLength();
			if (d < min) {
				min = d;
				target.startPoint.copy(_tampTarget.startPoint);
				target.endPoint.copy(_tampTarget.endPoint);
			}

			// AB -> V2V0
			_lineSegmentEdge.setPoints(triangle.c, triangle.a);
			this.closestLineToLineSegment(_lineSegmentEdge, _tampTarget);
			d = _tampTarget.getLength();
			if (d < min) {
				min = d;
				target.startPoint.copy(_tampTarget.startPoint);
				target.endPoint.copy(_tampTarget.endPoint);
			}

			return target;
		}
	}

}

const _vec3_1 = new Vector3();
const _vec3_2 = new Vector3();
const _mat4_1 = new Matrix4();

const _up = new Vector3(0, 1, 0);

const _ray = new Ray();

const _lineVec0 = new Vector3();

const _lineVec1 = new Vector3();
const _lineVec2 = new Vector3();
const _lineVec3 = new Vector3();

const _lineSegmentEdge = new LineSegment();

const _tampTarget = new LineSegment();