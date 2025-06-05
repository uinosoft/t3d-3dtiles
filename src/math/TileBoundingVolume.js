import { Sphere, Vector2, Vector3 } from 't3d';
import { TileOBB } from './TileOBB.js';
import { EllipsoidRegion } from './EllipsoidRegion.js';

export class TileBoundingVolume {

	constructor() {
		this.sphere = null;
		this.obb = null;
		this.region = null;
	}

	setOBBData(data, transform) {
		const obb = new TileOBB();

		obb.setFromCenterAndAxes(
			_vec3_4.set(data[0], data[1], data[2]),
			_vec3_1.set(data[3], data[4], data[5]),
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
			if (ray.intersectSphere(sphere, _vec3_1)) {
				sphereDistSq = sphere.containsPoint(ray.origin) ? 0 : ray.origin.distanceToSquared(_vec3_1);
			}
		}

		if (obb) {
			if (obb.intersectRay(ray, _vec3_1)) {
				obbDistSq = obb.containsPoint(ray.origin) ? 0 : ray.origin.distanceToSquared(_vec3_1);
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

const _vec3_1 = new Vector3();
const _vec3_2 = new Vector3();
const _vec3_3 = new Vector3();
const _vec3_4 = new Vector3();