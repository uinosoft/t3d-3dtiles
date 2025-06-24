import { Matrix4, Vector2, Vector3, MathUtils } from 't3d';
import { Ellipsoid } from './Ellipsoid.js';

export class EllipsoidRegion extends Ellipsoid {

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
			const midLat = MathUtils.mapLinear(0.5, 0, 1, latRange.x, latRange.y);
			const midLon = MathUtils.mapLinear(0.5, 0, 1, lonRange.x, lonRange.y);

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
		_invMatrix.setFromMatrix3(target.rotation).inverse();

		const points = this._getPoints(true);

		// get the center of the region
		target.box.makeEmpty();
		for (let i = 0, l = points.length; i < l; i++) {
			_center.copy(points[i]).applyMatrix4(_invMatrix);
			target.box.expandByPoint(_center);
		}
		target.box.getCenter(_center);
		_center.applyMatrix3(target.rotation);

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
		const { latRange, lonRange, heightRange } = this;

		const midLat = MathUtils.mapLinear(0.5, 0, 1, latRange.x, latRange.y);
		const midLon = MathUtils.mapLinear(0.5, 0, 1, lonRange.x, lonRange.y);

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
			const height = MathUtils.mapLinear(z, 0, 1, heightRange.x, heightRange.y);
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
const _center = new Vector3();
const _invMatrix = new Matrix4();

const PI = Math.PI;
const HALF_PI = PI / 2;

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