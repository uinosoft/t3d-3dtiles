import { Spherical, Vector3 } from 't3d';
import { swapToGeoFrame, latitudeToSphericalPhi } from './GeoUtils.js';

export class Ellipsoid {

	constructor(radius = new Vector3(1, 1, 1)) {
		this.radius = radius;
	}

	getCartographicToPosition(lat, lon, height, target) {
		// From Cesium function Ellipsoid.cartographicToCartesian
		// https://github.com/CesiumGS/cesium/blob/665ec32e813d5d6fe906ec3e87187f6c38ed5e49/packages/engine/Source/Core/Ellipsoid.js#L396
		this.getCartographicToNormal(lat, lon, _norm);

		const radius = this.radius;
		_vec3_1.copy(_norm);
		_vec3_1.x *= radius.x ** 2;
		_vec3_1.y *= radius.y ** 2;
		_vec3_1.z *= radius.z ** 2;

		const gamma = Math.sqrt(_norm.dot(_vec3_1));
		_vec3_1.multiplyScalar(1 / gamma);

		return target.copy(_vec3_1).addScaledVector(_norm, height);
	}

	getCartographicToNormal(lat, lon, target) {
		_spherical.set(1, latitudeToSphericalPhi(lat), lon);
		target.setFromSpherical(_spherical).normalize();

		// swap frame from the t3d.js frame to the geo coord frame
		swapToGeoFrame(target);
		return target;
	}

}

const _norm = new Vector3();
const _spherical = new Spherical();

const _vec3_1 = new Vector3();