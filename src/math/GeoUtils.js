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

export function swapToGeoFrame(target) {
	const { x, y, z } = target;
	target.x = z;
	target.y = x;
	target.z = y;
}

export function latitudeToSphericalPhi(latitude) {
	return -latitude + Math.PI / 2;
}