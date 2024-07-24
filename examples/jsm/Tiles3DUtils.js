import * as t3d from 't3d';

export function setTiles3DToOrigin(tiles3D, parent) {
	const sphere = new t3d.Sphere();
	tiles3D.getBoundingSphere(sphere);

	const position = sphere.center;

	const surfaceDirection = position.clone().normalize();
	const up = new t3d.Vector3(0, 1, 0);
	const rotationToNorthPole = rotationBetweenDirections(surfaceDirection, up);

	parent.quaternion.x = rotationToNorthPole.x;
	parent.quaternion.y = rotationToNorthPole.y;
	parent.quaternion.z = rotationToNorthPole.z;
	parent.quaternion.w = rotationToNorthPole.w;

	tiles3D.position.sub(position);
}

function rotationBetweenDirections(dir1, dir2) {
	const rotation = new t3d.Quaternion();
	const a = new t3d.Vector3();
	a.crossVectors(dir1, dir2);
	rotation.x = a.x;
	rotation.y = a.y;
	rotation.z = a.z;
	rotation.w = 1 + dir1.clone().dot(dir2);
	rotation.normalize();

	return rotation;
}