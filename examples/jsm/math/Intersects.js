import { Triangle, Matrix4, Sphere } from 't3d';
import { OBB } from 't3d-3dtiles';
import { LineSegment } from './LineSegment.js';

/**
 * Get the closest intersection point between a path and the 3D Tiles.
 * The path is defined by a line segment with a distance threshold,
 * so the actual area is a capsule.
 * @param {Object} tiles - The 3D Tiles Renderer.
 * @param {LineSegment} lineSegment - The line segment.
 * @param {Number} minDistance - The distance threshold.
 * @return {Object|null} - The intersection object or null.
 */
export function getIntersectWithPath(tiles, lineSegment, minDistance) {
	// Get obb in tiles root space
	const obb = _lineSegment_1.copy(lineSegment)
		// .applyMatrix4(tiles.group.worldMatrixInverse)
		.getOBB(minDistance, _obb_1);
	// TODO: obb.applyMatrix4 has bug, so don't use this until it's fixed
	// const obb = lineSegment.getOBB(minDistance, _obb_1);
	// obb.applyMatrix4(tiles.group.worldMatrixInverse);

	let intersectWithTile = null;

	tiles.activeTiles.forEach(tile => {
		if (!intersectsOBB(tile.cached.boundingVolume, obb)) return;

		// TODO: traverse scene
		const positions = tile.cached.geometry[0].attributes.a_Position.buffer.array;
		const indexes = tile.cached.geometry[0].index.buffer.array;

		// TODO: use mesh worldMatrix
		const toWorldMatrix = tile.cached.scene.children[0].worldMatrix;
		const toTileRootMatrix = _mat4_2.copy(toWorldMatrix).premultiply(_mat4_1);

		for (let i = 0; i < indexes.length; i = i + 3) {
			// TODO: more complex buffer structure, like non-indexed buffer

			const { a, b, c } = _triangle_1;

			const indexA = indexes[i], indexB = indexes[i + 1], indexC = indexes[i + 2];
			a.fromArray(positions, 3 * indexA);
			b.fromArray(positions, 3 * indexB);
			c.fromArray(positions, 3 * indexC);

			// Early out if obb doesn't intersect with the triangle's bounding sphere
			// Should compute triangle circumcenter to get a more accurate sphere?
			_sphere_1.setFromPoints(_triangle_1_array).applyMatrix4(toTileRootMatrix);
			if (!obb.intersectsSphere(_sphere_1)) continue;

			a.applyMatrix4(toWorldMatrix);
			b.applyMatrix4(toWorldMatrix);
			c.applyMatrix4(toWorldMatrix);

			_lineSegment_1.startPoint.set(0, 0, 0);
			_lineSegment_1.endPoint.set(0, 0, 0);
			lineSegment.closestLineToTriangle(_triangle_1, _lineSegment_1);

			// TODO: use distance squared
			const closestDis = _lineSegment_1.getLength();
			if (closestDis < minDistance) {
				if (!intersectWithTile || closestDis < intersectWithTile.distance) {
					intersectWithTile = {
						tile,
						distance: closestDis,
						startPoint: _lineSegment_1.startPoint.toArray(),
						endPoint: _lineSegment_1.endPoint.toArray()
					};
				}
			}
		}
	});

	return intersectWithTile;
}

const _obb_1 = new OBB();
const _mat4_1 = new Matrix4();
const _mat4_2 = new Matrix4();
const _lineSegment_1 = new LineSegment();
const _triangle_1 = new Triangle();
const _triangle_1_array = [_triangle_1.a, _triangle_1.b, _triangle_1.c];
const _sphere_1 = new Sphere();

function intersectsOBB(boundingVolume, obb) {
	const _sphere = boundingVolume.sphere;
	const _obb = boundingVolume.obb;

	// Early out if we don't hit this tile sphere
	if (_sphere && !obb.intersectsSphere(_sphere)) {
		return false;
	}

	// Early out if we don't this this tile box
	if (_obb && !_obb.intersectsOBB(obb)) {
		return false;
	}

	return true;
}