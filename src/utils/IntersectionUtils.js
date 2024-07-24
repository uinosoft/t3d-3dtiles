import { Matrix4, Ray, Vector3 } from 't3d';

export const raycastTraverse = (tile, tiles3D, ray, intersects, localRay = null) => {
	const { activeTiles } = tiles3D;
	const boundingVolume = tile.cached.boundingVolume;

	// reuse the ray when traversing the tree
	if (localRay === null) {
		localRay = _ray_1;
		_mat4_1.copy(tiles3D.worldMatrix).inverse();
		localRay.copy(ray).applyMatrix4(_mat4_1);
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
export const raycastTraverseFirstHit = (tile, tiles3D, ray, localRay = null) => {
	const { activeTiles } = tiles3D;

	// reuse the ray when traversing the tree
	if (localRay === null) {
		localRay = _ray_1;
		_mat4_1.copy(tiles3D.worldMatrix).inverse();
		localRay.copy(ray).applyMatrix4(_mat4_1);
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

		if (boundingVolume.intersectRay(localRay, _vec3_1)) {
			_vec3_1.applyMatrix4(tiles3D.worldMatrix);
			array.push({
				distance: _vec3_1.distanceToSquared(ray.origin),
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

export const distanceSort = (a, b) => {
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

const _mat4_1 = new Matrix4();
const _ray_1 = new Ray();
const _vec3_1 = new Vector3();

const _hitArray = [];