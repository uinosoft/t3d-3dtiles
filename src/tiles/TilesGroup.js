import { Object3D, Matrix4 } from 't3d';

// Specialization of "Group" that only updates world matrices of children if
// the transform has changed since the last update and ignores the "force"
// parameter under the assumption that the children tiles will not move.
const tempMat = new Matrix4();
export class TilesGroup extends Object3D {

	constructor(tilesRenderer) {
		super();
		this.isTilesGroup = true;
		this.name = 'TilesRenderer.TilesGroup';
		this.tilesRenderer = tilesRenderer;
		this.worldMatrixInverse = new Matrix4();
	}

	raycast(ray, intersects) {
		// returning "false" ends raycast traversal
		if (this.tilesRenderer.optimizeRaycast) {
			this.tilesRenderer.raycast(ray, intersects);
			return false;
		}

		return true;
	}

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
				this.worldMatrixInverse.copy(tempMat).invert();

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