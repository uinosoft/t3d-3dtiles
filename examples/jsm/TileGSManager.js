import { Vector3, Quaternion, Matrix4, Matrix3 } from 't3d';

const position = new Vector3();
const quaternion = new Quaternion();
const scale = new Vector3();
const mat4_1 = new Matrix4();
const mat4_2 = new Matrix4();
const mat3_1 = new Matrix3();
const mat3_2 = new Matrix3();

export class TileGSManager {

	constructor(GSMesh) {
		this._gsDataMap = new Map();
		this._gsMap = new Map();
		this.GSMesh = GSMesh;
	}

	updateTile(tile) {
		const internalData = this.GSMesh.internalData;
		tile.scene.traverse(child => {
			if (child.splatBuffer) {
				if (tile.visible) {
					internalData.vertexCount += child.splatBuffer.vertexCount;
					this._gsMap.set(child.uuid, child);
				} else {
					internalData.vertexCount -= child.splatBuffer.vertexCount;
					this._gsMap.delete(child.uuid);
				}
			}
		});

		if (internalData.vertexCount > internalData.positions.length / 3) {
			internalData.positions = new Float32Array(internalData.vertexCount * 3);
			internalData.colors = new Uint8Array(internalData.vertexCount * 4);
			internalData.covariances = new Float32Array(internalData.vertexCount * 6);
		}

		let offset = 0;
		this._gsMap.forEach((child, uuid) => {
			if (child.splatBuffer) {
				if (!this._gsDataMap.has(uuid)) {
					this._gsDataMap.set(uuid, storeGSToTileSpace(child.splatBuffer, child.worldMatrix));
				}
				const data = this._gsDataMap.get(uuid);
				internalData.positions.set(data.position, offset * 3);
				internalData.colors.set(data.color, offset * 4);
				internalData.covariances.set(data.covariances, offset * 6);
				offset += data.vertexCount;
			}
		});

		this.GSMesh.updateData();
	}

}

function storeGSToTileSpace(splatBuffer, matrix) {
	const data = {
		position: new Float32Array(splatBuffer.vertexCount * 3),
		color: new Uint8Array(splatBuffer.vertexCount * 4),
		covariances: new Float32Array(splatBuffer.vertexCount * 6),
		vertexCount: splatBuffer.vertexCount
	};
	let offset = 0;
	data.color.set(splatBuffer.colors, offset * 4);
	for (let i = 0; i < splatBuffer.vertexCount; i++) {
		position.fromArray(splatBuffer.positions, i * 3);
		quaternion.fromArray(splatBuffer.rotations, i * 4);
		scale.fromArray(splatBuffer.scales, i * 3);
		mat4_1.compose(position, quaternion, scale);
		mat4_2.copy(matrix);
		mat4_1.multiply(mat4_2);
		mat4_1.decompose(position, quaternion, scale);
		data.position.set(position.toArray(), offset * 3);

		mat3_1.set(
			scale.x, 0, 0,
			0, scale.y, 0,
			0, 0, scale.z
		);
		quaternion.set(quaternion.w, quaternion.x, quaternion.y, quaternion.z);
		mat4_1.makeRotationFromQuaternion(quaternion);
		mat3_2.setFromMatrix4(mat4_1);

		mat3_2.multiply(mat3_1);
		mat3_1.copy(mat3_2).transpose().premultiply(mat3_2);
		data.covariances[i * 6 + 0] = mat3_1.elements[0]; // xx
		data.covariances[i * 6 + 1] = mat3_1.elements[3]; // xy
		data.covariances[i * 6 + 2] = mat3_1.elements[6]; // xz
		data.covariances[i * 6 + 3] = mat3_1.elements[4]; // yy
		data.covariances[i * 6 + 4] = mat3_1.elements[7]; // yz
		data.covariances[i * 6 + 5] = mat3_1.elements[8]; // zz

		offset += 1;
	}
	return data;
}