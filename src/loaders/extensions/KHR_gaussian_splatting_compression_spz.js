
import { MathUtils } from 't3d';

export class KHR_gaussian_splatting_compression_spz {

	static getGeometry(bufferView, spzLoader) {
		return new Promise(function(resolve, reject) {
			const bufferViewTypedArray = new Uint8Array(bufferView);
			spzLoader(bufferViewTypedArray, {
				unpackOptions: { coordinateSystem: 'LDB' }
			}).then(gcloud => {
				const splatBuffer = KHR_gaussian_splatting_compression_spz.convertInternalDataToSplat(gcloud.numPoints, gcloud.positions, gcloud.rotations, gcloud.scales, gcloud.alphas, gcloud.colors);
				splatBuffer._isSplatBuffer = true;
				resolve(splatBuffer);
			}).catch(error => {
				reject(new Error('Failed to load SPZ: ' + error.message));
			});
		});
	}

	static convertInternalDataToSplat(vertexCount, positions, rotations, scales, alphas, colors) {
		const rotationsOut = new Float32Array(vertexCount * 4);
		const positionOut = new Float32Array(vertexCount * 3);
		for (let i = 0; i < vertexCount; i++) {
			const off = i * 4;
			const x = rotations[off + 0];
			const y = rotations[off + 1];
			const z = rotations[off + 2];
			const w = rotations[off + 3];

			rotationsOut[off + 0] = w;
			rotationsOut[off + 1] = x;
			rotationsOut[off + 2] = y;
			rotationsOut[off + 3] = z;
			const posOff = i * 3;
			positionOut[posOff + 0] = -positions[posOff + 0];
			positionOut[posOff + 1] = -positions[posOff + 1];
			positionOut[posOff + 2] = positions[posOff + 2];
		}

		const colorsOut = new Uint8Array(vertexCount * 4);
		for (let i = 0; i < vertexCount; i++) {
			colorsOut[i * 4 + 0] = MathUtils.clamp(colors[i * 3 + 0] * 255, 0, 255);
			colorsOut[i * 4 + 1] = MathUtils.clamp(colors[i * 3 + 1] * 255, 0, 255);
			colorsOut[i * 4 + 2] = MathUtils.clamp(colors[i * 3 + 2] * 255, 0, 255);
			colorsOut[i * 4 + 3] = MathUtils.clamp(alphas[i] * 255, 0, 255);
		}

		return {
			vertexCount,
			positions: positionOut,
			rotations: rotationsOut,
			scales,
			colors: colorsOut
		};
	}

}