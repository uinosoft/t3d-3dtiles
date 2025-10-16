import { KHR_gaussian_splatting_compression_spz } from './KHR_gaussian_splatting_compression_spz.js';
/**
 * KHR_gaussian_splatting extension
 * https://github.com/CesiumGS/glTF/blob/draft-splat-spz/extensions/2.0/Khronos/KHR_gaussian_splatting
 */
export class KHR_gaussian_splatting {

	static getGeometry(isCompressed, extensions, bufferViews, spzLoader) {
		if (isCompressed) {
			const bufferViewIndex = extensions.bufferView;
			const bufferView = bufferViews[bufferViewIndex];
			return KHR_gaussian_splatting_compression_spz.getGeometry(bufferView, spzLoader);
		} else {
			// TATO Uncompressed data processing
		}
	}

}

