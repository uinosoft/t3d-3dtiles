import { GLTFUtils } from 't3d/addons/loaders/glTF/GLTFUtils.js';
import { getUrlExtension } from '../../../core/renderer/utilities/urlExtension.js';

export class IndexParser {

	static parse(context, loader) {
		const { url, options } = context;
		const buffer = options.buffer;

		const _isGLB = isGLB(url);

		if (_isGLB) {
			const glbData = GLTFUtils.parseGLB(buffer);
			context.gltf = glbData.gltf;
			context.buffers = glbData.buffers;
		} else {
			context.gltf = buffer;
		}
	}

}

const isGLB = url => {
	return getUrlExtension(url) === 'glb';
};