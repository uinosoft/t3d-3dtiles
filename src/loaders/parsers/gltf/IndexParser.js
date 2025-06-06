import { GLTFUtils } from 't3d/addons/loaders/glTF/GLTFUtils.js';

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
	return url.substring(url.lastIndexOf('.') + 1) === 'glb';
};