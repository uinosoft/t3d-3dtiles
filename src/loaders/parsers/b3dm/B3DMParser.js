import { GLTFUtils } from 't3d/addons/loaders/glTF/GLTFUtils.js';

export class B3DMParser {

	static parse(context, loader) {
		const glbBytes = new Uint8Array(
			context.options.buffer,
			context.batchTableEnd,
			context.header.byteLength - context.batchTableEnd
		);

		const glbData = GLTFUtils.parseGLB(glbBytes.slice().buffer);

		context.gltf = glbData.gltf;
		context.buffers = glbData.buffers;
	}

}
