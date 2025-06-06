import { GLTFUtils } from 't3d/addons/loaders/glTF/GLTFUtils.js';

export class I3DMParser {

	static parse(context, loader) {
		const bodyBytes = new Uint8Array(
			context.options.buffer,
			context.batchTableEnd,
			context.header.byteLength - context.batchTableEnd
		);

		let promise = null;
		if (context.header.gltfFormat === 1) {
			promise = Promise.resolve(bodyBytes);
		} else {
			const externalUri = GLTFUtils.resolveURL(GLTFUtils.decodeText(bodyBytes), context.path);
			promise = loader.loadFile(externalUri, 'arraybuffer').then(buffer => new Uint8Array(buffer));
		}

		return promise.then(glbBytes => {
			const glbData = GLTFUtils.parseGLB(glbBytes.slice().buffer);

			context.gltf = glbData.gltf;
			context.buffers = glbData.buffers;
		});
	}

}
