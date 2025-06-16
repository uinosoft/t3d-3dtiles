import { getUrlExtension } from '../../utilities/urlExtension.js';

export class LoadParser {

	static parse(context, loader) {
		const extension = getUrlExtension(context.url);
		let pr = null;
		if (extension === 'gltf') {
			pr = loader.loadFile(context.url)
				.then(buffer => {
					context.options.buffer = buffer;
				});
		} else {
			pr = loader.loadFile(context.url, 'arraybuffer')
				.then(buffer => {
					context.options.buffer = buffer;
				});
		}
		return pr;
	}

}