import { Object3D } from 't3d';

export class CMPTRootParser {

	static parse(context, loader) {
		const { tiles, options, path } = context;

		const adjustmentTransform = options.adjustmentTransform;

		const promises = [];

		for (const i in tiles) {
			const { type, buffer } = tiles[i];

			const config = {
				fetchOptions: options.fetchOptions,
				path,
				buffer: buffer.slice().buffer
			};
			if (type === 'b3dm' || type === 'i3dm') {
				config.adjustmentTransform = adjustmentTransform;
			}

			const _loader = loader._loaders.get(type);
			if (_loader) {
				promises.push(_loader.load(`${path}/temp.${type}`, config));
			}
		}

		return Promise.all(promises).then(results => {
			const group = new Object3D();

			results.forEach(result => {
				group.add(result.root);
			});

			return {
				tiles: results,
				root: group
			};
		});
	}

}