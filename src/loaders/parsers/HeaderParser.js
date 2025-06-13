import { getUrlExtension } from '../../utils/urlExtension.js';

/**
 * This parser is used to parse the header of a 3D Tiles resource.
 * For 'b3dm', 'i3dm', 'pnts' and 'cmpt' formats.
 */
export class HeaderParser {

	static parse(context, loader) {
		const buffer = context.options.buffer;

		// TODO: this should be able to take a uint8array with an offset and length
		const dataView = new DataView(buffer);

		// 32-byte header for i3dm and 28-byte header for the others.

		// 4 bytes for the magic, can be 'b3dm', 'i3dm', 'pnts', 'cmpt' for now.
		// 'vctr' is not supported yet.

		const magic =
            String.fromCharCode(dataView.getUint8(0)) +
            String.fromCharCode(dataView.getUint8(1)) +
            String.fromCharCode(dataView.getUint8(2)) +
            String.fromCharCode(dataView.getUint8(3));

		const urlExtension = getUrlExtension(context.url);

		if (magic !== urlExtension) {
			throw `Not a ${urlExtension} type resource, with url ${context.url}!`;
		}

		// 4 bytes for the version number.

		const version = dataView.getUint32(4, true);

		if (version !== 1) {
			throw `${urlExtension} version must be 1, with url ${context.url}!`;
		}

		// 4 bytes for the byte length of the entire tile content.

		const byteLength = dataView.getUint32(8, true);

		if (byteLength !== buffer.byteLength) {
			throw `${urlExtension} data byte length check failed, with url ${context.url}!`;
		}

		// output the header information to the context and return
		// if the tile content is cmpt.

		if (urlExtension === 'cmpt') {
			// 4 bytes for the tiles length
			const tilesLength = dataView.getUint32(12, true);

			context.header = {
				magic,
				version,
				byteLength,
				tilesLength
			};

			return;
		}

		// 4 bytes for the byte length of the feature table JSON.

		const featureTableJSONByteLength = dataView.getUint32(12, true);

		// 4 bytes for the byte length of the feature table binary.

		const featureTableBinaryByteLength = dataView.getUint32(16, true);

		// 4 bytes for the byte length of the batch table JSON.

		const batchTableJSONByteLength = dataView.getUint32(20, true);

		// 4 bytes for the byte length of the batch table binary.

		const batchTableBinaryByteLength = dataView.getUint32(24, true);

		// 4 bytes for the gltf format if the tile content format is i3dm.

		let gltfFormat = null;
		if (urlExtension === 'i3dm') {
			gltfFormat = dataView.getUint32(28, true);
		}

		// output the header information to the context.

		context.header = {
			magic,
			version,
			byteLength,
			featureTableJSONByteLength,
			featureTableBinaryByteLength,
			batchTableJSONByteLength,
			batchTableBinaryByteLength,
			gltfFormat
		};

		return;
	}

}