import { BatchTable } from '../../utils/BatchTable.js';
import { FeatureTable } from '../../utils/FeatureTable.js';

/**
 * This parser is used to parse the feature table and batch table of a 3D Tiles resource.
 * For 'b3dm', 'i3dm' and 'pnts' formats.
 */
export class TableParser {

	static parse(context, loader) {
		const { header, options } = context;
		const buffer = options.buffer;

		const featureTableStart = header.magic === 'i3dm' ? 32 : 28;
		const featureTableEnd = featureTableStart + header.featureTableJSONByteLength + header.featureTableBinaryByteLength;
		const batchTableStart = featureTableEnd;
		const batchTableEnd = batchTableStart + header.batchTableJSONByteLength + header.batchTableBinaryByteLength;

		// parse the feature table

		const featureTableBuffer = buffer.slice(
			featureTableStart,
			featureTableEnd
		);

		const featureTable = new FeatureTable(
			featureTableBuffer,
			0,
			header.featureTableJSONByteLength,
			header.featureTableBinaryByteLength
		);

		// parse the batch table

		let batchSize;
		if (header.magic === 'b3dm') {
			batchSize = featureTable.getData('BATCH_LENGTH');
		} else if (header.magic === 'i3dm') {
			batchSize = featureTable.getData('INSTANCES_LENGTH');
		} else if (header.magic === 'pnts') {
			batchSize = featureTable.getData('BATCH_LENGTH') || featureTable.getData('POINTS_LENGTH');
		} else {
			throw `Unrecognized magic: ${header.magic}!`;
		}

		const batchTableBuffer = buffer.slice(
			batchTableStart,
			batchTableEnd
		);

		const batchTable = new BatchTable(
			batchTableBuffer,
			batchSize,
			0,
			header.batchTableJSONByteLength,
			header.batchTableBinaryByteLength
		);

		// output the tables to the context

		context.featureTable = featureTable;
		context.batchTable = batchTable;
		context.batchTableEnd = batchTableEnd;
	}

}