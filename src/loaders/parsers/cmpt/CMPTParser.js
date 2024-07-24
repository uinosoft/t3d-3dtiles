export class CMPTParser {

	static parse(context, loader) {
		const buffer = context.options.buffer;
		const tilesLength = context.header.tilesLength;

		const tiles = [];
		let offset = 16;

		for (let i = 0; i < tilesLength; i++) {
			const tileView = new DataView(buffer, offset, 12);

			const tileMagic =
				String.fromCharCode(tileView.getUint8(0)) +
				String.fromCharCode(tileView.getUint8(1)) +
				String.fromCharCode(tileView.getUint8(2)) +
				String.fromCharCode(tileView.getUint8(3));

			const tileVersion = tileView.getUint32(4, true);

			const byteLength = tileView.getUint32(8, true);

			const tileBuffer = new Uint8Array(buffer, offset, byteLength);

			tiles.push({
				type: tileMagic,
				buffer: tileBuffer,
				version: tileVersion
			});

			offset += byteLength;
		}

		context.tiles = tiles;
	}

}
