export class GLTFExtensionsPlugin {

	constructor(options) {
		options = {
			dracoLoader: null,
			ktxLoader: null,
			...options
		};

		this.tiles = null;

		this.dracoLoader = options.dracoLoader;
		this.ktxLoader = options.ktxLoader;
	}

	init(tiles) {
		this.tiles = tiles;

		if (this.dracoLoader) {
			tiles._loaders.get('b3dm').setDRACOLoader(this.dracoLoader);
			tiles._loaders.get('i3dm').setDRACOLoader(this.dracoLoader);
			tiles._loaders.get('cmpt').setDRACOLoader(this.dracoLoader);
			tiles._loaders.get('gltf').setDRACOLoader(this.dracoLoader);
		}

		if (this.ktxLoader) {
			this._loaders.get('b3dm').setKTX2Loader(this.ktxLoader);
			this._loaders.get('i3dm').setKTX2Loader(this.ktxLoader);
			this._loaders.get('cmpt').setKTX2Loader(this.ktxLoader);
			this._loaders.get('gltf').setKTX2Loader(this.ktxLoader);
		}
	}

}