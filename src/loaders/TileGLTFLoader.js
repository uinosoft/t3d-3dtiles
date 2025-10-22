import { GLTFLoader } from 't3d/addons/loaders/glTF/GLTFLoader.js';
import { KHR_gaussian_splatting } from './extensions/KHR_gaussian_splatting.js';
import { NodeParser } from './parsers/gltf/NodeParser.js';
import { PrimitiveParser } from './parsers/gltf/PrimitiveParser.js';
import { IndexParser } from './parsers/gltf/IndexParser.js';
import { AccessorParser } from './parsers/gltf/AccessorParser.js';
/**
 * TileGLTFLoader is a gltf loader
 * that extends the default gltf loader with additional parsers
 * to support the tile gltf format.
 */
export class TileGLTFLoader extends GLTFLoader {

	constructor(manager) {
		super(manager);
		this.replaceParser(IndexParser, 0);
		this.replaceParser(AccessorParser, 8);
		this.replaceParser(PrimitiveParser, 9);
		this.replaceParser(NodeParser, 10);
	    this.extensions.set('KHR_gaussian_splatting', KHR_gaussian_splatting);

		this.autoParseConfig.materials.push('KHR_gaussian_splatting');
	}

	setSPZLoader(spzLoader) {
		this._spzLoader = spzLoader;
		return this;
	}

	getSPZLoader() {
		return this._spzLoader;
	}

}