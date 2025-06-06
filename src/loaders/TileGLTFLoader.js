import { GLTFLoader } from 't3d/addons/loaders/glTF/GLTFLoader.js';

import { IndexParser } from './parsers/gltf/IndexParser.js';

/**
 * TileGLTFLoader is a gltf loader
 * that extends the default gltf loader with additional parsers
 * to support the tile gltf format.
 */
export class TileGLTFLoader extends GLTFLoader {

	constructor(manager) {
		super(manager);
		this.replaceParser(IndexParser, 0);
	}

}