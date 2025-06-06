import { GLTFLoader } from 't3d/addons/loaders/glTF/GLTFLoader.js';

import { ReferenceParser } from 't3d/addons/loaders/glTF/parsers/ReferenceParser.js';
import { Validator } from 't3d/addons/loaders/glTF/parsers/Validator.js';
import { BufferParser } from 't3d/addons/loaders/glTF/parsers/BufferParser.js';
import { BufferViewParser } from 't3d/addons/loaders/glTF/parsers/BufferViewParser.js';
import { ImageParser } from 't3d/addons/loaders/glTF/parsers/ImageParser.js';
import { TextureParser } from 't3d/addons/loaders/glTF/parsers/TextureParser.js';
import { MaterialParser } from 't3d/addons/loaders/glTF/parsers/MaterialParser.js';
import { AccessorParser } from 't3d/addons/loaders/glTF/parsers/AccessorParser.js';
import { PrimitiveParser } from 't3d/addons/loaders/glTF/parsers/PrimitiveParser.js';
import { NodeParser } from 't3d/addons/loaders/glTF/parsers/NodeParser.js';
import { SkinParser } from 't3d/addons/loaders/glTF/parsers/SkinParser.js';
import { SceneParser } from 't3d/addons/loaders/glTF/parsers/SceneParser.js';
import { AnimationParser } from 't3d/addons/loaders/glTF/parsers/AnimationParser.js';

import { HeaderParser } from './parsers/HeaderParser.js';
import { TableParser } from './parsers/TableParser.js';
import { B3DMParser } from './parsers/b3dm/B3DMParser.js';
import { B3DMRootParser } from './parsers/b3dm/B3DMRootParser.js';

import { KHR_techniques_webgl } from './extensions/KHR_techniques_webgl.js';

/**
 * B3DMLoader is a loader for the B3DM format.
 */
export class B3DMLoader extends GLTFLoader {

	constructor(manager) {
		super(manager, [
			HeaderParser, // insert HeaderParser
			TableParser, // insert TableParser
			B3DMParser, // insert B3DMParser
			ReferenceParser,
			Validator,
			BufferParser,
			BufferViewParser,
			ImageParser,
			TextureParser,
			MaterialParser,
			AccessorParser,
			PrimitiveParser,
			NodeParser,
			SkinParser,
			SceneParser,
			AnimationParser,
			B3DMRootParser // insert B3DMRootParser
		]);

		this.extensions.set('KHR_techniques_webgl', KHR_techniques_webgl);

		this.autoParseConfig.materials.push('KHR_techniques_webgl');
	}

}
