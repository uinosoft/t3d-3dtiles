import { GLTFLoader } from '../libs/glTF/GLTFLoader.js';

import { ReferenceParser } from '../libs/glTF/parsers/ReferenceParser.js';
import { Validator } from '../libs/glTF/parsers/Validator.js';
import { BufferParser } from '../libs/glTF/parsers/BufferParser.js';
import { BufferViewParser } from '../libs/glTF/parsers/BufferViewParser.js';
import { ImageParser } from '../libs/glTF/parsers/ImageParser.js';
import { TextureParser } from '../libs/glTF/parsers/TextureParser.js';
import { AccessorParser } from '../libs/glTF/parsers/AccessorParser.js';
import { PrimitiveParser } from '../libs/glTF/parsers/PrimitiveParser.js';
import { NodeParser } from '../libs/glTF/parsers/NodeParser.js';
import { SkinParser } from '../libs/glTF/parsers/SkinParser.js';
import { SceneParser } from '../libs/glTF/parsers/SceneParser.js';
import { AnimationParser } from '../libs/glTF/parsers/AnimationParser.js';

import { HeaderParser } from './parsers/HeaderParser.js';
import { TableParser } from './parsers/TableParser.js';
import { B3DMParser } from './parsers/b3dm/B3DMParser.js';
import { MaterialParser } from './parsers/b3dm/MaterialParser.js';
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
			MaterialParser, // replace MaterialParser
			AccessorParser,
			PrimitiveParser,
			NodeParser,
			SkinParser,
			SceneParser,
			AnimationParser,
			B3DMRootParser // insert B3DMRootParser
		]);

		this.extensions.set('KHR_techniques_webgl', KHR_techniques_webgl);
	}

}
