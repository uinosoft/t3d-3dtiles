import { GLTFLoader } from '../libs/glTF/GLTFLoader.js';

import { ReferenceParser } from '../libs/glTF/parsers/ReferenceParser.js';
import { Validator } from '../libs/glTF/parsers/Validator.js';
import { BufferParser } from '../libs/glTF/parsers/BufferParser.js';
import { BufferViewParser } from '../libs/glTF/parsers/BufferViewParser.js';
import { ImageParser } from '../libs/glTF/parsers/ImageParser.js';
import { TextureParser } from '../libs/glTF/parsers/TextureParser.js';
import { AccessorParser } from '../libs/glTF/parsers/AccessorParser.js';
import { NodeParser } from '../libs/glTF/parsers/NodeParser.js';
import { SkinParser } from '../libs/glTF/parsers/SkinParser.js';
import { SceneParser } from '../libs/glTF/parsers/SceneParser.js';
import { AnimationParser } from '../libs/glTF/parsers/AnimationParser.js';

import { HeaderParser } from './parsers/HeaderParser.js';
import { TableParser } from './parsers/TableParser.js';
import { I3DMParser } from './parsers/i3dm/I3DMParser.js';
import { MaterialParser } from './parsers/i3dm/MaterialParser.js';
import { PrimitiveParser } from './parsers/i3dm/PrimitiveParser.js';
import { I3DMRootParser } from './parsers/i3dm/I3DMRootParser.js';

/**
 * I3DMLoader is a loader for the I3DM format.
 */
export class I3DMLoader extends GLTFLoader {

	constructor(manager) {
		super(manager, [
			HeaderParser, // insert HeaderParser
			TableParser, // insert TableParser
			I3DMParser, // insert I3DMParser
			ReferenceParser,
			Validator,
			BufferParser,
			BufferViewParser,
			ImageParser,
			TextureParser,
			MaterialParser, // replace MaterialParser
			AccessorParser,
			PrimitiveParser, // replace PrimitiveParser
			NodeParser,
			SkinParser,
			SceneParser,
			AnimationParser,
			I3DMRootParser // insert I3DMSetupParser
		]);
	}

}
