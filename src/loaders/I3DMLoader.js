import { GLTFLoader } from 't3d/addons/loaders/glTF/GLTFLoader.js';

import { ReferenceParser } from 't3d/addons/loaders/glTF/parsers/ReferenceParser.js';
import { Validator } from 't3d/addons/loaders/glTF/parsers/Validator.js';
import { BufferParser } from 't3d/addons/loaders/glTF/parsers/BufferParser.js';
import { BufferViewParser } from 't3d/addons/loaders/glTF/parsers/BufferViewParser.js';
import { ImageParser } from 't3d/addons/loaders/glTF/parsers/ImageParser.js';
import { TextureParser } from 't3d/addons/loaders/glTF/parsers/TextureParser.js';
import { AccessorParser } from 't3d/addons/loaders/glTF/parsers/AccessorParser.js';
import { NodeParser } from 't3d/addons/loaders/glTF/parsers/NodeParser.js';
import { SkinParser } from 't3d/addons/loaders/glTF/parsers/SkinParser.js';
import { SceneParser } from 't3d/addons/loaders/glTF/parsers/SceneParser.js';
import { AnimationParser } from 't3d/addons/loaders/glTF/parsers/AnimationParser.js';

import { HeaderParser } from './parsers/HeaderParser.js';
import { TableParser } from './parsers/TableParser.js';
import { I3DMParser } from './parsers/i3dm/I3DMParser.js';
import { MaterialParser } from './parsers/i3dm/MaterialParser.js';
import { PrimitiveParser } from './parsers/i3dm/PrimitiveParser.js';
import { I3DMRootParser } from './parsers/i3dm/I3DMRootParser.js';

import { KHR_techniques_webgl_i } from './extensions/KHR_techniques_webgl_i.js';
import { KHR_materials_unlit_i } from './extensions/KHR_materials_unlit_i.js';
import { KHR_materials_pbrSpecularGlossiness_i } from './extensions/KHR_materials_pbrSpecularGlossiness_i.js';
import { KHR_materials_clearcoat_i } from './extensions/KHR_materials_clearcoat_i.js';

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

		this.extensions.set('KHR_techniques_webgl', KHR_techniques_webgl_i);

		this.extensions.set('KHR_materials_unlit', KHR_materials_unlit_i);
		this.extensions.set('KHR_materials_pbrSpecularGlossiness', KHR_materials_pbrSpecularGlossiness_i);
		this.extensions.set('KHR_materials_clearcoat', KHR_materials_clearcoat_i);

		this.autoParseConfig.materials.push('KHR_techniques_webgl');
	}

}
