import { GLTFLoader } from '../libs/glTF/GLTFLoader.js';

import { HeaderParser } from './parsers/HeaderParser.js';
import { TableParser } from './parsers/TableParser.js';
import { PNTSRootParser } from './parsers/pnts/PNTSRootParser.js';

/**
 * PNTSLoader is a loader for the PNTS format.
 */
export class PNTSLoader extends GLTFLoader {

	constructor(manager) {
		super(manager, [
			HeaderParser,
			TableParser,
			PNTSRootParser
		]);
	}

}
