import { GLTFLoader } from '../libs/glTF/GLTFLoader.js';

import { HeaderParser } from './parsers/HeaderParser.js';
import { CMPTParser } from './parsers/cmpt/CMPTParser.js';
import { CMPTRootParser } from './parsers/cmpt/CMPTRootParser.js';

import { B3DMLoader } from './B3DMLoader.js';
import { I3DMLoader } from './I3DMLoader.js';
import { PNTSLoader } from './PNTSLoader.js';

/**
 * CMPTLoader is a loader for the CMPT format.
 */
export class CMPTLoader extends GLTFLoader {

	constructor(manager) {
		super(manager, [
			HeaderParser,
			CMPTParser,
			CMPTRootParser
		]);

		const b3dmLoader = new B3DMLoader(manager);
		const i3dmLoader = new I3DMLoader(manager);
		const pntsLoader = new PNTSLoader(manager);

		this._loaders = new Map([
			['b3dm', b3dmLoader],
			['i3dm', i3dmLoader],
			['pnts', pntsLoader]
		]);
	}

	setDRACOLoader(dracoLoader) {
		for (const loader of this._loaders.values()) {
			loader.setDRACOLoader(dracoLoader);
		}
		return super.setDRACOLoader(dracoLoader);
	}

	setKTX2Loader(ktx2Loader) {
		for (const loader of this._loaders.values()) {
			loader.setKTX2Loader(ktx2Loader);
		}
		return super.setKTX2Loader(ktx2Loader);
	}

}
