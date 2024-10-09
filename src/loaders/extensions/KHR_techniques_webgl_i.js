import { KHR_techniques_webgl } from './KHR_techniques_webgl.js';
import { InstancedPBRMaterial } from '../../materials/InstancedPBRMaterial.js';

export class KHR_techniques_webgl_i extends KHR_techniques_webgl {

	static getMaterial() {
		return new InstancedPBRMaterial();
	}

}