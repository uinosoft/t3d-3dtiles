import { KHR_materials_clearcoat } from '../../libs/glTF/extensions/KHR_materials_clearcoat.js';
import { InstancedPBRMaterial } from '../../materials/InstancedPBRMaterial.js';

/**
 * Clearcoat Materials Extension
 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_materials_clearcoat
 */
export class KHR_materials_clearcoat_i extends KHR_materials_clearcoat {

	static getMaterial() {
		return new InstancedPBRMaterial();
	}

}