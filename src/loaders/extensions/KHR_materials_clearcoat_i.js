import { KHR_materials_clearcoat } from 't3d/addons/loaders/glTF/extensions/KHR_materials_clearcoat.js';
import { InstancedPBRMaterial } from '../../materials/InstancedPBRMaterial.js';

export class KHR_materials_clearcoat_i extends KHR_materials_clearcoat {

	static getMaterial() {
		return new InstancedPBRMaterial();
	}

}