import { KHR_materials_clearcoat } from 't3d/addons/loaders/glTF/extensions/KHR_materials_clearcoat.js';
import { InstancedPBRMaterial } from 't3d/addons/materials/InstancedMaterial.js';

export class KHR_materials_clearcoat_i extends KHR_materials_clearcoat {

	static getMaterial() {
		return new InstancedPBRMaterial();
	}

}