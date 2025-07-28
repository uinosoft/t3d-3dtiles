import { Color3 } from 't3d';
import { KHR_materials_pbrSpecularGlossiness } from 't3d/addons/loaders/glTF/extensions/KHR_materials_pbrSpecularGlossiness.js';
import { InstancedPBRMaterial } from 't3d/addons/materials/InstancedMaterial.js';

export class KHR_materials_pbrSpecularGlossiness_i extends KHR_materials_pbrSpecularGlossiness {

	static getMaterial() {
		const material = new InstancedPBRMaterial(); // fallback to InstancedPBRMaterial
		material.specular = new Color3(0x111111);
		return material;
	}

}