import { Color3 } from 't3d';
import { KHR_materials_pbrSpecularGlossiness } from '../../libs/glTF/extensions/KHR_materials_pbrSpecularGlossiness.js';
import { InstancedPBRMaterial } from '../../materials/InstancedPBRMaterial.js';

export class KHR_materials_pbrSpecularGlossiness_i extends KHR_materials_pbrSpecularGlossiness {

	static getMaterial() {
		const material = new InstancedPBRMaterial();
		material.specular = new Color3(0x111111);
		return material;
	}

}