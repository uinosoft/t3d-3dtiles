import { PBRMaterial } from 't3d';

/**
 * KHR_techniques_webgl extension
 * https://github.com/KhronosGroup/glTF/blob/main/extensions/2.0/Archived/KHR_techniques_webgl/README.md
 * This extension has been archived, so we only provide a basic implementation.
 */
export class KHR_techniques_webgl {

	static getMaterial() {
		return new PBRMaterial();
	}

	static parseParams(material, extension, textures) {
		const { values } = extension;
		const { u_diffuse } = values;

		if (u_diffuse) {
			material.diffuseMap = textures[u_diffuse.index];
			material.diffuseMapCoord = u_diffuse.texCoord || 0;
		}
	}

}