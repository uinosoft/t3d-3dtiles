import { BasicMaterial, ShaderLib, MATERIAL_TYPE } from 't3d';
import { instancingParsVert, instancingPositionVert, instancingNormalVert } from './InstancedShaderChunks.js';

export class InstancedBasicMaterial extends BasicMaterial {

	constructor() {
		super();
		this.type = MATERIAL_TYPE.SHADER;
		this.shaderName = 'TILE_I_BASIC';
		this.vertexShader = vertexShader;
		this.fragmentShader = ShaderLib.basic_frag;
		this.defines.USE_INSTANCING = true;
	}

}

InstancedBasicMaterial.prototype.isInstancedBasicMaterial = true;

let vertexShader = ShaderLib.basic_vert;

vertexShader = vertexShader.replace('#include <logdepthbuf_pars_vert>', `
    #include <logdepthbuf_pars_vert>
    ${instancingParsVert}
`);
vertexShader = vertexShader.replace('#include <pvm_vert>', `
    ${instancingPositionVert}
    #include <pvm_vert>
`);
vertexShader = vertexShader.replace('#include <normal_vert>', `
    ${instancingNormalVert}
    #include <normal_vert>
`);