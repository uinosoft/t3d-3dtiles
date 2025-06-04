import { MATERIAL_TYPE, ShaderLib } from 't3d';

// Adjusts the provided material to support fading in and out using a bayer pattern.
export function wrapFadeMaterial(material) {
	material.shaderName = `${material.shaderName || material.type}_fade`;

	material.defines.FEATURE_FADE = 0;

	material.uniforms.fadeIn = 0;
	material.uniforms.fadeOut = 0;

	const fragmentShader = material.fragmentShader ||
		(material.type === MATERIAL_TYPE.BASIC ?
			ShaderLib.basic_frag : ShaderLib.pbr_frag);

	material.type = MATERIAL_TYPE.SHADER;

	material.vertexShader = material.type === MATERIAL_TYPE.BASIC ?
		ShaderLib.basic_vert : ShaderLib.pbr_vert;
	material.fragmentShader = fragmentShader
		.replace(/void main\(/, value => /* glsl */`
			#if FEATURE_FADE

			// adapted from https://www.shadertoy.com/view/Mlt3z8
			float bayerDither2x2(vec2 v) {
				return mod(3.0 * v.y + 2.0 * v.x, 4.0);
			}

			float bayerDither4x4(vec2 v) {
				vec2 P1 = mod(v, 2.0);
				vec2 P2 = floor(0.5 * mod(v, 4.0));
				return 4.0 * bayerDither2x2(P1) + bayerDither2x2(P2);
			}

			uniform float fadeIn;
			uniform float fadeOut;

			#endif

			${value}
		`)
		.replace(/#include <end_frag>/, value => /* glsl */`
			${value}

			#if FEATURE_FADE

			float bayerValue = bayerDither4x4(floor(mod(gl_FragCoord.xy, 4.0)));
			float bayerBins = 16.0;
			float dither = (0.5 + bayerValue) / bayerBins;
			if (dither >= fadeIn) {
				discard;
			}

			if (dither < fadeOut) {
				discard;
			}

			#endif
		`);

	return material.uniforms;
}