import { Mesh, PlaneGeometry, ShaderMaterial } from 't3d';

export class PivotPointMesh extends Mesh {

	constructor() {
		super(new PlaneGeometry(0, 0), new PivotMaterial());
		this.renderOrder = Infinity;
	}

}

class PivotMaterial extends ShaderMaterial {

	constructor() {
		super(pivotShader);
		this.depthWrite = false;
		this.depthTest = false;
		this.transparent = true;
	}

}

const pivotShader = {
	name: 'PivotPoint',
	uniforms: {
		resolution: [512, 512],
		size: 15,
		thickness: 2,
		opacity: 1
	},
	vertexShader: /* glsl */`
		attribute vec3 a_Position; 
		attribute vec2 a_Uv;

		uniform mat4 u_ProjectionView;
		uniform mat4 u_Model;

		uniform float pixelRatio;
		uniform float size;
		uniform float thickness;
		uniform vec2 resolution;

		varying vec2 v_Uv;

		void main() {
			v_Uv = a_Uv;

			float aspect = resolution.x / resolution.y;
			vec2 offset = a_Uv * 2.0 - vec2(1.0);
			offset.y *= aspect;

			vec4 screenPoint = u_ProjectionView * u_Model * vec4(a_Position, 1.0);
			screenPoint.xy += offset * (size + thickness) * screenPoint.w / resolution.x;

			gl_Position = screenPoint;
		}
	`,
	fragmentShader: /* glsl */`
		uniform float size;
		uniform float thickness;
		uniform float opacity;

		varying vec2 v_Uv;

		void main() {
			float ht = 0.5 * thickness;
			float planeDim = size + thickness;
			float offset = (planeDim - ht - 2.0) / planeDim;
			float texelThickness = ht / planeDim;

			vec2 vec = v_Uv * 2.0 - vec2(1.0);
			float dist = abs(length(vec) - offset);
			float fw = fwidth(dist) * 0.5;
			float a = smoothstep(texelThickness - fw, texelThickness + fw, dist);

			gl_FragColor = vec4(1, 1, 1, opacity * (1.0 - a));
		}
	`
};