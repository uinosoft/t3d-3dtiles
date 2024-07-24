import { Mesh, Geometry, PointsMaterial, Attribute, Buffer, VERTEX_COLOR } from 't3d';

export class PNTSRootParser {

	static parse(context, loader) {
		const { featureTable } = context;

		const POINTS_LENGTH = featureTable.getData('POINTS_LENGTH');
		const POSITION = featureTable.getData('POSITION', POINTS_LENGTH, 'FLOAT', 'VEC3');
		const RGB = featureTable.getData('RGB', POINTS_LENGTH, 'UNSIGNED_BYTE', 'VEC3');
		const RGBA = featureTable.getData('RGBA', POINTS_LENGTH, 'UNSIGNED_BYTE', 'VEC4');

		// check unsupported features

		[
			// Global Properties
			'QUANTIZED_VOLUME_OFFSET',
			'QUANTIZED_VOLUME_SCALE',
			'CONSTANT_RGBA',
			'BATCH_LENGTH',

			// Per-point Properties
			'POSITION_QUANTIZED',
			'RGB565',
			'NORMAL',
			'NORMAL_OCT16P',
			'BATCH_ID'
		].forEach(feature => {
			if (feature in featureTable.header) {
				console.warn(`PNTSLoader: Unsupported FeatureTable feature "${feature}" detected.`);
			}
		});

		// generate root

		const geometry = new Geometry();
		geometry.addAttribute('a_Position', new Attribute(new Buffer(POSITION, 3), 3, 0, true));
		geometry.computeBoundingBox();
		geometry.computeBoundingSphere();

		const material = new PointsMaterial();
		material.size = 2;
		material.sizeAttenuation = false;

		if (RGB !== null) {
			geometry.addAttribute('a_Color', new Attribute(new Buffer(RGB, 3), 3, 0, true));
			material.vertexColors = VERTEX_COLOR.RGB;
		} else if (RGBA !== null) {
			geometry.addAttribute('a_Color', new Attribute(new Buffer(RGBA, 4), 4, 0, true));
			material.vertexColors = VERTEX_COLOR.RGBA;
		}

		const root = new Mesh(geometry, material);

		// output mesh to root

		context.root = root;

		// fix rtc center

		const rtcCenter = featureTable.getData('RTC_CENTER');
		if (rtcCenter) {
			root.position.x += rtcCenter[0];
			root.position.y += rtcCenter[1];
			root.position.z += rtcCenter[2];
		}
	}

}
