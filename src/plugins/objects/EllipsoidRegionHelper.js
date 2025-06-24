import { Mesh, Geometry, BoxGeometry, Attribute, Buffer, LineMaterial, Vector3, MathUtils } from 't3d';
import { EdgesBuilder } from 't3d/addons/geometries/builders/EdgesBuilder.js';
import { EllipsoidRegion } from '../../math/EllipsoidRegion.js';

const _pos = new Vector3();

function getRegionGeometry(ellipsoidRegion) {
	// retrieve the relevant fields
	const {
		latRange,
		lonRange,
		heightRange
	} = ellipsoidRegion;

	const { x: latStart, y: latEnd } = latRange;
	const { x: lonStart, y: lonEnd } = lonRange;
	const { x: heightStart, y: heightEnd } = heightRange;

	// get the attributes
	const geometry = new BoxGeometry(1, 1, 1, 32, 32);
	const { a_Position: position } = geometry.attributes;

	// perturb the position buffer into an ellipsoid region
	for (let i = 0, l = position.buffer.count; i < l; i++) {
		_pos.fromArray(position.buffer.array, i * 3);

		const lat = MathUtils.mapLinear(_pos.x, -0.5, 0.5, latStart, latEnd);
		const lon = MathUtils.mapLinear(_pos.y, -0.5, 0.5, lonStart, lonEnd);

		let height = heightStart;
		if (_pos.z < 0) {
			height = heightEnd;
		}
		ellipsoidRegion.getCartographicToPosition(lat, lon, height, _pos);
		_pos.toArray(position.buffer.array, i * 3);
	}

	return geometry;
}

export class EllipsoidRegionHelper extends Mesh {

	constructor(ellipsoidRegion = new EllipsoidRegion(), color = 0xffff00) {
		super(new Geometry(), new LineMaterial());

		this.ellipsoidRegion = ellipsoidRegion;

		this.material.diffuse.setHex(color);

		this.update();

		this.raycast = () => {}; // disable raycasting
	}

	update() {
		this.geometry.dispose();

		const regionGeometry = getRegionGeometry(this.ellipsoidRegion);
		const { positions } = EdgesBuilder.getGeometryData(
			regionGeometry.attributes.a_Position.buffer.array,
			regionGeometry.index.buffer.array,
			{ thresholdAngle: 80 }
		);
		const geometry = new Geometry();
		geometry.addAttribute('a_Position', new Attribute(new Buffer(new Float32Array(positions), 3)));
		// geometry.computeBoundingBox();
		// geometry.computeBoundingSphere();

		this.geometry = geometry;
		this.geometry.computeBoundingBox();
		this.geometry.computeBoundingSphere();
	}

	dispose() {
		this.geometry.dispose();
		this.material.dispose();
	}

}