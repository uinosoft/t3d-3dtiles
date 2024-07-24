import { BasicMaterial, Mesh, Geometry, Buffer, Attribute, DRAW_MODE } from 't3d';

export class OBBHelper extends Mesh {

	constructor(obb, color = 0xffff00) {
		const lineBoxVertices = new Float32Array([
			-0.5, -0.5, -0.5, 0.5, -0.5, -0.5,
			0.5, -0.5, -0.5, 0.5, 0.5, -0.5,
			0.5, 0.5, -0.5, -0.5, 0.5, -0.5,
			-0.5, 0.5, -0.5, -0.5, -0.5, -0.5,
			-0.5, -0.5, 0.5, 0.5, -0.5, 0.5,
			0.5, -0.5, 0.5, 0.5, 0.5, 0.5,
			0.5, 0.5, 0.5, -0.5, 0.5, 0.5,
			-0.5, 0.5, 0.5, -0.5, -0.5, 0.5,
			-0.5, -0.5, -0.5, -0.5, -0.5, 0.5,
			0.5, -0.5, -0.5, 0.5, -0.5, 0.5,
			0.5, 0.5, -0.5, 0.5, 0.5, 0.5,
			-0.5, 0.5, -0.5, -0.5, 0.5, 0.5
		]);
		const lineBoxGeometry = new Geometry();
		lineBoxGeometry.addAttribute('a_Position', new Attribute(new Buffer(lineBoxVertices, 3)));
		lineBoxGeometry.computeBoundingBox();
		lineBoxGeometry.computeBoundingSphere();

		const lineMaterial = new BasicMaterial();
		lineMaterial.drawMode = DRAW_MODE.LINES;

		super(lineBoxGeometry, lineMaterial);

		this.obb = obb;

		this.setColor(color);
	}

	setColor(color) {
		this.material.diffuse.setHex(color);
	}

	updateMatrix(force) {
		const { box, rotation } = this.obb;
		box.getCenter(this.position);
		box.getSize(this.scale);
		this.quaternion.setFromRotationMatrix(this.matrix.setFromMatrix3(rotation));

		super.updateMatrix(force);
	}

}