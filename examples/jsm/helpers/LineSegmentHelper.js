import { Object3D, SphereGeometry, BasicMaterial, Mesh, Geometry, Buffer, Attribute, DRAW_MODE, Vector3 } from 't3d';

export class LineSegmentHelper extends Object3D {

	constructor(lineSegment, pointRadius = 1, pointColor = 0x00ffff, lineColor = 0x00ffff) {
		super();

		// create start and end points

		const pointGeometry = new SphereGeometry(1);
		const pointMaterial = new BasicMaterial();

		this.start = new Mesh(pointGeometry, pointMaterial);
		this.start.scale.setScalar(pointRadius);
		this.add(this.start);

		this.end = new Mesh(pointGeometry, pointMaterial);
		this.end.scale.setScalar(pointRadius);
		this.add(this.end);

		// create line

		const lineArray = new Float32Array([0, 0, -0.5, 0, 0, 0.5]);
		const lineGeometry = new Geometry();
		lineGeometry.addAttribute('a_Position', new Attribute(new Buffer(lineArray, 3)));
		lineGeometry.computeBoundingBox();
		lineGeometry.computeBoundingSphere();

		const lineMaterial = new BasicMaterial();
		lineMaterial.drawMode = DRAW_MODE.LINES;

		this.line = new Mesh(lineGeometry, lineMaterial);
		this.line.raycast = () => {}; // disable raycasting
		this.add(this.line);

		//

		this.lineSegment = lineSegment;

		this.setPointColor(pointColor);
		this.setLineColor(lineColor);
	}

	setPointColor(color) {
		// both points share the same material
		// so we only need to update one of them
		this.start.material.diffuse.setHex(color);
	}

	setLineColor(color) {
		this.line.material.diffuse.setHex(color);
	}

	updateMatrix(force) {
		const { lineSegment } = this;
		const { startPoint, endPoint } = lineSegment;

		// TODO set center to group position

		this.start.position.copy(startPoint);
		this.end.position.copy(endPoint);

		lineSegment.getCenter(this.line.position);
		this.line.scale.z = lineSegment.getLength();
		_vec3_1.subVectors(endPoint, startPoint).normalize();
		this.line.quaternion.setFromUnitVectors(_zAxis, _vec3_1);

		super.updateMatrix(force);
	}

}

const _vec3_1 = new Vector3();
const _zAxis = new Vector3(0, 0, 1);