import { Vector3, Matrix4, Quaternion, Attribute, Buffer } from 't3d';

export class I3DMRootParser {

	static parse(context, loader) {
		const { featureTable, root, options } = context;

		const INSTANCES_LENGTH = featureTable.getData('INSTANCES_LENGTH');
		const POSITION = featureTable.getData('POSITION', INSTANCES_LENGTH, 'FLOAT', 'VEC3');
		const NORMAL_UP = featureTable.getData('NORMAL_UP', INSTANCES_LENGTH, 'FLOAT', 'VEC3');
		const NORMAL_RIGHT = featureTable.getData('NORMAL_RIGHT', INSTANCES_LENGTH, 'FLOAT', 'VEC3');
		const SCALE = featureTable.getData('SCALE', INSTANCES_LENGTH, 'FLOAT', 'SCALAR');
		const SCALE_NON_UNIFORM = featureTable.getData('SCALE_NON_UNIFORM', INSTANCES_LENGTH, 'FLOAT', 'VEC3');

		// check unsupported features

		[
			// Global Properties
			'QUANTIZED_VOLUME_OFFSET',
			'QUANTIZED_VOLUME_SCALE',
			'EAST_NORTH_UP',

			// Per-instance Properties
			'POSITION_QUANTIZED',
			'NORMAL_UP_OCT32P',
			'NORMAL_RIGHT_OCT32P'
		].forEach(feature => {
			if (feature in featureTable.header) {
				console.warn(`I3DMLoader: Unsupported FeatureTable feature "${feature}" detected.`);
			}
		});

		// set instance matrix for all geometries

		const averageVector = new Vector3();
		for (let i = 0; i < INSTANCES_LENGTH; i++) {
			averageVector.x += POSITION[i * 3 + 0] / INSTANCES_LENGTH;
			averageVector.y += POSITION[i * 3 + 1] / INSTANCES_LENGTH;
			averageVector.z += POSITION[i * 3 + 2] / INSTANCES_LENGTH;
		}

		const instances = [];

		root.traverse(child => {
			if (child.isMesh) {
				const { geometry } = child;
				geometry.instanceCount = INSTANCES_LENGTH;

				const instanceMatrix = new Attribute(new Buffer(new Float32Array(INSTANCES_LENGTH * 16), 16), 16);
				instanceMatrix.divisor = 1;
				geometry.addAttribute('instanceMatrix', instanceMatrix);

				// Center the instance around an average point to avoid jitter at large scales.
				// Transform the average vector by matrix world so we can account for any existing
				// transforms of the instanced mesh.
				child.updateMatrix(true);
				child.position.copy(averageVector).applyMatrix4(child.worldMatrix);

				instances.push(child);
			}
		});

		for (let i = 0; i < INSTANCES_LENGTH; i++) {
			// position
			tempPos.fromArray(POSITION, i * 3).sub(averageVector);

			// rotation
			if (NORMAL_UP) {
				tempUp.fromArray(NORMAL_UP, i * 3);
				tempRight.fromArray(NORMAL_RIGHT, i * 3);
				tempFwd.crossVectors(tempRight, tempUp).normalize();

				tempMat.set(
					tempRight.x, tempUp.x, tempFwd.x, 0,
					tempRight.y, tempUp.y, tempFwd.y, 0,
					tempRight.z, tempUp.z, tempFwd.z, 0,
					0, 0, 0, 1
				);

				tempQuat.setFromRotationMatrix(tempMat);
			} else {
				tempQuat.set(0, 0, 0, 1);
			}

			// scale
			if (SCALE) {
				tempSca.set(SCALE[i], SCALE[i], SCALE[i]);
			} else if (SCALE_NON_UNIFORM) {
				tempSca.fromArray(SCALE_NON_UNIFORM, i * 3);
			} else {
				tempSca.set(1, 1, 1);
			}

			// TODO instance matrix should be applied to model root
			tempMat.transform(tempPos, tempSca, tempQuat).multiply(options.adjustmentTransform);

			for (let j = 0, l = instances.length; j < l; j++) {
				const { geometry } = instances[j];
				const instanceArray = geometry.getAttribute('instanceMatrix').buffer.array;
				tempMat.toArray(instanceArray, i * 16);
				geometry.version++;
			}
		}

		// fix rtc center

		const rtcCenter = featureTable.getData('RTC_CENTER');
		if (rtcCenter) {
			root.position.x += rtcCenter[0];
			root.position.y += rtcCenter[1];
			root.position.z += rtcCenter[2];
		}
	}

}

const tempFwd = new Vector3();
const tempUp = new Vector3();
const tempRight = new Vector3();
const tempPos = new Vector3();
const tempQuat = new Quaternion();
const tempSca = new Vector3();
const tempMat = new Matrix4();