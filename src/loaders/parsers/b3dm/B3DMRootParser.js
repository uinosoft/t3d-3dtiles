export class B3DMRootParser {

	static parse(context, loader) {
		const { root, featureTable, options } = context;

		// fix rtc center

		const rtcCenter = featureTable.getData('RTC_CENTER');
		if (rtcCenter) {
			root.position.x += rtcCenter[0];
			root.position.y += rtcCenter[1];
			root.position.z += rtcCenter[2];
		}

		if (options.adjustmentTransform) {
			root.matrix.transform(root.position, root.scale, root.quaternion);
			root.matrix.multiply(options.adjustmentTransform);
			root.matrix.decompose(root.position, root.quaternion, root.scale);
		}
	}

}