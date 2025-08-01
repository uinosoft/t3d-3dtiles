import { PIXEL_TYPE, PIXEL_FORMAT } from 't3d';

export function safeTextureGetByteLength(tex) {
	let bytes = 0;

	if (tex.isTexture2D) {
		const { image, format, type, generateMipmaps } = tex;
		if (image) {
			bytes = getTextureByteLength(image.width, image.height, format, type);
			bytes *= generateMipmaps ? 1.3333 : 1;
		}
	} else if (tex.isTextureCube) {
		const { images, format, type, generateMipmaps } = tex;
		if (images) {
			images.forEach(image => {
				bytes += getTextureByteLength(image.width, image.height, format, type);
			});
			// TODO consider tex.mipmaps
			bytes *= generateMipmaps ? 1.3333 : 1;
		}
	}

	return bytes;
}

// Returns the estimated number of bytes used by the object
export function estimateBytesUsed(object) {
	const dedupeSet = new Set();

	let totalBytes = 0;
	object.traverse(c => {
		if (c.geometry) {
			for (const attrName in c.geometry.attributes) {
				const attr = c.geometry.getAttribute(attrName);
				const buffer = attr && attr.buffer;
				if (buffer && !dedupeSet.has(buffer)) {
					totalBytes += buffer.array.BYTES_PER_ELEMENT * buffer.count * buffer.stride;
					dedupeSet.add(buffer);
				}
			}
		}

		if (c.material) {
			const material = c.material;
			for (const key in material) {
				const value = material[key];
				if (value && value.isTexture && !dedupeSet.has(value)) {
					totalBytes += safeTextureGetByteLength(value);
					dedupeSet.add(value);
				}
			}
		}
	});


	return totalBytes;
}

function getTextureByteLength(width, height, format, type) {
	const typeInfo = getTextureTypeByteLength(type);
	switch (format) {
		case PIXEL_FORMAT.ALPHA:
			return width * height;
		case PIXEL_FORMAT.RED:
			return ((width * height) / typeInfo.components) * typeInfo.byteLength;
		case PIXEL_FORMAT.RED_INTEGER:
			return ((width * height) / typeInfo.components) * typeInfo.byteLength;
		case PIXEL_FORMAT.RG:
			return ((width * height * 2) / typeInfo.components) * typeInfo.byteLength;
		case PIXEL_FORMAT.RG_INTEGER:
			return ((width * height * 2) / typeInfo.components) * typeInfo.byteLength;
		case PIXEL_FORMAT.RGB:
			return ((width * height * 3) / typeInfo.components) * typeInfo.byteLength;
		case PIXEL_FORMAT.RGBA:
			return ((width * height * 4) / typeInfo.components) * typeInfo.byteLength;
		case PIXEL_FORMAT.RGBA_INTEGER:
			return ((width * height * 4) / typeInfo.components) * typeInfo.byteLength;
		// S3TC/DXT
		case PIXEL_FORMAT.RGB_S3TC_DXT1:
		case PIXEL_FORMAT.RGBA_S3TC_DXT1:
			return Math.floor((width + 3) / 4) * Math.floor((height + 3) / 4) * 8;
		case PIXEL_FORMAT.RGBA_S3TC_DXT3:
		case PIXEL_FORMAT.RGBA_S3TC_DXT5:
			return Math.floor((width + 3) / 4) * Math.floor((height + 3) / 4) * 16;
		// PVRTC
		case PIXEL_FORMAT.RGB_PVRTC_2BPPV1:
		case PIXEL_FORMAT.RGBA_PVRTC_2BPPV1:
			return (Math.max(width, 16) * Math.max(height, 8)) / 4;
		case PIXEL_FORMAT.RGB_PVRTC_4BPPV1:
		case PIXEL_FORMAT.RGBA_PVRTC_4BPPV1:
			return (Math.max(width, 8) * Math.max(height, 8)) / 2;
		// ETC
		case PIXEL_FORMAT.RGB_ETC1:
		case PIXEL_FORMAT.RGB_ETC2:
			return Math.floor((width + 3) / 4) * Math.floor((height + 3) / 4) * 8;
		case PIXEL_FORMAT.RGBA_ETC2_EAC:
			return Math.floor((width + 3) / 4) * Math.floor((height + 3) / 4) * 16;
		// ASTC
		case PIXEL_FORMAT.RGBA_ASTC_4x4:
			return Math.floor((width + 3) / 4) * Math.floor((height + 3) / 4) * 16;
		// BPTC
		case PIXEL_FORMAT.RGBA_BPTC:
			return Math.ceil(width / 4) * Math.ceil(height / 4) * 16;
	}

	console.warn(`Unable to determine texture byte length for ${format} format.`);

	// fallback: RGBA8
	return width * height * 4;
}

function getTextureTypeByteLength(type) {
	switch (type) {
		case PIXEL_TYPE.UNSIGNED_BYTE:
		case PIXEL_TYPE.BYTE:
			return { byteLength: 1, components: 1 };
		case PIXEL_TYPE.UNSIGNED_SHORT:
		case PIXEL_TYPE.SHORT:
		case PIXEL_TYPE.HALF_FLOAT:
			return { byteLength: 2, components: 1 };
		case PIXEL_TYPE.UNSIGNED_SHORT_4_4_4_4:
		case PIXEL_TYPE.UNSIGNED_SHORT_5_5_5_1:
			return { byteLength: 2, components: 4 };
		case PIXEL_TYPE.UNSIGNED_INT:
		case PIXEL_TYPE.INT:
		case PIXEL_TYPE.FLOAT:
			return { byteLength: 4, components: 1 };
		case PIXEL_TYPE.UNSIGNED_INT_24_8:
			return { byteLength: 4, components: 1 };
		default:
			return { byteLength: 1, components: 1 };
	}
}
