import { Sphere, Vector3 } from 't3d';
import { OBJECT_FRAME } from '../math/Ellipsoid.js';

export class ReorientationPlugin {

	constructor(options) {
		options = {
			up: '+z',
			recenter: true,

			lat: null,
			lon: null,
			height: 0,
			...options
		};

		this.tiles = null;

		this.up = options.up.toLowerCase().replace(/\s+/, '');
		this.lat = options.lat;
		this.lon = options.lon;
		this.height = options.height;
		this.recenter = options.recenter;
		this._callback = null;
	}

	init(tiles) {
		this.tiles = tiles;

		this._callback = () => {
			const { up, lat, lon, height, recenter } = this;

			if (lat !== null && lon !== null) {
				// if the latitude and longitude are provided then remove the position offset
				this.transformLatLonHeightToOrigin(lat, lon, height);
			} else {
				const { ellipsoid } = tiles;
				const minRadii = Math.min(ellipsoid.radius.x, ellipsoid.radius.y, ellipsoid.radius.z);
				tiles.getBoundingSphere(sphere);
				if (sphere.center.getLength() > minRadii * 0.5) {
					// otherwise see if this is possibly a tile set on the surface of the globe based on the positioning
					const cart = {};
					ellipsoid.getPositionToCartographic(sphere.center, cart);
					this.transformLatLonHeightToOrigin(cart.lat, cart.lon, cart.height);
					console.log(cart.lat, cart.lon, cart.height);
				} else {
					// lastly fall back to orienting the up direction to +Y
					tiles.euler.set(0, 0, 0);
					switch (up) {
						case 'x': case '+x':
							tiles.euler.z = Math.PI / 2;
							break;
						case '-x':
							tiles.euler.z = -Math.PI / 2;
							break;

						case 'y': case '+y':
							break;
						case '-y':
							tiles.euler.z = Math.PI;
							break;

						case 'z': case '+z':
							tiles.euler.x = -Math.PI / 2;
							break;
						case '-z':
							tiles.euler.x = Math.PI / 2;
							break;
					}

					tiles.position
						.copy(sphere.center)
						.applyEuler(tiles.euler)
						.multiplyScalar(-1);
				}
			}

			if (!recenter) {
				tiles.position.setScalar(0);
			}

			tiles.removeEventListener('load-tile-set', this._callback);
		};

		tiles.addEventListener('load-tile-set', this._callback);
	}

	transformLatLonHeightToOrigin(lat, lon, height = 0) {
		const tiles = this.tiles;
		const { ellipsoid } = tiles;

		// get ENU orientation (Z facing north and X facing west) and position
		ellipsoid.getRotationMatrixFromAzElRoll(lat, lon, 0, 0, 0, tiles.matrix, OBJECT_FRAME);
		ellipsoid.getCartographicToPosition(lat, lon, height, vec);

		// adjust the tiles matrix
		tiles.matrix
			.setPosition(vec)
			.invert()
			.decompose(tiles.position, tiles.quaternion, tiles.scale);
		tiles.updateMatrix();
	}

	dispose() {
		const tiles = this.tiles;

		tiles.position.setScalar(0);
		tiles.quaternion.identity();
		tiles.scale.set(1, 1, 1);

		tiles.addEventListener('load-tile-set', this._callback);
	}

}

const sphere = new Sphere();
const vec = new Vector3();