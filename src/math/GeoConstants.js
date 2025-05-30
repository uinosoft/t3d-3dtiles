import { Vector3 } from 't3d';
import { WGS84_RADIUS, WGS84_HEIGHT } from '../constants.js';
import { Ellipsoid } from './Ellipsoid.js';

export const WGS84_ELLIPSOID = new Ellipsoid(new Vector3(WGS84_RADIUS, WGS84_RADIUS, WGS84_HEIGHT));
WGS84_ELLIPSOID.name = 'WGS84 Earth';