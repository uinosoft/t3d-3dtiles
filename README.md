# t3d-3dtiles

[![NPM Package][npm]][npm-url]

A [3D Tiles](https://www.ogc.org/standard/3dtiles/) extension for [t3d.js](https://github.com/uinosoft/t3d.js).

[Examples](https://uinosoft.github.io/t3d-3dtiles/examples/)

## Licensing

This project is licensed under the BSD 3-Clause License. See the [LICENSE](./LICENSE) file for details.

### Dependencies

- [NASA-AMMOS/3DTilesRendererJS](https://github.com/NASA-AMMOS/3DTilesRendererJS): Licensed under the Apache License 2.0. Copyright © 2020 California Institute of Technology. See [LICENSE](./src/core/LICENSE) for details. The `core` module is directly referenced from `NASA-AMMOS/3DTilesRendererJS`, while other parts have been modified to work with the `t3d.js` rendering engine while maintaining the original functionality.

## Quick Start

Most of the interfaces and usage are similar to [NASA-AMMOS/3DTilesRendererJS](https://github.com/NASA-AMMOS/3DTilesRendererJS/blob/master/README.md). 

Here is the basic usage of t3d-3dtiles:

````javascript
import { TilesRenderer } from 't3d-3dtiles';

// Create tiles with tileset URI
const tiles = new TilesRenderer('./path/to/tileset.json');

// Add tiles.group to scene
scene.add(tiles.group);

// Add camera to tiles, you can add multiple cameras
tiles.addCamera(camera);

// Set screen size for tiles
tiles.setResolution(width, height);

function loop(count) {
    requestAnimationFrame(loop);

    ...
    	
    tiles.update(); // Update tiles every frame

    forwardRenderer.render(scene, camera);
}
requestAnimationFrame(loop);
````

[npm]: https://img.shields.io/npm/v/t3d-3dtiles
[npm-url]: https://www.npmjs.com/package/t3d-3dtiles