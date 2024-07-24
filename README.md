# t3d-3dtiles

[![NPM Package][npm]][npm-url]

An extension for [3D Tiles](https://www.ogc.org/standard/3dtiles/) based on [t3d.js](https://github.com/uinosoft/t3d.js), requires t3d.js version `v0.2.7` or later.
Inspired by [NASSA-AMMOS / 3DTilesRendererJS](https://github.com/NASA-AMMOS/3DTilesRendererJS).

[Examples](https://uinosoft.github.io/t3d-3dtiles/examples/)

## Usage

Here is the basic usage of t3d-3dtiles:

````javascript
// Create tiles3D with tileset URI
const tilesetURI = "./path/to/tileset.json";
const tiles3D = new Tiles3D(tilesetURI);

// Add tiles3D to scene
scene.add(tiles3D);

// Add camera to tiles3D, you can add multiple cameras
tiles3D.addCamera(camera);

// Set screen size for tiles3D
tiles3D.resize(width, height);

function loop(count) {
    requestAnimationFrame(loop);

    ...
    	
    tiles3D.update(); // Update tiles3D every frame

    forwardRenderer.render(scene, camera);
}
requestAnimationFrame(loop);
````

[npm]: https://img.shields.io/npm/v/t3d-3dtiles
[npm-url]: https://www.npmjs.com/package/t3d-3dtiles