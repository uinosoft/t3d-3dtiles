const CesiumIon = {};

CesiumIon.fetchAssetJson = async function(assetId, accessToken) {
	const url = new URL(`https://api.cesium.com/v1/assets/${assetId}/endpoint`);
	url.searchParams.append('access_token', accessToken);

	const res = await fetch(url, { mode: 'cors' });

	if (!res.ok) throw `${res.status} : ${res.statusText}`;

	return await res.json();
};

CesiumIon.initTiles3D = function(tiles3D, assetJson) {
	const assetVersion = new URL(assetJson.url).searchParams.get('v');

	tiles3D.fetchOptions.mode = 'cors';
	tiles3D.fetchOptions.headers = {};
	tiles3D.fetchOptions.headers.Authorization = `Bearer ${assetJson.accessToken}`;

	tiles3D.preprocessURL = uri => {
		uri = new URL(uri);
		if (/^http/.test(uri.protocol)) {
			uri.searchParams.append('v', assetVersion);
		}
		return uri.toString();
	};
};

export default CesiumIon;