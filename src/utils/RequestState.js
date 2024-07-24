/**
 * State of the request.
 *
 * @enum {Number}
 */
const RequestState = {
	UNLOADED: 0,

	LOADING: 1,

	PARSING: 2,

	LOADED: 3,

	FAILED: 4
};
export default Object.freeze(RequestState);
