export const instancingParsVert = `
    #ifdef USE_INSTANCING
        attribute mat4 instanceMatrix;
    #endif
`;

export const instancingPositionVert = `
    #ifdef USE_INSTANCING
        transformed = (instanceMatrix * vec4(transformed, 1.0)).xyz;
    #endif
`;

export const instancingNormalVert = `
    #ifdef USE_INSTANCING
        #ifdef USE_INSTANCING
            objectNormal = (transposeMat4(inverseMat4(instanceMatrix)) * vec4(objectNormal, 0.0)).xyz;
        #endif

        #ifdef USE_TANGENT
            objectTangent = (transposeMat4(inverseMat4(instanceMatrix)) * vec4(objectTangent, 0.0)).xyz;
        #endif
    #endif
`;