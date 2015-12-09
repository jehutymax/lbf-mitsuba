//
// Created by Rafael Campos on 11/30/15.
//

#ifndef MITSUBA_LBF_SAMPLER_CONFIG_H
#define MITSUBA_LBF_SAMPLER_CONFIG_H

namespace lbf {

#include <stdio.h>

typedef float SampleElem;

// To enable LBF, set the following to 1.
// 0 eliminates LBF-related code from Mitsuba.
#define ENABLE_LBF 1

//#define SAMPLER_API __declspec(dllexport)

// Avoids samples on pixel boundaries.
#define SAMPLE_EPSILON 0.01

    enum EOffset {
        EXcoord = 0,
        EYcoord = 1,
        EColor1 = 0,
        EColor2 = 1,
        EColor3 = 2,
        EWorldX = 0,
        EWorldY = 1,
        EWorldZ = 2,
        ENormX = 3,
        ENormY = 4,
        ENormZ = 5,
        ETexture1X = 6,
        ETexture1Y = 7,
        ETexture1Z = 8,
        ETexture2X = 9,
        ETexture2Y = 10,
        ETexture2Z = 11,
        EVisibility = 12,
    };

}

#endif //MITSUBA_LBF_SAMPLER_CONFIG_H
