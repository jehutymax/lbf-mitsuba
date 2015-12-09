//
// Created by Rafael Campos on 11/30/15.
//

#ifndef MITSUBA_LBF_FEATURESAMPLER_H
#define MITSUBA_LBF_FEATURESAMPLER_H

#include <stddef.h>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <string>
#include <LBF/config.h>
#include "lbf_config.h"

namespace lbf {
    class FeatureSampler {
    public:
        static void initialize(size_t width, size_t height, size_t samplesPerPixel);
        static void processData(char* sceneName);
        static bool shouldFeatureBeSaved(size_t x, size_t y, size_t k);
        static int getPixelCount(size_t x, size_t y);

        static size_t getWidth();
        static size_t getHeight();
        static size_t getSamplesPerPixel();
        static size_t getNumberOfSamples();

        static SampleElem getFeature(size_t x, size_t y, size_t k, EOffset offset);
        static SampleElem getFeature(size_t index, bool isTexture2);

        static void setWidth(size_t width);
        static void setHeight(size_t height);
        static void setSamplesPerPixel(size_t samplesPerPixel);

        static void setPosition(size_t x, size_t y, size_t k, SampleElem position, EOffset offset);
        static void setPosition(size_t index, SampleElem position, bool saveSample,
                                int n, size_t x, size_t y, EOffset offset);

        static void setColor(size_t x, size_t y, size_t k, SampleElem color, EOffset offset);
        static void setColor(size_t index, SampleElem color, bool saveSample,
                             int n, size_t x, size_t y, EOffset offset);

        static void setFeature(size_t x, size_t y, size_t k, SampleElem feature, EOffset offset);
        static void setFeature(size_t index, SampleElem feature, bool saveSample,
                               int n, size_t x, size_t y, EOffset offset);


    private:
        static size_t getIndex(size_t x, size_t y, size_t &k, bool &saveSample);
        static void processTexture2Data();
        static void generateBufferIndex(int &bufferIndex, float &bufferNormFactor, int n);
        static void saveRenderTime(char* sceneName, char* inputFolder, char* name);

        static int m_numOfBuffers;
        static int* m_featureInd;
        static int* m_pixelSampleCount;
        static size_t m_width;
        static size_t m_height;
        static size_t m_spp;
        static size_t m_numOfSamples;
        static SampleElem* m_pixelData;
        static SampleElem* m_texture2Data;
        static SampleElem* m_sampleMean;
        static SampleElem* m_varData;
        static SampleElem** m_pixelBuffer;
    };
}

#endif //MITSUBA_LBF_FEATURESAMPLER_H
