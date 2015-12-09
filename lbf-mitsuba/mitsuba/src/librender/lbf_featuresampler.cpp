//
// Created by Rafael Campos on 11/30/15.
//

#include <mitsuba/render/lbf_featuresampler.h>

int lbf::FeatureSampler::m_numOfBuffers;
int* lbf::FeatureSampler::m_featureInd;
int* lbf::FeatureSampler::m_pixelSampleCount;
size_t lbf::FeatureSampler::m_width;
size_t lbf::FeatureSampler::m_height;
size_t lbf::FeatureSampler::m_spp;
size_t lbf::FeatureSampler::m_numOfSamples;
lbf::SampleElem* lbf::FeatureSampler::m_pixelData;
lbf::SampleElem* lbf::FeatureSampler::m_texture2Data;
lbf::SampleElem* lbf::FeatureSampler::m_sampleMean;
lbf::SampleElem* lbf::FeatureSampler::m_varData;
lbf::SampleElem** lbf::FeatureSampler::m_pixelBuffer;

void lbf::FeatureSampler::initialize(size_t width, size_t height, size_t samplesPerPixel)
{
    m_width = width;
    m_height = height;
    m_spp = samplesPerPixel;
    m_numOfSamples = width * height;

    // initialization
    m_pixelData = new lbf::SampleElem[SAMPLE_LENGTH * m_numOfSamples];
    m_sampleMean = new lbf::SampleElem[SAMPLE_LENGTH * m_numOfSamples];
    m_varData = new lbf::SampleElem[SAMPLE_LENGTH * m_numOfSamples];
    m_texture2Data = new lbf::SampleElem[NUM_OF_TEXTURE_2 * m_numOfSamples * m_spp];
    memset(m_pixelData, 0, SAMPLE_LENGTH * m_numOfSamples * sizeof(lbf::SampleElem));
    memset(m_sampleMean, 0, SAMPLE_LENGTH * m_numOfSamples * sizeof(lbf::SampleElem));
    memset(m_varData, 0, SAMPLE_LENGTH * m_numOfSamples * sizeof(lbf::SampleElem));
    memset(m_texture2Data, 0, NUM_OF_TEXTURE_2 * m_numOfSamples * m_spp * sizeof(lbf::SampleElem));

    m_numOfBuffers = 2;
    m_pixelBuffer = new lbf::SampleElem*[m_numOfBuffers];
    for (int i = 0; i < m_numOfBuffers; ++i) {
        m_pixelBuffer[i] = new lbf::SampleElem[SAMPLE_LENGTH * m_numOfSamples];
        memset(m_pixelBuffer[i], 0, SAMPLE_LENGTH * m_numOfSamples * sizeof(lbf::SampleElem));
    }

    // number of samples processed for each pixel
    m_pixelSampleCount = new int[m_numOfSamples];
    m_featureInd = new int[m_numOfSamples];
    memset(m_pixelSampleCount, 0, m_numOfSamples * sizeof(int));
    memset(m_featureInd, 0, m_numOfSamples * sizeof(int));
}


void lbf::FeatureSampler::generateBufferIndex(int &bufferIndex, float &bufferNormFactor, int n) {
    bufferIndex = ((n - 1) >= (m_spp / m_numOfBuffers));
    bufferNormFactor = ((n - 1) % (m_spp / m_numOfBuffers)) + 1;
    assert(bufferNormFactor > 0 && bufferNormFactor <= (m_spp / m_numOfBuffers));
    assert(bufferIndex >= 0 && bufferIndex < m_numOfBuffers);
}


void lbf::FeatureSampler::processTexture2Data() {
    assert(m_numOfBuffers >= 1);
    int imageSize = m_width * m_height;

    for (size_t i = 0; i < m_numOfSamples; ++i) {
        int pixelTex2Index = SAMPLE_LENGTH * i + FEATURE + ETexture2X;
        memset(&m_pixelData[pixelTex2Index], 0, NUM_OF_TEXTURE_2 * sizeof(lbf::SampleElem));
        memset(&m_sampleMean[pixelTex2Index], 0, NUM_OF_TEXTURE_2 * sizeof(lbf::SampleElem));
        memset(&m_varData[pixelTex2Index], 0, NUM_OF_TEXTURE_2 * sizeof(lbf::SampleElem));

        // Update mean and variance::
        for (int q = 0; q < m_spp; ++q) {
            // todo: check this statement for correctness
            int tex2index = NUM_OF_TEXTURE_2 * m_spp * i + (q * NUM_OF_TEXTURE_2);
            for (int k = 0; k < NUM_OF_TEXTURE_2; ++k) {
                lbf::SampleElem meanDelta = m_texture2Data[tex2index + k] - m_pixelData[pixelTex2Index + k];
                m_pixelData[pixelTex2Index + k] += meanDelta / (q + 1);

                int bufInd;
                float bufNormVector;
                generateBufferIndex(bufInd, bufNormVector, q + 1);
                int bufferTex2Ind = i + (k + TEXTURE_2_X) * imageSize;
                if(bufNormVector == 1.0f)
                {
                    m_pixelBuffer[bufInd][bufferTex2Ind] = 0;
                }
                meanDelta = m_texture2Data[tex2index + k] - m_pixelBuffer[bufInd][bufferTex2Ind];
                m_pixelBuffer[bufInd][bufferTex2Ind] += meanDelta / bufNormVector;

                lbf::SampleElem varDelta = m_texture2Data[tex2index + k] - m_sampleMean[pixelTex2Index + k];
                m_sampleMean[pixelTex2Index + k] += varDelta / (q + 1);
                m_varData[pixelTex2Index + k] += varDelta * (m_texture2Data[tex2index + k] - m_sampleMean[pixelTex2Index + k]);
            }
        }
    }
}


void lbf::FeatureSampler::saveRenderTime(char *sceneName, char *inputFolder, char *name) {
    std::string filename = std::string(sceneName);
    int lastSlash = filename.find_last_of("\\");
    int lastBackslash = filename.find_last_of("/");
    lastSlash = (lastSlash > lastBackslash) ? lastSlash : lastBackslash;
    if (lastSlash != std::string::npos) {
        std::string tempFile = filename.substr(0, lastSlash);
        std::string tempName = filename.substr(lastSlash + 1, std::string::npos);
        strcpy(inputFolder, tempFile.c_str());
        strcpy(name, tempName.c_str());
    }
    else {
        sprintf(inputFolder, ".");
        strcpy(name, filename.c_str());
    }

//    timer.stop();
    char tempFilename[1000];
    sprintf(tempFilename, "%s\\%s_timing.txt", inputFolder, name);
    FILE* fp = fopen(tempFilename, "wt");

    if(!fp)
    {
        fprintf(stderr, "Could not open file %s\n", tempFilename);
        getchar();
        exit(-1);
    }

//    fprintf(fp, "Render Time: %f sec\n", timer.getMilliseconds());
//    fclose(fp);
//    timer.reset();
//    timer.start();
}


void lbf::FeatureSampler::processData(char *filename) {
    // record time
    char inputFolder[BUFFER_SIZE];
    char sceneName[BUFFER_SIZE];
    saveRenderTime(filename, inputFolder, sceneName);

    // update texture2 values
    processTexture2Data();
    delete[] m_texture2Data;
    delete[] m_sampleMean;
    delete[] m_pixelSampleCount;
    delete[] m_featureInd;
//    renderTimer.stop();

    // run neural network and filter
    //todo: add the call to the neural network
}

bool lbf::FeatureSampler::shouldFeatureBeSaved(size_t x, size_t y,
                                               size_t k) {
    if(x < 0 || x >= m_width || y < 0 || y>= m_height)
        return false;
    return (m_featureInd[x + y * m_width] == k);
}


// the getters and the setters
lbf::SampleElem lbf::FeatureSampler::getFeature(size_t x, size_t y,
                                           size_t k, EOffset offset) {
    bool saveSample;
    size_t index = getIndex(x, y, k, saveSample);
    if (index != -1) {
        bool isTexture2 = (offset >= ETexture2X && offset <= ETexture2Z);
        size_t featureOffset;
        if (!isTexture2) {
            featureOffset = SAMPLE_LENGTH * index + FEATURE + offset;
        } else {
            featureOffset = NUM_OF_TEXTURE_2 * index * m_spp + k * NUM_OF_TEXTURE_2 + (offset - ETexture2X);
        }
        return getFeature(featureOffset, isTexture2);
    }

    return 0;
}


lbf::SampleElem lbf::FeatureSampler::getFeature(size_t index, bool isTexture2) {
    if(!isTexture2) {
        assert(index < SAMPLE_LENGTH * m_numOfSamples);
        return m_pixelData[index];
    } else {
        assert(index < NUM_OF_TEXTURE_2 * m_numOfSamples * m_spp);
        return m_texture2Data[index];
    }
}


size_t lbf::FeatureSampler::getWidth() {
    return m_width;
}


size_t lbf::FeatureSampler::getHeight() {
    return m_height;
}


size_t lbf::FeatureSampler::getSamplesPerPixel() {
    return m_spp;
}


size_t lbf::FeatureSampler::getNumberOfSamples() {
    assert(m_numOfSamples == m_width * m_height * m_spp);
    return m_numOfSamples;
}


size_t lbf::FeatureSampler::getIndex(size_t x, size_t y,
                                     size_t &k, bool &saveSample) {
    saveSample = false;

    if(x >= 0 && x < m_width && y >= 0 && y < m_height && m_height >= 0 && m_width >= 0 && k >= 0 && k < m_spp)
    {
        size_t index = x + y * m_width;
        m_pixelSampleCount[index] = MAX(k+1, m_pixelSampleCount[index]);
        assert(index < m_width * m_height);
        return index;
    }

    return -1;
}


int lbf::FeatureSampler::getPixelCount(size_t x, size_t y) {
    if(x < 0 || x >= m_width || y < 0 || y >= m_height) {
        return 0;
    }

    size_t index = x + y * m_width;
    assert(index >= 0 && index < m_numOfSamples);
    return m_pixelSampleCount[index];
}


void lbf::FeatureSampler::setPosition(size_t x, size_t y,
                                      size_t k, lbf::SampleElem position, EOffset offset) {
    bool saveSample;
    size_t index = getIndex(x, y, k, saveSample);

    if(index != -1) {
        setPosition(SAMPLE_LENGTH * index + POSITION + offset, position, saveSample, m_pixelSampleCount[index], x, y, offset);
    }
}


void lbf::FeatureSampler::setPosition(size_t index, lbf::SampleElem position, bool saveSample, int n,
                                      size_t x, size_t y,
                                      EOffset offset) {
    assert(index < SAMPLE_LENGTH * m_numOfSamples);
    assert(n != 0);

    // Update mean and variance
    lbf::SampleElem meanDelta = position - m_pixelData[index];
    m_pixelData[index] += meanDelta / n;

    int bufInd;
    float bufNormFactor;
    generateBufferIndex(bufInd, bufNormFactor, n);
    size_t bufPosInd = x + y * m_width + (POSITION + offset) * m_width * m_height;
    meanDelta = position - m_pixelBuffer[bufInd][bufPosInd];
    m_pixelBuffer[bufInd][bufPosInd] += meanDelta / bufNormFactor;

    lbf::SampleElem varDelta = position - m_sampleMean[index];
    m_sampleMean[index] += varDelta / n;
    m_varData[index] += varDelta * (position - m_sampleMean[index]);
}


void lbf::FeatureSampler::setColor(size_t x, size_t y,
                                   size_t k, lbf::SampleElem color, EOffset offset) {
    bool saveSample;
    size_t index = getIndex(x, y, k, saveSample);

    if(index != -1) {

        setColor(SAMPLE_LENGTH * index + COLOR + offset, color, saveSample, m_pixelSampleCount[index], x, y, offset);

        int colorIndex = m_spp * NUM_OF_COLORS * index + k*NUM_OF_COLORS + offset;
        assert(colorIndex < (NUM_OF_COLORS * m_width * m_height * m_spp));
    }

}

void lbf::FeatureSampler::setColor(size_t index, lbf::SampleElem color, bool saveSample, int n,
                                   size_t x, size_t y, EOffset offset) {

    assert(index < SAMPLE_LENGTH * m_numOfSamples);
    assert(n != 0);

    // Update mean and variance
    lbf::SampleElem meanDelta = color - m_pixelData[index];
    m_pixelData[index] += meanDelta / n;

    int bufInd;
    float bufNormFactor;
    generateBufferIndex(bufInd, bufNormFactor, n);
    size_t bufColorInd = x + y * m_width + (COLOR + offset) * m_width * m_height;
    meanDelta = color - m_pixelBuffer[bufInd][bufColorInd];
    m_pixelBuffer[bufInd][bufColorInd] += meanDelta / bufNormFactor;

    lbf::SampleElem varDelta = color - m_sampleMean[index];
    m_sampleMean[index] += varDelta / n;
    m_varData[index] += varDelta * (color - m_sampleMean[index]);
}

void lbf::FeatureSampler::setFeature(size_t x, size_t y,
                                     size_t k, lbf::SampleElem feature, EOffset offset) {
    bool saveSample;
    size_t index = getIndex(x, y, k, saveSample);

    if(index != -1) {

        if(offset == ETexture1X) {
            m_featureInd[index]++;
        }

        bool isTexture2 = false;
        if(offset >= ETexture2X && offset <= ETexture2Z) {
            isTexture2 = true;
            m_texture2Data[NUM_OF_TEXTURE_2 * m_spp * index + k * NUM_OF_TEXTURE_2 + (offset - ETexture2X)] = feature;
        }

        setFeature(SAMPLE_LENGTH * index + FEATURE + offset, feature, saveSample, m_pixelSampleCount[index], x, y, offset);

    }

}

void lbf::FeatureSampler::setFeature(size_t index, lbf::SampleElem feature, bool saveSample, int n,
                                     size_t x, size_t y,
                                     EOffset offset) {
    assert(index < SAMPLE_LENGTH * m_numOfSamples);
    assert(n != 0);

    // Update mean and variance
    lbf::SampleElem meanDelta = feature - m_pixelData[index];
    m_pixelData[index] += meanDelta / n;

    int bufInd;
    float bufNormFactor;
    generateBufferIndex(bufInd, bufNormFactor, n);
    size_t bufFeatInd = x + y * m_width + (FEATURE + offset) * m_width * m_height;
    meanDelta = feature - m_pixelBuffer[bufInd][bufFeatInd];
    m_pixelBuffer[bufInd][bufFeatInd] += meanDelta / bufNormFactor;

    lbf::SampleElem varDelta = feature - m_sampleMean[index];
    m_sampleMean[index] += varDelta / n;
    m_varData[index] += varDelta * (feature - m_sampleMean[index]);
}


void lbf::FeatureSampler::setWidth(size_t width) {
    m_width = width;
    m_numOfSamples = m_width * m_height;
}


void lbf::FeatureSampler::setHeight(size_t height) {
    m_height = height;
    m_numOfSamples = m_width * m_height;
}


void lbf::FeatureSampler::setSamplesPerPixel(size_t samplesPerPixel) {
    m_spp = samplesPerPixel;
}






