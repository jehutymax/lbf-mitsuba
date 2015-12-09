/* 
 * File:   multisampler.h
 * Author: Fabrice Rousselle
 *
 * Created on 31. mars 2013, 18:48
 */
#include <stdafx.h>
#include "multisampler.h"

#include <vector>
#include <limits>
#include <algorithm>

#include "montecarlo.h"
#include "camera.h"

#include "stratified.h"
#include "lowdiscrepancy.h"

extern int pbrtSamplesPerPixel;

MultiSampler::MultiSampler(int xstart, int xend, int ystart, int yend,
    int spp, float sopen, float sclose, float threshold, int nIterations,
    int sppInit, const Film *film, bool use_ld_samples)
    : Sampler(xstart, xend, ystart, yend, spp, sopen, sclose),
      use_ld_samples(use_ld_samples),
      xPixelCount(film->xResolution),
      yPixelCount(film->yResolution),
      _nIterations(nIterations) {

    sppInitReq = sppInit;
    samplesBuf = NULL;
	numOfBuffers = 2;
	assert(numOfBuffers == 2);
	assert(xPixelCount > 0);
	assert(yPixelCount > 0);
    scrambling.resize(numOfBuffers);
    initBase(NULL);
    
}


// Constructor used to create sub-samplers during the initialization phase. The
// sub-samplers are the one that do the actual job, the 'main' sampler is only
// used to create these.
MultiSampler::MultiSampler(const MultiSampler *parent, int xstart,
    int xend, int ystart, int yend, int pass)
    : Sampler(parent->xPixelStart, parent->xPixelEnd, parent->yPixelStart,
      parent->yPixelEnd, parent->samplesPerPixel, parent->shutterOpen,
      parent->shutterClose),
      use_ld_samples(parent->use_ld_samples),
      xPixelCount(parent->xPixelCount),
      yPixelCount(parent->yPixelCount),
      _nIterations(parent->_nIterations),
      film(parent->film) {
    _xPos = xstart;
    _yPos = ystart;
    _xStartSub = xstart; _xEndSub = xend;
    _yStartSub = ystart; _yEndSub = yend;
    sppInitReq = parent->sppInitReq;
    samplesBuf = NULL;
	numOfBuffers = 2;
    scrambling.resize(numOfBuffers);
    initBase(parent, pass);
}


// Constructor for sub-samplers during the adaptive phase. The sub-samplers do
// the actual jobs, while the "main" sampler only distributes the workload.
// This version of the sampler uses a sampling map to drive the sampling.
MultiSampler::MultiSampler(const MultiSampler *parent, int xstart,
    int xend, int ystart, int yend,
    MultiScramblingInfo *scrambling)
    : Sampler(parent->xPixelStart, parent->xPixelEnd, parent->yPixelStart,
      parent->yPixelEnd, parent->samplesPerPixel, parent->shutterOpen,
      parent->shutterClose),
      use_ld_samples(parent->use_ld_samples),
      xPixelCount(parent->xPixelCount),
      yPixelCount(parent->yPixelCount),
      _nIterations(parent->_nIterations),
      film(parent->film) {
    // Update the sampler state
    adaptive = true;
    samplerInit = NULL;
    samplesBuf = NULL;
	numOfBuffers = 2;
    this->scrambling.resize(numOfBuffers);
    this->scrambling[0] = scrambling;
    
    // Each sub-sampler is tasked to sample a specific tile of the whole image.
    _xPos = xstart;
    _yPos = ystart;
    _xStartSub = xstart; _xEndSub = xend;
    _yStartSub = ystart; _yEndSub = yend;
    
    _isMainSampler = false;
    
    if (parent == NULL) {
        Severe("oups");
    }
}


MultiSampler::~MultiSampler() {
    if (samplerInit != NULL) delete samplerInit;
    if (_isMainSampler) {
        for (size_t i = 0; i < scrambling.size(); i++) {
            delete [] scrambling[i];
        }
    }
    if (samplesBuf != NULL)
        delete [] samplesBuf;
}


void MultiSampler::initBase(const MultiSampler * parent, int pass) {
    // The MultiSampler distributes the samples over two passes.
    sppInit = samplesPerPixel / scrambling.size();
    if (_nIterations > 0) {
        sppInit = min(Floor2Int(sppInitReq / scrambling.size()), sppInit);
    }
    
    adaptive = false;

    // Construct the sampler for the initialization phase
    if (use_ld_samples) {
        samplerInit = NULL;
    }
    else {
        samplerInit = new RandomSampler(_xStartSub, _xEndSub, _yStartSub, _yEndSub,
            sppInit, shutterOpen, shutterClose);
    }

    // Compute the total number of pixels to be generated
    int nPix = xPixelCount * yPixelCount;
    int nSamplesInit = (scrambling.size() * sppInit) * nPix;
    int nSamplesAdapt = samplesPerPixel * nPix - nSamplesInit;
    pixelsToSampleTotal = Ceil2Int(float(nSamplesAdapt) / samplesPerPixel);

    if (parent != NULL) {
        _xPos = _xStartSub;
        _yPos = _yStartSub;
        _isMainSampler = false;
        scrambling[0] = parent->scrambling[pass];
    }
    else {
        _xPos = xPixelStart;
        _yPos = yPixelStart;
        _isMainSampler = true;
        int nPixInit = (xPixelEnd-xPixelStart) * (yPixelEnd-yPixelStart);
        for (size_t i = 0; i < scrambling.size(); i++) {
            scrambling[i] = new MultiScramblingInfo[nPixInit];
        }
      
    }
}


Sampler *MultiSampler::GetSubSampler(int num, int count) {
    // The MultiSampler performs the sampling over two passes, one for each
    // destination buffer. The first half of the samplers will correspond to the
    // first pass, while the second half corresponds to the second pass. To
    // enable this, we simply modify the given 'num' and 'count' to cycle twice.

    count = count / scrambling.size();
    int pass = num / count;
    num = num % count;
    if (!adaptive) {
        int x0, x1, y0, y1;
        ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
        if (x0 == x1 || y0 == y1) return NULL;
        return new MultiSampler(this, x0, x1, y0, y1, pass);
    }
    else {
        // Compute this job's tile
        int x0, x1, y0, y1;
        ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
        
        // Ensure we don't go outside the sampling map bounds
        x0 = max(x0, 0); x1 = min(x1, xPixelCount);
        y0 = max(y0, 0); y1 = min(y1, yPixelCount);
        
        // Use the appropriate sampling map
        return new MultiSampler(this, x0, x1, y0, y1, scrambling[pass]);
    }
}

int MultiSampler::GetMoreSamplesMapLD(Sample *samples, RNG &rng) {

    // Nothing to do for degenerate patch
    if (_xStartSub == _xEndSub || _yStartSub == _yEndSub)
        return 0;
    
    // During the initialization phase, we drawn samples one by one
    if (!adaptive) {
        // Move to the next line of the tile if the current is done
        if (_xPos == _xEndSub) {
            _xPos = _xStartSub;
            _yPos++;
        }
        
        // Stop if we processed all lines of the tile
        if (_yPos == _yEndSub)
            return 0;
        
        // Allocate the "samples buffer" needed for low-discrepancy sampling
        if (samplesBuf == NULL)
            samplesBuf = new float[LDPixelSampleFloatsNeeded(samples, sppInit)];
        
        // Draw the samples
        MyLDPixelSample(_xPos, _yPos, shutterOpen, shutterClose, sppInit, samples, rng, scrambling[0]);
        _xPos++;
        
        return sppInit;
    }
    
    // Allocate the "samples buffer" needed for low-discrepancy sampling
    if (samplesBuf == NULL)
        samplesBuf = new float[LDPixelSampleFloatsNeeded(samples, samplesPerPixel)];
    
    // Go over the tile until we get some samples
    for ( ; _yPos < _yEndSub; _yPos++) {
        for ( ; _xPos < _xEndSub; _xPos++) {
            // Get requested sample count for current pixel
            int pix = _xPos + _yPos * xPixelCount;
    
            // Update the rounding error
            int nSamples = nPixelSamples;//min(samplesPerPixel, nSamples);
            
            // Skip this pixel if there's no sample to draw
            if (nSamples <= 0)
                continue;
            
            // Draw the samples
            MyLDPixelSample(_xPos, _yPos, shutterOpen, shutterClose, nSamples, samples, rng, scrambling[0]);

            // Move to next pixel and return sample count
            _xPos++;
            return nSamples;
        }
        _xPos = _xStartSub;
    }
    
    return 0;
}

Sampler *CreateMultiSampler(const ParamSet &params,
                                const Film *film, const Camera *camera) {
    int spp = -1;
	if(pbrtSamplesPerPixel != 0) {
		spp = pbrtSamplesPerPixel;
	} else {
		spp = params.FindOneInt("pixelsamples", 8);
		pbrtSamplesPerPixel = spp;
	}
	assert(spp > 0);
    float th = params.FindOneFloat("threshold", std::numeric_limits<float>::infinity());
    int nIterations = 1;
    int sppInit = spp;

    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    
    // Low-discrepancy
    bool use_ld_samples = true;
	assert(use_ld_samples);
	assert(nIterations == 1);
    return new MultiSampler(xstart, xend, ystart, yend, spp, camera->shutterOpen,
        camera->shutterClose, th, nIterations, sppInit, film, use_ld_samples);
}


void MultiSampler::MyLDPixelSample(int xPos, int yPos, float shutterOpen,
    float shutterClose, int nPixelSamples, Sample *samples, RNG &rng,
    MultiScramblingInfo *scramblingArray) {
	
    // Prepare temporary array pointers for low-discrepancy camera samples
    float *buf = samplesBuf;
    float *imageSamples = buf; buf += 2 * nPixelSamples;
    float *lensSamples = buf;  buf += 2 * nPixelSamples;
    float *timeSamples = buf;  buf += nPixelSamples;

    // Prepare temporary array pointers for low-discrepancy integrator samples
    uint32_t count1D = samples[0].n1D.size();
    uint32_t count2D = samples[0].n2D.size();
    const uint32_t *n1D = count1D > 0 ? &samples[0].n1D[0] : NULL;
    const uint32_t *n2D = count2D > 0 ? &samples[0].n2D[0] : NULL;
    float **oneDSamples = ALLOCA(float *, count1D);
    float **twoDSamples = ALLOCA(float *, count2D);
    for (uint32_t i = 0; i < count1D; ++i) {
        oneDSamples[i] = buf;
        buf += n1D[i] * nPixelSamples;
    }
    for (uint32_t i = 0; i < count2D; ++i) {
        twoDSamples[i] = buf;
        buf += 2 * n2D[i] * nPixelSamples;
    }
    
    // Get reference to the current pixel's scrambling info
    int pix = (xPos-xPixelStart) + (yPos-yPixelStart) * (xPixelEnd-xPixelStart);
    MultiScramblingInfo &scrambling = scramblingArray[pix];
    
    // Define scrambling seeds on first call
    if (scrambling._nGenerated == 0) {
        int nSeeds = 5 + count1D + 2 * count2D;
        scrambling._seeds.resize(nSeeds);
        for (int i = 0; i < nSeeds; i++)
            scrambling._seeds[i] = rng.RandomUInt();
        scrambling._image = &scrambling._seeds[0];
        scrambling._lens = scrambling._image + 2;
        scrambling._time = scrambling._lens  + 2;
        scrambling._oneD = scrambling._time  + 1;
        scrambling._twoD = scrambling._oneD  + count1D;
    }
    
    // Generate low-discrepancy pixel samples
    MyLDShuffleScrambled2D(1, nPixelSamples, scrambling._nGenerated, imageSamples, rng, scrambling._image);
    MyLDShuffleScrambled2D(1, nPixelSamples, scrambling._nGenerated, lensSamples,  rng, scrambling._lens);
    MyLDShuffleScrambled1D(1, nPixelSamples, scrambling._nGenerated, timeSamples,  rng, *scrambling._time);
    for (uint32_t i = 0; i < count1D; ++i)
        MyLDShuffleScrambled1D(n1D[i], nPixelSamples, scrambling._nGenerated, oneDSamples[i], rng, scrambling._oneD[i]);
    for (uint32_t i = 0; i < count2D; ++i)
        MyLDShuffleScrambled2D(n2D[i], nPixelSamples, scrambling._nGenerated, twoDSamples[i], rng, &scrambling._twoD[2*i]);

    // Initialize _samples_ with computed sample values
    for (int i = 0; i < nPixelSamples; ++i) {
        samples[i].imageX = xPos + imageSamples[2*i];
        samples[i].imageY = yPos + imageSamples[2*i+1];
		if(((size_t) samples[i].imageX) == xPos + 1) {
			samples[i].imageX -= SAMPLE_EPSILON;
		} 
		if(((size_t) samples[i].imageY) == yPos + 1) {
			samples[i].imageY -= SAMPLE_EPSILON;
		}
        samples[i].time = Lerp(timeSamples[i], shutterOpen, shutterClose);
        samples[i].lensU = lensSamples[2*i];
        samples[i].lensV = lensSamples[2*i+1];
        // Copy integrator samples into _samples[i]_
        for (uint32_t j = 0; j < count1D; ++j) {
            int startSamp = n1D[j] * i;
            for (uint32_t k = 0; k < n1D[j]; ++k)
                samples[i].oneD[j][k] = oneDSamples[j][startSamp+k];
        }
        for (uint32_t j = 0; j < count2D; ++j) {
            int startSamp = 2 * n2D[j] * i;
            for (uint32_t k = 0; k < 2*n2D[j]; ++k)
                samples[i].twoD[j][k] = twoDSamples[j][startSamp+k];
        }
       
    }
    
    scrambling._nGenerated += nPixelSamples;
}

