
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// For checking memory leaks
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif

// main/pbrt.cpp*
#include "stdafx.h"
#include "api.h"
#include "probes.h"
#include "parser.h"
#include "parallel.h"
#include "../LBF/Globals.h"
#include "../core/Timer.h"

extern char sceneName[BUFFER_SIZE]; 
extern int pbrtSamplesPerPixel; 
extern Timer timer;

// main program
int main(int argc, char *argv[]) {

	// Record total time
    timer.Start();

    Options options;
    vector<string> filenames;
	pbrtSamplesPerPixel = 0;

    // Process command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--ncores")) options.nCores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--outfile")) options.imageFile = argv[++i];
        else if (!strcmp(argv[i], "--quick")) options.quickRender = true;
        else if (!strcmp(argv[i], "--quiet")) options.quiet = true;
        else if (!strcmp(argv[i], "--verbose")) options.verbose = true;
		else if (!strcmp(argv[i], "--spp")) pbrtSamplesPerPixel = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            printf("usage: pbrt [--ncores n] [--outfile filename] [--quick] [--quiet] "
                   "[--verbose] [--help] <filename.pbrt> ...\n");
            return 0;
        }
        else {
			filenames.push_back(argv[i]);
		}
    }

    // Print welcome banner
    if (!options.quiet) {
        printf("pbrt version %s of %s at %s [Detected %d core(s)]\n",
               PBRT_VERSION, __DATE__, __TIME__, NumSystemCores());
        printf("Copyright (c)1998-2012 Matt Pharr and Greg Humphreys.\n");
        printf("CS613 - Machine Learning - Rafael Campos\n");
        printf("Submitted as part of the Learning Based Filtering Project\n");
        fflush(stdout);
    }
    pbrtInit(options);
    // Process scene description
    PBRT_STARTED_PARSING();
    if (filenames.size() == 0) {
        // Parse scene from standard input
        ParseFile("-");
    } else {

        // Parse scene from input files
        for (u_int i = 0; i < filenames.size(); i++) {
			if(i == 0) {
			
				size_t lastslash = filenames[i].find_last_of("\\"); 
				if (lastslash != std::string::npos){ 
					std::string tempFile = filenames[i].substr(0, lastslash + 1);
					strcpy(sceneName, tempFile.c_str()); 	
				}
				strcpy(sceneName, filenames[i].c_str());

			}
            if (!ParseFile(filenames[i])) {
                Error("Couldn't open scene file \"%s\"", filenames[i].c_str());
			}
		}
    }
    pbrtCleanup();

	timer.Stop();
	float runtime = timer.Time();
	
	char tmpFileName[1000];
	_splitpath(filenames[0].c_str(), NULL, NULL, tmpFileName, NULL);
	sprintf(tmpFileName, "%s_timing.txt", sceneName);

	FILE* fp;
	fopen_s(&fp, tmpFileName, "at");

	if(!fp)
	{
		fprintf(stderr, "Could not open file %s\n", tmpFileName);
		getchar();
		exit(-1);
	}

	fprintf(fp, "Total Time: %f sec", runtime);
	fclose(fp);

	printf("Total Time: %f\n", runtime);
   
	// Dump memory leaks
	_CrtDumpMemoryLeaks();

	return 0;
}


