//
//  OFDMEngine.h
//  OFDM
//
//  Created by Stephen MacKinnon on 1/17/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//

#ifndef __OFDM__OFDMEngine__
#define __OFDM__OFDMEngine__

#include <iostream>
#include "ofdm_params.h"
#include <math.h>
#include <vector>
#include <complex>
#define _USE_MATH_DEFINES

using namespace std;

class OFDMEngine {
public:
    OFDMEngine();
    vector<double> Modulate( unsigned char* data, int iDataLength );
    void Demodulate( std::vector<double> *data );
    double FrameDetect( std::vector<double>* data );
    void FFTTest();
    vector<double> GenerateHeader();

//private:
    void normalize( double &val );
    vector<double> filter( vector<double> &b, double a, vector<double> &x );
};

#endif /* defined(__OFDM__OFDMEngine__) */
