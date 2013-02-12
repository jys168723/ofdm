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
    OFDMEngine( int diffRef0, int diffRef1 ); // Temporary ctor for debugging
    
    vector<double> Modulate( unsigned char* data, int iDataLength );
    void Demodulate( std::vector<double> &symbRx, bool bLastFrame, int iUnpad );
    double FrameDetect( std::vector<double>* data );
    vector<double> GenerateHeader();
    
    vector<double> GetDataRx() { return m_dataRx; }
    vector<double> GetPhase() { return m_phase; }

private:
    vector<double> filter( vector<double> &b, double a, vector<double> &x );
    
    vector<double> m_dataRx,
                   m_phase;
    // Temp diff ref members
    int m_iDiffRef0, m_iDiffRef1;
};

#endif /* defined(__OFDM__OFDMEngine__) */
