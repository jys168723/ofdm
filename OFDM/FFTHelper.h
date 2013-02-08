//
//  FFTHelper.h
//  OFDM
//
//  Created by Stephen MacKinnon on 2/7/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//

#ifndef OFDM_FFTHelper_h
#define OFDM_FFTHelper_h

#include <complex>                  // Must include complex before fftw3 so fftw3 treats
#include "api/fftw3.h"

using namespace std;

namespace FFTHelper {
    template <typename T>
    extern vector<complex<T>> fft_complex_1d( vector<complex<T>> &in, bool isIfft ) {
        int n= in.size();
        
        fftw_complex *fftOut= (fftw_complex*)fftw_malloc( n*sizeof(fftw_complex) );
        
        fftw_plan fftPlan;
        if ( isIfft )
            fftPlan= fftw_plan_dft_1d( n, (fftw_complex*)&in[0], fftOut, FFTW_BACKWARD, FFTW_ESTIMATE );
        else
            fftPlan= fftw_plan_dft_1d( n, (fftw_complex*)&in[0], fftOut, FFTW_FORWARD, FFTW_ESTIMATE );
        
        fftw_execute(fftPlan);
        
        vector<complex<T>> outVector( n, 0 );
        for( uint i=0; i<n; ++i )
            outVector[i]= complex<T>( fftOut[i][0]/(float)n, fftOut[i][1]/(float)n );
        
        return outVector;
    }
    
    template <typename T>
    extern vector<vector<complex<T>>> fft_complex_2d( vector<vector<complex<T>>> &in, bool isIfft ) {
        int nx= (int)in.size();
        int ny= (int)in[0].size();
        
        vector<vector<complex<T>>> out( nx, vector<complex<T>>(ny, 0) );
        
        for( uint i=0; i<nx; ++i ) {
            vector<complex<T>> fftResult= fft_complex_1d( in[i], isIfft );
            
            for( uint j=0; j<ny; ++j )
                out[i][j]= fftResult[j];
        }
        
        return out;
    }
}

#endif
