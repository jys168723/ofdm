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
    extern vector<complex<T>> fft_complex_1d( vector<complex<T>> &in, bool isIfft, bool bScale ) {
        int n= static_cast<int>( in.size() );
        
        fftw_complex *fftOut= (fftw_complex*)fftw_malloc( n*sizeof(fftw_complex) );
        
        fftw_plan fftPlan;
        if ( isIfft )
            fftPlan= fftw_plan_dft_1d( n, (fftw_complex*)&in[0], fftOut, FFTW_BACKWARD, FFTW_ESTIMATE );
        else
            fftPlan= fftw_plan_dft_1d( n, (fftw_complex*)&in[0], fftOut, FFTW_FORWARD, FFTW_ESTIMATE );
        
        fftw_execute(fftPlan);
        
        vector<complex<T>> outVector( n, 0 );
        for( uint i=0; i<n; ++i ) {
            if( bScale )
                outVector[i]= complex<T>( fftOut[i][0]/(float)n, fftOut[i][1]/(float)n );
            else
                outVector[i]= complex<T>( fftOut[i][0], fftOut[i][1] );
        }
            
        
        return outVector;
    }
    
    // Takes a column-wise 2D fft of in. Identical to Matlab fft
    template <typename T>
    extern vector<vector<complex<T>>> fft_complex_2d( vector<vector<complex<T>>> &in, bool isIfft, bool bScale ) {
        int nx= (int)in.size();
        int ny= (int)in[0].size();
        
        vector<vector<complex<T>>> out( nx, vector<complex<T>>(ny, 0) );
        
        for( uint iCol=0; iCol<ny; ++iCol ) {
            vector<complex<T>> fftInput( nx, 0 );
            for( uint iRow=0; iRow<nx; ++iRow )
                fftInput[iRow]= in[iRow][iCol];
            
            vector<complex<T>> fftResult= fft_complex_1d( fftInput, isIfft, bScale );
            
            for( uint iRow=0; iRow<nx; ++iRow )
                out[iRow][iCol]= fftResult[iRow];
        }
        
        return out;
    }
}

#endif
