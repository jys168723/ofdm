//
//  MatrixHelper.h
//  OFDM
//
//  Created by Stephen MacKinnon on 1/21/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//

#ifndef __OFDM__MatrixHelper__
#define __OFDM__MatrixHelper__

#include <iostream>
#include <vector>
#include <complex>
#include "expection_definitions.h"

using namespace std;

namespace MatrixHelper {
    
    // Function accepts an MxN complex vector in and returns its NxM conjugate transpose
    template <typename T>
    extern vector<vector<complex<T>>> conjTranspose2d(vector<vector<complex<T>>> &in) {
        // Initialize output vector
        vector<vector<complex<T>>> out( in[0].size(), vector<complex<T>>(in.size()) );
        
        for( uint iRow=0; iRow<out.size(); ++iRow ) {
            for( uint iCol=0; iCol<out[0].size(); ++iCol )
                out[iRow][iCol]= conj(in[iCol][iRow]);
        }
        return out;
    } // end conjTranspose2d()
    
    // Function accepts an MxN vector in and returns its NxM transpose
    template <typename T>
    extern vector<vector<T>> transpose2d(vector<vector<T>> &in) {
        // Initialize output vector
        vector<vector<T>> out;
        out.resize(in[0].size());
        for( uint i=0; i<out.size(); ++i )
            out[i].resize(in.size());
        
        for( uint iRow=0; iRow<out.size(); ++iRow ) {
            for( uint iCol=0; iCol<out[0].size(); ++iCol )
                out[iRow][iCol]= in[iCol][iRow];
        }
        return out;
    } // end transpose2d()
    
    // Function accepts a 2d vector and returns its MxN reshape
    // Vector must contain M*N elements
    template <typename T>
    extern vector<vector<T>> reshape2d(vector<vector<T>> &in, uint m, uint n) {
        if( in.size()*in[0].size() != m*n ) {
            throw EXCEPTION_ARRAY_MISMATCH;
            return NULL;
        }
        
        // Initialize output vector
        vector<vector<T>> out;
        out.resize(m);
        for( uint i=0; i<m; ++i )
            out[i].resize(n);
        
        for( uint iRow=0; iRow<out.size(); ++iRow ) {
            for( uint iCol=0; iCol<out[0].size(); ++iCol )
                out[iRow][iCol]= in[iCol][iRow];
        }
        return out;
    } // end reshape2d()

} // end namespace MatrixHelper

#endif /* defined(__OFDM__MatrixHelper__) */
