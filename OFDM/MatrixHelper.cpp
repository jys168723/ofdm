//
//  MatrixHelper.cpp
//  OFDM
//
//  Created by Stephen MacKinnon on 1/21/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//

#include "MatrixHelper.h"

namespace MatrixHelper {
template <typename T>
static vector<vector<complex<T>>> conjTranspose2d(vector<vector<complex<T>>> &in) {
    // Initialize output vector
    vector<vector<complex<T>>> out;
    out.resize(in[0].size());
    for( uint i=0; i<out.size(); ++i )
        out[i].resize(in.size());
    
    for( uint iRow=0; iRow<out.size(); ++iRow ) {
        for( uint iCol=0; iCol<out[0].size(); ++iCol )
            out[iRow][iCol]= conj(in[iCol][iRow]);
    }
    return out;
}
}