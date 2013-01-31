//
//  MathHelper.h
//  OFDM
//
//  Created by Stephen MacKinnon on 1/29/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//

#ifndef OFDM_MathHelper_h
#define OFDM_MathHelper_h

#include <math.h>

using namespace std;

namespace MathHelper {
    
    // Function mimics Matlab's round() function, which rounds up numbers whose
    // decimal part is >= 0.5
    extern double Round( double val ) {
        double decimal= val - floor(val);
        double whole= floor(val);
        if( decimal >= 0.5 )
            whole+= 1;

        return whole;
    }
    
    // Float implementation of round
    extern float Round( float val ) {
        return static_cast<float>( Round(static_cast<double>(val)) );
    }
}


#endif
