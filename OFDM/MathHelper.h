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
    inline double round( double val ) {
        double decimal= val - floor(val);
        double whole= floor(val);
        if( decimal >= 0.5 )
            whole+= 1;

        return whole;
    } // end function round()


    // Function mimics Matlab's rem() function, which is equivalent to the modulo (%)
    // operator, except for when the numerator is negative
    inline int rem( int numerator, int denom ) {
        int result= numerator % denom;
        if( numerator < 0 )
            result*= -1;
        return result;
    } // end function rem()
    
    
    // Function sets values whos absolute value is < 0.000001 to 0
    inline void normalize( double &val ) {
        if( abs(val) < 0.000001 )
            val= 0;
    } // end function normalize()
}


#endif
