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
#include "MathHelper.h"

using namespace std;

namespace MatrixHelper {
    
    // Function accepts an MxN complex vector in and returns its NxM conjugate transpose
    template <typename T>
    inline vector<vector<complex<T>>> conj_transpose( vector<vector<complex<T>>> &in ) {
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
    extern vector<vector<T>> transpose( vector<vector<T>> &in ) {
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
    extern vector<vector<T>> reshape( vector<vector<T>> &in, uint m, uint n ) {
        if( in.size()*in[0].size() != m*n ) {
            throw EXCEPTION_ARRAY_DIMENSION_MISMATCH;
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
    
    
    // Function returns a 1d representation of a 2d vector. Reshape is column-wise
    // by default, but can be row-wise by passing false to arg isColumnwise
    template <typename T>
    extern vector<double> flatten( vector<vector<T>> &in, bool isColumnwise= true ) {
        vector<T> out( in.size()*in[0].size(), 0 );
        
        uint uNumIter= 0;
        // Column-wise
        if( isColumnwise ) {
            for( uint iCol=0; iCol<in[0].size(); ++iCol ) {
                for( uint iRow=0; iRow<in.size(); ++iRow ) {
                    out[uNumIter]= in[iRow][iCol];
                    ++uNumIter;
                }
            }
        }
        // Row-wise
        else {
            for( uint iRow=0; iRow<in.size(); ++iRow ) {
                for( uint iCol=0; iCol<in[0].size(); ++iCol ) {
                    out[uNumIter]= in[iRow][iCol];
                    ++uNumIter;
                }
            }
        }
        
        return out;
    } // end function reshape2dTo1d()
    
    
    // Function accepts a 1d vector and returns its MxN column-wise reshape
    template <typename T>
    extern vector<vector<T>> reshape( vector<T> &in, uint m, uint n ) {
        if( in.size() != m*n ) {
            throw EXCEPTION_ARRAY_DIMENSION_MISMATCH;
            return vector<vector<T>>();
        }
        
        uint uNumIter= 0;
        // Initialize output vector
        vector<vector<T>> out( m, vector<T>(n, 0));
        for( uint iCol=0; iCol<n; ++iCol ) {
            for( uint iRow=0; iRow<m; ++iRow ) {
                out[iRow][iCol]= in[uNumIter];
                ++uNumIter;
            }
        }
        return out;
    } // end function reshape()
    
    
    // Function prints out the elements in a 1d vector
    template <typename T>
    extern void print_vector( vector<T> &in, bool bNewLines= false ) {
        cout<<endl;
        for( uint i=0; i<in.size(); ++i ) {
            if( bNewLines )
                cout<<endl<<i<<": ";
            cout<<in[i]<<"  ";
        }
    } // end function print_vector()
    
    
    // Function prints out the elements in a 2d vector to the console
    template <typename T>
    extern void print_vector( vector<vector<T>> &in ) {
        for( uint iRow=0; iRow<in.size(); ++iRow ) {
            cout<<endl<<iRow<<": ";
            for( uint iCol=0; iCol<in[0].size(); ++iCol )
                cout<<in[iRow][iCol]<<", ";
        }
    } // end function print_vector()

    
    // Function normalizes the real and imaginary components of a complex number
    template <typename T>
    inline void normalize( complex<T> &val ) {
        double real= val.real();
        double imag= val.imag();
        MathHelper::normalize(real);
        MathHelper::normalize(imag);
        val= complex<double>( real, imag );
    } // end function normalize()
    

    // Function normalizes each element in a vector of complex numbers
    template <typename T>
    inline void nomalize( vector<complex<T>> &in ) {
        for( uint i=0; i<in.size(); ++i ) {
            complex<T> value= in[i];
            normalizeComplex(value);
            in[i]= value;
        }
    } // end function normalize()
    
    
    // Function normalizes each element in a 2D vector of complex numbers
    template <typename T>
    inline void nomalize( vector<vector<complex<T>>> &in ) {
        for( uint iRow=0; iRow<in.size(); ++iRow ) {
            for( uint iCol=0; iCol<in[0].size(); ++iCol ) {
                complex<T> value= in[iRow][iCol];
                normalize(value);
                in[iRow][iCol]= value;
            }
        }
    } // end function normalize()
    
    
    // Returns a vector containing the difference and approximate derivative of an input vector
    template <typename T>
    inline vector<T> diff( vector<T> &in ) {
        vector<T> out( in.size()-1, 0 );
        for( uint i=0; i<out.size(); ++i ) {
            out[i]= in[i+1] - in[i];
        }
        return out;
    } // end function diff()
    
    
    // Returns a vector containing the column-wise difference and approximate derivative of an input vector
    template <typename T>
    inline vector<vector<T>> diff( vector<vector<T>> &in ) {
        vector<vector<T>> out(in.size()-1, vector<T>(in[0].size(), 0) );
        for( uint iCol=0; iCol<out[0].size(); ++iCol ) {
            for( uint iRow=0; iRow<out.size(); ++iRow ) {
                out[iRow][iCol]= in[iRow+1][iCol]-in[iRow][iCol];
            }
        }
        return out;
    } // end function diff()
    
} // end namespace MatrixHelper

#endif /* defined(__OFDM__MatrixHelper__) */
