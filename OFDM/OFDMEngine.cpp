//
//  OFDMEngine.cpp
//  OFDM
//
//  Created by Stephen MacKinnon on 1/17/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//

#include "OFDMEngine.h"
#include "boost/multi_array.hpp"
#include <complex>                  // Must include complex before fftw3 so fftw3 treats
#include "api/fftw3.h"              // fftw3_complex as c++ complex
#include <stdlib.h>
#include <time.h>
#include <exception>
#include "expection_definitions.h"
#include "MatrixHelper.h"

using namespace std;

OFDMEngine::OFDMEngine() {
    
}

std::vector<double> OFDMEngine::Modulate( unsigned char *pData, long lDataLength ) {
    // Symbols per carrier for this frame
    uint uCarrierSymbCount= ceil( lDataLength/CARRIER_COUNT );
    
    // 2D array
    vector<vector<unsigned char>> dataMatrix( uCarrierSymbCount+1, vector<unsigned char>(CARRIER_COUNT, 0) );
    
    // Populate 2D array with data. Serial to parallel
    uint uNumIter= 0;
    for( uint iRow= 0; iRow< dataMatrix.size(); ++iRow ) {
        for( uint iCol= 0; iCol<dataMatrix[0].size(); ++iCol ) {
            // Populate first row with 0s initially,
            if( iRow == 0 )
                dataMatrix[iRow][iCol]= 0;
            // If we run out of data, populate remainder of array with 0s
            else if( uNumIter-dataMatrix[0].size()+1 > lDataLength )
                dataMatrix[iRow][iCol]= 0;
            else
                dataMatrix[iRow][iCol]= static_cast<unsigned char>( pData[uNumIter-dataMatrix[0].size()] );
            ++uNumIter;
        }
    }
    
//    // DEBUGGING: print out matrix
//    cout<<"2D Array";
//    // MatrixHelper::print2dMatrix(dataMatrix);
//    for( uint iRow= 0; iRow<dataMatrix.size(); ++iRow ) {
//        cout<<endl;
//        for( uint iCol= 0; iCol<dataMatrix[0].size(); ++iCol ) {
//            cout<<static_cast<int>( dataMatrix[iRow][iCol] )<<"   ";
//        }
//    }
    
    ////////////////////////////////////////////////
    //
    //      Differential Encoding
    //
    ////////////////////////////////////////////////
    
    // Add diff ref row
    ++uCarrierSymbCount;
    
    // Seed for random number generator
    srand( static_cast<uint>(time(NULL)) );
    
    // Store diff ref in first row of matrix
    for( uint iCol=0; iCol<dataMatrix[0].size(); ++iCol ) {
        dataMatrix[0][iCol]= ceil( ( rand()/static_cast<float>(RAND_MAX) )*pow(2, SYMB_SIZE)+0.5f );
    }
    
    // DEBUGGING
    dataMatrix[0][0]=2;
    dataMatrix[0][1]=4;
    //dataMatrix[0][2]=4;
    //dataMatrix[0][3]=4;
    
    // Differential encoding using diff ref
    for( uint iRow= 1; iRow<dataMatrix.size(); ++iRow ) {
        for( uint iCol= 0; iCol<dataMatrix[0].size(); ++iCol ) {
            dataMatrix[iRow][iCol]= (dataMatrix[iRow][iCol]+dataMatrix[iRow-1][iCol]) % static_cast<uint>(pow(2, SYMB_SIZE));
        }
    }
    
    // DEBUGGING: print out matrix
    cout<<"2D Array: After differential encoding";
    for( uint iRow= 0; iRow< uCarrierSymbCount; ++iRow ) {
        cout<<endl;
        for( uint iCol= 0; iCol< CARRIER_COUNT; ++iCol ) {
            cout<<static_cast<float>(dataMatrix[iRow][iCol])<<"   ";
        }
    }
    
    ////////////////////////////////////////////////
    //
    //      PSK Modulation
    //
    ////////////////////////////////////////////////
    
    //complex_float_array complexMatrix(boost::extents[uCarrierSymbCount][CARRIER_COUNT]);
    vector<vector<complex<float>>> complexMatrix( uCarrierSymbCount, vector<complex<float>>( CARRIER_COUNT, 0));
    //MatrixHelper::init2dMatrix(complexMatrix, uCarrierSymbCount, CARRIER_COUNT);

    for( uint iRow=0; iRow<complexMatrix.size(); ++iRow ) {
        for( uint iCol= 0; iCol<complexMatrix[0].size(); ++iCol ) {
            complexMatrix[iRow][iCol]= polar( 1.0f, static_cast<float>( dataMatrix[iRow][iCol]*(2.0f*M_PI/pow(2.0f, SYMB_SIZE))) );
        }
    }
    
    // DEBUGGING: print out matrix
//    cout<<"After PSK Modulation: "<<endl;
//    for( uint iRow=0; iRow<uCarrierSymbCount; ++iRow ) {
//        cout<<endl;
//        for( uint iCol=0; iCol<CARRIER_COUNT; ++iCol ) {// CC = 3
//            cout<<static_cast<int>(complexMatrix[iRow][iCol].real())<<" + "<<static_cast<int>(complexMatrix[iRow][iCol].imag())<<"i"<<"   ";
//        }
//    }
    
    ///////////////////////////////////////////////////////////
    //
    //      Assign IFFT bins to carriers and imaged carriers
    //
    ///////////////////////////////////////////////////////////
    
    // Initialize 2D spectrum_tx array
    vector< vector<complex<double>> > spectrumTx( uCarrierSymbCount, vector<complex<double>>(IFFT_SIZE, 0));
    
    cout<<endl<<"Carriers:"<<endl;
    for( uint i=0; i<CARRIER_COUNT; ++i )
        cout<<CARRIERS[i]<<", ";
    cout<<endl<<"Conj carriers:"<<endl;
    for( uint i=0; i<CARRIER_COUNT; ++i )
        cout<<CONJ_CARRIERS[i]<<", ";
    
    // Matlab: spectrum_tx(:,carriers) = complex_matrix;
    for( uint iRow=0; iRow<uCarrierSymbCount; ++iRow ) {
        for( uint iCol=0; iCol<CARRIER_COUNT; ++iCol ) {
            spectrumTx[iRow][ CARRIERS[iCol]-1 ]= complexMatrix[iRow][iCol];
        }
    }
    // Matlab: spectrum_tx(:,conj_carriers) = conj(complex_matrix);
    for( uint iRow=0; iRow<uCarrierSymbCount; ++iRow ) {
       for( uint iCol=0; iCol<CARRIER_COUNT; ++iCol ) {
            spectrumTx[iRow][ CONJ_CARRIERS[iCol]-1 ]= conj( complexMatrix[iRow][iCol] );
        }
    }
    
    // Perform IFFT on conjugate transpose of spectrumTx
    vector< vector<complex<double>> > spectrumTx_transp= MatrixHelper::conjTranspose2d<double>(spectrumTx);
    
    // DEBUGGING: print out matrix
    cout<<endl<<"Spectrum tx' pre ifft: "<<endl;
    cout<<"Dimensions: "<<spectrumTx_transp.size()<<" x "<<spectrumTx_transp[0].size()<<endl;
    for( uint iRow=0; iRow<spectrumTx_transp.size(); ++iRow ) {
        cout<<endl;
        for( uint iCol=0; iCol<spectrumTx_transp[0].size(); ++iCol ) {// CC = 3
            spectrumTx_transp[iRow][iCol]= complex<double>( floor(spectrumTx_transp[iRow][iCol].real()), floor(spectrumTx_transp[iRow][iCol].imag()) );
            cout<<static_cast<int>(spectrumTx_transp[iRow][iCol].real())<<" + "<<static_cast<int>(spectrumTx_transp[iRow][iCol].imag())<<"i"<<"   ";
        }
    }
    
    //complex<double> *in= &spectrumTx_transp[0][0];
    double *out= (double*)fftw_malloc(sizeof(double)*2*spectrumTx_transp.size()*spectrumTx_transp[0].size());
    
//    cout<<"IN: "<<endl;
//    uNumIter= 0;
//    for( uint iRow=0; iRow<spectrumTx_transp.size(); ++iRow ) {
//        cout<<endl;
//        for( uint iCol=0; iCol<spectrumTx_transp[0].size(); ++iCol ) {
//            cout<<in[uNumIter].real()<<" + "<<in[uNumIter].imag()<<"i  ";
//            ++uNumIter;
//        }
//    }
    
    fftw_plan p= fftw_plan_dft_c2r_2d((int)spectrumTx_transp.size(), (int)spectrumTx_transp[0].size(), (fftw_complex*)&spectrumTx_transp[0][0], out, FFTW_ESTIMATE);
    
    fftw_execute(p);
    
    fftw_destroy_plan(p);
    
    // DEBUGGING: print IFFT results
    cout<<"IFFT result: ";
    for( uint iRow=0; iRow<uCarrierSymbCount; ++iRow ) {
        cout<<endl;
        for( uint iCol=0; iCol<IFFT_SIZE; ++iCol ) {
            cout<<out[iRow*iCol]<<"    ";
        }
    }
    
    ///////////////////////////////////////////////////////////
    //
    //      Add a periodic guard time
    //
    ///////////////////////////////////////////////////////////
    
    // MATLAB: end_symb = IFFT_SIZE

    vector< vector<double> > signalTx;
    signalTx.resize(uCarrierSymbCount);
    for( uint i=0; i<uCarrierSymbCount; ++i )
        signalTx[i].resize(IFFT_SIZE+GUARD_TIME);
    
    // Populate signalTx
    for( uint iRow=0; iRow<uCarrierSymbCount; ++iRow ) {
        for( uint iCol=0; iCol<IFFT_SIZE+GUARD_TIME; ++iCol ) {
            if( iCol < GUARD_TIME )
                signalTx[iRow][iCol]= out[iRow*(IFFT_SIZE-GUARD_TIME+iCol)];
            else
                signalTx[iRow][iCol]= out[iRow*(iCol-GUARD_TIME)];
        }
    }
    
    // Transpose signalTx
    signalTx= MatrixHelper::transpose2d(signalTx);
    
    // Reshape along columns
    // i.e.
    // 1 2 3 ---> 1 4 2 5 3 6
    // 4 5 6
    vector<double> output;
    output.resize(signalTx.size()*signalTx[0].size());
    for( uint iRow=0; iRow<signalTx.size(); ++iRow ) {
        for( uint iCol=0; iCol<signalTx[0].size(); ++iCol ) {
            output[iRow*iCol]= signalTx[iRow][iCol];
        }
    }
    return output;
}

void OFDMEngine::Demodulate( unsigned char *data, long lDataLength ) {
    
}