//
//  OFDMEngine.cpp
//  OFDM
//
//  Created by Stephen MacKinnon on 1/17/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//

#include "OFDMEngine.h"
#include <complex>                  // Must include complex before fftw3 so fftw3 treats
#include "api/fftw3.h"              // fftw3_complex as c++ complex
#include <stdlib.h>
#include <time.h>
#include "MatrixHelper.h"

using namespace std;

OFDMEngine::OFDMEngine() {
    
}

std::vector<double> OFDMEngine::Modulate( unsigned char *pData, long lDataLength ) {
    cout<<"Carriers: "<<endl;
    for( uint i=0; i<CARRIER_COUNT; ++i )
        cout<<CARRIERS[i]<<", ";
    cout<<"Conj Carriers: "<<endl;
    for( uint i=0; i<CARRIER_COUNT; ++i )
        cout<<CONJ_CARRIERS[i]<<", ";
    
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
    
    // DEBUGGING - manually setting diff ref to match Matlab's
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
    
    vector<vector<complex<float>>> complexMatrix( uCarrierSymbCount, vector<complex<float>>( CARRIER_COUNT, 0));

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
        for( uint iCol=0; iCol<spectrumTx_transp[0].size(); ++iCol ) {
            cout<<static_cast<int>(spectrumTx_transp[iRow][iCol].real())<<" + "<<static_cast<int>(spectrumTx_transp[iRow][iCol].imag())<<"i"<<"   ";
        }
    }
    
    
    ///////////////////////////////////////////////////////////
    //
    //      FFT/IFFT
    //
    ///////////////////////////////////////////////////////////
    
    int nx= static_cast<int>( spectrumTx_transp.size() );
    int nc= static_cast<int>( spectrumTx_transp[0].size() );
    int ny= (nc*2)-1;
    
    fftw_complex *in= (fftw_complex*)fftw_malloc(nx*nc*sizeof(fftw_complex));//&spectrumTx_transp[0][0];
    double *out= (double*)fftw_malloc(nx*ny*sizeof(double));
    fftw_complex *out2= (fftw_complex*)fftw_malloc(nx*nc*sizeof(fftw_complex));
    // Populate in
    cout<<"\n\nInput:";
    for( uint i=0; i<nx; ++i ) {
        cout<<endl;
        for( uint j=0; j<nc; ++j ) {
            in[i*nc+j][0]= floor( spectrumTx_transp[i][j].real() );
            in[i*nc+j][1]= floor( spectrumTx_transp[i][j].imag() );
            
            // Print out input
            cout<<in[i*nc+j][0]<<" + "<<in[i*nc+j][1]<<", ";
        }
    }
    
    // Make c2r plan
    fftw_plan ifftPlan= fftw_plan_dft_c2r_2d(nx, ny, in, out, FFTW_ESTIMATE);
    
    // Execute c2r plan
    fftw_execute(ifftPlan);
    
    // Print IFFT result
//    cout<<"\nIFFT result: ";
//    for( uint i=0; i<nx; ++i ) {
//        for( uint j=0; j<ny; ++j ) {
//            cout<<out[i*ny+j]/(double)(nx*ny)<<", ";
//        }
//    }
    
    // Make r2c plan
    fftw_plan fftPlan= fftw_plan_dft_r2c_2d(nx, ny, out, out2, FFTW_ESTIMATE);
    
    // Execute r2c plan
    fftw_execute(fftPlan);
    
    // Print FFT results
    cout<<"\nFFT result: ";
    for( uint i=0; i<nx; ++i ) {
        cout<<endl;
        for( uint j=0; j<nc; ++j ) {
            double real= out2[i*nc+j][0]/(nx*ny);
            double imag= out2[i*nc+j][1]/(nx*ny);
            normalize(real);
            normalize(imag);
            cout<<(int)real<<" + "<<(int)imag<<", ";
        }
    }
    
    fftw_destroy_plan(ifftPlan);
    fftw_destroy_plan(fftPlan);
    
    
    ///////////////////////////////////////////////////////////
    //
    //      Add a periodic guard time
    //
    ///////////////////////////////////////////////////////////
    
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

void OFDMEngine::Demodulate( std::vector<double> *data, long lDataLength ) {
    uint uSymbPeriod= IFFT_SIZE + GUARD_TIME;
    
    // Reshape the linear time waveform into FFT segments
    uint uNumCol= floor( data->size()/(float)uSymbPeriod );
    vector<vector<double>> symbRxMatrix(uSymbPeriod, vector<double>(uNumCol, 0));
    
    // TODO: Verify whether this should be a row or column-wise reshape
    uint uNumIter= 0;
    for( uint iRow=0; iRow<symbRxMatrix.size(); ++iRow ) {
        for( uint iCol=0; iCol<symbRxMatrix[0].size(); ++iCol ) {
            symbRxMatrix[iRow][iCol]= (*data)[uNumIter];
        }
    }
    // Remove the periodic guard time
}

void OFDMEngine::FFTTest() {
    // r2c - input size = n real numbers, output size = n/2+1 complex
    // allocate enough memory in realIn for n/2+1= 5 complex numbers
    int n= 8;
    int nc= (n/2)+1;
    double* realIn= (double*)fftw_malloc(sizeof(double)*n);
    fftw_complex* complexOut= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nc);
    
    // Populate real input
    for( uint i=0; i<n; ++i )
        realIn[i]= i;    
    
    // Create fftw plans
    fftw_plan fftPlan= fftw_plan_dft_r2c_1d(n, realIn, complexOut, FFTW_ESTIMATE);
    
    // Verification
    cout<<"Input:"<<endl;
    for( uint i=0; i<8; ++i )
        cout<<realIn[i]<<", ";
    cout<<endl;
    
    // Execute FFT
    fftw_execute(fftPlan);
    
    // Print result
    // 0th and n/2-th elements of the complex output are purely real
    cout<<"FFT Result"<<endl;
    for( uint i=0; i<nc; ++i ) {
        cout<<complexOut[i][0]<<" + "<<complexOut[i][1]<<"i ";
    }
    cout<<endl;
    
    // Allocate IFFT in
    //double* realOut= (double*)fftw_malloc(sizeof(fftw_complex)*nc);
    double* in2= (double*)(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nc);
    
    fftw_plan ifftPlan= fftw_plan_dft_c2r_1d(n, (fftw_complex*)complexOut, (double*)in2, FFTW_ESTIMATE);
    
    fftw_execute(ifftPlan);
    
    // Print IFFT results
    cout<<"IFFT Results"<<endl;
    for( uint i=0; i<n; ++i )
        cout<<in2[i]/(double)(n)<<", ";
    
    fftw_free(realIn);
    fftw_free(complexOut);
    fftw_free(in2);
    
}

void OFDMEngine::normalize( double &val ) {
    if( val < 0.000001 && val > -0.000001 )
        val= 0;
}