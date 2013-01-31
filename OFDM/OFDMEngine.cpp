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
#include <algorithm>
#include "MatrixHelper.h"
#include "MathHelper.h"

using namespace std;

OFDMEngine::OFDMEngine() {
    
}

// Function takes in byte data, modulates it, and returns a pointer to the
// modulated, double-precision data
vector<double> OFDMEngine::Modulate( unsigned char *pData, int iDataLength ) {
    // Print out carriers
    cout<<"Carriers: "<<endl;
    for( uint i=0; i<CARRIER_COUNT; ++i )
        cout<<CARRIERS[i]<<", ";
    cout<<"Conj Carriers: "<<endl;
    for( uint i=0; i<CARRIER_COUNT; ++i )
        cout<<CONJ_CARRIERS[i]<<", ";
    
    // Symbols per carrier for this frame
    uint uCarrierSymbCount= ceil( iDataLength/CARRIER_COUNT );
    
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
            else if( uNumIter-dataMatrix[0].size()+1 > iDataLength )
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
    dataMatrix[0][0]=4;
    dataMatrix[0][1]=1;
    //dataMatrix[0][2]=4;
    //dataMatrix[0][3]=4;
    
    // Differential encoding using diff ref
    for( uint iRow= 1; iRow<dataMatrix.size(); ++iRow ) {
        for( uint iCol= 0; iCol<dataMatrix[0].size(); ++iCol ) {
            dataMatrix[iRow][iCol]= (dataMatrix[iRow][iCol]+dataMatrix[iRow-1][iCol]) % static_cast<uint>(pow(2, SYMB_SIZE));
        }
    }
    
    // DEBUGGING: print out matrix
//    cout<<"2D Array: After differential encoding";
//    for( uint iRow= 0; iRow< uCarrierSymbCount; ++iRow ) {
//        cout<<endl;
//        for( uint iCol= 0; iCol< CARRIER_COUNT; ++iCol ) {
//            cout<<static_cast<float>(dataMatrix[iRow][iCol])<<"   ";
//        }
//    }
    
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
    
    // DEBUGGING: Print out matrix
    //cout<<"After assign IFFT bins:";
    for( uint iRow=0; iRow<spectrumTx.size(); ++iRow ) {
        //cout<<endl;
        //cout<<"Row "<<iRow<<": ";
        for( uint iCol=0; iCol<spectrumTx[0].size(); ++iCol ) {
            double real= spectrumTx[iRow][iCol].real();
            double imag= spectrumTx[iRow][iCol].imag();
            normalize(real);
            normalize(imag);
            spectrumTx[iRow][iCol]= complex<double>(real, imag);
            //cout<<real<<" + "<<imag<<"i, ";
        }
    }
    
    // Perform IFFT on conjugate transpose of spectrumTx
    vector< vector<complex<double>> > spectrumTx_transp= MatrixHelper::conjTranspose2d<double>(spectrumTx);
    
    // DEBUGGING: print out matrix
//    cout<<endl<<"Spectrum tx' pre ifft: "<<endl;
//    cout<<"Dimensions: "<<spectrumTx_transp.size()<<" x "<<spectrumTx_transp[0].size()<<endl;
//    for( uint iRow=0; iRow<spectrumTx_transp.size(); ++iRow ) {
//        cout<<endl;
//        cout<<"Row "<<iRow<<": ";
//        for( uint iCol=0; iCol<spectrumTx_transp[0].size(); ++iCol ) {
//            cout<<spectrumTx_transp[iRow][iCol].real()<<" + "<<spectrumTx_transp[iRow][iCol].imag()<<", ";
//        }
//    }
    
    
    ///////////////////////////////////////////////////////////
    //
    //      IFFT
    //
    ///////////////////////////////////////////////////////////
    
    int nx= static_cast<int>( spectrumTx_transp.size() );
    int nc= static_cast<int>( spectrumTx_transp[0].size() );

    fftw_complex *in;
    in= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nx*nc);
    
    // Perform IFFT in-place
    fftw_plan p= fftw_plan_dft_2d(nx, nc, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    // Populate in
    //cout<<"\n\nInput:";
    for( uint i=0; i<nx; ++i ) {
        //cout<<endl;
        //cout<<"Row "<<i<<": ";
        for( uint j=0; j<nc; ++j ) {
            in[i*nc+j][0]= spectrumTx_transp[i][j].real();
            in[i*nc+j][1]= spectrumTx_transp[i][j].imag();
            
            // Print out input
            //cout<<i*nc+j<<": "<<in[i*nc+j][0]<<" + "<<in[i*nc+j][1]<<", ";
        }
    }

    // Execute IFFT
    fftw_execute(p);
    
    // signal_tx = real(in)
    vector<vector<double>> signal_tx( nx, vector<double>(nc, 0) );
    for( uint i=0; i<nx; ++i ) {
        for( uint j=0; j<nc; ++j ) {
            signal_tx[i][j]= in[i*nc+j][0];
        }
    }
    
    // DEBUGGING: Verify IFFT by taking FFT
//    fftw_complex *out;
//    out= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nx*nc);
//    fftw_plan p2 = fftw_plan_dft_2d(nx, nc, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//    fftw_execute(p2);
//   
//    cout<<"Output:"<<endl;
//    for( uint i=0; i<nx; ++i ) {
//        cout<<endl;
//        cout<<"Row "<<i<<": ";
//        for( uint j=0; j<nc; ++j ) {
//            double real= out[i*nc+j][0]/(nx*nc);
//            double imag= out[i*nc+j][1]/(nx*nc);
//            normalize(real);
//            normalize(imag);
//            cout<<real<<" + "<<imag<<"i, ";
//            //cout<<out[i*nc+j][0]/(nx*nc)<<" + "<<out[i*nc+j][1]<<", ";
//        }
//    }
    
    // Transpose signal_tx
    signal_tx= MatrixHelper::transpose2d<double>( signal_tx );
    
    ///////////////////////////////////////////////////////////
    //
    //      Add a periodic guard time
    //
    ///////////////////////////////////////////////////////////
    
    // Print matrix
//    cout<<"Matrix without guard time:";
//    for( uint i=0; i<signal_tx.size(); ++i ) {
//        cout<<endl<<"Row "<<i<<": ";
//        for( uint j=0; j<signal_tx[0].size(); ++j ) {
//            cout<<signal_tx[i][j]<<", ";
//        }
//    }
    
    // Prepend columns (last_col - GUARD_TIME : last_col) to signal_tx
    // TODO: There is probably a better way to do this
    vector<vector<double>> prepend_matrix( signal_tx.size(), vector<double>(GUARD_TIME, 0) );
    for( uint i=0; i<prepend_matrix.size(); ++i ) {
        for( uint j=0; j<prepend_matrix[0].size(); ++j )
            prepend_matrix[i][j]= signal_tx[i][signal_tx[0].size()-GUARD_TIME+j];
    }
    
    for( uint i=0; i<signal_tx.size(); ++i ) {
        signal_tx[i].reserve( signal_tx[i].size()+GUARD_TIME );
        for( uint j=0; j<GUARD_TIME; ++j )
            signal_tx[i].insert( signal_tx[i].begin()+j, prepend_matrix[i][j] );
    }
    
    signal_tx= MatrixHelper::transpose2d<double>( signal_tx );
    
    // Print matrix
    cout<<"Signal tx transp:";
    for( uint i=0; i<signal_tx.size(); ++i ) {
        cout<<endl<<"Row "<<i<<": ";
        for( uint j=0; j<signal_tx[0].size(); ++j ) {
            cout<<signal_tx[i][j]<<", ";
        }
    }
        
    vector<double> signal_tx_1d( signal_tx.size()*signal_tx[0].size(), 0 );
    // Populate signal_tx_1d column-wise
    uNumIter= 0;
    for( uint j=0; j<signal_tx[0].size(); ++j ) {
        for( uint i=0; i<signal_tx.size(); ++i ) {
            signal_tx_1d[uNumIter]=signal_tx[i][j];
            ++uNumIter;
        }
    }
    
//    cout<<"Signal tx 1d"<<endl;
//    for( uint i=0; i<signal_tx_1d.size(); ++i )
//        cout<<signal_tx_1d[i]<<", ";
    
    return signal_tx_1d;
}

void OFDMEngine::Demodulate( std::vector<double> *symbRx ) {
   
    uint uSymbPeriod= IFFT_SIZE + GUARD_TIME;
    
    // Reshape the linear time waveform into FFT segments
    uint uNumCol= floor( symbRx->size()/(float)uSymbPeriod );
    vector<vector<double>> symbRxMatrix(uSymbPeriod, vector<double>(uNumCol, 0));

    uint uNumIter= 0;
    for( uint iCol=0; iCol<symbRxMatrix[0].size(); ++iCol ) {
        for( uint iRow=0; iRow<symbRxMatrix.size(); ++iRow ) {
            symbRxMatrix[iRow][iCol]= (*symbRx)[uNumIter];
            ++uNumIter;
        }
    }
    
    // DEBUGGING: Print out matrix
    cout<<"SymbRxMatrix:";
    for( uint iRow=0; iRow<symbRxMatrix.size(); ++iRow ) {
        cout<<endl;
        for( uint iCol=0; iCol<symbRxMatrix[0].size(); ++iCol ) {
            cout<<symbRxMatrix[iRow][iCol]<<", ";
        }
    }
    
    // Remove the periodic guard time
    for( uint i=0; i<GUARD_TIME; ++i )
        symbRxMatrix.erase( symbRxMatrix.begin() );
    
    // DEBUGGING: Print out matrix
    cout<<"SymbRxMatrix:";
    for( uint iRow=0; iRow<symbRxMatrix.size(); ++iRow ) {
        cout<<endl;
        for( uint iCol=0; iCol<symbRxMatrix[0].size(); ++iCol ) {
            cout<<symbRxMatrix[iRow][iCol]<<", ";
        }
    }
    
    // Take FFT of the received time wave to obtain data spectrum
    int nc= static_cast<int>( symbRxMatrix[0].size()/2.0f + 1 );
    int ny= static_cast<int>( symbRxMatrix[0].size() );
    int nx= static_cast<int>( symbRxMatrix.size() );
    fftw_complex *fftResult= (fftw_complex*)fftw_malloc( nc*nx*sizeof(fftw_complex) );
    double *fftInput= (double*)&symbRxMatrix;
    
    fftw_plan fftPlan= fftw_plan_dft_r2c_2d( nx, ny, fftInput, fftResult, FFTW_ESTIMATE );
    fftw_execute(fftPlan);
    
    // Print FFT results
    cout<<"\nFFT result: ";
    for( uint i=0; i<nx; ++i ) {
        cout<<endl;
        for( uint j=0; j<ny; ++j ) {
            double real= fftResult[i*nc+j][0]/(nx*ny);
            double imag= fftResult[i*nc+j][1]/(nx*ny);
            normalize(real);
            normalize(imag);
            cout<<real<<" + "<<imag<<", ";
        }
    }
    
    //    fftw_complex *out2= (fftw_complex*)fftw_malloc(nx*nc*sizeof(fftw_complex));
    //    // Make r2c plan
    //    fftw_plan fftPlan= fftw_plan_dft_r2c_2d(nx, ny, out, out2, FFTW_ESTIMATE);
    //
    //    // Execute r2c plan
    //    fftw_execute(fftPlan);
    //
    //    // Print FFT results
    //    cout<<"\nFFT result: ";
    //    for( uint i=0; i<nx; ++i ) {
    //        cout<<endl;
    //        for( uint j=0; j<nc; ++j ) {
    //            double real= out2[i*nc+j][0]/(nx*ny);
    //            double imag= out2[i*nc+j][1]/(nx*ny);
    //            normalize(real);
    //            normalize(imag);
    //            cout<<(int)real<<" + "<<(int)imag<<", ";
    //        }
    //    }
    //    fftw_destroy_plan(fftPlan);
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

double OFDMEngine::FrameDetect( std::vector<double>* data ) {
    // Take abs of data
    vector<double> signal( data->size(), 0 );
    for( uint i=0; i<data->size(); ++i )
        signal[i]= abs( (*data)[i] );
    
    // Sampled version of the signal
    vector<double> sampSignal( ceil(signal.size()/(float)ENVELOPE), 0 );
    for( uint i=0; i<sampSignal.size(); ++i )
        sampSignal[i]= signal[i*ENVELOPE];
    
    vector<double> ones( round(SYMB_PERIOD/static_cast<double>(ENVELOPE)), 1 );
    vector<double> movSum= filter( ones, 1, sampSignal );
    movSum.erase( movSum.begin(), movSum.begin()+round(SYMB_PERIOD/static_cast<double>(ENVELOPE)) );
    
    double minElement= *min_element( movSum.begin(), movSum.end() );
    double apprx= minElement*ENVELOPE+SYMB_PERIOD;

    // Move back by approximately 110% of the symbol period to start searching
    uint idxStart= round( apprx-1.1*SYMB_PERIOD );
    
    // Look into the narrowed down window
    // Matlab: mov_sum = filter(ones(1,symb_period),1,signal(idx_start:round(apprx+symb_period/3)));
    ones.resize( SYMB_PERIOD );
    //generate( ones.begin(), ones.end(), 1 ); // Populates ones with ones
    movSum.resize( round(apprx+SYMB_PERIOD/3.0f) );
    signal.erase( signal.begin(), signal.begin()+idxStart );
    signal.erase( signal.begin()+round(apprx+SYMB_PERIOD/3.0f), signal.end() );
    movSum= filter( ones, 1, signal );
    movSum.erase( movSum.begin(), movSum.begin()+SYMB_PERIOD );
    
    double nullSig= *min_element( movSum.begin(), movSum.end() );
    //double startSymb= min( )
    return 0;
} // end OFDMEngine::FrameDetect()


// Function generates a header and trailer (exact copy of the header)
vector<double> OFDMEngine::GenerateHeader() {
    vector<double> header( 2*(IFFT_SIZE+GUARD_TIME), 0 );
    double f= 0.25;
    for( uint i=0; i<IFFT_SIZE+GUARD_TIME; ++i ) {
        header[i]= sin(i*2*M_PI*f);
    }
    
    f= f/(M_PI*2.0/3.0);
    for( uint i= IFFT_SIZE+GUARD_TIME; i<2*(IFFT_SIZE+GUARD_TIME); ++i ) {
        header[i]= sin(i*2*M_PI*f);
    }
    
    return header;
} // end OFDMEngine::GenerateHeader()


// Function "normalizes" low numbers to 0
void OFDMEngine::normalize( double &val ) {
    if( abs(val) < 0.000001 )
        val= 0;
} // end OFDMEngine::normalize()


///////////////////////////////////////////////////////////////////////////////////////
//
// Implementation of matlab's filter() function, which is described
// by the difference equation:
//
// a(1)y(n) = b(1)x(n) + b(2)x(n-1) + ... + b(Nb)x(n-Nb+1)
// - a(2)y(n-1) - ... - a(Na)y(n-Na+1)
//
// For our purposes, we only need a to be a scalar, therefore
// y(n) = [b(1)x(n) + b(2)x(n-1) + ... + b(Nb)x(n-Nb+1)]/a
//
// Note: matlab vectors are 1-indexed
///////////////////////////////////////////////////////////////////////////////////////

vector<double> OFDMEngine::filter( vector<double> &b, double a, vector<double> &x ) {
    vector<double> y( x.size(), 0 );
    // 20 - 5 + 1
    for( int iy=0; iy<y.size(); ++iy ) {
        double val= 0;
        for( int ib=0; ib<b.size(); ++ib ) {
            if( iy-ib >= 0 )
                val= val + b[ib]*x[iy-ib];
        }
        y[iy]= val / a;
    }
    return y;
} // end OFDMEngine::filter()
