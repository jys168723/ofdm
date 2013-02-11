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
#include "FFTHelper.h"

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
    dataMatrix[0][0]=1;
    dataMatrix[0][1]=1;
    //dataMatrix[0][2]=4;
    //dataMatrix[0][3]=4;
    
    // Differential encoding using diff ref
    for( uint iRow= 1; iRow<dataMatrix.size(); ++iRow ) {
        for( uint iCol= 0; iCol<dataMatrix[0].size(); ++iCol ) {
            dataMatrix[iRow][iCol]= (dataMatrix[iRow][iCol]+dataMatrix[iRow-1][iCol]) % static_cast<uint>(pow(2, SYMB_SIZE));
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
    vector< vector<complex<double>> > spectrumTx_transp= MatrixHelper::conj_transpose<double>(spectrumTx);
    
    ///////////////////////////////////////////////////////////
    //
    //      IFFT
    //
    ///////////////////////////////////////////////////////////
    
    vector<vector<complex<double>>> signal_tx_complex= FFTHelper::fft_complex_2d(spectrumTx_transp, true, true);
    
    //MatrixHelper::print2dVector(signal_tx_complex);
    
    signal_tx_complex= MatrixHelper::conj_transpose(signal_tx_complex);
    
    // Only need real part
    vector<vector<double>> signal_tx( signal_tx_complex.size(), vector<double>(signal_tx_complex[0].size(), 0) );
    for( uint i=0; i<signal_tx.size(); ++i ) {
        for( uint j=0; j<signal_tx[0].size(); ++j ) {
            double real= signal_tx_complex[i][j].real();
            MathHelper::normalize(real);
            signal_tx[i][j]= real;
        }
    }
    
    ///////////////////////////////////////////////////////////
    //
    //      Add a periodic guard time
    //
    ///////////////////////////////////////////////////////////
    
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
    
    signal_tx= MatrixHelper::transpose<double>( signal_tx );

    vector<double> signal_tx_1d= MatrixHelper::flatten(signal_tx);

    return signal_tx_1d;
} // end function Modulate()


void OFDMEngine::Demodulate( std::vector<double> &symbRx ) {
    
    uint uSymbPeriod= IFFT_SIZE + GUARD_TIME;
    
    // Reshape the linear time waveform into FFT segments
    int iNumCol= floor( symbRx.size()/(float)uSymbPeriod );
    vector<vector<double>> symbRxMatrix= MatrixHelper::reshape(symbRx, uSymbPeriod, iNumCol);
    
    // Remove the periodic guard time
    symbRxMatrix.erase( symbRxMatrix.begin(), symbRxMatrix.begin()+GUARD_TIME );
    
    // Create a matrix of complex data for FFT input. Set all imaginary values to 0
    vector<vector<complex<double>>> fftInputMatrix( symbRxMatrix.size(), vector<complex<double>>( symbRxMatrix[0].size(), 0) );
    for( uint iRow=0; iRow<fftInputMatrix.size(); ++iRow ) {
        for( uint iCol=0; iCol<fftInputMatrix[0].size(); ++iCol ) {
            fftInputMatrix[iRow][iCol]= complex<double>( symbRxMatrix[iRow][iCol], 0 );
        }
    }
    
    vector<vector<complex<double>>> fftResult= FFTHelper::fft_complex_2d(fftInputMatrix, false, false);
    
    MatrixHelper::nomalize(fftResult);   // TODO: Do we want to normalize (low values to 0)?
    fftResult= MatrixHelper::conj_transpose(fftResult);

    // Extract columns of data on IFFT bins of all carriers only
    vector<vector<complex<double>>> rx_spectrum_matrix( fftResult.size(), vector<complex<double>>(CARRIER_COUNT, 0) );
    for( uint iCarrier=0; iCarrier<CARRIER_COUNT; ++iCarrier ) {
        for( uint iRow=0; iRow<rx_spectrum_matrix.size(); ++iRow ) {
            rx_spectrum_matrix[iRow][iCarrier]= fftResult[iRow][CARRIERS[iCarrier]-1];
        }
    }
    
    cout<<"\nrx_spectum_matrix:\n";
    MatrixHelper::print_vector(rx_spectrum_matrix);
    
    ///////////////////////////////////////////////////////////
    //
    //  PSK Demodulation
    //
    ///////////////////////////////////////////////////////////
    
    vector<vector<double>> rx_phase( rx_spectrum_matrix.size(), vector<double>(rx_spectrum_matrix[0].size(), 0) );
    for( uint iRow=0; iRow<rx_phase.size(); ++iRow ) {
        for( uint iCol=0; iCol<rx_phase[0].size(); ++iCol ) {
            // TODO: We might not want to round here...
            int phaseRad= round( arg( rx_spectrum_matrix[iRow][iCol] )*(180.0f/(float)M_PI) );
            // Make negative phases positive
            rx_phase[iRow][iCol]= MathHelper::rem( static_cast<int>(phaseRad)+360, 360 );
        }
    }
    
    vector<vector<double>> decoded_phase= MatrixHelper::diff(rx_phase);
    
    // Make negative phases positive
    for( uint iRow=0; iRow<decoded_phase.size(); ++iRow ) {
        for( uint iCol=0; iCol<decoded_phase[0].size(); ++iCol ) {
            decoded_phase[iRow][iCol]= MathHelper::rem( decoded_phase[iRow][iCol]+360, 360 );
        }
    }
    
    MatrixHelper::print_vector(decoded_phase);
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
    vector<double> header( HEAD_LEN, 0 );
    double f= 0.25;
    for( uint i=0; i<HEAD_LEN; ++i ) {
        header[i]= sin( i*2*M_PI*f );
    }
    
    f= f / ( M_PI*2.0/3.0 );
    for( uint i= 0; i<HEAD_LEN; ++i ) {
        header[i]= header[i] + sin( i*2*M_PI*f );
    }
    
    return header;
} // end OFDMEngine::GenerateHeader()


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