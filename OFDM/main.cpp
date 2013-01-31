//
//  main.cpp
//  OFDM
//
//  Created by Stephen MacKinnon on 1/15/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//

#include <algorithm>
#include <iostream>
#include <fstream> 
#include "OFDMEngine.h"

using namespace std;

// Forward declarations
unsigned char* readDataFromFile( const char* filename );
bool writeDataToFile( unsigned char* data, const char *filename );
double variance( vector<double> &in );

// Hacky for now, will be unnecessary when we establish a
// fixed buffer size
int iDataLength= 64;

int main(int argc, const char * argv[])
{
    
    OFDMEngine* pEngine= new OFDMEngine();
    unsigned char data[64];
     
    // Populate data
    // NOTE: I'm not sure if our input data is going to be a 1d array
    //       or a matrix, if it is a matrix, there is some additional
    //       reshaping that must take place (see w, h in matlab code)
    cout<<"Raw input data:"<<endl;
    for( uint i=0; i<64; ++i ) {
        data[i]= static_cast<unsigned char>(i+1);
        cout<<i<<": "<<static_cast<float>(data[i])<<endl;
    }
    unsigned char* pData= &data[0];
    
    cout<<"Carriers:"<<endl;
    for( uint i=0; i<CARRIER_COUNT; ++i )
        cout<<CARRIERS[i]<<" ";
    cout<<endl<<"Conj Carriers:"<<endl;
    for( uint i=0; i<CARRIER_COUNT; ++i )
        cout<<CONJ_CARRIERS[i]<<" ";
    
//    vector<double> demodData( 64, 0 );
//    cout<<"input data:";
//    for( uint i=0; i<demodData.size(); ++i ) {
//        demodData[i]= i+1;
//        cout<<demodData[i]<<endl;
//    }
//    
//    pEngine->Demodulate( &demodData );
    
    ////////////////////////////////////////////////////////////////////////////////
    //
    //      OFDM Transmitter
    //
    ////////////////////////////////////////////////////////////////////////////////

    // Generate header
    vector<double> header= pEngine->GenerateHeader();
    // Frame guard of 0s
    vector<double> frameGuard( SYMB_PERIOD, 0 );
    // Time wave tx will store all data to be transmitted
    vector<double> timeWaveTx;
    
    double framePower= 0;
    
    long lModulatedData= 0;
    while( lModulatedData < iDataLength ) {
        //long lFrameLen= 16;
        int iFrameLen= min( static_cast<int>(SYMB_PER_FRAME*CARRIER_COUNT), iDataLength );
        vector<double> timeSignalTx= pEngine->Modulate( &pData[lModulatedData], iFrameLen );
        
        // Append timeSignalTx and frameGuard to timeWaveTx
        timeWaveTx.reserve( timeSignalTx.size()+frameGuard.size() );
        timeWaveTx.insert( timeWaveTx.end(), timeSignalTx.begin(), timeSignalTx.end() );
        timeWaveTx.insert( timeWaveTx.end(), frameGuard.begin(), frameGuard.end() );
        
        // Calculate frame power
        framePower= variance( timeWaveTx );
        
        lModulatedData+= iFrameLen;
    }
    
    // Multiply header/trailer by power
    for( uint i=0; i<header.size(); ++i )
        header[i]*= framePower;
    
    // Add header/trailer to beginning and end of timeWaveTx
    timeWaveTx.reserve( header.size()*2 );
    timeWaveTx.insert( timeWaveTx.begin(), header.begin(), header.end() );
    timeWaveTx.insert( timeWaveTx.end(), header.begin(), header.end() );
    
    
    ////////////////////////////////////////////////////////////////////////////////
    //
    //      OFDM Receiver
    //
    ////////////////////////////////////////////////////////////////////////////////
    
    vector<double> timeWave;
    vector<double> timeWaveRx( timeWaveTx );
    
    uint uUnpad= 0,
         uStartX= 0,
         uEndX= (uint)timeWaveRx.size();
    bool bLastFrame= false;
    
    if( iDataLength % CARRIER_COUNT != 0 )
        uUnpad= CARRIER_COUNT - (iDataLength % CARRIER_COUNT);
    
    uint uNumFrame= ceil( iDataLength * (WORD_SIZE / static_cast<float>(SYMB_SIZE)) / (SYMB_PER_FRAME*CARRIER_COUNT) );
    
    vector<double>::iterator timeWaveRxIterator= timeWaveRx.begin();
    
    for( uint i=0; i<uNumFrame; ++i ) {
        uint uReserveSize= 0;
        if( i == 1 )
            uReserveSize= min( uEndX, (HEAD_LEN+SYMB_PERIOD * ( (SYMB_PER_FRAME + 1) / 2 + 1) ) );
        else
            uReserveSize= min( uEndX, ( (uStartX-1) + ( SYMB_PERIOD * ((SYMB_PER_FRAME + 1) / 2 + 1)) ) );
        
        // Populate timeWave with uReserveSize elements from timeWaveTx
        timeWave.reserve(uReserveSize);
        timeWave.insert(timeWave.end(), timeWaveRxIterator, timeWaveRxIterator+uReserveSize );
        //timeWaveRxIterator+= uReserveSize;
        
        // Detect the data frame that only contains the useful information
        // NOTE: We need to use OFDMFrameDetect() for this, but for now I
        //       know where the frame begins because we have a perfect
        //       simulated communication channel
        int iFrameStart= distance( timeWaveRx.begin(), timeWaveRxIterator );
        int iFrameEnd;
        if( i==uNumFrame ) {
            bLastFrame= true;
            iFrameEnd= min( (double)uEndX, static_cast<double>( (iFrameStart-1) + SYMB_PERIOD*(1+ceil(remainder(1, CARRIER_COUNT*SYMB_PER_FRAME)/(double)CARRIER_COUNT)) ) );
        } else
            iFrameEnd= min( iFrameStart-1+(SYMB_PER_FRAME+1)*SYMB_PERIOD, uEndX );
        
        // MATLAB: time_wave = time_wave_rx(frame_start:frame_end);
        timeWave.resize( iFrameEnd-iFrameStart );
        timeWaveRx.insert( timeWaveRx.begin(), timeWaveRx.begin()+iFrameStart, timeWaveRx.begin()+iFrameEnd );
        
        uStartX= iFrameEnd-SYMB_PERIOD;
        
        // Demodulate the received time signal
        
    }
    
    return 0;
}

// Function reads all data from a file into a buffer
unsigned char* readDataFromFile( const char* filename ) {
    FILE *file;
    
	// Open file
	file = fopen(filename, "r");
	if (!file)
	{
		fprintf(stderr, "Unable to open file %s", filename);
		return false;
	}
	
	//Get file length
	fseek(file, 0, SEEK_END);
	iDataLength= ftell(file);
	fseek(file, 0, SEEK_SET);
    
	//Allocate memory
	unsigned char* data=(unsigned char *)malloc(iDataLength+1);
	if (!data)
	{
		fprintf(stderr, "Memory error!");
        fclose(file);
		return nullptr;
	}
    
	//Read file contents into buffer
	fread(data, iDataLength, 1, file);
	fclose(file);
    
    return data;
    
} // end Transmitter::ReadDataFromFile()


bool writeDataToFile( unsigned char* data, const char *filename ) {
    FILE* file;
    
    file= fopen(filename, "wb");
    if (!file)
	{
		fprintf(stderr, "Unable to open file %s", filename);
		return false;
	}
    
    fwrite(data, sizeof(unsigned char), iDataLength, file);
    fclose(file);
    
    return true;
}

// Function calculates the variance of a 1d vector of doubles
double variance( vector<double> &in ) {
    double total= 0;
    for( uint i=0; i<in.size(); ++i ) {
        total+= in[i];
    }
    double avg= total / static_cast<double>(in.size());
    
    double var= 0;
    for( uint i=0; i<in.size(); ++i ) {
        var+= pow( avg-in[i], 2 );
    }
    return var / static_cast<double>(in.size());
}
