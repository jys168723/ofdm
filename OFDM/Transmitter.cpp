//
//  Transmitter.cpp
//  OFDM
//
//  Created by Stephen MacKinnon on 1/16/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//

#include "Transmitter.h"

using namespace std;

Transmitter::Transmitter() :
    m_pData() { }


// Function reads all data from a file into a buffer
bool Transmitter::ReadDataFromFile( const char* filename ) {
    FILE *file;
    
	//Open file
	file = fopen(filename, "r");
	if (!file)
	{
		fprintf(stderr, "Unable to open file %s", filename);
		return false;
	}
	
	//Get file length
	fseek(file, 0, SEEK_END);
	m_lDataLength=ftell(file);
	fseek(file, 0, SEEK_SET);
    
	//Allocate memory
	m_pData=(unsigned char *)malloc(m_lDataLength+1);
	if (!m_pData)
	{
		fprintf(stderr, "Memory error!");
        fclose(file);
		return false;
	}
    
	//Read file contents into buffer
	fread(m_pData, m_lDataLength, 1, file);
	fclose(file);
    
    return true;
    
} // end Transmitter::ReadDataFromFile()


bool Transmitter::WriteDataToFile( const char *filename ) {
    FILE* file;
    
    file= fopen(filename, "wb");
    if (!file)
	{
		fprintf(stderr, "Unable to open file %s", filename);
		return false;
	}
    
    fwrite(m_pData, sizeof(unsigned char), m_lDataLength, file);
    fclose(file);
    
    return true;
}


void Transmitter::GenerateHeader() {
    
}

void Transmitter::TransmitData() {
    
    // TODO: Generate header/trailer
    
    char frameGuard[SYMB_PERIOD];
    for( unsigned int i= 0; i<SYMB_PERIOD; ++i ) {
        frameGuard[i]= '0';
    }
    
    unsigned int symb_per_carrier= ceil( m_lDataLength/CARRIER_COUNT );
    
    if( symb_per_carrier > SYMB_PER_FRAME ) { // We have multiple frames
        
        
    }
}

void Transmitter::Modulate( char* pData, unsigned int uDataLen, unsigned int uIfftSize, unsigned int uCarriers, unsigned int uConjCarriers, unsigned int uCarrierCount, unsigned int uSymbSize, unsigned int uGuardTime ) {
//    unsigned int uCarrierSymbCount= ceil( static_cast<float>(uDataLen)/static_cast<float>(uCarrierCount) );
//    char cModulatedData[uDataLen];
//    for( unsigned int i= 0; i<uDataLen; ++i ) {
//        cModulatedData[i]= pData[i];
//    }
//    if( uDataLen/uCarrierCount != uCarrierSymbCount ) {
//        char cZeros
//    }
}
