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

// Hacky for now, will be unnecessary when we establish a
// fixed buffer size
long lDataLength= 16u;

int main(int argc, const char * argv[])
{
    
    OFDMEngine* pEngine= new OFDMEngine();
    unsigned char data[16];
    cout<<"Raw data:"<<endl;
    for( uint i=0; i<16; ++i ) {
        data[i]= static_cast<unsigned char>(i+1);
        cout<<i<<": "<<static_cast<float>(data[i])<<endl;
    }
    unsigned char* pData= &data[0];
    
    long lModulatedData= 0;
    while( lModulatedData < lDataLength ) {
        long lFrameLen= 16;//min( lDataLength-lModulatedData, static_cast<long>(SYMB_PER_FRAME*CARRIER_COUNT) );
        pEngine->Modulate( &pData[lModulatedData], lFrameLen );
        lModulatedData+= lFrameLen;
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
	lDataLength= ftell(file);
	fseek(file, 0, SEEK_SET);
    
	//Allocate memory
	unsigned char* data=(unsigned char *)malloc(lDataLength+1);
	if (!data)
	{
		fprintf(stderr, "Memory error!");
        fclose(file);
		return nullptr;
	}
    
	//Read file contents into buffer
	fread(data, lDataLength, 1, file);
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
    
    fwrite(data, sizeof(unsigned char), lDataLength, file);
    fclose(file);
    
    return true;
}
