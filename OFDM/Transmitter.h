//
//  Transmitter.h
//  OFDM
//
//  Created by Stephen MacKinnon on 1/16/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//

#ifndef __OFDM__Transmitter__
#define __OFDM__Transmitter__

#include <fstream>
#include "ofdm_params.h"

class Transmitter {
public:
    Transmitter();
    
    bool ReadDataFromFile( const char* filename );
    bool WriteDataToFile( const char* filename );
    void GenerateHeader( );
    void Modulate( char* pData, unsigned int uDataLen, unsigned int uIfftSize, unsigned int uCarriers, unsigned int uConjCarriers, unsigned int uCarrierCount, unsigned int uSymbSize, unsigned int uGuardTime );
    void TransmitData( );
    
private:
    unsigned char* m_pData;
    unsigned long m_lDataLength;
    
};

#endif /* defined(__OFDM__Transmitter__) */
