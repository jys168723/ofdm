//
//  ofdm_params.h
//  OFDM
//
//  Created by Stephen MacKinnon on 1/16/13.
//  Copyright (c) 2013 Stephen MacKinnon. All rights reserved.
//
#pragma once

#ifndef OFDM_ofdm_params_h
#define OFDM_ofdm_params_h

#include <math.h>

typedef unsigned int uint;

static const uint IFFT_SIZE= 16;//1024;          // IFFT size must be a power of 2
static const uint CARRIER_COUNT= 2;//128;       // Carrier count must be < ifft_size/2 - 2
static const uint SYMB_SIZE= 2;             // Bits per symbol
static const uint WORD_SIZE= 8;             // Bits per word of source data
static const uint GUARD_TIME= IFFT_SIZE/4;  // Length of guard interval for each symbol period
static const uint SYMB_PER_FRAME=           // Number of symbols per carrier in each frame for transmission
    ceil(8192/CARRIER_COUNT);
static const uint SYMB_PERIOD=              // Length of one symbol period including guard time
    IFFT_SIZE + GUARD_TIME;
static const uint HEAD_LEN= SYMB_PERIOD*8;  // Length of header and trailer of the transmitted data
static const uint ENVELOPE_SIZE=            // Size of the envelope detector
    ceil(SYMB_PERIOD/256)+1;

// More elegant way of doing this?
inline static const uint calculateCarrierSpacing() {
    uint uCarrierSpacing= 0u;
    
    while( CARRIER_COUNT*uCarrierSpacing <= (IFFT_SIZE/2 - 2) ) ++uCarrierSpacing;
    --uCarrierSpacing;
    
    return uCarrierSpacing;
}

static const uint CARRIER_SPACING= calculateCarrierSpacing();  // Spacing for carriers distributed in IFFT bins

inline static const uint* calculateCarriers() {
    uint* pCarriers= (uint*)malloc(sizeof(uint)*CARRIER_COUNT);
    
    const uint midFreq= IFFT_SIZE/4.0f;
    const uint firstCarrier= midFreq-floor( ((float)CARRIER_COUNT-1.0f)*(float)CARRIER_SPACING/2.0f );
    //const float LAST_CARRIER= MID_FREQ+floor( ((float)CARRIER_COUNT-1.0f)*(float)CARRIER_SPACING/2.0f );
    
    for( uint i= 0; i<CARRIER_COUNT; ++i ) pCarriers[i]= firstCarrier + (CARRIER_SPACING*i);
    
    return pCarriers;
}

static const uint* CARRIERS= calculateCarriers();

inline static const uint* calculateConjCarriers() {
    //uint conjCarriers[CARRIER_COUNT];
    uint* pConjCarriers= (uint*)malloc(sizeof(uint)*CARRIER_COUNT);
    
    for( uint i= 0; i<CARRIER_COUNT; ++i ) pConjCarriers[i]= IFFT_SIZE-CARRIERS[i]+2.0f;
    
    return pConjCarriers;
}

static const uint* CONJ_CARRIERS= calculateConjCarriers();

#endif
