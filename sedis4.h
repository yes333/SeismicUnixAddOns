#ifndef SEDIS4_H
#define SEDIS4_H

//=============================================================================
//
// header file for SEDIS IV (Rev. 1.05.01) output data
// written by: Xuelin Qiu, 2006/09/24
// modified by Sanyu Ye, 2012/07
//
//=============================================================================


typedef struct
 {
    char SEC;
    char MIN;
    char HOUR;
    char DAY;		// Day of week. Not used.
    char MDAY;
    char MONTH;
    short YEAR;		// 2 bytes
 } TimeDate;            // 8 bytes

typedef struct {
    char DataHeader[12];	// Data ID "SeismicData" 0 
    char HeaderSize;		// size of header 80 bytes
    char ConfigWord;		// see below
    char ChannelBitMap;		// see below
    //unsigned short BlockSamples;	// Samples per one minute Block
    unsigned char BlockSamples1;	// lower byte. Samples per one minute Block
    unsigned char BlockSamples2;	// upper byte. Samples per one minute Block
    char SampleBytes;		// Bytes per one sample
    TimeDate SampleTime;	// time of first sample in the block
    unsigned short Bat;		// Battery voltage = 50/1024*Bat in Volts
    unsigned short Temp;	// temperature = ( 5000/1024*Temp - 600)/10 °C    //30 bytes
    char Rev;			// software revision number 2 means 1.05
    char Board;			// Serial number of ADC board.
    char Reserved1;
  
    unsigned char NumberSV;	// number of usable satellites in Drift measurement
    short Drift;		// Drift / 32,768 = Drift correction in msec.
    TimeDate SedisTime;
    TimeDate GPSTime;
    char Reserved2;

    char Compass;		// Switcher. If 0xff then data below not valid, if 0 then position or if 1 then compass ASCII string.
    char PositionStr[23];	// ASCII string compass or position see below
    char Reserved3[3];
} Sedis_Header;			// size 80 bytes

#endif