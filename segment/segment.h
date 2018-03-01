#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI (3.1415926535897932384626433832795)

#define SAMPLING_FREQUENCY		20000 
#define DOWNRATE				50			// Downsampling factor for intensity
#define MAX_SIG_LENGTH			200000

#define NUMBER_CHANNEL  128      		        /* maxmimum number of filters */
#define MINCF			50
#define MAXCF			8000

#define BW_CORRECTION       1.019      			/* ERB bandwidth correction 4th order */

struct gammaTone
{
	float cf, bw;
	float midEarCoeff;
	int delay;
	float p[4];
	float q[4];
};

// Global variables
int SigLength;
float Input[MAX_SIG_LENGTH];

float *InitIntensity[NUMBER_CHANNEL], *TSmoothIntensity[NUMBER_CHANNEL], *TFSmoothIntensity[NUMBER_CHANNEL], *DiffTFSmoothIntensity[NUMBER_CHANNEL];

float Passband, Stopband, Ripple; 

gammaTone fChan[NUMBER_CHANNEL];


int ReadInput(char *filename);

// Auditory periphery
void AudiPeriph();

float HzToERBRate(float Hz);

float ERBRateToHz(float ERBRate);
	
void gammaToneFilter(float *input, float *output, gammaTone fChan, int sigLength);


// Smoothing
void ExtractEv();

void Smooth(float tScale);


// Functions for generating a lowpass filter
void kaiserPara(float delta, float transBw, int &fLength, float &beta);
	
void kaiserLowPass(float *filter, int fLength, float beta, float wn);

float bessi0(float x);

#define MAX_FRAME		10000
#define MAX_SEGMENT		10000

int Mask[MAX_FRAME][NUMBER_CHANNEL];
float Ev[MAX_FRAME][NUMBER_CHANNEL];
float diffEv[MAX_FRAME][NUMBER_CHANNEL];

struct chanOnset
{
	int pos, fpos, next;
};

struct chanOffset
{
	int pos, next;
	float weight;
};

chanOnset onS[MAX_SEGMENT][NUMBER_CHANNEL];
chanOffset offS[MAX_SEGMENT][NUMBER_CHANNEL];

struct segment
{
	int onTime[NUMBER_CHANNEL];
	int offTime[NUMBER_CHANNEL];

	int sChan, eChan;
};

segment eSeg[MAX_SEGMENT], onFront[MAX_SEGMENT], offFront[MAX_SEGMENT];

int numChanOnset[NUMBER_CHANNEL];
int numChanOffset[NUMBER_CHANNEL];

int numSegment;
int numOnFront;
int numOffFront;

float Theta[NUMBER_CHANNEL];

void readEnvelope(char *filename);
void freqSmooth(int freqScale);
void edgeDetect();
void getOnset(float theta);
void getOffset(float theta);

