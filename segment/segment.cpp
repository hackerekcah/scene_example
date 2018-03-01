#include "segment.h"

void initialTheta()
{
	FILE *fp;
	fp=fopen("c:/testing/gamma/info.dat", "r");

	for(int n=0; n<NUMBER_CHANNEL; n++)
		fscanf(fp, "%f\n", &Theta[n]);

	fclose(fp);
}

int ReadInput(char *filename)
{
	FILE *fp;
	int sigLength=0;

	if ((fp = fopen(filename, "r")) == NULL){
		printf("Cannot open input file!\n");
		exit(0);
	}

	while (!feof(fp))
	{
		float f;
		fscanf(fp, "%f\n", &f);
		Input[sigLength]=f;
		sigLength++;
	}
	fclose(fp);

	return(sigLength);
}

void AudiPeriph()
{	
	float lowerERB, upperERB, spaceERB;
	float cf;
	int chan;

	lowerERB = HzToERBRate(MINCF);
	upperERB = HzToERBRate(MAXCF);
  
	spaceERB = (NUMBER_CHANNEL > 1) ? (upperERB-lowerERB)/(NUMBER_CHANNEL-1) : 0;
	
	for (chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		cf = ERBRateToHz(lowerERB + chan*spaceERB);
		fChan[chan].cf = cf;
		fChan[chan].bw = 24.7*(cf*0.00437 + 1.0) * BW_CORRECTION;
		fChan[chan].delay = 1.5/PI/fChan[chan].bw*float(SAMPLING_FREQUENCY);		
	}
}

float HzToERBRate(float Hz)
{
	return( 21.4*log10(Hz*0.00437 + 1.0) );
}

float ERBRateToHz(float ERBRate)
{
	return( (pow(10, ERBRate/21.4) - 1) / 0.00437 );
}

void gammaToneFilter(float *input, float *output, gammaTone fChan, int sigLength)
{
	float dt, twoPiT, gain, z;
	float f1, f2;
	float x[4], y[4];
	int i;
	int n;
	
	dt= 1/float(SAMPLING_FREQUENCY);

	twoPiT=2 * PI * dt;
	//gain = fChan.midEarCoeff * pow(twoPiT * fChan.bw, 4) / 3.0;
	gain = pow(twoPiT * fChan.bw, 4) / 3.0;
	z = exp(-twoPiT * fChan.bw);

	f1 = cos(fChan.cf * twoPiT) * z;
	f2 = sin(fChan.cf * twoPiT) * z;

	for (i=0; i<4; i++)
	{
		fChan.p[i] = 0;
		fChan.q[i] = 0;
	}

	for (n=0; n<sigLength; n++)
	{
		output[n] = fChan.p[3] * gain;
		for (i=0; i<4; i++)
		{
			x[i] = f1*fChan.p[i] - f2*fChan.q[i];
			y[i] = f2*fChan.p[i] + f1*fChan.q[i];
		}

		fChan.p[0] = input[n] * f1 + x[0];
		fChan.q[0] = input[n] * f2 + y[0];
		
		fChan.p[1] = fChan.p[0] + x[1];
		fChan.q[1] = fChan.q[0] + y[1];
		
		fChan.p[2] = fChan.p[1] + x[1] + x[2];
		fChan.q[2] = fChan.q[1] + y[1] + y[2];
		
		fChan.p[3] = fChan.p[2] + x[1] + 2*x[2] + x[3];
		fChan.q[3] = fChan.q[2] + y[1] + 2*y[2] + y[3];
	}
}

void ExtractEv()
{
	float *gOut;

	float *filter, beta;
	int chan, fLength, m; 
	int n, tim;
	
	gOut = new float[SigLength];

	kaiserPara(Ripple, float(Stopband - Passband) / SAMPLING_FREQUENCY, fLength, beta);
	filter = new float[fLength+1];
	kaiserLowPass(filter, fLength, beta, float(Passband + Stopband) / SAMPLING_FREQUENCY);

	int shift=fLength/2;
	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		gammaToneFilter(Input, gOut, fChan[chan], SigLength);

		for(n=0; n<SigLength/DOWNRATE; n++)
		{
			InitIntensity[chan][n] = 0;
		
			for(m=0; m<=fLength; m++)
			{
				tim = n*DOWNRATE + shift - m;
				if ( (tim >= 0) && (tim < SigLength) ) 
					InitIntensity[chan][n] += fabs(gOut[tim]) * filter[m];
			}

			InitIntensity[chan][n] = 10*log10(fabs(InitIntensity[chan][n])+1);
		}
	}
	
	delete [] filter;
	delete [] gOut;
}

void Smooth(float tScale)
{
	float *filter, beta;
	int chan, fLength, m; 
	int n, tim;

	Passband = tScale;
	Stopband = tScale+10;
	Ripple = 0.01;

	kaiserPara(Ripple, float(Stopband - Passband) / (SAMPLING_FREQUENCY/DOWNRATE), fLength, beta);
	filter = new float[fLength+1];
	kaiserLowPass(filter, fLength, beta, float(Passband + Stopband) / (SAMPLING_FREQUENCY/DOWNRATE) );

	int shift=fLength/2;
	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		for(n=0; n<SigLength; n++)
		{
			TSmoothIntensity[chan][n] = 0;
		
			for(m=0; m<=fLength; m++)
			{
				tim = n+shift-m;
				if ( (tim >= 0) && (tim < SigLength) ) 
					TSmoothIntensity[chan][n] += InitIntensity[chan][tim] * filter[m];
			}
		}
	}
	
	delete [] filter;
}

void kaiserPara(float delta, float transBw, int &fLength, float &beta)
{
	float a, len;
	
	a= -20 * log10(delta);

	if (a <= 21) beta = 0;
	else if (a<= 50) beta = 0.5842 * pow(1.0*(a-21), 0.4) + 0.07889 * (a-21);
	else beta = 0.1102 * (a - 8.7);

	len = (a - 7.95) / 14.36 / transBw;
	fLength = int(len);
	if ((len - fLength) < 0.5) fLength++;
	else fLength+=2;

	if (fLength%2 != 0) fLength++;
}
	
void kaiserLowPass(float *filter, int fLength, float beta, float wn)
{
	int tim, step;
	float k, sum;

	for (tim=0; tim<=fLength; tim++)
	{
		k = 2*tim/float(fLength) - 1;
		filter[tim] = bessi0( beta*sqrt(1- k*k)) / bessi0( beta );
	}

	sum=0;
	for (tim=0; tim<=fLength; tim++)
	{
		step = tim - fLength/2;
		if (step !=0) filter[tim] *= sin(wn * PI * step) / PI / step;
		else filter[tim] *= wn;

		sum += filter[tim];
	}
}

void kaiserBandPass(float *filter, int fLength, float beta, float wc, float wn)
{
	int tim;

	kaiserLowPass(filter, fLength, beta, wn);

	for (tim=0; tim<=fLength; tim++)
		filter[tim] *= 2 * cos((tim-fLength/2) * wc * PI);
}

float bessi0(float x)
{
	float ax,ans;
	float y;

	ax = fabs(x);
	if (ax < 3.75)
	{
		y = x/3.75;
		y *= y;
		ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y*0.45813e-2)))));
	}
	else
	{
		y = 3.75 / ax;
		ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
	}
	return ans;
}

void freqSmooth(int freqScale)
{
	float c1=1/sqrt(2*PI)/freqScale;
	float c2=0.5/freqScale/freqScale;
	
	float *h;
	h=new float[8*freqScale+1];

	for(int n=0; n<=8*freqScale; n++)
		h[n]=c1*exp(-pow(n-4*freqScale, 2)*c2);

	for(int f=0; f<SigLength; f++)
	{
		for(int chan=0; chan<NUMBER_CHANNEL; chan++)
		{
			Ev[f][chan]=0;
			diffEv[f][chan]=0;
			for(int n=0; n<=8*freqScale; n++)
			{
				int p=float(chan+4*freqScale-n);
				if( (p>=0) && (p<NUMBER_CHANNEL) )
				{
					Ev[f][chan] += h[n]*TSmoothIntensity[p][f];

					if(f<(SigLength-1)) diffEv[f][chan] += h[n]*(TSmoothIntensity[p][f+1]-TSmoothIntensity[p][f]);
					else diffEv[f][chan]=0;
				}
			}
		}
	}

	delete h;
}

void getOFSet(int chan, float theta1, float theta2)
{
	int n, numOn=0, numOff=0;
			
	for(n=1; n<(SigLength-2); n++)
	{
		if ( (diffEv[n][chan]>diffEv[n-1][chan]) && (diffEv[n][chan]>=diffEv[n+1][chan]) && (diffEv[n][chan]>theta1))
		{
			onS[numOn][chan].pos=n;
			if (numOn>0) onS[numOn-1][chan].next=n;

			numOn++;
		}

		if ( (diffEv[n][chan]<=diffEv[n-1][chan]) && (diffEv[n][chan]<diffEv[n+1][chan]) && (diffEv[n][chan]<theta2))
		{
			offS[numOff][chan].pos=n;
			offS[numOff][chan].weight=diffEv[n][chan];

			if (numOff>0) offS[numOff-1][chan].next=n;

			numOff++;
		}
	}

	numChanOnset[chan]=numOn;
	if (numOn>0) onS[numOn-1][chan].next=SigLength-2;

	numChanOffset[chan]=numOff;
	if (numOff>0) offS[numOff-1][chan].next=SigLength-2;
}			               

void edgeDetect(float paraEdge)
{
	int chan, n;

	float mean=0;
	float std=0;

	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		for(n=40; n<(SigLength-40); n++)
		{
			float diff = diffEv[n][chan];
			mean += diff;
			std += pow(diff, 2);
		}
	}

	mean /= float((SigLength-80)*NUMBER_CHANNEL);
	std /= float((SigLength-80)*NUMBER_CHANNEL);

	std = sqrt(std - pow(mean, 2));

	for(chan=0; chan<NUMBER_CHANNEL; chan++)
		getOFSet(chan, mean + paraEdge*std, 100);
}
    
void precisePosition(int dis)
{
	int n=0, chan;
	int mark[NUMBER_CHANNEL];

	while (n<numSegment)
	{
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			mark[chan]=0;

		for(chan=eSeg[n].sChan; chan<=eSeg[n].eChan; chan++)
		{
			int index=-1, mDis=dis, m;

			for(m=0; m<numChanOnset[chan]; m++)
			{
				if ( abs(onS[m][chan].pos-eSeg[n].onTime[chan]) < mDis ){
					mDis=abs(onS[m][chan].pos-eSeg[n].onTime[chan]);
					index=m;
				}
			}
			
			if (index>=0){ eSeg[n].onTime[chan]=onS[index][chan].pos; mark[chan]=1;	}
		
			index=-1;
			mDis=dis;
		
			for(m=0; m<numChanOffset[chan]; m++)
			{
				if ( (abs(offS[m][chan].pos-eSeg[n].offTime[chan]) < mDis) && (offS[m][chan].pos>eSeg[n].onTime[chan]) ){
					mDis=abs(offS[m][chan].pos-eSeg[n].offTime[chan]);
					index=m;
				}
			}

			if (index>=0) eSeg[n].offTime[chan]=offS[index][chan].pos;
		}
		
/*		while ( (eSeg[n].sChan<eSeg[n].eChan) && (mark[eSeg[n].sChan]==0) ){ eSeg[n].sChan++; }
		
		while ( (eSeg[n].eChan>eSeg[n].sChan) && (mark[eSeg[n].eChan]==0) ){ eSeg[n].eChan--; }

		if( (eSeg[n].sChan+3)>=eSeg[n].eChan)
		{
			numSegment--; 
			for(int k=n; k<numSegment; k++)
				eSeg[k]=eSeg[k+1];
		}

		else */n++;
	}
}

int findMatchOnset(int chan1, int chan2, int pos, int dis)
{	
	int mDis=dis, index=-1, m;

	for(m=0; m<numChanOnset[chan1]; m++)
		if ( abs(onS[m][chan1].pos-pos) < mDis )
		{
			mDis=abs(onS[m][chan1].pos-pos);
			index=m;
		}

	if (index>=0)
	{
		for(m=0; m<numChanOnset[chan2]; m++)
			if ( abs(onS[m][chan2].pos-onS[index][chan1].pos) < mDis )
			{
				index = -1;
				break;
			}
	}
	
	return(index);
}

int findMatchOffset(int chan1, int chan2, int pos, int dis)
{	
	int mDis=dis, index=-1, m;

	for(m=0; m<numChanOffset[chan1]; m++){
		if ( abs(offS[m][chan1].pos-pos) < mDis ){
			mDis=abs(offS[m][chan1].pos-pos);
			index=m;
		}
	}

	if (index>=0)
	{
		for(m=0; m<numChanOffset[chan2]; m++)
			if ( abs(offS[m][chan2].pos-offS[index][chan1].pos) < mDis )
			{
				index = -1;
				break;
			}
	}
	
	return(index);
}

/*float crossChannel(int chan1, int chan2, int sp, int ep)
{
	float m1=0;
	float m2=0;

	float f=0;

	for(int n=sp; n<=ep; n++)
	{
		m1+=Input[n][chan1];
		m2+=Input[n][chan2];
	}

	m1/=float(ep-sp+1);
	m2/=float(ep-sp+1);

	for(n=sp; n<=ep; n++)
	{
		f+=pow(Input[n][chan1]-m1-Input[n][chan2]+m2, 2);
	}

	f/=float(ep-sp+1);

	return(sqrt(f));
}*/

float crossChannel(int chan1, int chan2, int sp, int ep)
{
	int n;
	//sp=sp+1;
	if (ep<sp) return(1);

	float *f1, *f2;
	f1=new float[ep-sp+1];
	f2=new float[ep-sp+1];

	for(n=0; n<=(ep-sp); n++){
		f1[n]=Ev[n+sp][chan1];
		f2[n]=Ev[n+sp][chan2];
		//f1[n]=Ev[n+sp][chan1]-Ev[n+sp-1][chan1];
		//f2[n]=Ev[n+sp][chan2]-Ev[n+sp-1][chan2];
	}

	float m1=0, m2=0;
	for(n=0; n<=(ep-sp); n++){
		m1+=f1[n]; m2+=f2[n];
	}

	m1/=float(ep-sp+1);
	m2/=float(ep-sp+1);

	for(n=0; n<=(ep-sp); n++){
		f1[n]-=m1;	f2[n]-=m2;
	}

	float s1=0, s2=0, s=0;
	for(n=0; n<=(ep-sp); n++){	
		s+=f1[n]*f2[n];
		s1+=f1[n]*f1[n];
		s2+=f2[n]*f2[n];
	}

	s/=sqrt(s1*s2+1e-300);

	delete f1;
	delete f2;

	return(-s);
}

int onSetExpand(int chan1, int chan2, int n, int dis, float theta)
{
	int sp=eSeg[n].onTime[chan2];
	int ep=eSeg[n].offTime[chan2];

	int index=-1;
	if (crossChannel(chan1, chan2, sp, ep)<theta)
	{
		index=findMatchOnset(chan1, chan2, sp, dis);

		if(index>=0)
		{
			eSeg[n].onTime[chan1]=onS[index][chan1].pos;
					
			index=-1;
			int mDis=dis;
			for(int m=0; m<numChanOffset[chan1]; m++)
			{
				if ( (offS[m][chan1].pos>eSeg[n].onTime[chan1]) && (abs(offS[m][chan1].pos-eSeg[n].offTime[chan2])<mDis) )
				{
					mDis=abs(offS[m][chan1].pos-eSeg[n].offTime[chan2]);
					index=m;
				}
			}

			if(index>=0) eSeg[n].offTime[chan1]=offS[index][chan1].pos;
		}
	}
	
	return(index);
}
	
void mergeSegment(int n, int m)
{
	int chan, sChan, eChan;
	
	if (eSeg[n].sChan>eSeg[m].sChan)
	{
		sChan=eSeg[n].sChan;

		for(chan=eSeg[m].sChan; chan<eSeg[n].sChan; chan++)
		{
			eSeg[n].onTime[chan]=eSeg[m].onTime[chan];
			eSeg[n].offTime[chan]=eSeg[m].offTime[chan];
		}
		eSeg[n].sChan=eSeg[m].sChan;
	}
	else sChan=eSeg[m].sChan;
	
	if (eSeg[n].eChan<eSeg[m].eChan)
	{
		eChan=eSeg[n].eChan;

		for(chan=eSeg[n].eChan+1; chan<=eSeg[m].eChan; chan++)
		{
			eSeg[n].onTime[chan]=eSeg[m].onTime[chan];
			eSeg[n].offTime[chan]=eSeg[m].offTime[chan];
		}
		eSeg[n].eChan=eSeg[m].eChan;
	}
	else eChan=eSeg[m].eChan;

	for(chan=sChan; chan<=eChan; chan++)
	{
		if (eSeg[n].offTime[chan]<eSeg[m].offTime[chan]) eSeg[n].offTime[chan]=eSeg[m].offTime[chan];
	}

	numSegment--;
	for(int k=m; k<numSegment; k++)
		eSeg[k]=eSeg[k+1];
}

void freqExpand(int dis, float theta)
{
	int n, m, chan;

	for(n=0; n<numSegment; n++)
	{
		int index=0;
		while ( (eSeg[n].sChan>0) && (index>=0) )
		{
			index=onSetExpand(eSeg[n].sChan-1, eSeg[n].sChan, n, dis, theta);
			if (index>=0) eSeg[n].sChan--;
		}

		index=0;
		while ( (eSeg[n].eChan<(NUMBER_CHANNEL-1)) && (index>=0) )
		{
			index=onSetExpand(eSeg[n].eChan+1, eSeg[n].eChan, n, dis, theta);
			if (index>=0) eSeg[n].eChan++;
		}
	}	           
     
	n=0;
	while (n<(numSegment-1))
	{
		m=n+1;

	    while (m<numSegment)
		{
			int sChan = (eSeg[n].sChan>eSeg[m].sChan) ? eSeg[n].sChan:eSeg[m].sChan;
			int eChan = (eSeg[n].eChan<eSeg[m].eChan) ? eSeg[n].eChan:eSeg[m].eChan;

			for(chan=sChan; chan<=eChan; chan++)
			{
				if (eSeg[n].onTime[chan]==eSeg[m].onTime[chan])
				{
					mergeSegment(n, m);
					m=n;
					break;
				}
			}		    
			m++;
		}
		n++;
	}
}

void removeOnset()
{
	int n, chan;
	for(chan=0; chan<NUMBER_CHANNEL; chan++)
		for(int frame=0; frame<SigLength; frame++)
			Mask[frame][chan]=0;

	for(n=0; n<numSegment; n++)
		for(chan=eSeg[n].sChan; chan<=eSeg[n].eChan; chan++)
		{
			int sp=eSeg[n].onTime[chan];
			int ep=eSeg[n].offTime[chan];

			for(int frame=sp; frame<=ep; frame++)
				Mask[frame][chan]=1;
		}

	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		n=0;
		while (n<numChanOnset[chan])
		{
			int pos=onS[n][chan].pos;
			if (Mask[pos][chan]==1)
			{
				numChanOnset[chan]--;

				for(int k=n; k<numChanOnset[chan]; k++)
					onS[k][chan]=onS[k+1][chan];
			}
			else n++;
		}
	}
}

void chanOFsetMatch()
{
	int n, chan, m;

	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		m=0;
		for(n=0; n<numChanOnset[chan]; n++)
		{
			while( (m<numChanOffset[chan]) && (offS[m][chan].pos<onS[n][chan].pos) ){ m++; }
			int sIn=m;

			while( (m<numChanOffset[chan]) && (offS[m][chan].pos<onS[n][chan].next) ){ m++; }
			int eIn=m-1;

			if ((eIn<sIn) || (sIn==numChanOffset[chan]) ) onS[n][chan].fpos=onS[n][chan].next-1;
			
			else
			{
				onS[n][chan].fpos=offS[sIn][chan].pos;
				float weight=offS[sIn][chan].weight;

				for(m=sIn+1; m<=eIn; m++)
				{
					if (offS[m][chan].weight<weight) 
					{
						onS[n][chan].fpos=offS[m][chan].pos;
						weight=offS[m][chan].weight;
					}
				}
			}
		}
	}			
}

void onsetFront(float theta, int dis)
{
	int n, chan, mark[2][MAX_SEGMENT];
	numOnFront=0;
	
	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
	    if(chan>0)
		{
			for(n=0; n<numChanOnset[chan-1]; n++)
				mark[0][n]=mark[1][n];
		}

		for(n=0; n<numChanOnset[chan]; n++)
		{
			int sp=onS[n][chan].pos;
			int ep=onS[n][chan].fpos;
			int m;

			int index=-1;
			if( (chan>0) && (crossChannel(chan-1, chan, sp, ep)<theta) )
				index=findMatchOnset(chan-1, chan, sp, dis);

			if (index>=0) m=mark[0][index];
			else{ m=numOnFront;	onFront[m].sChan=chan; numOnFront++; }

			mark[1][n]=m; onFront[m].eChan=chan;
			onFront[m].onTime[chan]=sp; onFront[m].offTime[chan]=ep;
		}    
	}

	n=0;
	while (n<numOnFront)
	{
		if(onFront[n].eChan<=(onFront[n].sChan+3))
		{
			numOnFront--;
			for(int k=n; k<numOnFront; k++)
				onFront[k]=onFront[k+1];
		}
		else n++;
	}
}

void offsetFront(int dis)
{
	int n, chan, mark[2][MAX_SEGMENT];
	numOffFront=0;
	
	for(chan=0; chan<NUMBER_CHANNEL; chan++)
	{
	    if(chan>0)
		{
			for(n=0; n<numChanOffset[chan-1]; n++)
				mark[0][n]=mark[1][n];
		}

		for(n=0; n<numChanOffset[chan]; n++)
		{
			int ep=offS[n][chan].pos;
			int m;

			int index=-1;
			if(chan>0) index=findMatchOffset(chan-1, chan, ep, dis);

			if (index>=0) m=mark[0][index];
			else{ 
				m=numOffFront; offFront[m].sChan=chan; numOffFront++; }

			mark[1][n]=m; offFront[m].eChan=chan;
			offFront[m].offTime[chan]=ep;
		}    
	}

	n=0;
	while(n<numOffFront)
	{
		if(offFront[n].eChan<=offFront[n].sChan)
		{
			numOffFront--;
			for(int k=n; k<numOffFront; k++)
				offFront[k]=offFront[k+1];
		}
		else n++;
	}
}

void linearFit(int *mark, int n)
{
	int chan;
	int sp=eSeg[n].sChan, ep=eSeg[n].eChan;

	while ((sp<ep) && (mark[sp]==1)) { sp++; }
	while ((sp<ep) && (mark[ep]==1)) { ep--; }

	if(sp==ep)
	{
		for(chan=eSeg[n].sChan; chan<=eSeg[n].eChan; chan++)
			eSeg[n].offTime[chan]=eSeg[n].offTime[sp];
		return;
	}

	chan=sp; 
	int state=0;
	int p1, p2;
	while (chan<ep)
	{
		if( (mark[chan]==0) && (state==0) ) p1=chan;
		else if( (mark[chan]==0) && (state==1) )
		{
			state=0;
			p2=chan;
			float d=float(eSeg[n].offTime[p2]-eSeg[n].offTime[p1])/float(p2-p1);

			for(chan=p1+1; chan<=p2-1; chan++)
				eSeg[n].offTime[chan]=eSeg[n].offTime[p1]+int(float(chan-p1)*d);
			
			p1=p2;
			chan=p2;
		}
		else if(mark[chan]==1){ state=1; }

		chan++;
	}

	for(chan=sp-1; chan>=eSeg[n].sChan; chan--)
		eSeg[n].offTime[chan]=eSeg[n].offTime[sp];

	for(chan=ep+1; chan<=eSeg[n].eChan; chan++)
		eSeg[n].offTime[chan]=eSeg[n].offTime[ep];
}

void matchOFFront()
{
	int mark[NUMBER_CHANNEL], chan, n;

	for(n=numSegment; n<(numSegment+numOnFront); n++)
	{
		eSeg[n]=onFront[n-numSegment];
	
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			mark[chan]=0;

		for(chan=eSeg[n].sChan; chan<=eSeg[n].eChan; chan++)
			mark[chan]=1;

		int unMatchChan=eSeg[n].eChan-eSeg[n].sChan+1;
		while (unMatchChan>0)
		{
			int overlap=0;
			int index=-1;
			for(int m=0; m<numOffFront; m++)
			{
				int sChan = (eSeg[n].sChan>offFront[m].sChan) ? eSeg[n].sChan:offFront[m].sChan;
				int eChan = (eSeg[n].eChan<offFront[m].eChan) ? eSeg[n].eChan:offFront[m].eChan;

				int count=0;
				for(chan=sChan; chan<=eChan; chan++)
				{
					if( (eSeg[n].offTime[chan]==offFront[m].offTime[chan]) && (mark[chan]>0) )count++;
				}
				
				if(count>overlap){ overlap=count, index=m; }
			}
			
			if(index>=0)
			{
				int sChan = (eSeg[n].sChan>offFront[index].sChan) ? eSeg[n].sChan:offFront[index].sChan;
				int eChan = (eSeg[n].eChan<offFront[index].eChan) ? eSeg[n].eChan:offFront[index].eChan;

				for(int chan=sChan; chan<=eChan; chan++)
					if (mark[chan]>0){
						eSeg[n].offTime[chan]=offFront[index].offTime[chan];
						mark[chan]=0;
						unMatchChan--;
					}
			}
			else
			{
				linearFit(mark, n);
				unMatchChan=0;
			}
		}
	}
	
	numSegment += numOnFront;
}

void segmentGenerate(float freqScale, int type, float paraEdge, float paraCChan, int dis)
{
	freqSmooth(freqScale);

	edgeDetect(paraEdge);

	precisePosition(dis);

	freqExpand(dis, paraCChan);

	if(type==1) 
	{
		removeOnset();

		chanOFsetMatch();

		onsetFront(paraCChan, dis);
	
		offsetFront(dis);

		matchOFFront();
	}
}

int main(int argc, char *argv[])
{	
	// Initialize data
	if (argc != 3) 
	{
		printf("Ussage: Segment inputFile outputFile\n");
		exit(0);
	}

	SigLength = ReadInput(argv[1]);
	printf("%d samples\n", SigLength);

	// Auditory Periphery
	int chan;
	for (chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		InitIntensity[chan] = new float[SigLength/DOWNRATE];
		TSmoothIntensity[chan] = new float[SigLength/DOWNRATE];
	}
			
	Passband = 30;	Stopband = 60;	Ripple = 0.01;
	AudiPeriph();
	ExtractEv();

	SigLength /= DOWNRATE;

	// Initialize parameters

	float timeScale[3], freqScale[3], typeScale[3], paraEdge[3], paraCChan[3];

	timeScale[0]=4; timeScale[1]=14; timeScale[2]=14; 
	typeScale[0]=1; typeScale[1]=1; typeScale[2]=0;
	freqScale[0]=6; freqScale[1]=6; freqScale[2]=0.5;
	paraEdge[0]=1; paraEdge[1]=1; paraEdge[2]=1;
	paraCChan[0]=-0.95; paraCChan[1]=-0.95; paraCChan[2]=-0.85;


	// Smooth and generate Segments
	int frame, scale;
	//initialTheta();

	numSegment=0;
	int oldTScale=0;
	for(scale=0; scale<3; scale++)
	{
		Smooth(timeScale[scale]);
		
		segmentGenerate(freqScale[scale], typeScale[scale], paraEdge[scale], paraCChan[scale], 8);
	}
	
	SigLength /= 4;

	for(chan=0; chan<NUMBER_CHANNEL; chan++)
		for(frame=0; frame<SigLength; frame++)
			Mask[frame][chan]=0;

	for(int n=0; n<numSegment; n++)
		for(chan=eSeg[n].sChan; chan<=eSeg[n].eChan; chan++)
		{
			int sp=floor(float(eSeg[n].onTime[chan]-2)/4);
			int ep=ceil(float(eSeg[n].offTime[chan]+2)/4);

			if (sp<0) sp=0;
			if (ep>=(SigLength)) ep=SigLength-1;
			if (sp>0) sp--;
			if (ep<(SigLength-1)) ep++;

			for(frame=sp; frame<=ep; frame++)
			{
				if(Mask[frame][chan]==0) Mask[frame][chan]=n+1;
				else if(eSeg[Mask[frame][chan]-1].onTime[chan]<eSeg[n].onTime[chan])
					Mask[frame][chan]=n+1;
			}
		}

	FILE *fp;
	fp=fopen(argv[2], "w");
	for(frame=0; frame<SigLength; frame++)
	{
		for(chan=0; chan<NUMBER_CHANNEL; chan++)
			fprintf(fp, " %3d", Mask[frame][chan]);

		fprintf(fp, "\n");
	}	
	fclose(fp);
				
	for (chan=0; chan<NUMBER_CHANNEL; chan++)
	{
		delete [] InitIntensity[chan];
		delete [] TSmoothIntensity[chan];
	}
	return 0;

}
