function fea = extract_reverb_fea(audiofile)

% [sig, fs] = audioread('reverb_features/ReverbAudio/SuburbanGarage_reverb.wav');
[sig, fs] = audioread(audiofile);

% if two channel, extract the first channel
sig = sig(:,1);

% draw spectrogram
% m for mel scale as frequency axis unit
% w for output wav on top of spectrogram
% spgrambw (sig, fs ,'mw');
% soundsc(sig,fs);

%% STEP1: STFT, log spectrogram in dB
%ms
win_len = round(0.064 * fs);
%overlap
overlp = round(3/4 * win_len);

% By default, number of fft is greater than win_len, by some order of 2
% column of stft is half of size of fft + 1
stft = spectrogram(sig, win_len, overlp);

mag_spec = abs(stft);

log_spec = 10*log10(mag_spec + eps);


%% STEP2: apply melfilterbank

% melfilterbank fft size should be the same with stft fft size
nfft = 2*(size(stft,1) -1);

%generate triangular filterbank, 26 bands in [0,8k]hz range
% 26 filterbank
% 0-8k hz
% mod
%   h: specify fl/fh in hz
%   g: plot
melfb=melbankm(26,nfft,fs,0,8000,'hg');

% output of mel log spectrogram
s_mel = melfb * log_spec;

%% STEP3: find peaks in each freqeuncy band(channel)
for numChan=1:size(s_mel,1)
    % peak value > 10dB, peak distance > 5
    [pks{numChan},locs{numChan}] = findpeaks(s_mel(numChan,:),'MinPeakHeight',10, 'MinPeakDistance', 5);
end

%% STEP4: Calculate slopes for each peak in each channel

MAX_FRAME_INDEX = size(s_mel,2);

% for each channel
for numChan=1:size(s_mel,1)
    
    %vector of peaks location
    pkLocs = locs{numChan};
    
     slopes{numChan}=[];
    % for each peaks
    for pkIdx = 1:length(pkLocs)
        if isempty(pkLocs)
            slopes{numChan} = [];
        else
            % do a linear regression and get slopes
            startloc = pkLocs(pkIdx);
            endloc = min(pkLocs(pkIdx)+5,MAX_FRAME_INDEX);
            Y = s_mel(numChan,startloc:endloc)';
            X = [ones(length(Y),1) (1:length(Y))'];
            
            % Y=XB -> B=X\Y, slope=B(2)
            B = X\Y;
            
            %slopes >0 should be wrong?
            if B(2) < 0
                slopes{numChan} = [slopes{numChan}, B(2)]; 
            end
        end
        
    end
end

%% STEP5: For each channel, calculate Slope means over time
slope_means=[];
for numChan = 1:length(slopes)
    slope_means(numChan) = mean(slopes{numChan});
end

% think NaN as 0
slope_means(isnan(slope_means)) = 0;

%% STEP6: Slope mean over bands, Skewness, BassRatio, TrebleRatio 
mean_over_bands = mean(slope_means);
skew = skewness(slope_means);

low = [2 3];
mid = [12 13];
high = [24 25];
%bass ratio: low/mid
BR= sum(slope_means(low)) / (sum(slope_means(mid)) - eps);
%for significantly large BR
BR = min(BR,100);

%treble ratio
TR = sum(slope_means(mid)) / (sum(slope_means(high)) - eps);
% for significantly large TR
TR = min(TR,100);

%% all features concat
fea = [slope_means mean_over_bands skew BR TR];

% plot to show the feature
% plot(fea);
end
