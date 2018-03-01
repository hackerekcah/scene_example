addpath('cochleagram');
file_name = 'outsidefield_reverb';
[s, fs] = audioread(['reverb_features/ReverbAudio/', file_name, '.wav']);
sig = s(:,1);

%gammatone filter
r = gammatone(sig, 128, [50,11050], fs);

%generate 
a = cochleagram(r, fs*0.02);
cochplot(a, [50,11050], 10);