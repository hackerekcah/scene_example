file_name = 'a002_0_10';

s = audioread(['audio_lda_mis/audio/', file_name, '.wav']);
s = s(:,1);
s = resample(s,20000, 44100);
save(['segment/',file_name, '.dat'], 's', '-ascii');

seg_res = load(['segment/', file_name, '.out']);
imagesc(seg_res');
set(gca, 'YDir', 'normal'); 