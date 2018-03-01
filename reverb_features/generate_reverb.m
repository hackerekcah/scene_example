function [reverb, fs_ir] = generate_reverb(clean_file, ir_file)

[clean, fs_ori] = audioread(clean_file);

[ir, fs_ir] = audioread(ir_file);

sig = resample(clean, fs_ir, fs_ori);

reverb = conv(sig,ir);

% normalize to [-1,1]
reverb = reverb ./ max(abs(reverb));

% soundsc(sig, fs_ir);
% soundsc(ir, fs_ir);
% soundsc(reverb, fs_ir);
end