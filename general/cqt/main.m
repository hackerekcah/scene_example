meta_map = create_metamap('../meta.txt');

cqt_files = dir('audio/*.h5');
for file = cqt_files'
    file_path = strcat('audio/',file.name);
    
    [pathstr, name, ext] = fileparts(file_path);
    %lookup meta_map to get label
    label = meta_map(fullfile(pathstr,name));
    
    cqt_spectrogram(file_path, label);
    break;
end
