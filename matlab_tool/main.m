scenes = {'resi', 'bus', 'park', 'forest', 'train', 'tram', 'general'};

for scene = scenes
    scene_str = scene{1};

    rootdir = '/home/songhongwei/data_home/scene_example/';
    sceneroot = strcat(rootdir, scene_str, '/');
    metafile = strcat(sceneroot, 'tool/', scene_str, '.txt');
    cqt_dirs = strcat(sceneroot, 'cqt/audio/*.h5');
    cqt_audio_root = strcat(sceneroot, 'cqt/audio/');

    meta_map = create_metamap(metafile);

    cqt_files = dir(cqt_dirs);
    for file = cqt_files'
        file_path = strcat(cqt_audio_root,file.name);

        [pathstr, name, ext] = fileparts(file_path);
        %lookup meta_map to get label
        label = meta_map(fullfile('audio/',name));

        cqt_spectrogram(file_path, label);
        %break;
    end
    
end


