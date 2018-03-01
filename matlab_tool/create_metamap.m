function meta_map = create_metamap(meta_file)

    meta_map = containers.Map;

    fid = fopen(meta_file);
    line = fgetl(fid);
    
    while ischar(line)
        str_arr = split(line);
        
        % keys of map must be char type, conver string to char array
        meta_map(char(str_arr(1))) = char(str_arr(2));
        line = fgetl(fid);
    end
    
    fclose(fid);
end
