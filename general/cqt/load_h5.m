function fea = load_h5(h5file_name)
    data_h5 = hdf5read(h5file_name,'/cqt');
    %h5disp(h5file_name)

    fea = zeros(291,334);
    
    %convert h5array to matlab array
    for i = 1: size(data_h5, 1)
        fea(:,i) = data_h5(i).Data;
    end

end