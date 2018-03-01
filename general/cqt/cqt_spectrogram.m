function cqt_spectrogram(file_path, label)
    figure;
    [path_,file_name,ext_] = fileparts(file_path);
    
    cqt_matrix = load_h5(file_path);

    imagesc(cqt_matrix);
    xlabel('frames');
    ylabel('frequency');
    colorbar;
    set(gca, 'YDir', 'normal'); 
    title(strcat(file_name,',',label), 'interpreter', 'None');
   
    fig_path = strcat('png/',file_name,'.png');
    
    saveas(gcf, fig_path);

end