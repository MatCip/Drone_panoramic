function [images,n_images,fov] = read_images(sdirectory,ext)
    fov = 66;
    
    ext_str = strcat('/*.',ext);
    files = dir([sdirectory ext_str]);
    n_images = length(files);
    images = cell(n_images,1);
    
    
    for i = 1:n_images
        filename =files(i).name;
        
        numb=sscanf(filename,'i%d');
        
        if(~isempty(numb))
        im_name=strcat(sdirectory,'/',filename);
        images{i} = imread(im_name);
        
        end
    end
end