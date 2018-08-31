clear;
clc;
close all;


%% load images
clc;disp('Load Images...')

dir_input = 'input';
subfolders=dir(dir_input);
img_struct=struct ;
cont=1;


for i=1:length(subfolders)
    
current_folder_name = subfolders(i).name;
angle=sscanf(current_folder_name,'angle_%d');

if(~isempty(angle))
        dir=strcat('input/',current_folder_name);
       

        %you must indicate the correct extension of the images!!
        [images,n_images,fov] = read_images(dir,'jpg'); 
        
        %% merge
        if(n_images>=2)
        image_fi = merge(images,n_images,0,'cil_projection');
      
        
        else
         image_fi=images{1};
          
        end
        name_pan=strcat('panoramic_angle_',num2str(angle),'.jpg');
        imwrite(image_fi,name_pan);
        img_struct(cont).panoramic=image_fi;
        img_struct(cont).angle=angle;
        cont=cont+1;
        

end
end


%% vertical merge
dim=size(img_struct);
if(dim(2)~=1)
    im_final=verticale_merge(img_struct);
    imwrite(im_final,'final_merge.jpg');
end


