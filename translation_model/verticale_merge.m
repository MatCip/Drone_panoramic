function im_final=verticale_merge(img_struct) 


dim=size(img_struct);
n_imgs=dim(2);
imgs_new=cell(n_imgs,1);

%% padding the images 
img_struct=padding(img_struct,n_imgs); 

%% rotating all the image of 90°
for i=1:n_imgs
imgs_new{i}=img_struct(i).panoramic;
end


for i=1:n_imgs
    imgs_new{i}=imrotate(imgs_new{i},90);
end
[imgs_new,images_c_gray]= project_col_image_C(imgs_new,33);


% %% only for windows user 
% 
% 
% for i=1:n_imgs
% [image, descriptors, locs]=sift(images_c_gray{i});
%  img_sift(i).descriptors=descriptors;
%  img_sift(i).locs=locs;
% end
% for i=1:n_imgs-1
%    [match,num]= match2(images_c_gray{i},img_sift(i).descriptors, img_sift(i).locs,images_c_gray{i+1},img_sift(i+1).descriptors, img_sift(i+1).locs,0.5);
% end

%%  merging and re-rotating 
im_final=merge(imgs_new,n_imgs,1);
im_final=imrotate(im_final,-90);


end


%% padding function

function img_struct_new = padding (img_struct,n_imgs)
img_struct_new=img_struct;
dimensions=zeros(n_imgs,3);
for i=1:n_imgs
     dim=size(img_struct(i).panoramic);
     dimensions(i,:)=dim;
     
end
[M_y,~]=min(dimensions(:,1));
[M_x,~]=min(dimensions(:,2));

for i=1:n_imgs
     %dim=size(img_struct(i).panoramic);
     %image_new=uint8(zeros(M_y,M_x,3));
     
     image_new = img_struct(i).panoramic(1:M_y,1:M_x,:);
     %for j=1:dim(1)
     %    for k=1:dim(2)
     %         image_new(j,k,1)=img_struct(i).panoramic(j,k,1);
     %         image_new(j,k,2)=img_struct(i).panoramic(j,k,2);
     %         image_new(j,k,3)=img_struct(i).panoramic(j,k,3);
     %    end
     %end
     %figure
     %imshow(image_new);
     %disp(size(image_new))
     img_struct_new(i).panoramic = image_new;
     
end



end