%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT
% images      : cell containing the set of images that has to be merged,
% n_im        : (optional) select how many images have to be merged,
% str_project : (optional) select the projection of the images to use.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
% image_fi    : merged image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function image_fi = merge(images,n_im,cond,str_project)
%% select the number or the image
    n_images = 0;
    images_gray=cell(n_im,1);
    switch nargin
        case 1
            n_images = size(images,1);
        case 3
            n_images = n_im;
            images_new = cell(n_images);
            for i=1:n_images
                images_new{i} = images{i};
                images_gray{i}=rgb2gray(images{i});
            end
            images = images_new;
            
        case 4
            n_images = n_im;
            images_new = cell(n_images);
            for i=1:n_images
                images_new{i} = images{i};
               
            end
            % project the images to the cilinder
             fov=70;
             images = images_new;
             [images,images_gray]=project_col_image_C(images,fov);
            
           
    end
    
    %% run vlfeat toolbox
    run('vlfeat/toolbox/vl_setup')
    
    %% SIFT
    clc;disp('Calculation of the SIFT...')
    frames = cell(n_images,1);
    descriptors = cell(n_images,1);

    for i=1:n_images
        [frames{i},descriptors{i}] = vl_sift(single(images_gray{i}));
    end
    
    %% find matches
    clc;disp('Calculation of the Matches...')
    matches = cell(n_images,1);
    scores = cell(n_images,1);

    for i=1:n_images-1
        [matches{i}, scores{i}] = vl_ubcmatch(descriptors{i+1}, descriptors{i},0.5);
    end
    
    %% RANSAC computation
    clc;disp('RANSAC...')
    
    M = cell(n_images-1,1);
    T = cell(n_images-1,1);
    for i=1:n_images-1
        [M{i},T{i}] =RANSAC(matches{i},frames{i+1},frames{i},3);
    end
    
    %% parameters
    dim1 = size(images{1},1);
    dim2 = size(images{1},2);

    %% get the dimensions of the final images and the corners of all images
    dim1_fi = round(n_images/2)*dim1;
    dim2_fi = n_images*dim2;
    
    im_fi = uint8(zeros(dim1_fi,dim2_fi,3));
    shift = round(dim1_fi/2-dim1/2);
    
    %creation of the slides
    diapositiva = cell(n_images,1);
    im_fi1 = im_fi;
    im_fi1((1+shift):(dim1+shift),1:dim2,:) = images{1}(1:dim1,1:dim2,:);
    diapositiva{1} = im_fi1;
    x = 1:dim2;
    y = 1:dim1;
    [X,Y] = meshgrid(x,y);
    indices = zeros(dim1,dim2,2);
    indices(:,:,1) = X; indices(:,:,2) = Y;
    corn1=zeros(n_images,2);
    corn2=zeros(n_images,2);
    corn3=zeros(n_images,2);
    center=zeros(n_images,2);
    corn1(1,:) = [1,1+shift];
    corn2(1,:) = [dim2,1+shift];
    corn3(1,:) = [dim2,dim1+shift];
    center(1,:) = round((corn1(1,:)+corn3(1,:))/2);
    mask_images = cell(n_images,1);
%     
 
    mask_temp = uint8(ones(dim1_fi,dim2_fi,3)*2);
    mask_temp(corn1(1,2):corn3(1,2),corn1(1,1):corn2(1,1),:) = 255;
    mask_images{1} = mask_temp;
    cont=1;
    for j=2:n_images
        imag_temp = uint8(zeros(dim1_fi,dim2_fi,3));
        mask_temp = uint8(ones(dim1_fi,dim2_fi,3)*j*2);
        for m=1:dim2
            for n=1:dim1
                p = [indices(n,m,1),indices(n,m,2)]';
                p_new = M{j-1}*p + T{j-1}; p_new = round(p_new);

                if n==1 && m==1
                    corn1(j,1)=p_new(1);
                    corn1(j,2)=p_new(2)+shift;
                end
                
                if n==1 && m==dim2
                    corn2(j,1)=p_new(1);
                    corn2(j,2)=p_new(2)+shift;
                end
                
                if n==dim1 && m==dim2
                    corn3(j,1)=p_new(1);
                    corn3(j,2)=p_new(2)+shift;
                    
                end

                imag_temp(p_new(2)+shift,p_new(1),:) = images{j}(n,m,:);
                mask_temp(p_new(2)+shift,p_new(1),:) = 255;
                indices(n,m,:) = p_new;

            end
        end
       
        center(j,:) = round((corn1(j,:)+corn3(j,:))/2);
        diapositiva{j} = imag_temp;
        mask_images{j} = mask_temp;
        figure('Name','Slide Images','NumberTitle','off');
        imshow(diapositiva{j}); hold on;
        plot(corn1(j,1),corn1(j,2),'r*');
        plot(corn2(j,1),corn2(j,2),'r*');
        plot(corn3(j,1),corn3(j,2),'r*');
        plot(center(j,1),center(j,2),'b*');
    end
    
    %% overlap mask
    mask_overlap = cell(n_images-1,1);
    for i=2:n_images
        mask_overlap{i}=double(mask_images{i})-double(mask_images{i-1});
      
    for n=1:dim1_fi
        for m=1:dim2_fi
         if( mask_overlap{i}(n,m,1)==0 &&  mask_overlap{i}(n,m,2)==0 && mask_overlap{i}(n,m,3)==0)
            %compute the distance from the center
            d1=sqrt((n-center(i-1,2))^2+(m-center(i-1,1))^2);
            d2=sqrt((n-center(i,2))^2+(m-center(i,1))^2);
            im_fi(n,m,:)=uint8((d2/(d1+d2))*double(diapositiva{i-1}(n,m,:))+(d1/(d1+d2))*double(diapositiva{i}(n,m,:)));
            
         end
         if(  mask_overlap{i}(n,m,1)==255-(i-1)*2 &&  mask_overlap{i}(n,m,2)==255-(i-1)*2 && mask_overlap{i}(n,m,3)==255-(i-1)*2)
         im_fi(n,m,:)=diapositiva{i}(n,m,:);
         end
         if( mask_overlap{i}(n,m,1)<0 &&  mask_overlap{i}(n,m,2)<0 && mask_overlap{i}(n,m,3)<0)
         im_fi(n,m,:)=diapositiva{i-1}(n,m,:);
         end
        end
    end

     figure('Name','Temp Merge','NumberTitle','off')
     imshow(im_fi);  
    end
    
    %% resize
    bound_down=corn3(n_images,2);
    bound_right=corn2(n_images,1);
    image_fi=zeros(bound_down-(1+shift),bound_right,3);
    i1=1;
    i2=1;
    
    for y=(1+shift):bound_down
        for x=1:bound_right 
            
            image_fi(i1,i2,:)=im_fi(y,x,:);
            i2=i2+1;
           
 
        end
        i2=1;
        i1=i1+1;
    end
    image_fi=(uint8(image_fi));
    imshow(image_fi);

end
