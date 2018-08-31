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
            images = images_new;
            
            % project the images to the cilinder
            fov = 20;
            [images, images_gray] = project_col_image_C(images,fov);
           
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

    for i=1:n_images
        if i<n_images
            [matches{i}, scores{i}] = vl_ubcmatch(descriptors{i}, descriptors{i+1});
        end
        if i==n_images
            [matches{i}, scores{i}] = vl_ubcmatch(descriptors{i}, descriptors{1});
        end
    end
    
    %% RANSAC computation
    clc;disp('RANSAC...')

    T = zeros(n_images-1,2);
    for i=1:(n_images-1)
        T(i,:) = RANSAC(frames{i},frames{i+1},matches{i});
        if(cond==1)
         T(i,2)=- T(i,2);
        end
        
      
    end
    
    %% parameters
    dim1 = size(images{1},1);
    dim2 = size(images{1},2);
    bit_depth = 8;

    %% get the dimensions of the final images and the corners of all images
    T_max_y = 0;
    T_min_y = 0;
    corn1 = zeros(n_images,2);
    corn1(1,:) = [1,1]; %[x,y]
    for i=1:n_images-1
        T_temp_y = sum(T(1:i,2));
        corn1(i+1,:) = [corn1(i,1)+T(i,1)-1,corn1(i,2)+T(i,2)];
        if T_temp_y > T_max_y
            T_max_y = T_temp_y;
        end
        if T_temp_y < T_min_y
            T_min_y = T_temp_y;
        end
    end
    T_max_y = T_max_y+1;
    T_min_y = T_min_y+1;

    corn1(:,2) = corn1(:,2)-T_min_y+1;
    corn2 = zeros(n_images,2);
    corn3 = zeros(n_images,2);
    corn4 = zeros(n_images,2);
    for i=1:n_images
        corn2(i,:) = corn1(i,:) + [dim2-1,0];
        corn3(i,:) = corn1(i,:) + [0,dim1-1];
        corn4(i,:) = corn1(i,:) + [dim2-1,dim1-1];
    end

    dim1_fi = T_max_y-T_min_y+dim1;
    dim2_fi = sum(T(:,1)-1)+dim2;

    %% computation of the masks
    mask_images = cell(n_images,1);
    for i=1:n_images
        m = zeros(dim1_fi,dim2_fi,3);
        m(corn1(i,2):corn3(i,2),corn1(i,1):corn2(i,1),:) = ones(dim1,dim2,3);
        mask_images{i} = logical(m);
    end

    mask_overlap = cell(n_images-1,1);
    for i=2:n_images
        m1 = double(mask_images{i});
        m_sum = zeros(size(m1));
        m_temp = zeros(size(m1));
        m_temp2 = zeros(size(m1));
        for j=1:(i-1)
            m2 = double(mask_images{j});
            m_sum = m_temp + m2;
            m_temp = double(m_sum>0);

            m_sum = m_temp + m1;
            m_temp2 = double(m_sum>1);
        end
        mask_overlap{i-1} = m_temp2;
        
        
    end

    %% corner of the overlap regions
    corn1_over = zeros(n_images-1,2);
    corn2_over = zeros(n_images-1,2);
    corn3_over = zeros(n_images-1,2);
    overlap_size = zeros(n_images-1,1);
    for i=1:(n_images-1)
        mask = mask_overlap{i}(:,:,1);
        mask_row = ceil(mean(mask));
        mask_col = ceil(mean(mask'));
        c = find(mask_row~=0);
        corn1_over(i,1) = c(1);
        corn3_over(i,1) = corn1_over(i,1);
        corn2_over(i,1) = c(end);
        c = find(mask_col~=0);
        corn1_over(i,2) = c(1);
        corn2_over(i,2) = corn1_over(i,2);
        corn3_over(i,2) = c(end);
        overlap_size(i) = corn2_over(i,1)-corn1_over(i,1)+1;
    end

    %% merge masks
    mask_merge = cell(n_images-1,1);
    beta = 0.4;
    corn1_merge = zeros(n_images-1,2);
    corn2_merge = zeros(n_images-1,2);
    corn3_merge = zeros(n_images-1,2);
    for i=1:n_images-1
        mask = logical(zeros(dim1_fi,dim2_fi,3));
        x_center = floor((corn1_over(i,1)+corn2_over(i,1))/2);
        y1 = corn1_over(i,2); 
        y2 = corn3_over(i,2);
        x1 = floor(x_center-(overlap_size(i)/2)*beta);
        x2 = floor(x_center+(overlap_size(i)/2)*beta);
        mask(y1:y2,x1:x2,:) = 1;
        mask_merge{i} = mask;
        corn1_merge(i,1) = x1;
        corn1_merge(i,2) = y1;
        corn2_merge(i,1) = x2;
        corn2_merge(i,2) = y1;
        corn3_merge(i,1) = x1;
        corn3_merge(i,2) = y2;
        overlap_size(i) = corn2_merge(i,1)-corn1_merge(i,1)+1;
    end
    
    mask_images_no_overlap = cell(n_images,1);
    mask_images_no_overlap{1} = logical(mask_images{1}-mask_merge{1});
    for i=2:(n_images-1)
        mask_images_no_overlap{i} = logical(mask_images{i}-mask_merge{i-1});
    end
    mask_images_no_overlap{n_images} = logical(mask_images{n_images}-mask_merge{n_images-1});
    
    %% generate transformed images
    images_tr = cell(n_images,1);
    for i=1:n_images
        im = uint8(zeros(dim1_fi,dim2_fi,3));
        im(mask_images{i}) = images{i};
        images_tr{i} = im;
    end

    %% merge
   disp('Merging...')
    image_fi = images_tr{1};
    
    
    for i=1:(n_images-1)
        disp(strcat('(im',num2str(i),',im',num2str(i+1),')'))
        im2 = images_tr{i+1};
        im_layer1 = uint8(zeros(size(im2)));
        im_layer2 = uint8(zeros(size(im2)));

        %in the overlap region
        for y=corn1_merge(i,2):corn3_merge(i,2)
            j=1;
            for x=corn1_merge(i,1):corn2_merge(i,1)
                alpha = j/overlap_size(i);
                im_layer1(y,x,:) = (1-alpha)*image_fi(y,x,:);
                im_layer2(y,x,:) = alpha*im2(y,x,:);
                image_fi(y,x,:) = im_layer1(y,x,:) + im_layer2(y,x,:);
                j=j+1;
                
            end
            
            for x=corn2_merge(i,1)+1:corn2(i+1,1)
                image_fi(y,x,:) = im2(y,x,:);
            end
        end
    end
    
    %% plots
    figure
    imshow(image_fi); hold on;
    plot(corn1(:,1),corn1(:,2),'g*');
    plot(corn2(:,1),corn2(:,2),'g*');
    plot(corn3(:,1),corn3(:,2),'g*');
    plot(corn1_over(:,1),corn1_over(:,2),'r*');
    plot(corn2_over(:,1),corn2_over(:,2),'r*');
    plot(corn3_over(:,1),corn3_over(:,2),'r*');
    plot(corn1_merge(:,1),corn1_merge(:,2),'b*');
    plot(corn2_merge(:,1),corn2_merge(:,2),'b*');
    plot(corn3_merge(:,1),corn3_merge(:,2),'b*');
    legend('image region corner','','','overlap region corner','','','merge region corner')
    
    %% erease the black borders
    y1 = max(corn1(:,2));
    y2 = min(corn3(:,2));
    image_fi = image_fi(y1:y2,:,:);

end
