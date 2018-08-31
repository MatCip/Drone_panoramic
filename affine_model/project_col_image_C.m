function [images_c, images_c_gray] = project_col_image_C(images,fov)
   disp('Projection...')
    n_images = size(images,1);
    images_c = cell(n_images,1);
    images_c_gray = cell(n_images,1);
    dim1 = size(images{1},1);
    dim2 = size(images{1},2);

    for i=1:n_images
        image = images{i};
        image_c = uint8(zeros(dim1,dim2,3));
       
        image_c(:,:,1) = projectIC(image(:,:,1),fov/2);
        image_c(:,:,2) = projectIC(image(:,:,2),fov/2);
        image_c(:,:,3) = projectIC(image(:,:,3),fov/2);
        images_c{i} = image_c;
        images_c_gray{i} = projectIC_gray(image,fov/2);
    end
end

%%
%projects the images on a cylindrical surface
% image : planar image
% imageC : image projected on the cylindrical surface
% angle: half FOV of the camera 
function [imageC] = projectIC(image,angle)

    ig = image;
    [h w] = size(ig);
    imageC = uint8(zeros(h,w));

    alpha = angle/180*pi;
    d = (w/2)/tan(alpha);
    r = d/cos(alpha);



    for x = -w/2+1:w/2
        for y = -h/2+1:h/2

           x1 = d * tan(x/r);
           y1 = y * (d/r) /cos(x/r);

           if x1 < w/2 && x1 > -w/2+1 && y1 < h/2 && y1 > -h/2+1 
                imageC(y+(h/2), x+(w/2) ) = ig(round(y1+(h/2)),round(x1+(w/2)));
            end
        end
    end
end

function [imageC] = projectIC_gray(image,angle)

    ig = rgb2gray(image);
    [h w] = size(ig);
    imageC = uint8(zeros(h,w));

    alpha = angle/180*pi;
    d = (w/2)/tan(alpha);
    r = d/cos(alpha);



    for x = -w/2+1:w/2
        for y = -h/2+1:h/2

           x1 = d * tan(x/r);
           y1 = y * (d/r) /cos(x/r);

           if x1 < w/2 && x1 > -w/2+1 && y1 < h/2 && y1 > -h/2+1 
                imageC(y+(h/2), x+(w/2) ) = ig(round(y1+(h/2)),round(x1+(w/2)));
            end
        end
    end
end


