function imgs = load_images(img_name, images_folder) 
    %% load images
    if nargin < 2
        images_folder = '../images/';
    end
    target_images_folder = strcat(images_folder, 'target/');
    gsource_images_folder = strcat(images_folder, 'gray_source/');
    csource_images_folder = strcat(images_folder, 'color_source/');


    target_image_name = strcat(target_images_folder, img_name);
    gsource_image_name = strcat(gsource_images_folder, img_name);
    csource_image_name = strcat(csource_images_folder, img_name);
    target_image = im2double(imread(target_image_name));
    gsource_image = im2double(imread(gsource_image_name));
    csource_image = im2double(imread(csource_image_name));
    imgs.target_image = target_image;
    imgs.gsource_image = gsource_image;
    imgs.csource_image = csource_image;
end
