function success = save_image(img, img_name, SAVE, output_folder)
    if nargin < 3
        SAVE = true;
        output_folder = '../output/';
    elseif nargin < 4
        output_folder = '../output/';
    end
    
    success = false;
    if SAVE
        imwrite(img, strcat(output_folder, img_name));
        success = true;
    end
end