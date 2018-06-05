%Resize images to reduce cost.

img_name = 'beach2';
img_format = '.jpg';

im = imread(['./../data/' img_name img_format]);

im_rs = imresize(im, 0.5);
imshow(im_rs);

imwrite(im_rs, ['./../data/' img_name '_rs' img_format]);