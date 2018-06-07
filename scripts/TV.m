%Total variation test:

clear all;
close all;

%%

c_im = imread('C:\Users\Saulo\Documents\GitHub\ImageColorization\results\beach2.jpg');

lab_in = rgb2lab(c_im);

lab_out = zeros(size(lab_in));
lab_out(:,:,1) = lab_in(:,:,1);

%% Data:
for i = 2:3

    c_channel = lab_in(:,:,i);
    figure; imshow(c_channel, []);
    title('Input channel');

    %data term
    d = reshape(c_channel, size(c_channel, 1)*size(c_channel, 2), 1);

    %Sparse Hessian matrix
    r = [5 -2 zeros(1, length(d) - 2)];
    c = r';
    A = sptoeplitz(c, r);

    %Sparse linear system solution
    x = A \ d;

    lab_out(:,:,i) = reshape(x, size(c_channel, 1), size(c_channel, 2));
    figure; imshow(lab_out(:,:,i),[]);
    title('Output channel');
end

out_rgb = lab2rgb(lab_out);

figure; imshow(out_rgb); title('Regularized image');