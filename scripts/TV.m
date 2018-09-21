%Total variation test:

% clear all;
close all;

%%

c_im = target.rgb;
figure; imshow(c_im);
w = weights;
% w = w/(max(w) - min(w));
% c_im = imread('C:\Users\Saulo\Documents\GitHub\ImageColorization\data\kodim09.png');

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
    W = sparse(1:length(w), 1:length(w), w);
    A = sptoeplitz(c, r);
    A = A + W;
    
    %Sparse linear system solution
    x = A \ (w.*d);

    lab_out(:,:,i) = reshape(x, size(c_channel, 1), size(c_channel, 2));
    figure; imshow(lab_out(:,:,i),[]);
    title('Output channel');
end

out_rgb = lab2rgb(lab_out);

figure; imshow(out_rgb); title('Regularized image');