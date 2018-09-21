source.image = imread('./../data/beach1.jpg');

sz = size(source.image);
rsf = 30;

sz_rs = [sz(1) - mod(sz(1),rsf) sz(2) - mod(sz(2),rsf) sz(3)];

L = rgb2lab(source.image(1:sz_rs(1), 1:sz_rs(2),:));
% imshow(L)

Ls = imresize(L, 1/rsf)
Lsa = imresize(Ls, rsf);

Lr = zeros(sz_rs(1), sz_rs(2), 3);

Lr(:,:,1) = L(:,:,1);
Lr(:,:,2:3) = Lsa(:,:,2:3);
im = lab2rgb(Lr);
imshow(im)

figure;
imshow(source.image)
