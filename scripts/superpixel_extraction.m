Nsp = 500;

imc = imread('./../data/chamomile3.jpg');
img = rgb2gray(imc);

[Lc, Nc] = superpixels(imc, Nsp);
[Lg, Ng] = superpixels(img, Nsp);

figure;
BMc = boundarymask(Lc);
BMg = boundarymask(Lg);
imshow([imoverlay(imc, BMc, 'r') imoverlay(img, BMg, 'r')]);
