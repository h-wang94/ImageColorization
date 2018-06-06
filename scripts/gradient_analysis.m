%Gradient analysis
close all;

im = rgb2lab(imread('C:\Users\Saulo\Documents\GitHub\ImageColorization\data\eagle1.jpg'));
Gmag = zeros(size(im));
Gdir = zeros(size(im));

for i = 1:3
    [Gmag(:,:,i), Gdir(:,:,i)] = imgradient(im(:,:,i), 'prewitt');
    
    Gmagmm = (Gmag(:,:,i) - min(min(Gmag(:,:,i))))/(max(max(Gmag(:,:,i))) - min(min(Gmag(:,:,i))));

    figure; imshow(Gmagmm.*Gdir(:,:,i),[]); 
    colormap hot; colorbar; caxis([-180 180]);
    title(['Natural: channel ' num2str(i)]);
end

%%

im = rgb2lab(imread('C:\Users\Saulo\Documents\GitHub\ImageColorization\results\beach2.jpg'));
Gmag = zeros(size(im));
Gdir = zeros(size(im));

for i = 1:3
    [Gmag(:,:,i), Gdir(:,:,i)] = imgradient(im(:,:,i), 'sobel');
    
    Gmagmm = (Gmag(:,:,i) - min(min(Gmag(:,:,i))))/(max(max(Gmag(:,:,i))) - min(min(Gmag(:,:,i))));

    figure; imshow(Gmagmm.*Gdir(:,:,i),[]); 
	colormap hot; colorbar; caxis([-180 180]);
    title(['Artificial: channel ' num2str(i)]);
end
