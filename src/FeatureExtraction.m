function [FeatVectors, FeatLens]  = FeatureExtraction(img_gray, ftsParams, samples)
%Compute the feature vector for each sample pixel of the input image (img_gray)
%
%Current feature list:
%1: Pixel luminance             
%2: Window standard deviation   
%3: Gabor filter banks          
%4: DCT window
%5: DFT window
%6: Dense SIFT
%7: Haralick texture descriptor

if nargin < 3
    %Whole image (for target)
    idxs = 1:size(img_gray,1)*size(img_gray,2);
else
    idxs = sub2ind(size(img_gray), samples(1,:), samples(2,:));
end

%Parameters input:
activeFeats = ftsParams.features;
vectorizeFeats = ftsParams.vectorize;

stdWS = ftsParams.stdWS;
gbParams.wl = ftsParams.gbWl;
gbParams.oris = ftsParams.gbOr;
cosineWS = ftsParams.dctWS;
fourierWS = ftsParams.dftWS;
siftf.patchsize = ftsParams.siftPS;
siftf.gridspacing = ftsParams.siftGs;
haraWS = ftsParams.haraWS;

%% Features:
%TODO: mudar para alocação estática

%LPF before feature extraction to remove noise.
% img_gray = imfilter(img_gray,fspecial('gaussian',5,0.5),'same','replicate');

FeatVectors = [];
FeatLens = [];
if (activeFeats(1))
%Luminance of pixel
  FeatVectors = [FeatVectors;
                 minmaxNormalization(img_gray(idxs), vectorizeFeats(1))];
  
  FeatLens = [FeatLens 1];
end

if (activeFeats(2))
%Std dev neighborhood
  stds = sd_neighborhood(img_gray, stdWS);
  FeatVectors = [FeatVectors;
                 minmaxNormalization(stds(idxs), vectorizeFeats(2))];

  FeatLens = [FeatLens 1];
end

if (activeFeats(3))
%Gabor filter bank:
    gaborBank = gabor(gbParams.wl, gbParams.oris);
    [gaborMag, ~] = imgaborfilt(img_gray, gaborBank);

    nFilters = size(gaborMag, 3);
    for i = 1:nFilters
        gabors = gaborMag(:,:,i);
%         imshow(gabors);
        FeatVectors = [FeatVectors;
                       minmaxNormalization(gabors(idxs), vectorizeFeats(3))];
    end
  
  FeatLens = [FeatLens nFilters];
end

if (activeFeats(4))
%Discrete Cosine Transform
    dcts = WindowFeature(img_gray, 'dct', cosineWS);
    
    FeatVectors = [FeatVectors;
                   minmaxNormalization(dcts(:,idxs), vectorizeFeats(4))];
                 
  FeatLens = [FeatLens cosineWS*cosineWS];
end

if (activeFeats(5))
%Discrete Fourier Transform window
    dfts = WindowFeature(img_gray, 'dft', fourierWS);
    
    FeatVectors = [FeatVectors;
                   minmaxNormalization(dfts(:,idxs), vectorizeFeats(5))];
                 
  FeatLens = [FeatLens fourierWS*fourierWS];
end

if (activeFeats(6))
%Dense SIFT
    pad_frame = 3;
    
    % Zero padding to compensate for size change.
    padded_image = padarray(img_gray, [pad_frame, pad_frame]);
    sift_v = dense_sift(padded_image, siftf.patchsize, siftf.gridspacing);
    
    FeatVectors = [FeatVectors;
                   minmaxNormalization(sift_v(:,idxs), vectorizeFeats(6))];
                 
  FeatLens = [FeatLens 128];
end

if (activeFeats(7))
%Haralick Texture Features
  htf = WindowFeature(img_gray, 'haralick', haraWS);
  
  htf = htf(find(~isnan(sum(htf,2))),:);
  
  FeatVectors = [FeatVectors;
                 minmaxNormalization(htf(:,idxs), vectorizeFeats(7))];
   
  FeatLens = [FeatLens size(htf,1)];
end

if (activeFeats(8))
%DAISY descriptor
    dzy = compute_daisy(img_gray);
    disp('Testing...');
end

if (true)
    figure(200);
    for i = 1:size(FeatVectors,1)
        imshow(reshape(FeatVectors(i,:), size(img_gray, 1), size(img_gray, 2)));
        title(['Feature ' num2str(i)]);
    end
end



end

function Y = minmaxNormalization(X, vec)
%Computes the minmax of a set of samples
% The result is to put every sample inside a hypercube of lenght one in
% each dimension.
    dim = size(X, 1);
    
    if (dim == 1)
        Y = (X - min(X))/(max(X) - min(X));
    else
        if (vec)
            Y = zeros(size(X));
            for i = 1:dim
                Y(i,:) = minmaxNormalization(X(i,:), false);
            end
        else
            Y = (X - min(min(X)))/(max(max(X)) - min(min(X)));
        end
    end
end

