function [FeatVectors, featWeights]  = FeatureExtraction(img_gray, activeFeats, samples)
%Compute the feature vector for each sample pixel of the input image (img_gray)
%Normalize each feature so they have same influence.
%-scalar: x' = (x - minx) / (maxx - minx)
%-vector:
%
%Current feature list:
%1: Pixel luminance             [0-1]
%2: Window standard deviation   [0-0.5]
%3: Gabor filter banks          [
%4: DCT window                  [

if nargin < 3
    idxs = 1:size(img_gray,1)*size(img_gray,2);
else
    idxs = sub2ind(size(img_gray), samples(1,:), samples(2,:));
end

%% Features:
%TODO: mudar para alocação estática

FeatVectors = [];
featWeights = [];

if (activeFeats(1))
%Luminance of pixel
    FeatVectors = [FeatVectors;
                   img_gray(idxs)];               
    featWeights = [featWeights;
                   1];
end

if (activeFeats(2))
%Std dev neighborhood
    stds = sd_neighborhood(img_gray, 5);
    FeatVectors = [FeatVectors;
                   stds(idxs)];
    featWeights = [featWeights;
                   1];
end

if (activeFeats(3))
%Gabor filter bank:
    gaborBank = gabor(2.^(1:1), 0:-30:-150);
    [gaborMag, ~] = imgaborfilt(img_gray, gaborBank);

    nFilters = size(gaborMag, 3);
    for i = 1:nFilters
        aux = gaborMag(:,:,i);
        FeatVectors = [FeatVectors;
                       aux(idxs)];
    end
    featWeights = [featWeights;
                   ones(nFilters,1)/nFilters];
end

if (activeFeats(4))
%Discrete Cosine Transform
    kWS = 7;
    f_len = kWS*kWS;
    dcts = WindowFeature(img_gray, 'dct');
%     dcts = IsotropicScaling(dcts);
    
    FeatVectors = [FeatVectors;
                   dcts(:,idxs)];
    featWeights = [featWeights;
                   ones(f_len,1)/f_len];
end

if (activeFeats(5))
%Discrete Fourier Transform window
    kWS = 11;
    f_len = kWS*kWS;
    dfts = WindowFeature(img_gray, 'dft', kWS);
%     dfts = IsotropicScaling(dfts);
    
    FeatVectors = [FeatVectors;
                   dfts(:,idxs)];
    featWeights = [featWeights;
                   ones(f_len,1)/f_len];
end

if (activeFeats(6))
%Dense SIFT
    patchsize=8;
    gridspacing=1;
    
    pad_frame = 3;
    padded_image = padarray(img_gray, [pad_frame, pad_frame]);
    
    sift_v =dense_sift(padded_image, patchsize, gridspacing);
    
    FeatVectors = [FeatVectors;
                   sift_v(:,idxs)];
    featWeights = [featWeights;
                   ones(128,1)/128];
end

if (activeFeats(7))
%DAISY descriptor
    dzy = compute_daisy(img_gray);
    disp('Testing...');
end

%% Normalize Features

%Each column of coeff contains coefficients for one principal component
%and the columns are in descending order of component variance. 
% By default, pca centers the data and uses the singular value decomposition (SVD) algorithm.
coeff = pca(FeatVectors');

featWeights = featWeights/sum(featWeights);

end

