function [FeatVectors, featWeights]  = FeatureExtraction(img_gray, featActive, samples)
%Compute the feature vector for each sample pixel of the input image (img_gray)
%Also determines the weight for each feature.
%Current feature list:
%1: Pixel luminance             [0-1]
%2: Window standard deviation   [?]
%3: Gabor filter banks          [

if nargin < 3
    idxs = 1:size(img_gray,1)*size(img_gray,2);
else
    idxs = sub2ind(size(img_gray), samples(1,:), samples(2,:));
end

%% Features:
%TODO: mudar para alocação estática

FeatVectors = [];
featWeights = [];

if (featActive(1))
%Luminance
    FeatVectors = [FeatVectors;
                   img_gray(idxs)];
    featWeights = [featWeights;
                   1];
end

if (featActive(2))
%Std dev neighborhood
    stds = sd_neighborhood(img_gray, 5);
    FeatVectors = [FeatVectors;
                   stds(idxs)];
    featWeights = [featWeights;
                   1];
end

% %LPF/HPF
% h = fspecial('log', 3);
% imf = NormalizeImage(imfilter(img_gray,h));
% feat_vect = [feat_vect;
%              imf(idxs)];
% h = fspecial('gaussian', 7);
% imf = NormalizeImage(imfilter(img_gray,h));
% feat_vect = [feat_vect;
%              imf(idxs)];

if (featActive(3))
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

%Normalize weights
featWeights = featWeights/sum(featWeights);
end

