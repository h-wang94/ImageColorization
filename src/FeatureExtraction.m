function [FeatVectors, featWeights]  = FeatureExtraction(img_gray, activeFeats, samples)
%Compute the feature vector for each sample pixel of the input image (img_gray)
%
%Current feature list:
%1: Pixel luminance             
%2: Window standard deviation   
%3: Gabor filter banks          
%4: DCT window
%5: DFT window
%6: Dense SIFT

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
                   minmaxNormalization(img_gray(idxs), false)];               
    featWeights = [featWeights;
                   1];
end

% if (activeFeats())
% %Mean luminance w
%     
% end

if (activeFeats(2))
%Std dev neighborhood
    stds = sd_neighborhood(img_gray, 5);
    FeatVectors = [FeatVectors;
                   minmaxNormalization(stds(idxs), false)];
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
                       minmaxNormalization(aux(idxs), false)];
    end
    featWeights = [featWeights;
                   ones(nFilters,1)/nFilters];
end

if (activeFeats(4))
%Discrete Cosine Transform
    kWS = 7;
    f_len = kWS*kWS;
    dcts = WindowFeature(img_gray, 'dct', kWS);
    
    FeatVectors = [FeatVectors;
                   minmaxNormalization(dcts(:,idxs), false)];
    featWeights = [featWeights;
                   ones(f_len,1)/f_len];
end

if (activeFeats(5))
%Discrete Fourier Transform window
    kWS = 7;
    f_len = kWS*kWS;
    dfts = WindowFeature(img_gray, 'dft', kWS);
    
    FeatVectors = [FeatVectors;
                   minmaxNormalization(dfts(:,idxs), false)];
    featWeights = [featWeights;
                   ones(f_len,1)/f_len];
end

if (activeFeats(6))
%Dense SIFT
    patchsize=8;
    gridspacing=1;
    desclen = 128;
    pad_frame = 3;
    
    % Zero padding to compensate for size change.
    padded_image = padarray(img_gray, [pad_frame, pad_frame]);
    sift_v = dense_sift(padded_image, patchsize, gridspacing);
    
    FeatVectors = [FeatVectors;
                   minmaxNormalization(sift_v(:,idxs), false)];
    featWeights = [featWeights;
                   ones(desclen,1)/desclen];
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
% coeff = pca(FeatVectors');

featWeights = featWeights/sum(featWeights);
% featWeights = ones(length(featWeights),1);

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

