function [FeatVectors, FeatLens]  = FeatureExtraction(img_gray, ftsParams, samples)
%Compute the feature vector for each sample pixel of the input image (img_gray)
%
%Current feature list:
%low:     (1)Pixel intensity             
%medium:  (2)Derivatives: Gradient directions weighted by magnitude.
%         (3)Frequency: Gabor filter banks          
%         (4)Invariants: Dense SIFT
%high:    (5)Object: Saliency (after superpixel extraction)

%DEPRECATED:
%STD window
%DCT window
%DFT window

if nargin < 3
    %Whole image (for target)
    idxs = 1:size(img_gray,1)*size(img_gray,2);
else
    idxs = sub2ind(size(img_gray), samples(1,:), samples(2,:));
end

%Features Parameters:
gbParams.wl = ftsParams.gbWl;
gbParams.oris = ftsParams.gbOr;
siftf.patchsize = ftsParams.siftPS;
siftf.gridspacing = ftsParams.siftGs;

%LPF before feature extraction to remove noise.
% img_gray = imfilter(img_gray,fspecial('gaussian',5,0.5),'same','replicate');

            
FeatVectors = [];
FeatLens = [];

%% Intensity of pixel
FeatVectors = [FeatVectors;
               minmaxNormalization(img_gray(idxs), [])];

FeatLens = [FeatLens 1];

%% Gradient (Magnitude x Direction)
[Gmag, Gdir] = imgradient(img_gray,'intermediate');

FeatVectors = [FeatVectors;
               minmaxNormalization(Gmag(idxs), []);
               Gdir(idxs)];

FeatLens = [FeatLens 1 1];

%% Gabor filter bank:
%Test> 19/09:
wavelengthMin = 4/sqrt(2);
wavelengthMax = hypot(size(img_gray,1),size(img_gray,2));
n = floor(log2(wavelengthMax/wavelengthMin));
wavelength = 2.^(0:(n-2-1)) * wavelengthMin;
deltaTheta = 15;
orientation = 0:deltaTheta:(180-deltaTheta);
gbParams.wl = wavelength;
gbParams.oris = orientation;

gaborBank = gabor(gbParams.wl, gbParams.oris);
[gaborMag, ~] = imgaborfilt(img_gray, gaborBank);

n_filters = size(gaborMag, 3);
for i = 1:n_filters
  gabors = gaborMag(:,:,i);
  FeatVectors = [FeatVectors;
                 minmaxNormalization(gabors(idxs), true)];
end

FeatLens = [FeatLens n_filters];

%% Dense SIFT
pad_frame = 3;

% Zero padding to compensate for size change.
padded_image = padarray(img_gray, [pad_frame, pad_frame]);
sift_v = dense_sift(padded_image, siftf.patchsize, siftf.gridspacing);

FeatVectors = [FeatVectors;
               minmaxNormalization(sift_v(:,idxs), true)];

FeatLens = [FeatLens 128];

%%   
if (false)
    figure(200);
    for i = 1:size(FeatVectors,1)
        imshow(reshape(FeatVectors(i,:), size(img_gray, 1), size(img_gray, 2)),[]);
%         imwrite(reshape(FeatVectors(i,:), size(img_gray, 1), size(img_gray, 2)), ...
%           ['./../results/feats_temp/' num2str(i) '_r.png'], 'png');
        title(['Feature ' num2str(i)]);
        pause(0.1);
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