function [FeatVectors, FeatLens]  = FeatureExtraction(img_gray, ftsParams, samples)
%Compute the feature vector for each sample pixel of the input image (img_gray)
%
%Current feature list:
%low:     (1)Pixel intensity             
%medium:  (2)Derivatives: Gradient magnitude
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

%Feature Extraction Parameters:
activeFeats = ftsParams.features;
vectorizeFeats = ftsParams.vectorize;

gbParams.wl = ftsParams.gbWl;
gbParams.oris = ftsParams.gbOr;
siftf.patchsize = ftsParams.siftPS;
siftf.gridspacing = ftsParams.siftGs;

%LPF before feature extraction to remove noise.
% img_gray = imfilter(img_gray,fspecial('gaussian',5,0.5),'same','replicate');

%% Features:
            
FeatVectors = [];
FeatLens = [];

if (activeFeats(1))
%Intensity of pixel
  FeatVectors = [FeatVectors;
                 minmaxNormalization(img_gray(idxs), vectorizeFeats(1))];
  
  FeatLens = [FeatLens 1];
end

if (activeFeats(2))
%Gradient  
  [Gmag, ~] = imgradient(img_gray,'sobel');
  
  FeatVectors = [FeatVectors;
                 minmaxNormalization(Gmag(idxs), vectorizeFeats(8))];
  
  FeatLens = [FeatLens 1];
end

if (activeFeats(3))
%Gabor filter bank:
  gaborBank = gabor(gbParams.wl, gbParams.oris);
  [gaborMag, ~] = imgaborfilt(img_gray, gaborBank);

  n_filters = size(gaborMag, 3);
  for i = 1:n_filters
    gabors = gaborMag(:,:,i);
    FeatVectors = [FeatVectors;
                   minmaxNormalization(gabors(idxs), vectorizeFeats(3))];
  end
  
  FeatLens = [FeatLens n_filters];
end

if (activeFeats(4))
%Dense SIFT
  pad_frame = 3;

  % Zero padding to compensate for size change.
  padded_image = padarray(img_gray, [pad_frame, pad_frame]);
  sift_v = dense_sift(padded_image, siftf.patchsize, siftf.gridspacing);

  FeatVectors = [FeatVectors;
                 minmaxNormalization(sift_v(:,idxs), vectorizeFeats(6))];

  FeatLens = [FeatLens 128];
end
  
if (false)
    figure(200);
    for i = 1:size(FeatVectors,1)
        imshow(reshape(FeatVectors(i,:), size(img_gray, 1), size(img_gray, 2)));
%         imwrite(reshape(FeatVectors(i,:), size(img_gray, 1), size(img_gray, 2)), ...
%           ['./../results/feats_temp/' num2str(i) '_r.png'], 'png');
        title(['Feature ' num2str(i)]);
        pause;
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