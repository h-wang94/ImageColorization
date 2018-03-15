function feat_vect = FeatureExtraction(img_gray, feature_bool, samples)
%Computes the features indicated by feature_bool over the idxs of img_gray.
%Current feature list:
%1: Pixel luminance             [0-1]
%2: Window standard deviation   [?]

if nargin < 3
    idxs = 1:size(img_gray,1)*size(img_gray,2);
else
    idxs = sub2ind(size(img_gray), samples(2,:), samples(1,:));
end
    

%% Features:
%TODO: mudar para alocação estática

feat_vect = [];

%Luminance
feat_vect = [feat_vect;
             img_gray(idxs)];

%Std dev neighborhood
stds = sd_neighborhood(img_gray, 5);
feat_vect = [feat_vect;
             stds(idxs)];
   
%Filter test:
h = fspecial('log', 3);
imf = NormalizeImage(imfilter(img_gray,h));
feat_vect = [feat_vect;
             imf(idxs)];


h = fspecial('gaussian', 7);
imf = NormalizeImage(imfilter(img_gray,h));
feat_vect = [feat_vect;
             imf(idxs)];



%LBP:

end

