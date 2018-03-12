function feat_vect = FeatureExtraction(img_gray, feature_bool, samples)
%Computes the features indicated by feature_bool over the idxs of img_gray.
%Feature list:
%1: Pixel luminance
%2: Window standard deviation

if nargin < 3
    idxs = 1:size(img_gray,1)*size(img_gray,2);
else
    idxs = sub2ind(size(img_gray), samples(2,:), samples(1,:));
end
    

%% Features:
%TODO: mudar para estatico
feat_vect = [];

%Luminance
feat_vect = [feat_vect;
             img_gray(idxs)];

%Std dev neighborhood
stds = sd_neighborhood(img_gray, 5);
feat_vect = [feat_vect;
             stds(idxs)];

end

