function [labels, scores] = PredictSuperpixelsClassesKNN(nb_classes, nb_dists, Kfs, nClasses, mcCost)
%Predict superpixel classes using the idea from MATLAB predict function:
%https://www.mathworks.com/help/stats/classificationknn.predict.html
% Chosen class maximizes the posterior probability

%Determine classes of Kfs neighborhood of each pixel
nb_classes = nb_classes(:,1:Kfs);

%Generate weights based on distance in feature space
fsdist_w = 1./nb_dists(:,1:Kfs);
fsdist_w = fsdist_w ./ repmat(sum(fsdist_w, 2), [1, Kfs]);

%Calculate posterior probability (class given observation)
pP = zeros(size(nb_classes, 1), nClasses);
for sp_i = 1:size(pP, 1)
  for cl_i = 1:nClasses
    clw = (nb_classes(sp_i,:) == cl_i ).*fsdist_w(sp_i,:);
    pP(sp_i, cl_i) = sum(clw);%/sum(fsdist_w(sp_i,:));
  end
end
scores = pP;

clCost = pP*mcCost;
[values, labels] = min(clCost, [], 2);

%% Tolerance cost differences
tol = clCost - repmat(values, 1, nClasses) < 1e-1;
tol = sum(tol,2);
%Mark label as not sure
labels(find(tol > 1)) = -1;

end

