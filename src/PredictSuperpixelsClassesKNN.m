function [labelsScores, labelsCosts, doubts, clScores, clCost] = ...
  PredictSuperpixelsClassesKNN(nb_classes, nb_dists, Kfs, nClasses, mcCost)
%Predict superpixel classes using the idea from MATLAB predict function:
%https://www.mathworks.com/help/stats/classificationknn.predict.html
% Chosen class is the one that minimizes the cost based on posterior
% probability.
%-Kfs < 0, normal with K = Kfs.
%-Kfs = 0, all superpixels.
%-Kfs > 0, equality with K = Kfs.
%-mcCost: If present, labels minimize cost. If absent, maximizes posterior. 

if (Kfs < 0)
  Kfs = abs(Kfs);
  
  %Samples the structures using only the closest.
  nb_classes = nb_classes(:,1:Kfs);
  nb_dists = nb_dists(:,1:Kfs);
elseif (Kfs == 0)
  %Do nothing.
  Kfs = size(nb_classes, 2);
else
  eq_dists = [];
  n_p_c = [];
  for spi = 1:nClasses
    %Transposes because MATLAB is col-major and we are working with rows.
    t = nb_dists';
    t = t((nb_classes' == spi));
    t = reshape(t, length(t)/size(nb_dists,1), size(nb_dists,1))';
    eq_dists = [eq_dists t(:,1:min([Kfs size(t, 2)]))];
    n_p_c = [n_p_c min([Kfs size(t, 2)])];
  end
  %nb_dists receives the closest K from each class.
  nb_dists = eq_dists;
  
  nb_classes = [];
  %Not the smartest solution. But adapts to already running code.
  for i = 1:nClasses
    nb_classes = [nb_classes i*ones(size(nb_dists,1), n_p_c(i))];
  end
end

%Generate weights based on distance in feature space
fsdist_w = 1./nb_dists.^4;
fsdist_w = fsdist_w ./ repmat(sum(fsdist_w, 2), [1, size(fsdist_w,2)]);

%Calculate posterior probability (class given observation)
pP = zeros(size(nb_classes, 1), nClasses);
for sp_i = 1:size(pP, 1)
  for cl_i = 1:nClasses
    clw = (nb_classes(sp_i,:) == cl_i).*fsdist_w(sp_i,:);
    pP(sp_i, cl_i) = sum(clw); %Already normalized.
  end
end

%Scores and costs assignemnt
clScores = pP;
clCost = pP*mcCost;

%Labeling
[scores, labelsScores] = max(clScores, [], 2);
[~, labelsCosts] = min(clCost, [], 2);

%% Scores Label doubt:

tol = (repmat(scores, 1, nClasses) - clScores < 0.1);
tol = sum(tol,2);
%Mark label as not sure
doubts = (tol > 1);

end

