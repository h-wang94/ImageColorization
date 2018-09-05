function [labels, scores, clCost] = PredictSuperpixelsClassesBalancedKNN(nb_classes, nb_dists, Kfs, ...
  nClasses, mcCost, DOUBT)
%Predict superpixel classes using the idea from MATLAB predict function:
%https://www.mathworks.com/help/stats/classificationknn.predict.html
% Chosen class maximizes the posterior probability

eq_dists = [];
for spi = 1:nClasses
  %Transposes because MATLAB is col-major and we are working with rows.
  t = nb_dists';
	t = t((nb_classes' == spi));
	t = reshape(t, length(t)/size(nb_dists,1), size(nb_dists,1))';
	eq_dists = [eq_dists t(:,1:Kfs)];
end

%Generate weights based on distance in feature space
fsdist_w = 1./eq_dists.^4;
fsdist_w = fsdist_w ./ repmat(sum(fsdist_w, 2), [1, nClasses*Kfs]);

%Calculate posterior probability (class given observation)
pP = zeros(size(nb_classes, 1), nClasses);
for sp_i = 1:size(pP, 1)
  for cl_i = 1:nClasses
    clw = fsdist_w(sp_i,(1+(cl_i-1)*Kfs):cl_i*Kfs);
    pP(sp_i, cl_i) = sum(clw);%/sum(fsdist_w(sp_i,:)); 
  end
end
scores = pP;

clCost = pP*mcCost;
[values, labels] = min(clCost, [], 2);

%% Tolerance cost differences
if (DOUBT)
  tol = clCost - repmat(values, 1, nClasses) < 1e-1;
  tol = sum(tol,2);
  %Mark label as not sure
  labels(find(tol > 1)) = -1;
end
end

