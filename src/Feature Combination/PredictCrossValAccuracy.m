function [pAcc, testLabels] = PredictCrossValAccuracy(ldata, tdata, kK, mcCost)
%Computes the Accuracy of the customized Predict function on a cross
%validation setting.

[neighbor_idxs, neighbor_dists] = knnsearch(ldata(:,1:(end-1)), tdata(:,1:(end-1)), ...
  'K', kK);
llabels = ldata(:,end);
neighbor_classes = llabels(neighbor_idxs);

[testLabels, ~] = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, kK, ...
  length(mcCost), mcCost);

pAcc = sum(testLabels == tdata(:,end))/sum(testLabels ~= -1);

end

