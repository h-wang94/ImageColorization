function [] = FeatureCombinationSearch(source, target, mc_cost, nClusters)
%Test the combinations of available features to determine the best feature
%set from a visual perspective

%Feature subspace combination test
K = 1:15;
idx_range = cumsum(target.fvl);

for c = 100:(2^(length(FP.features)-1) - 1)
  act_feats = flip(dec2bin(c,length(target.fvl)));

  %Generate feature subset:
  act_idxs = [];
  if(str2num(act_feats(1)))
    act_idxs = [act_idxs 1];
  end
  for f = 2:(length(target.fvl)-1)
    if(str2num(act_feats(f)))
      act_idxs = [act_idxs (idx_range(f-1)+1):(idx_range(f))];
    end
  end

  %Call fitcknn for the current subset of features:
  cvloss_w = zeros(numel(K),1);
  for k = 1:numel(K)
    knn_w = fitcknn(source.fv_sp(act_idxs,:)', source.sp_clusters',...
        'NumNeighbors', K(k), 'CrossVal', 'On', ...
        'Cost', mc_cost, 'BreakTies', 'nearest');
    cvloss_w(k) = kfoldLoss(knn_w);
  end
  [~, minKw_idx] = min(cvloss_w);
  %Plot
  fh = figure; plot(K, cvloss_w);
  xlabel('Number of nearest neighbors');
  ylabel('10 fold classification error');
  title(['Subset #' num2str(c) ' Active: ' act_feats]);
  saveas(fh, ['./../results/' num2str(c) '_' act_feats '_plot.png'], 'png');
  close(fh);

  [neighbor_idxs, neighbor_dists] = knnsearch(source.fv_sp(act_idxs,:)', target.fv_sp(act_idxs,:)', ...
    'K', source.nSuperpixels); 
  neighbor_classes = source.sp_clusters(neighbor_idxs);

  labels = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, minKw_idx, nClusters, ...
    mc_cost);
  labels_m = modeTies(neighbor_classes(:,1:minKw_idx));

  out_1 = lab2rgb(CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
    neighbor_classes, labels));
  out_2 = lab2rgb(CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
    neighbor_classes, labels_m));

  imwrite([out_1 out_2], ['./../results/' num2str(c) '_' act_feats '.png'], 'png');
end

end