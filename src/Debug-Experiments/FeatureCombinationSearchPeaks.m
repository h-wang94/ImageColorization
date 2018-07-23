function [] = FeatureCombinationSearchPeaks(source, samples, target, mc_cost, nClusters)
%Test the combinations of available features to determine the best feature
%set from a visual perspective

%Feature subspace combination test
K = 1:15;
idx_range = cumsum(target.fvl);
STATS = 3;

for c = 1:(2^(length(target.fvl)) - 1)
  act_feats = flip(dec2bin(c,length(target.fvl)));

  %Generate feature subset:
  act_idxs = [];
  if(str2num(act_feats(1)))
    for s = 1:STATS
      act_idxs = [act_idxs 1 + (s-1)*idx_range(end)];
    end
  end
  for f = 2:(length(target.fvl))
    if(str2num(act_feats(f)))
      for s = 1:STATS    
        act_idxs = [act_idxs ...
          ((idx_range(f-1)+1):(idx_range(f))) + (s-1)*idx_range(end)];
      end
    end
  end
  act_idxs = sort(act_idxs);
  
  %Call fitcknn for the current subset of features:
  kPeaks = zeros(numel(K),1);
  for k = 1:numel(K)
    [~, Np] = SourceNearestColors(K(k), source, source.fv_sp(act_idxs,:), samples);
    kPeaks(k) = Np;
  end
  [~, minKw_idx] = min(kPeaks);
  %Plot
  fh = figure; plot(K, kPeaks);
  xlabel('Number of nearest neighbors');
  ylabel('Number of peaks on Distance Metric');
  title(['Subset #' num2str(c) ' Active: ' act_feats]);
  saveas(fh, ['./../results/' num2str(c) '_' act_feats '_plot.png'], 'png');
  close(fh);

  %% Colorization
  
  %Matching:
  [neighbor_idxs, neighbor_dists] = knnsearch(source.fv_sp(act_idxs,:)', target.fv_sp(act_idxs,:)');
  out_mt = lab2rgb(CopyClosestSuperpixelAvgColor(source, target, neighbor_idxs));
  
  %Classification:
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

  imwrite(out_mt, ['./../results/' num2str(c) '_' act_feats '_m.png'], 'png');
  imwrite([out_1 out_2], ['./../results/' num2str(c) '_' act_feats '.png'], 'png');
end

end