function [distsSPComb, distsMedians] = FeatureCombinationSearch(source, samples, target, ... 
  mc_cost, nClusters, K, type, outColor)
%Test the combinations of available features to determine the best feature
%set from a visual perspective

%TODO: ajustar feature extraction
nSTATS = 3;

NfeatsCombinations = 2^(length(target.fvl)) - 1;

distsSPComb = zeros(NfeatsCombinations, length(source.validSuperpixels));
distsMedians = zeros(2, NfeatsCombinations);
for c = 1:NfeatsCombinations
  act_feats = flip(dec2bin(c,length(target.fvl)));

  %Generate features subset:
  act_idxs = FeatureSubset(act_feats, target.fvl, nSTATS);
  
  switch type
    case 'peaks'
      kPeaks = zeros(numel(K),1);
      for k = 1:numel(K)
        [~, Np] = SourceNearestColors(K(k), source, source.fv_sp(act_idxs,:), samples);
        kPeaks(k) = Np;
      end
      [~, minK_idx] = min(kPeaks);
      
      %Plot
      fh = figure; plot(K, kPeaks);
      xlabel('Number of nearest neighbors');
      ylabel('Number of peaks on Distance Metric');
      title(['Subset #' num2str(c) ' Active: ' act_feats]);
      saveas(fh, ['./../results/' num2str(c) '_' act_feats '_plot.png'], 'png');
      close(fh);
      
    case 'medianColor'
      distsMedians(:,c) = ClusterDistNormMedian([source.fv_sp(act_idxs,:)' source.sp_clusters'], ...
        mc_cost, 'pdist');
      
      %Compute table of distances
%       for k = 1:numel(K)
%         [Np] = SourceSPNNColorsMedianDist(K(k), source, source.fv_sp(act_idxs,:), samples);
%         distsSPComb(c,:) = Np;
%       end
      
    case 'randSubspace'
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
      
    otherwise
      disp('Invalid Feature Combination Search type');
  end

  %% Colorization

  if (outColor)
    %Matching:
    [neighbor_idxs, neighbor_dists] = knnsearch(source.fv_sp(act_idxs,:)', target.fv_sp(act_idxs,:)');
    out_mt = lab2rgb(CopyClosestSuperpixelAvgColor(source, target, neighbor_idxs));

    %Classification:
    [neighbor_idxs, neighbor_dists] = knnsearch(source.fv_sp(act_idxs,:)', target.fv_sp(act_idxs,:)', ...
      'K', source.nSuperpixels); 
    neighbor_classes = source.sp_clusters(neighbor_idxs);

    labels = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, minK_idx, nClusters, ...
      mc_cost);
    labels_m = modeTies(neighbor_classes(:,1:minK_idx));

    out_1 = lab2rgb(CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
      neighbor_classes, labels));
    out_2 = lab2rgb(CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
      neighbor_classes, labels_m));

    imwrite(out_mt, ['./../results/' num2str(c) '_' act_feats '_m.png'], 'png');
    imwrite([out_1 out_2], ['./../results/' num2str(c) '_' act_feats '.png'], 'png');
  end
end

end