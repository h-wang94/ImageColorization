function [distsSPComb, distsMedians] = FeatureCombinationSearch(source, target, samples, fvLen, ... 
  mc_cost, nClusters, K, metric, outColor)
%Test the combinations of available features to determine the best feature
%set from a visual perspective
%distMedians receives the ...

nSTATS = 3; %TODO: ajustar feature extraction
kPeakThresh = 25;

NfeatsCombinations = 2^(length(fvLen)) - 1;

%Accel
CD = squareform(pdist(source.sp_chrom'));

distsSPComb = zeros(NfeatsCombinations, length(source.validSuperpixels));
distsMedians = zeros(2, NfeatsCombinations);
for c = 1:NfeatsCombinations
  act_feats = flip(dec2bin(c,length(fvLen)));

  %Generate features subset:
  act_idxs = FeatureSubset(act_feats, fvLen, nSTATS);
  
  switch metric
    case 'peaks'
      %Computes the number of peaks based on median distances in color
      %space from k-NN on feature space.
      FD = squareform(pdist(source.fv_sp(act_idxs,:)'));
      [mcd, dmc] = SourceSPNNColorsDists(K, [], source.sp_chrom, ...
        source.validSuperpixels, source.lin_sp, samples, FD, CD);
    
      distsMedians(1,c) = sum(mcd > kPeakThresh);
      distsMedians(2,c) = sum(dmc > kPeakThresh);
      
    case 'cluster'
      %Computes clustering metrics (intra/inter-class distances).
      distsMedians(:,c) = ClusterDistNormMedian([source.fv_sp(act_idxs,:)' source.sp_clusters'], ...
        mc_cost, 'pdist');

    case 'randSubspace'
      %Call fitcknn for the current subset of features:
%       cvloss_w = zeros(numel(K),1);
%       for k = 1:numel(K)
%         knn_w = fitcknn(source.fv_sp(act_idxs,:)', source.sp_clusters',...
%             'NumNeighbors', K(k), 'CrossVal', 'On', ...
%             'Cost', mc_cost, 'BreakTies', 'nearest');
%         cvloss_w(k) = kfoldLoss(knn_w);
%       end
%       [~, minKw_idx] = min(cvloss_w);
%       %Plot
%       fh = figure; plot(K, cvloss_w);
%       xlabel('Number of nearest neighbors');
%       ylabel('10 fold classification error');
%       title(['Subset #' num2str(c) ' Active: ' act_feats]);
%       saveas(fh, ['./../results/' num2str(c) '_' act_feats '_plot.png'], 'png');
%       close(fh);
      
    otherwise
      disp('Invalid Feature Combination Search type');
  end
end

  %% Colorization
if (outColor)
  switch metric
    case 'peaks'
      eval = distsMedians(2,:);
    case 'cluster'
      eval = distsMedians(1,:)./distsMedians(2,:);
  end
  %Find best combination
  [~,c] = min(eval(5:end));
  act_feats = flip(dec2bin(c+4,length(fvLen)));
  act_idxs = FeatureSubset(act_feats, fvLen, nSTATS);

  %Classification and colorization:
  [neighbor_idxs, neighbor_dists] = knnsearch(source.fv_sp(act_idxs,:)', target.fv_sp(act_idxs,:)', ...
    'K', source.nSuperpixels); 
  neighbor_classes = source.sp_clusters(neighbor_idxs);
  labels = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, K, nClusters, ...
    mc_cost);
  rgb_out = lab2rgb(CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
    neighbor_classes, labels));
  imwrite(rgb_out, ['./../results/' num2str(c) '_' act_feats 'K' num2str(K) '.png'], 'png');
end

end