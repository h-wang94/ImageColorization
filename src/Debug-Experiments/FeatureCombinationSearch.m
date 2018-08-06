function [distsSPComb, distsMedians] = FeatureCombinationSearch(source, target, samples, fvLen, ... 
  mc_cost, nClusters, K, metric, outColor)
%Test the combinations of available features to determine the best feature
%set from a visual perspective
%distMedians receives the ...

STATS = 1; %TODO: ajustar feature extraction

NfeatsCombinations = 2^(length(fvLen)) - 1;

%Accel
cd = pdist(source.sp_chrom');
CD = squareform(cd);

distsSPComb = zeros(NfeatsCombinations, length(source.validSuperpixels));
distsMedians = zeros(2, NfeatsCombinations);
for c = 1:NfeatsCombinations
  act_feats = flip(dec2bin(c,length(fvLen)));

  %Generate features subset:
  act_idxs = FeatureSubset(act_feats, fvLen, STATS);
  
  switch metric
    case 'meanDist'
      FD = squareform(pdist(source.fv_sp(act_idxs,:)'));
      [mcd, dmc] = SourceSPNNColorsDists(K, [], source.sp_chrom, ...
        source.validSuperpixels, source.lin_sp, samples, FD, CD);

      distsMedians(1,c) = mean(mcd);
      distsMedians(2,c) = mean(dmc);

    case 'peaksDist'
      %Computes the number of peaks based on median distances in color
      %space from k-NN on feature space.
      FD = squareform(pdist(source.fv_sp(act_idxs,:)'));
      [mcd, dmc] = SourceSPNNColorsDists(K, [], source.sp_chrom, ...
        source.validSuperpixels, source.lin_sp, samples, FD, CD);

      distsMedians(1,c) = sum(mcd > cd(floor(length(cd)/4)));
      distsMedians(2,c) = sum(dmc > cd(floor(length(cd)/4)));

    case 'cluster'
      %Computes clustering metrics (intra/inter-class distances).
      distsMedians(:,c) = ClusterDistNormMedian([source.fv_sp(act_idxs,:)' source.sp_clusters'], ...
        mc_cost, 'pdist');

    case 'leaveOneOut'
      %Computes the Leave One Out cross validation over source image
      [neighbor_idxs, neighbor_dists] = knnsearch(source.fv_sp(act_idxs,:)', source.fv_sp(act_idxs,:)', ...
        'K', source.nSuperpixels); % Return all distances for further reference.
      neighbor_classes = source.sp_clusters(neighbor_idxs(:,2:end));
      [labels, ~] = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists(:,2:end), K, nClusters, ...
        mc_cost, false);

      colorNN = zeros(2,length(source.validSuperpixels));
      for i = 1:length(source.validSuperpixels)
        colorNN(:,i) = median(source.sp_chrom(:, neighbor_idxs(i,2:K+1)),2);
      end
      distsMedians(1,c) = mean( (labels ~= source.sp_clusters').* ...
        sqrt(sum((colorNN - source.sp_chrom).^2))' );
      
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
    case 'peaksDist'
      eval = distsMedians(2,:);
      [~,c] = min(eval);
    case 'cluster'
      eval = distsMedians(1,:)./distsMedians(2,:);
      [~,c] = min(eval(5:end));
      c = c + 4;
    case 'leaveOneOut'
      eval = distsMedians(1,:);
      [~,c] = min(eval);
  end
  %Find best combination
  act_feats = flip(dec2bin(c,length(fvLen)));
  act_idxs = FeatureSubset(act_feats, fvLen, STATS);

  %Classification and colorization:
  [neighbor_idxs, neighbor_dists] = knnsearch(source.fv_sp(act_idxs,:)', target.fv_sp(act_idxs,:)', ...
    'K', source.nSuperpixels); 
  neighbor_classes = source.sp_clusters(neighbor_idxs);
  labels = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, K, nClusters, ...
    mc_cost, true);
  rgb_out = lab2rgb(CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
    neighbor_classes, labels));
  imwrite(rgb_out, ['./../results/' num2str(c) '_' act_feats 'K' num2str(K) metric '.png'], 'png');
end

end