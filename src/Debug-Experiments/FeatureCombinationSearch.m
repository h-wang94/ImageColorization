function [distsSPComb, distsMedians] = FeatureCombinationSearch(source, target, samples, fvLen, ... 
  mc_cost, nClusters, K, metric, outName)
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
      neighbor_idxs = neighbor_idxs(:,2:end);
      neighbor_dists = neighbor_dists(:,2:end);
      neighbor_classes = source.sp_clusters(neighbor_idxs);
      [labels, ~] = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, K, nClusters, ...
        mc_cost, false);

      %Median of colors of assigned class.
      colorNN = zeros(2,length(source.validSuperpixels));
      for i = 1:length(source.validSuperpixels)
        class_instantes = find(neighbor_classes(i,:) == labels(i));
        colorNN(:,i) = source.sp_chrom(:, neighbor_idxs(i, class_instantes(1)));
      end
      distsMedians(1,c) = mean( (labels ~= source.sp_clusters').* ...
        sqrt(sum((colorNN - source.sp_chrom).^2))' );
      
    case 'randSubspace'
      
    otherwise
      disp('Invalid Feature Combination Search type');
  end
end

  %% Colorization
if (~isempty(outName))
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
  imwrite(rgb_out, ['./../results/' outName num2str(c) '_' 'K' num2str(K) metric '.png'], 'png');
end

end