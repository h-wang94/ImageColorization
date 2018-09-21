function fitness=LeaveOneEval(pop, source, K, mc_cost)
  
  fitness = zeros(size(pop,1),1)';
  for i = 1:size(pop,1)
    %Pairwise distance in Feature subset
    featSubset = source.fv_sp(pop(i,:),:);
    
    [neighbor_idxs, neighbor_dists] = knnsearch(featSubset', featSubset', ...
      'K', source.nSuperpixels); % Return all distances for further reference.
    neighbor_idxs = neighbor_idxs(:,2:end);
    neighbor_dists = neighbor_dists(:,2:end);
    neighbor_classes = source.sp_clusters(neighbor_idxs);
    [labels, ~] = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, K, nClusters, ...
      mc_cost, false);

    %Median of colors of assigned classes.
    colorNN = zeros(2,length(source.validSuperpixels));
    for spi = 1:length(source.validSuperpixels)
      class_instantes = find(neighbor_classes(spi,:) == labels(spi));
      colorNN(:,spi) = source.sp_chrom(:, neighbor_idxs(spi, class_instantes(1)));
    end
    
    fitness(i) = -mean( (labels ~= source.sp_clusters').* ...
      sqrt(sum((colorNN - source.sp_chrom).^2))' );
  end

end