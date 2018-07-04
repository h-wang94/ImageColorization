while(true)
  rn = randperm(216);
  rn = rn(1:50);
  
  [neighbor_idxs, neighbor_dists] = knnsearch(source.fv_sp(rn,:)', target.fv_sp(rn,:)', ...
      'K', source.nSuperpixels); % Return all distances for further reference.
    neighbor_classes = source.sp_clusters(neighbor_idxs);
  
  labels = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, 3, IP.nClusters, ...
      mc_cost);
  
  figure; imshow(lab2rgb(CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
      neighbor_classes, labels)));
    title('Predict with costs');
    
  pause;  
end