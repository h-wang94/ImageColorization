function fitness=ColorEval(pop, source, samples, colorDists, K, peakThresh)
  
  %Pairwise distance in Chrominance space.
  CD = colorDists;
  
  fitness = zeros(size(pop,1),1)';
  for i = 1:size(pop,1)
    %Pairwise distance in Feature subset
    featSubset = source.fv_sp(pop(i,:),:);
    FDsub = squareform(pdist(featSubset'));

    [~, dmc] = SourceSPNNColorsDists(K, [], source.sp_chrom, ...
      source.validSuperpixels, source.lin_sp, samples, FDsub, CD);
    
    fitness(i) = -sum(dmc > peakThresh);
  end
end