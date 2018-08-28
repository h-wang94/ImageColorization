function distance = FeaturesDistances(spfv1, spfM2)
  nFeats = 4;
  nStats = 2;  
  nSP = size(spfM2,1);
  fLen = length(spfv1);
  %Function pdist2 tries to compute element-wise first.
  assert(nSP ~= 1); 

  minmax = @(x) (x - min(x))/(max(x) - min(x)); 
  
  %Number of histogram bins
  nBins = (fLen - nStats*(nFeats+1))/nFeats;

  %Indexes (scalar x histogram features)
  idx = 1;
  scalar_idxs = zeros(1, fLen, 'logical');
  for fi = 1:(nFeats+1)
    scalar_idxs(idx:idx+1) = true;
    idx = idx + nStats + nBins;
  end

  %Compute distances of scalar features
  spfM1 = repmat(spfv1, nSP, 1);
  scalar_dists = (spfM1(:,scalar_idxs) - spfM2(:,scalar_idxs)).^2;
  
  %Compute histogram distances
  spfH1 = spfv1(:,~scalar_idxs);
  spfH2 = spfM2(:,~scalar_idxs);
  hists_dists = zeros(nSP,nFeats);
  for fi = 1:nFeats
    hists_dists(:,fi) = match_distance(spfH1(1+(fi-1)*nBins:fi*nBins), spfH2(:,1+(fi-1)*nBins:fi*nBins));
    hists_dists(:,fi) = minmax(hists_dists(:,fi));
  end
  
  distance = mean(scalar_dists,2) + mean(hists_dists,2);
end