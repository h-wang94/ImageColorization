function distance = FeaturesDistances(spfv1, spfv2)
  nSP = size(spfv2,1);
  if (nSP == 1)
    error('O MATLAB não documentou que eu precisaria fazer isso.');
  end;

  nBins = (length(spfv1) - 2 - 2*4)/4;

  scalar_dists = zeros(nSP,4);
  hists_dists = zeros(nSP,4);
  idx = 1;
  for fi = 1:4
    scalar_dists(fi) = norm(repmat(spfv1(idx:idx+1),nSP,1) - ...
      spfv2(:,idx:idx+1));
    hists_dists(fi) = match_distance(spfv1(idx+2:idx+2+nBins), spfv2(:,idx+2:idx+2+nBins));
    idx = idx + 2 + nBins;
  end
  
  distance = mean(scalar_dists) + mean(hists_dists);
end

% , 'Distance', @FeaturesDistances