function [dist2medianNN] = SourceSPNNColorsMedianDist(Knn, source, src_fv_sp, samples)
%TODO

nVSP = length(source.validSuperpixels);

%Compute superpixel color by averaging
sp_ab = zeros(2, nVSP);
for spi = 1:length(sp_ab)
  vsi = source.validSuperpixels(spi);
  vsi_idxs = source.lin_sp == vsi;

  sp_ab(:,spi) = mean(samples.ab(:,vsi_idxs),2);
end

%Pairwise distance in Chrominance and Feature spaces.
% cdist = pdist(sp_ab');
fdist = pdist(src_fv_sp');
% CD = squareform(cdist);
FD = squareform(fdist);

%Color median of nearest superpixel
dist2medianNN = zeros(1, nVSP);
for spi = 1:nVSP
  [~, nn] = sort(FD(:,spi));
  median_color = median(sp_ab(:,nn(2:Knn+1)),2);
  dist2medianNN(spi) = norm(median_color - sp_ab(:,spi));
end

end

