function [medianColDist, dist2MedianCol] = SourceSPNNColorsDistMedian(Knn, src_fv_sp, ...
  validSPs, SPIdxsLin,  samples)
%Computes the COLOR space distances from each superpixel to its nearest
%neighbors in FEATURE space.

nVSP = length(validSPs);

%Compute superpixel color by averaging (TODO)
sp_ab = zeros(2, nVSP);
for spi = 1:length(sp_ab)
  vsi = validSPs(spi);
  vsi_idxs = (SPIdxsLin == vsi);

  sp_ab(:,spi) = median(samples.ab(:,vsi_idxs),2);
  
  %DEBUG> Colors inside same superpixel
%   sum(vsi_idxs)
%   test = samples.ab(:,vsi_idxs);
%   subplot(1,2,1); scatter(test(1,:), test(2,:), '.'); hold on;
%   subplot(1,2,2); scatter(sp_ab(1,spi), sp_ab(2,spi), '.'); hold on;
%   pause(0.01)
%   drawnow;
end

%Pairwise distance in Chrominance and Feature spaces.
fdist = pdist(src_fv_sp');
FD = squareform(fdist);
cdist = pdist(sp_ab');
CD = squareform(cdist);

%Median of color distance for all spixels
medianColDist = zeros(1, nVSP);

%Distance to color median for all spixels
dist2MedianCol = zeros(1, nVSP);

for spi = 1:nVSP
  [~, nn] = sort(FD(:,spi));
  median_color = median(sp_ab(:,nn(2:Knn+1)),2);

  dist2MedianCol(spi) = norm(median_color - sp_ab(:,spi));
  medianColDist(spi) = median(CD(nn(2:Knn+1),spi));
end

end

