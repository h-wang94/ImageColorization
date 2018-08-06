function [medianColDist, dist2MedianCol] = SourceSPNNColorsDists(K, src_fv_sp, srcSPChrom, ...
  validSPs, SPIdxsLin,  samples, FD, CD)
%Computes the COLOR space distances from each superpixel to its nearest
%neighbors in FEATURE space.

nVSP = length(validSPs);

%Pairwise distance in Chrominance and Feature spaces.
if (nargin < 8)
  FD = squareform(pdist(src_fv_sp'));
  CD = squareform(pdist(srcSPChrom'));
end
  
%Median of color distance for all spixels
medianColDist = zeros(1, nVSP);
%Distance to color median for all spixels
dist2MedianCol = zeros(1, nVSP);
for spi = 1:nVSP
  [~, nn] = sort(FD(:,spi));
  
  median_color = median(srcSPChrom(:,nn(2:K+1)),2);
  dist2MedianCol(spi) = norm(median_color - srcSPChrom(:,spi));
  
  medianColDist(spi) = median(CD(nn(2:K+1),spi));
end

% %Test
% stem(dist2MedianCol, 'markerSize', 3); drawnow;
end

