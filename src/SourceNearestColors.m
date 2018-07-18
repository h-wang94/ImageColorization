function [NearestColorMedian, nPeaks] = SourceNearestColors(Knn, source, src_fv_sp, samples)
%TODO

nVSP = length(source.validSuperpixels);

%Compute average color of superpixels
sp_ab = zeros(2, nVSP);
for spi = 1:length(sp_ab)
  vsi = source.validSuperpixels(spi);
  vsi_idxs = source.lin_sp == vsi;

  sp_ab(:,spi) = mean(samples.ab(:,vsi_idxs),2);
end

%Pairwise distance in Chrominance and Feature spaces.
cdist = pdist(sp_ab');
fdist = pdist(src_fv_sp');
CD = squareform(cdist);
FD = squareform(fdist);

%Average the kNN color distance for all spixels
NearestColorMedian = zeros(1, nVSP);
for spi = 1:nVSP
  [~, nn] = sort(FD(:,spi));
  NearestColorMedian(spi) = median(CD(nn(2:Knn+1),spi));
end

%Compute number of peaks
nPeaks = sum(NearestColorMedian > 40);

end

