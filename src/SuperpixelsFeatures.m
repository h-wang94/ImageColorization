function [tgt_spfv, src_spfv] = SuperpixelsFeatures(source, samples, target)
%

%Histogram length:
nBins = 10;


%Number of superpixels of each image
src_nSP = length(source.validSuperpixels);

%Statitics of vector features:
t_fv = zeros(4, size(target.fv,2));
s_fv = zeros(4, size(samples.fv,2));

t_fv(1,:) = target.fv(1,:);
t_fv(2,:) = target.fv(2,:);
t_fv(3,:) = std(target.fv(3:110,:));
t_fv(4,:) = median(target.fv(111:end,:));

s_fv(1,:) = samples.fv(1,:);
s_fv(2,:) = samples.fv(2,:);
s_fv(3,:) = std(samples.fv(3:110, :));
s_fv(4,:) = median(samples.fv(111:end,:));

tgt_spfv = fillSPFV(t_fv, target.lin_sp, target.nSuperpixels, nBins);
src_spfv = fillSPFV(s_fv, source.lin_sp, length(source.validSuperpixels), nBins);

end

function sp_fv = fillSPFV(fv, linSPIdxs, nSP, nBins)
  sp_fv_len = (2 + nBins)*4;
  binsEdges = 0:(1/nBins):1;

  sp_fv = zeros(sp_fv_len,nSP);
  for i = 1:nSP
    mask = (linSPIdxs == i);
    
    idx = 1;
    for fi = 1:4
      vals = fv(fi,mask);
      %Mean
      sp_fv(idx,i) = mean(vals); 
      idx = idx + 1;
      %Std
      sp_fv(idx,i) = std(vals);
      idx = idx + 1;
      %Normalized histogram
      sp_fv(idx:(idx+nBins-1),i) = histcounts(vals,binsEdges) / length(vals);
      idx = idx + nBins;
    end
  end
  %MinMax normalization of each feature (for further distance computation)
  idx = 1;
  for fi = 1:4
    sp_fv(idx,:) = minmaxNormalization(sp_fv(idx,:));
    idx = idx + 1;
    sp_fv(idx,:) = minmaxNormalization(sp_fv(idx,:));
    idx = idx + nBins + 1;
  end
end

function Y = minmaxNormalization(X)
  Y = (X - min(X))/(max(X) - min(X));
end
