function [tgt_spfv, src_spfv] = SuperpixelsFeatures(source, samples, target)
%Computes the features of superpixels based on its pixel elements.

nBins = 10;  %TODO: Parameter!
nFeats = length(target.fvl);
ftsBounds = cumsum(target.fvl);
ftsBounds = [0 ftsBounds];

minmax = @(x) (x - min(x))/(max(x) - min(x));

%Feature Aggregation 1:
%- Scalars: remain the same.
%- Vectors: apply statistic.
agg = { @(x) x, ...
        @(x) x, ...
        @(x) minmax(std(x)), ...
        @(x) minmax(median(x))};

t_fv = zeros(nFeats, size(target.fv,2));
s_fv = zeros(nFeats, size(samples.fv,2));
for fi = 1:nFeats
  t_fv(fi,:) = agg{fi}(target.fv(ftsBounds(fi)+1:ftsBounds(fi+1),:));
  s_fv(fi,:) = agg{fi}(samples.fv(ftsBounds(fi)+1:ftsBounds(fi+1),:));
end

tgt_spfv = fillSPFV(t_fv, target.lin_sp, target.nSuperpixels, nBins, nFeats);
src_spfv = fillSPFV(s_fv, source.lin_sp, source.nSuperpixels, nBins, nFeats);

end

function sp_fv = fillSPFV(fv, linSPIdxs, nSP, nBins, nFeats)
  minmax = @(x) (x - min(x))/(max(x) - min(x));
  
  %Statistics to be computed for each feature
  binsEdges = 0:(1/nBins):1;
  statsLen = [1 1 nBins];
  statsBounds = [0 cumsum(statsLen)];
  stats = { @(x) mean(x), ...
            @(x) std(x), ...
            @(x) histcounts(x, binsEdges) / length(x)};
  
  
  %Compute the Superpixel Features
  sp_fv = zeros(sum(statsLen)*nFeats, nSP);
  for spi = 1:nSP
    mask = (linSPIdxs == spi);
    
    for fi = 1:nFeats
      vals = fv(fi,mask);
      
      for si = 1:length(stats)
        sp_fv((fi-1)*sum(statsLen) + (statsBounds(si)+1:statsBounds(si+1)), spi) = stats{si}(vals);
      end
    end
  end
  
  %MinMax normalization of each scalar feature (for balanced distance
  %computation).
  idx = 1;
  for fi = 1:nFeats
    sp_fv(idx,:) = minmax(sp_fv(idx,:));
    sp_fv(idx + 1,:) = minmax(sp_fv(idx + 1,:));
    idx = idx + nBins + 2;
  end
end
