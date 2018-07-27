function [fvSubsetIdxs] = FeatureSubset(activeFeats, featsLens, nSTATS)
%Given a set of features and the active indexes, returns the active subset.

fvSubsetIdxs = [];
feats_end_idx = cumsum(featsLens);


if(ischar(activeFeats))
activeFeatsNum = zeros(1,length(activeFeats));
  for i = 1:length(activeFeats)
    activeFeatsNum(i) = str2num(activeFeats(i));
  end
else
  activeFeatsNum = activeFeats;
end

if(activeFeatsNum(1))
  for s = 1:nSTATS
    fvSubsetIdxs = [fvSubsetIdxs 1 + (s-1)*feats_end_idx(end)];
  end
end
for f = 2:(length(feats_end_idx))
  if(activeFeatsNum(f))
    for s = 1:nSTATS    
      fvSubsetIdxs = [fvSubsetIdxs ...
        ((feats_end_idx(f-1)+1):(feats_end_idx(f))) + (s-1)*feats_end_idx(end)];
    end
  end
end

fvSubsetIdxs = sort(fvSubsetIdxs);

end

