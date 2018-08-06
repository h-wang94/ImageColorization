function [tgt_spfv, source_spfv] = SuperpixelFeatures(source, samples, target, SP_STATS)
%Builds the feature vector of superpixels by using statistics.

f_len = size(target.fv,1);

tgt_nSP = target.nSuperpixels;
src_nSP = length(source.validSuperpixels);

%List of Statistics:
stats = {@(x) mean(x,2), ...
         @(x) std(x,[],2), ...
         @(x) median(x,2)};

tgt_spfv = zeros(sum(SP_STATS)*f_len, tgt_nSP);
source_spfv = zeros(sum(SP_STATS)*f_len, src_nSP);
st_count = 0;
for sti = 1:length(SP_STATS)
  if (SP_STATS(sti))
    for i = 1:tgt_nSP
      tgt_spfv(1+(st_count*f_len):((st_count+1)*f_len),i) = ...
        stats{sti}(target.fv(:, target.lin_sp == i));
    end
    for i = 1:src_nSP
      vi = source.validSuperpixels(i);
      source_spfv(1+(st_count*f_len):((st_count+1)*f_len),i) = ...
        stats{sti}(samples.fv(:, source.lin_sp == vi));
    end
    st_count = st_count + 1;
  end
end

end