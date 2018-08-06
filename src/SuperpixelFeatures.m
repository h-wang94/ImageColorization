function [tgt_spfv, source_spfv] = SuperpixelFeatures(source, samples, target, SP_STATS)
%Builds the feature vector of superpixels by <insert way to concentrate information>

f_len = size(target.fv,1);

tgt_nSP = target.nSuperpixels;
src_nSP = length(source.validSuperpixels);

tgt_spfv = zeros(sum(SP_STATS)*f_len, tgt_nSP);
source_spfv = zeros(sum(SP_STATS)*f_len, src_nSP);

for i = 1:tgt_nSP
%   test = target.fv(:, target.lin_sp == i);
%   for d = 1:f_len
%     plot(test(d,:));
%     title(['SP ' num2str(i) ', f ' num2str(d)]);
%     ylim([0.0 1.0]);
%     drawnow;
% %     pause;
%    end

  tgt_spfv(1:f_len,i) = mean(target.fv(:, target.lin_sp == i), 2);
%   tgt_spfv(f_len+1:2*f_len,i) = std(target.fv(:, target.lin_sp == i), [], 2);
%   tgt_spfv(2*f_len+1:3*f_len,i) = median(target.fv(:, target.lin_sp == i), 2);    
end
for i = 1:src_nSP
  vi = source.validSuperpixels(i);
  source_spfv(1:f_len,i) = mean(samples.fv(:, source.lin_sp == vi), 2);
%   source_spfv(f_len+1:2*f_len,i) = std(samples.fv(:, source.lin_sp == vi), [], 2);
%   source_spfv(2*f_len+1:3*f_len,i) = median(samples.fv(:, source.lin_sp == vi), 2);
end

end