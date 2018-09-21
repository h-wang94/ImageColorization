function [samples_fv_r, target_fv_r, new_fvl] = DimensionalityReduction(samples_fv, target_fv, ...
  fvLens, REDUCTION)
%Changes the feature vector of both target and source according to a linear
%transformation (PCA/LDA).

DR_FEATURES = [false, false, false, false, false, true, false];
new_fvl = fvLens;
blockSize = sum(fvLens);
nSTATS = 3;

%Divide the STATS blocks:
for i = 1:nSTATS
  target_blocks{i} = target_fv(((i-1)*blockSize+1):i*blockSize,:);
  sample_blocks{i} = samples_fv(((i-1)*blockSize+1):i*blockSize,:);
end

for i = 1:length(DR_FEATURES)
  if (DR_FEATURES(i))
    single_act = zeros(1,length(DR_FEATURES));
    single_act(i) = 1;
    
    feat_sub = FeatureSubset(single_act, fvLens, nSTATS);
    [coeff, score, latent] = pca(samples_fv(feat_sub,:)');

    %Keep only the dimensions that have variances at least 1% of the
    %maximum. (Completes for block division).
    aboveThresVars = sum(latent/latent(1) > 1e-2);
    if (mod(aboveThresVars,nSTATS) ~= 0)
      aboveThresVars = aboveThresVars + (nSTATS - mod(aboveThresVars,nSTATS));
    end
    PC_coeff = coeff(:, 1:aboveThresVars);

    %Compute transformation
    target_transf = (target_fv(feat_sub,:)' * PC_coeff )';
    samples_transf = (samples_fv(feat_sub,:)' * PC_coeff )';
    
    %Replace the old block with the mapped one.
    for b = 1:nSTATS
      %Repeat for each block
      
      
    end
    blocks = feat_sub(2:end) - feat_sub(1:end-1);
    end_idxs = find(blocks > 1);
        
    %Update the feature length to reduced one.
    new_fvl(i) = aboveThresVars/nSTATS;    
  else
    %Receives the full block
%     target_fv_r = 
%     samples_fv_r = 
  end
  
  %Concatenate the blocks
  
  
end

% switch REDUCTION
%   case 1
%   %PCA
%   [coeff, score, latent] = pca(samples_fv');
% 
%   PC_coeff = coeff(:, latent/latent(1) > 1e-2);
% %         PC_coeff = coeff(:, 1:2);
% 
%   target_fv_r = (target_fv' * PC_coeff )';
%   samples_fv_r = (samples_fv' * PC_coeff )';
% 
%   case 2
% %     %LDA
% %     [fv_r, W] = FDA(samples_fv, clusters.idxs);%(samples.lin_idxs));
% % 
% %     target_fv_r = W'*target_fv;
% %     samples_fv_r = W'*samples_fv;
% 
%   otherwise
%     disp('Unrecognized dimensionality reduction technique');
% end

end