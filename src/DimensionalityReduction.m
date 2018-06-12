function [samples_fv_r, target_fv_r] = DimensionalityReduction(samples_fv, target_fv, REDUCTION)
%Changes the feature vector of both target and source according to a linear
%transformation (PCA/LDA).

    switch REDUCTION
        case 1
        %PCA
        [coeff, score, latent] = pca(samples_fv');

        PC_coeff = coeff(:, latent/latent(1) > 1e-2);
%         PC_coeff = coeff(:, 1:2);

        target_fv_r = (target_fv' * PC_coeff )';
        samples_fv_r = (samples_fv' * PC_coeff )';
        
        case 2
        %LDA
        [fv_r, W] = FDA(samples_fv, clusters.idxs(samples.lin_idxs));

        target_fv_r = W'*target_fv;
        samples_fv_r = W'*samples_fv;
        
        otherwise
            disp('Unrecognized dimensionality reduction technique');
    end

end