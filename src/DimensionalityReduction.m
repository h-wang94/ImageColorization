function DimensionalityReduction(samples, target, REDUCTION, PLOT)
%Changes the feature vector of both target and source according to a linear
%transformation (PCA/LDA).

    switch REDUCTION
        case 1
        %PCA
        [coeff, score, latent] = pca(samples.fv');

        PC_coeff = coeff(:, latent/latent(1) > 1e-2);
%         PC_coeff = coeff(:, 1:2);

        if(PLOT)
            samples.fv_pc = ( samples.fv' * PC_coeff )';
            
            figure; title('Source: Labeled samples in PC space'); hold on;
            for i = 1:IP.nClusters
                instances = find(samples.clusters == i);
                scatter(samples.fv_pc(1,instances), samples.fv_pc(2,instances),'.');
            end
        end

        target.fv = (target.fv' * PC_coeff )';
        samples.fv = (samples.fv' * PC_coeff )';
        
        case 2
        %LDA
        [fv_r, W] = FDA(samples.fv, clusters.idxs(samples.lin_idxs));

        if(PLOT)
            figure; title('Source: Labeled samples in LD space'); hold on;
            for i = 1:IP.nClusters
                instances = find(samples.clusters == i);
                scatter(fv_r(1,instances), fv_r(2,instances),'.');
            end
        end

        target.fv = W'*target.fv;
        samples.fv = W'*samples.fv;
        
        otherwise
            disp('Unrecognized dimensionality reduction technique');
    end

end