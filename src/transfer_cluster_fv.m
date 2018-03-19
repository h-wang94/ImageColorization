function transferred = transfer_cluster_fv(samples, target, clusters)
%
    target_fv = target.fv;
    transferred = zeros(size(target.image));
    transferred(:, :, 1) = target.luminance;
    
    [c_y, c_x] = size(target.luminance);
    occur = [];
    for j = 1:c_x
        for i = 1:c_y
            idx = (i-1) + (j-1)*c_y + 1;
            best_idx = compute_best_match(target_fv(:, idx), samples.fv);
            
            occur = [occur clusters.idxs(best_idx)];
%             transferred(i, j, 2:3) = normrnd(clusters.centroids(clusters.idxs(best_idx),:),...
%                                              clusters.stds(clusters.idxs(best_idx))); 
            transferred(i, j, 2:3) = clusters.centroids(clusters.idxs(best_idx),:);
        end
    end
    
    transferred(:,:,1) = transferred(:,:,1)*100;
end