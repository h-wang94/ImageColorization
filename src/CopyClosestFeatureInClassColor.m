function [lab_out, tiesIdx, candidates] = CopyClosestFeatureInClassColor(samples, target, clusters, K)
%Performs a classification on each pixel and assigns either the color of 
%the the class centroid or the color of theclosest pixel in feature space
%from the class it was assigned during classification.
kK = 11;

%Indexing
l_idxs = sub2ind(samples.sourceSize, samples.idxs(1,:), samples.idxs(2,:));
centroid_idx = clusters.idxs(l_idxs);

COLORS_CENTROID = false;
if (COLORS_CENTROID)
    colors = clusters.centroids(centroid_idx,:)';
else
    colors = samples.ab;
end

lab_out = zeros(size(target.image));
lab_out(:,:,1) = target.luminance*100;

%% Transfer
%Weights test:
% ws = repmat(samples.fv_w, 1, size(samples.fv, 2));
% wt = repmat(target.fv_w, 1, size(target.fv, 2));
% [bsft, dists] = knnsearch((ws.*samples.fv)', (wt.*target.fv)', 'K', kK);

%best samples for each target
[bsft, dists] = knnsearch(samples.fv', target.fv', 'K', K);

%Find the cluster of each of the closest samples
tiesIdx = [];
candidates = [];
for i = 1:(size(target.image,1)*size(target.image,2))
    %Index conversion
    [r, c] = ind2sub(size(target.luminance),i);
    
    %closest samples indexes
    csi = samples.lin_idxs(bsft(i,:),:);
    class_hipts = clusters.idxs(csi);
    %if there is no majority, returns the first which is the closest.
    class_majority = mode(class_hipts);
    
    [majority_instances, ~] = find(class_hipts == class_majority);
    %% Analysis
    if (length(majority_instances) <= kK/2)
        tiesIdx = [tiesIdx [r c]'];
    end
    
    candidates = [candidates csi];
    
    %% Color assignment
    [~, color_idxs] = intersect(samples.lin_idxs, csi(majority_instances));
    lab_out(r, c, 2:3) = mean(colors(:, color_idxs), 2);

    % In case find returns duplicates, use the first index.
%     color_idx = find(csi(majority_instances(1)) == samples.lin_idxs);
%     lab_out(r, c, 2:3) = colors(:, color_idx(1));

end

end

