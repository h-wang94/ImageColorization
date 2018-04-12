
function [lab_out, tiesIdx] = CopyClosestFeatureInClassColor(samples, target, clusters)
%Assigns for each pixel the color of the closest pixel in feature space
%from the closest class (majority voting).

if nargin < 3
    colors = samples.ab;
else
    %Indexing
	l_idxs = sub2ind(samples.sourceSize, samples.idxs(1,:), samples.idxs(2,:));
    centroid_idx = clusters.idxs(l_idxs);
%     colors = clusters.centroids(centroid_idx,:)';
    
    colors = samples.ab;
end

lab_out = zeros(size(target.image));
lab_out(:,:,1) = target.luminance;

[r, c] = size(target.luminance);
tiesIdx = [];
for j = 1:c
    for i = 1:r
        idx = (i-1) + (j-1)*r + 1;
        [best_idxs, ties] = BestMatchesFS(target.fv(:, idx), target.fv_w, samples.fv);
        
        class_hip = centroid_idx(best_idxs);
        %if there is no majority, returns the first which is the closest.
        class_majority = mode(class_hip); 
        
        [major_instances, ~] = find(class_hip == class_majority);
        lab_out(i, j, 2:3) = colors(:, best_idxs(major_instances(1)));
        
        %% Debug:
        if (ties > 1)
            tiesIdx = [tiesIdx [i j]'];
        end
    end
end

lab_out(:,:,1) = lab_out(:,:,1)*100;


end

