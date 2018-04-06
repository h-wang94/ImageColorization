function [lab_out, tiesIdx] = CopyClosestFeatureColor(samples, target, clusters)
%Color transfer technique that assigns for each target pixel the color the
%source pixel with smallest distance in feature space.

if nargin < 3
    colors = samples.ab;
else
    %Classification: colors are assigned according to class centroid.
	l_idxs = sub2ind(samples.sourceSize, samples.idxs(1,:), samples.idxs(2,:));
    centroid_idx = clusters.idxs(l_idxs);
    colors = clusters.centroids(centroid_idx,:)';
end

lab_out = zeros(size(target.image));
lab_out(:,:,1) = target.luminance;

[r, c] = size(target.luminance);
tiesIdx = [];
for j = 1:c
    for i = 1:r
        idx = (i-1) + (j-1)*r + 1;
        [best_idxs, ties] = BestMatchesFS(target.fv(:, idx), target.fv_w, samples.fv);

        lab_out(i, j, 2:3) = colors(:, best_idxs(1));
        
        %% Debug:
        if (ties > 1)
            tiesIdx = [tiesIdx [i j]'];
        end
    end
end

lab_out(:,:,1) = lab_out(:,:,1)*100;

end

