function clusters = ColorClustering(lab_img, nClusters, PLOT)
%TODO: Speed up !
ab = double(lab_img);
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab, nrows*ncols, 3); 

%% K-means clustering

% repeat the clustering 3 times to avoid local minima
[cluster_idx, C, ~, D] = kmeans(ab, nClusters, 'Distance', 'sqEuclidean', ...
                          'Replicates', 5);

if (PLOT)
    figure;
    subplot(1,2,1);
    hold on;
    for i = 1:nClusters
        idx = find(cluster_idx == i);
        scatter3(ab(idx,1), ab(idx,2), ab(idx,3), '.');
    end
    hold off;
    title('Color clusters');
end
   
%% Cluster statistics

minD = min(D');

stds = zeros(1, nClusters);
sizes = zeros(1, nClusters);

for i = 1:nClusters
    idx = cluster_idx == i;
    
    stds(i) = std(minD(idx));
    sizes(i) = sum(idx);
end

%% Output

clusters.idxs = cluster_idx;
clusters.centroids = C;
clusters.stds = stds;
clusters.cardin = sizes;

%% Compressed Recolorization

if (PLOT) 
    ab_out = zeros(nrows*ncols, 2);

    for i = 1:length(ab_out)
        ab_out(i,:) = normrnd(C(cluster_idx(i),2:3), 0);
    end

    % ab_out = reshape(ab_out, sz(1:2));
    im_out(:,:,1) = lab_img(:,:,1);
    im_out(:,:,2) = reshape(ab_out(:,1), [nrows, ncols]);
    im_out(:,:,3) = reshape(ab_out(:,2), [nrows, ncols]);

    subplot(1,2,2);
    imshow(lab2rgb(im_out));
    title(['Image colorized from ' num2str(nClusters) ' clusters (Centroid)']);
end

end

