function [tgt_lab_recol] = ClassScoreSpatialRelabeling(source, target, nClusters, Kfs, Kis, ...
  neighbors_list, tgt_lab, clusters)
%Relabeling based on classification scores and spatial coherence. 
%TODO: Only works in mode 3 so far.
%Parameters after (...) should be removed

%% Feature space class scores:
clScoresF = zeros(target.nSuperpixels, nClusters);
for i = 1:size(clScoresF,1)
  classes = source.sp_clusters(neighbors_list(i,:));
  clScoresF(i,:) = hist(classes, 1:nClusters);
end
clScoresF = clScoresF ./ repmat(sum(clScoresF, 2), 1, nClusters);
[delta_idx, ~] = max(clScoresF,[],2);
clScoresF_delta = (clScoresF == repmat(delta_idx, 1, nClusters));

%% Image space class scores:
[neighbor_sps, ctrd_dists] = ...
  knnsearch(target.sp_centroids', target.sp_centroids', 'K', Kis+1);
neighbor_sps = neighbor_sps(:,2:end);
ctrd_dists = 1./ctrd_dists(:,2:end);
ctrd_dists = ctrd_dists ./ repmat(sum(ctrd_dists, 2), 1, Kis);

% clScoresI = zeros(target.nSuperpixels, nClusters);
% for i = 1:size(clScoresI,1)
% %   neighbor_scores = clScoresF_delta(neighbor_sps(i,:)',:);
%   neighbor_scores = clScoresF(neighbor_sps(i,:)',:);
%   w = ctrd_dists(i,:);
%   clScoresI(i,:) = w*neighbor_scores;
% end
% clScoresI = clScoresI ./ repmat(sum(clScoresI, 2), 1, nClusters);

%% Combination
kLambda = 1.5;
clScores = clScoresF;
for rel_it = 1:3
    clScoresI = zeros(target.nSuperpixels, nClusters);
    for i = 1:size(clScoresI,1)
    %   neighbor_scores = clScoresF_delta(neighbor_sps(i,:)',:);
      neighbor_scores = clScores(neighbor_sps(i,:)',:);
      w = ctrd_dists(i,:);
      clScoresI(i,:) = w*neighbor_scores;
    end
    clScoresI = clScoresI ./ repmat(sum(clScoresI, 2), 1, nClusters);
    
  clScores = clScores + (kLambda / rel_it)*clScoresI;
  clScores = clScores ./ repmat(sum(clScores, 2), 1, nClusters);
end

%% Output
%Class changes
[~, clf] = max(clScoresF, [], 2);
[~, cl] = max(clScores, [], 2);

%Color changes
tgt_lab_recol = tgt_lab;
figure; imshow(lab2rgb(tgt_lab_recol));
title('Before relabeling');

mod = find(clf ~= cl);
for i = 1:length(mod)
  mod_idx = mod(i);
%   classes = source.sp_clusters(neighbors_list(mod_idx,:));
%   [~, majority_instances] = find(classes == cl(mod_idx));
% 
%   if (majority_instances == [])
%     
%   end
%   
%   %Matching superpixels ROI masks
  tgt_mask = (target.sp==mod_idx);
%   src_mask = (source.sp==neighbors_list(mod_idx, majority_instances(1)));
% 
%   %Prototype color transfer (Superpixel average)
%   for c = 2:3
%       mask_c = source.lab(:,:,c).*src_mask;
%       avg_sp = sum(sum(mask_c))/length(find(src_mask));
%       tgt_lab2(:,:,c) = tgt_lab2(:,:,c) + avg_sp*tgt_mask;
%   end
  for c = 2:3
    tgt_lab_recol(:,:,c) = tgt_lab_recol(:,:,c).*~tgt_mask + clusters.centroids(cl(mod_idx),c)*tgt_mask;
  end
end

figure;imshow(lab2rgb(tgt_lab_recol));
title('After relabeling');


end

