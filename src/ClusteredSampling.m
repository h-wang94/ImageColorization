function [samples, clusters] = ClusteredSampling(lab_img, clusters, nClusters, nSamples)
%Sampling function based on density of clusters.
%Samples are taken from each cluster in quantities proportional to cluster
%standard deviation.
%Returns approximately nSamples.

%Samples per cluster
spc = (clusters.stds/sum(clusters.stds))*nSamples;
% spc = (clusters.cardin/sum(clusters.cardin))*nSamples;
spc = round(spc);
spc = spc + (spc == 0);

%Randomize the indexes to be taken as samples from each cluster. 
RandIdx = rand(nClusters, max(spc));
C = repmat(clusters.cardin', 1, max(spc));
SamplesIdx = ceil(RandIdx.*C);

%Converts cluster index to image index
lin_idxs = [];
sampleCluster = [];
for i = 1:nClusters
    samplesCi = SamplesIdx(i, 1:spc(i));
    
    [src_idxs, ~] = find(clusters.idxs == i);
    lin_idxs = [lin_idxs; src_idxs(samplesCi)];
    sampleCluster = [sampleCluster; i*ones(spc(i),1)]; 
end

%% Return values
sz = size(lab_img);
[r, c] = ind2sub(sz(1:2), lin_idxs);
samples.idxs = [r'; c'];
samples.lin_idxs = lin_idxs;
samples.clusters = sampleCluster;

a = lab_img(:,:,2);
b = lab_img(:,:,3);
samples.ab = [a(lin_idxs)'; b(lin_idxs)'] ;

end