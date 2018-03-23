function samples = ClusteredSampling(lab_img, nClusters, nSamples)
%Sampling function based on density of clusters.
%Samples are taken from each cluster in quantities proportional to cluster
%standard deviation.
%Returns approximately nSamples.

clusters = ColorClustering(lab_img, nClusters, false);
if (nClusters ~= length(clusters.cardin))
    disp('Number of clusters is inconsistent');
end

%Samples per cluster
spc = (clusters.stds/sum(clusters.stds))*nSamples;
% spc = (clusters.cardin/sum(clusters.cardin))*nSamples;
spc = round(spc);
spc = spc + (spc == 0);

%Randomize the index to be taken as samples from each cluster. 
RandIdx = rand(nClusters, max(spc));
C = repmat(clusters.cardin', 1, max(spc));
SamplesIdx = round(RandIdx.*C);

%Converts cluster index to image index
lin_idxs = [];
for i = 1:nClusters
    samplesCi = SamplesIdx(i, 1:spc(i));
    
    [src_idxs, ~] = find(clusters.idxs == i);
    lin_idxs = [lin_idxs; src_idxs(samplesCi)];
end

%assigns the return values
sz = size(lab_img);
[r, c] = ind2sub(sz(1:2), lin_idxs);
samples.idxs = [r'; c'];

a = lab_img(:,:,2);
b = lab_img(:,:,3);
samples.ab = [a(lin_idxs)'; b(lin_idxs)'] ;

end