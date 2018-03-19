function samples = ClusteredSampling(lab_img, clusters , nClusters, nSamples)
%TODO> WIP

%TODO: incluir a chamada a colorclustering internamente.
if (nClusters ~= length(clusters.cardin))
    disp('Number of clusters is inconsistent');
end

%Samples per cluster
spc = (clusters.cardin/norm(clusters.cardin))*nSamples;

for i = 1:nClusters
    samps = rand(1, spc(i));
    
    [~, find_idxs] = find(clusters.idxs == i);
    
    t = find_idxs(samps)
end

% samples.idxs =
% samples.ab = 

end