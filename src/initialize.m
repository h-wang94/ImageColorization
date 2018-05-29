%% -----------------------------------------------------------------
% Colorization by example prototype.
% Author: Saulo Pereira
%
%TODO:
%-Remover contagem de ties quando acabar a fase de debug.
%-------------------------------------------------------------------

% clc
clear all; close all;

%% Setup
%TODO: arquivo
GRAPH =             true;
SAVE =              false;
SAMPLE_METHOD =     2;  %0 = brute-force, 1 = jittered-sampling, 2 = clustered-sampling
COL_METHOD =        1;  %0 = "regression", 1 = "classification"

%Parameters: 
nSamples =          2^10;
nClusters =         5;
features =          [true true false false false false false];

src_name = 'beach1_r.jpg';
tgt_name = 'beach2_r.jpg';

%% Input data (source and target images)

[source.image, target.image] = LoadImages(src_name, tgt_name, '../data/');

%% Color space conversion
source.lab = rgb2lab(source.image);

if (GRAPH)
    abs = reshape(source.lab(:,:,2:3), size(source.lab,1)*size(source.lab,2), 2);
    figure(1); scatter(abs(:,1), abs(:,2), '.'); hold on
    title('Source Lab chrominance distribution');
    
    drawnow;
end

%% Map source luminance to target luminance
tgt_lab = rgb2lab(cat(3, target.image, target.image, target.image));
target.luminance = tgt_lab(:,:,1)/100;

% target.luminance = target.image;
source.luminance = luminance_remap(source.lab, target.luminance, src_name == tgt_name);
% source.luminance = source.lab(:,:,1)/100;

%% Clustering
% Performs the clustering for sampling and/or classification.

if (SAMPLE_METHOD == 2 || COL_METHOD == 1)
    disp('Source color clustering'); tic;
    clusters = ColorClustering(source.lab, nClusters, GRAPH);
    
    if (nClusters ~= length(clusters.cardin))
        disp('Number of clusters is inconsistent');
    end
    
    toc;
end

%% Source sampling
disp('Source image sampling'); tic;

switch SAMPLE_METHOD
    case 0
    %No sampling:
    [samples.idxs, samples.ab] = FullSampling(source.lab);

    case 1
    %Jittered sampling:
    [samples.idxs, samples_ab] = JitterSampleIndexes(source.lab, nSamples);
    samples.idxs = [samples.idxs(2,:); samples.idxs(1,:)];
    samples.lin_idxs = sub2ind(size(source.luminance), samples.idxs(1,:), samples.idxs(2,:))';
    samples.ab = samples_ab(2:3,:);
    
    case 2
    %Clustered sampling:
    samples = ClusteredSampling(source.lab, clusters, nClusters, nSamples);
    
    otherwise
    disp('Invalid SAMPLE_METHOD');
end
samples.sourceSize = size(source.luminance);

toc;

if (GRAPH)
    figure(2); imshow(source.image); title('Samples from source'); hold on;
    %Invert coordinates because it is a plot over an image.
    scatter(samples.idxs(2,:), samples.idxs(1,:), '.r'); hold off;
    
    figure(1); title('Lab chrominance distribution (total x sampled)');
    scatter(samples.ab(1,:), samples.ab(2,:), 6, 'r');

    drawnow;
end

%% Feature extraction:
disp('Feature extraction'); tic;

[target.fv, target.fv_w] = FeatureExtraction(target.luminance, features);
[samples.fv, samples.fv_w] = FeatureExtraction(source.luminance, features, samples.idxs);

toc;

%Feature space analysis
if (false)
    d_ab = pdist(samples.ab');
%     D_ab = squareform(D_ab);
%     d_fv = pdist((samples.fv.*repmat(samples.fv_w,1,length(samples.idxs)))' );
    d_fv = pdist(samples.fv');
    
    figure; scatter(d_fv, d_ab, '.');
    title('Relationship between distances in Feature and Color spaces');
    xlabel('Feature distance (not scaled)');
    ylabel('Color distance');
end
if(GRAPH && SAMPLE_METHOD == 2)
    figure(10); title('Classes in feature space'); hold on; 
    figure(11); imshow(source.luminance); title('Samples by class'); hold on;
    for i = 1:nClusters
        instances = find(samples.clusters == i);
        figure(10); scatter(samples.fv(1,instances), samples.fv(2,instances),'.');
        figure(11); scatter(samples.idxs(2,instances), samples.idxs(1,instances),'.');
    end
    figure(10); hold off;
	figure(11); hold off;
end

%% Color transfer:
disp('Color transfer'); tic

switch COL_METHOD
    case 0
    [tgt_lab, tiesIdx] = CopyClosestFeatureColor(source, samples, target, true);
    
    case 1
    [tgt_lab, tiesIdx] = CopyClosestFeatureInClassColor(samples, target, clusters);
    
    otherwise
    disp('Invalid COL_METHOD');
end
toc;
%% Color space reconversion
tgt_rgb = lab2rgb(tgt_lab);

%% Show results
figure;
imshow(tgt_rgb); hold on;
if(~isempty(tiesIdx))
    scatter(tiesIdx(2,:), tiesIdx(1,:), '.k'); hold off;
end
title('Colorized result (ties marked)');
% figure;
% imshow(tgt_rgb);

if (false)
    figure; title('Source x Target Lab chrominance distribution'); hold on;
    abs = reshape(source.lab(:,:,2:3), size(source.lab,1)*size(source.lab,2), 2);
    scatter(abs(:,1), abs(:,2), '.'); 
    abs = reshape(tgt_lab(:,:,2:3), size(tgt_lab,1)*size(tgt_lab,2), 2);
    scatter(abs(:,1), abs(:,2), 6, 'g'); hold off;
end

%% save images
% success = save_image(tgt_color, new_name, SAVE);
