%% -----------------------------------------------------------------
% Exemplar-based colorization algorithm
% Author: Saulo Pereira
%
%-------------------------------------------------------------------

% clc
clear all; close all;

%% Input
%TODO: argin
[IP, OO] = InputAlgorithmParameters('default');

%Figure list:
% figs.
fColorDist = 50;
fLabelsFS = 51;
fLabelsImage = 52;
fAnalysisInput = 53;
fCandidatesImage = 54;
%fCandidatesFS = 55;

%% Input data (source and target images)
[source.image, target.image] = LoadImages(IP.sourceFile, IP.targetFile, IP.dataFolder);

%% Color space conversion
source.lab = rgb2lab(source.image);

if (OO.PLOT)
    abs = reshape(source.lab(:,:,2:3), size(source.lab,1)*size(source.lab,2), 2);
    figure(fColorDist); scatter(abs(:,1), abs(:,2), '.'); hold on
    title('Source Lab chrominance distribution');
    
    drawnow;
end

%% Map source luminance to target luminance
tgt_lab = rgb2lab(cat(3, target.image, target.image, target.image));
target.luminance = tgt_lab(:,:,1)/100;

% target.luminance = target.image;
source.luminance = luminance_remap(source.lab, target.luminance, IP.sourceFile == IP.targetFile);
% source.luminance = source.lab(:,:,1)/100;

%% Clustering
% Performs the clustering for sampling and/or classification.

if (IP.SAMPLE_METHOD == 2 || IP.COL_METHOD == 1)
    disp('Source color clustering'); tic;
    clusters = ColorClustering(source.lab, IP.nClusters, OO.PLOT);
    
    if (IP.nClusters ~= length(clusters.cardin))
        disp('Number of clusters is inconsistent');
    end
    
    toc;
end

%% Source sampling
disp('Source image sampling'); tic;

switch IP.SAMPLE_METHOD
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
    samples = ClusteredSampling(source.lab, clusters, IP.nClusters, IP.nSamples);
    
    otherwise
    disp('Invalid SAMPLE_METHOD');
end
samples.sourceSize = size(source.luminance);

toc;

if (OO.PLOT)
    figure; imshow(source.image); title('Samples from source'); hold on;
    %Invert coordinates because it is a plot over an image.
    scatter(samples.idxs(2,:), samples.idxs(1,:), '.r'); hold off;
    
    figure(fColorDist); title('Lab chrominance distribution (total x sampled)');
    scatter(samples.ab(1,:), samples.ab(2,:), 6, 'r');

    drawnow;
end

%% Feature extraction:
disp('Feature extraction'); tic;

[target.fv, target.fv_w] = FeatureExtraction(target.luminance, IP.features);
[samples.fv, samples.fv_w] = FeatureExtraction(source.luminance, IP.features, samples.idxs);

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
if(OO.ANALYSIS && IP.SAMPLE_METHOD == 2)
    figure(fLabelsFS); title('Source: Labeled samples in feature space'); hold on; 
    figure(fLabelsImage); imshow(source.luminance); title('Source: Labeled samples over image'); hold on;
    for i = 1:IP.nClusters
        instances = find(samples.clusters == i);
        figure(fLabelsFS); scatter(samples.fv(1,instances), samples.fv(2,instances),'.');
        figure(fLabelsImage); scatter(samples.idxs(2,instances), samples.idxs(1,instances),'.');
    end
    figure(fLabelsFS); hold off;
	figure(fLabelsImage); hold off;
end

if (OO.PLOT)
    figure; title('Target: Feature space distribution'); hold on;
    scatter(target.fv(1,:), target.fv(2,:), '.k'); hold off;
end

drawnow;

%% Principal components (TEST):
% disp('Extract principal components');
if(false)
    %PCA
    [coeff, score, latent] = pca(samples.fv');

    PC_coeff = coeff(:, latent/latent(1) > 1e-2);
    % PC_coeff = coeff(:, 1:2);

    samples.fv_pc = ( samples.fv' * PC_coeff )';
        
    fLabelsRDFS = 100;
    if(OO.PLOT)
        figure(fLabelsRDFS); title('Source: Labeled samples in PC space'); hold on;
        for i = 1:IP.nClusters
            instances = find(samples.clusters == i);
            scatter(samples.fv_pc(1,instances), samples.fv_pc(2,instances),'.');
        end
        figure(fLabelsRDFS); hold off;
    end
    
    target.fv = (target.fv' * PC_coeff )';
    samples.fv = (samples.fv' * PC_coeff )';
end
%LDA
if (false)
    [fv_r, W] = FDA(samples.fv, clusters.idxs(samples.lin_idxs));

    if(true)
        figure(fLabelsRDFS); title('Source: Labeled samples in LDA space'); hold on;
        for i = 1:IP.nClusters
            instances = find(samples.clusters == i);
            scatter(fv_r(1,instances), fv_r(2,instances),'.');
        end
        figure(fLabelsRDFS); hold off;
    end

    target.fv = W'*target.fv;
    samples.fv = W'*samples.fv;
end
    
%% Color transfer:
disp('Color transfer'); tic

switch IP.COL_METHOD
    case 0
    [tgt_lab, tiesIdx] = CopyClosestFeatureColor(source, samples, target, true);
    
    case 1
    [tgt_lab, tiesIdx, cddt_list] = CopyClosestFeatureInClassColor(samples, target, clusters);
    
    otherwise
    disp('Invalid COL_METHOD');
end
toc;
%% Color space reconversion
tgt_rgb = lab2rgb(tgt_lab);

%% Show results

if (OO.PLOT)
    figure; imshow(tgt_rgb); title('Colorized result (ties marked)'); hold on;
    if(~isempty(tiesIdx))
        scatter(tiesIdx(2,:), tiesIdx(1,:), '.k'); hold off;
    end
end

%% Analysis

%Generate candidate source image
figure(fCandidatesImage);
imshow(source.image); title('Source candidates');

%Generate cursor input image
fig = figure(fAnalysisInput); 
imshow(tgt_rgb); title('Colorized result (index)');
datacursormode on;
dcm_obj = datacursormode(fig);

%Arguments for Analysis Tool
AnalysisArguments.sourceSize = size(source.luminance);
% AnalysisArguments.sourceImage = source.image;
AnalysisArguments.cddt_list = cddt_list;
AnalysisArguments.targetSize = size(target.luminance);
AnalysisArguments.targetFS = target.fv;
AnalysisArguments.fCandidatesImage = fCandidatesImage;
AnalysisArguments.fCandidatesFS = fLabelsFS;
set(0,'userdata',AnalysisArguments);

%Overwrite update function
set(dcm_obj, 'UpdateFcn', @MyAnalysisTool)

%% Save Output images

% imwrite()
