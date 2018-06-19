%% -----------------------------------------------------------------
% Exemplar-based colorization algorithm
% Author: Saulo Pereira
%-------------------------------------------------------------------

% clc
clear all; close all;

%% Input parameters
%TODO: default to argin
[IP, FP, OO] = InputAlgorithmParameters('default');
figs = GenerateFiguresList;

%% Input data (source and target images)
[source.image, target.image] = LoadImages(IP.sourceFile, IP.targetFile, IP.dataFolder);

%>Flow control
SUPERPIXEL = IP.SUPERPIXEL;
COLOR_CLUSTERING = IP.SAMPLE_METHOD == 2 || IP.COL_METHOD == 1 || ...
  IP.COL_METHOD == 3 || IP.COL_METHOD == 4;

%% Color space conversion
source.lab = rgb2lab(source.image);

if (OO.PLOT)
  abs = reshape(source.lab(:,:,2:3), size(source.lab,1)*size(source.lab,2), 2);
  figure(figs.ColorDist); scatter(abs(:,1), abs(:,2), '.'); hold on
  title('Source Lab chrominance distribution');

  drawnow;
end

%% Luminance Remapping (source to target)
tgt_lab = rgb2lab(cat(3, target.image, target.image, target.image));
target.luminance = tgt_lab(:,:,1)/100;

source.luminance = luminance_remap(source.lab, target.luminance, IP.sourceFile == IP.targetFile);

%% Superpixel extraction
if (SUPERPIXEL)
  disp('Superpixel extraction'); tic;

  %TODO: superpixels should be a struct inside both source and target.
  [source.sp, source.lin_sp, source.sp_centroids, source.nSuperpixels] = ...
    SuperpixelExtraction(source.luminance, IP.nSuperpixels);
  [target.sp, target.lin_sp, target.sp_centroids, target.nSuperpixels] = ...
    SuperpixelExtraction(target.image, IP.nSuperpixels);

  toc;

  if (OO.PLOT)
    figure(figs.TargetSP); imshow(imoverlay(target.image, boundarymask(target.sp, 4), 'w')); 
    title('Target superpixels');
    figure(figs.SourceSP); imshow(imoverlay(source.image, boundarymask(source.sp, 4), 'w')); 
    title('Source superpixels');
  end
end

%% Color Clustering (Automatic Labeling)
% Performs the clustering for sampling and/or classification.

if (COLOR_CLUSTERING)
  disp('Source color clustering'); tic;
  clusters = ColorClustering(source.lab, IP.nClusters, OO.PLOT);

  toc;
end

if (SUPERPIXEL && COLOR_CLUSTERING)
  disp('Superpixel labeling'); tic;
  source.sp_clusters = zeros(1, max(source.lin_sp));

  for i = 1:length(source.sp_clusters)
      source.sp_clusters(i) = mode(clusters.idxs(source.lin_sp == i));
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
  [samples.idxs, samples_ab] = JitterSampleIndexes(source.lab, IP.nSamples);
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

  figure(figs.ColorDist); title('Lab chrominance distribution (total x sampled)');
  scatter(samples.ab(1,:), samples.ab(2,:), 6, 'r');

  drawnow;
end

%% Feature extraction
disp('Feature extraction'); tic;

target.fv = FeatureExtraction(target.luminance, FP);
samples.fv = FeatureExtraction(source.luminance, FP, samples.idxs);

toc;

if (SUPERPIXEL)
    disp('Superpixel feature averaging'); tic;
    [target.fv_sp, source.fv_sp] = SuperpixelFeatures(source, samples, target);
    
    target = rmfield(target, 'fv');
    samples = rmfield(samples, 'fv');
    
    toc;
end

% Principal components:
if (IP.DIM_RED)
  disp('Dimensionality Reduction on Feature Space'); tic;

  if (~SUPERPIXEL)
      [samples.fv, target.fv] = DimensionalityReduction(samples.fv, target.fv, IP.DIM_RED);
  else
      [source.fv_sp, target.fv_sp] = ...
        DimensionalityReduction(source.fv_sp, target.fv_sp, IP.DIM_RED);
  end
  toc;

%   if(OO.PLOT)
%     samples.fv_pc = ( samples.fv' * PC_coeff )';
% 
%     figure; title('Source: Labeled samples in PC space'); hold on;
%     for i = 1:IP.nClusters
%         instances = find(samples.clusters == i);
%         scatter(samples.fv_pc(1,instances), samples.fv_pc(2,instances),'.');
%     end
%   end
end

%% Feature space analysis

if(OO.ANALYSIS && COLOR_CLUSTERING)
  figure(figs.LabelsFS); title('Source: Labeled samples in feature space'); hold on; 
  figure(figs.LabelsImage); imshow(source.luminance); 
  title('Source: Labeled samples over image'); hold on;

  switch IP.COL_METHOD
    case 1
    for i = 1:IP.nClusters
        instances = find(samples.clusters == i);
        figure(figs.LabelsFS); scatter(samples.fv(1,instances), samples.fv(2,instances),'.');
        figure(figs.LabelsImage); scatter(samples.idxs(2,instances), samples.idxs(1,instances),'.');
    end

    case 3
    for i = 1:IP.nClusters
      sp_instances = find(source.sp_clusters == i);
      for j = 1:IP.nSuperpixels
%             TODO: with centroids
%             instances = find(source.lin_sp == sp_instances(j));
      end    

    end
  end
  figure(figs.LabelsFS); hold off;
  figure(figs.LabelsImage); hold off; 

  drawnow;
end

if (OO.ANALYSIS && OO.PLOT)
  figure; title('Target: Feature space distribution'); hold on;
  if (~SUPERPIXEL) 
      scatter(target.fv(1,:), target.fv(2,:), '.k'); hold off;
  else
      scatter(target.fv_sp(1,:), target.fv_sp(2,:), '.k'); hold off;
  end
end

%% Matching / Classification

switch IP.COL_METHOD
  case 0
  [neighbor_idx, neighbor_dists] = knnsearch(samples.fv', target.fv'); 
  case 1
  [neighbor_idx, neighbor_dists] = knnsearch(samples.fv', target.fv', 'K', IP.Kfs);
  labels = samples.lin_idxs(neighbor_idx,:);
  labels = clusters.idxs(labels);
  labels = reshape(labels, size(neighbor_idx,1), size(neighbor_idx,2));
  case 2
  [neighbor_idx, neighbor_dists] = knnsearch(source.fv_sp', target.fv_sp');
  case 3
  [neighbor_idx, neighbor_dists] = knnsearch(source.fv_sp', target.fv_sp', 'K', IP.Kfs);
  labels = source.sp_clusters(neighbor_idx);
  case 4
  %Same as 3
    
  otherwise
  disp('AHH');
  
end

%% Relabeling
%TODO

%% Color transfer:
%TODO: separate color transfer from classification (labeling)
% better control with two different steps.
% avoid code repetition.

disp('Color transfer'); tic

switch IP.COL_METHOD
  case 0
  [tgt_lab, best_dists] = CopyClosestFeatureColor(samples, target);

  case 1
  [tgt_lab, tiesIdx, neighbors_list] = ...
    CopyClosestFeatureInClassColor(samples, target, clusters, IP.Kfs);

  case 2
  [tgt_lab, neighbors_list] = CopyClosestSuperpixelAvgColor(source, target);

  case 3
  [tgt_lab, neighbors_list] = CopyClosestSuperpixelFromClassAvgColor(source, target, IP.Kfs);

  case 4
  [tgt_lab, neighbors_list] = CopyClosestSuperpixelFromClassScribble(source, target, IP.Kfs);

  otherwise
  disp('Invalid COL_METHOD');
end

toc;

%% TEST> Class scores:
%Commit
ClassScoreSpatialRelabeling(source, target, IP.nClusters, IP.Kfs, IP.Kis, ...
  neighbors_list, tgt_lab, clusters);

%% Color space reconversion
target.rgb = lab2rgb(tgt_lab);

%% Show results
figure; imshow(target.rgb);
% error('Parei propositalmente aqui');

%% Save Output images
if (OO.SAVE)
    imwrite(target.rgb, ['./../results/' 'CM' num2str(IP.COL_METHOD) '_' IP.targetFile], 'png');
end

%% Result Analysis (TODO: create function)
MatchingAnalysis(IP.COL_METHOD, figs, source, target, neighbors_list);
