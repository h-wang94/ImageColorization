%% -----------------------------------------------------------------
% Exemplar-based colorization algorithm
% Author: Saulo Pereira
%-------------------------------------------------------------------
input_file = 'default';
clearvars -except batch_out;

%% Input parameters
[IP, FP, OO] = InputAlgorithmParameters(input_file);
figs = GenerateFiguresList;

%% Input data (source and target images)
[source.image, target.image] = LoadImages(IP.sourceFile, IP.targetFile, IP.dataFolder);

%% Color space conversion
source.lab = rgb2lab(source.image);

if (OO.PLOT)
  ShowColorDistribution(source.image, source.lab);
end

%% Luminance Remapping (source to target)
tgt_lab = rgb2lab(cat(3, target.image, target.image, target.image));
target.luminance = tgt_lab(:,:,1)/100;

source.luminance = luminance_remap(source.lab, target.luminance, IP.sourceFile == IP.targetFile);

%% Color Clustering (Automatic Labeling)
% Performs the clustering for sampling and/or classification.
if (IP.COLOR_CLUSTERING)
  disp('Source color clustering'); tic;
  clusters = ColorClustering(source.lab, IP.nClusters, IP.CL_CHANNELS, OO.PLOT);

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
if (OO.PLOT && ~IP.SUPERPIXEL)
  figure; imshow(source.image); title('Samples from source'); hold on;
  %Invert coordinates because it is a plot over an image.
  scatter(samples.idxs(2,:), samples.idxs(1,:), '.r'); hold off;

  figure(figs.ColorDist); title('Lab chrominance distribution (total x sampled)');
  scatter(samples.ab(1,:), samples.ab(2,:), 6, 'r');

  drawnow;
end

%% Feature extraction
dataName = IP.sourceFile(1:end-2-3);

try
  disp('Feature loading (try)'); tic;

  load(['./../temp/' dataName '_full']);
  toc;
catch
  disp('Feature extraction'); tic;

  [target_fv, target_fvl] = FeatureExtraction(target.luminance, FP);
  [samples_fv, samples_fvl] = FeatureExtraction(source.luminance, FP);
  toc;

  save(['./../temp/' dataName '_full'], 'target_fv', 'samples_fv', 'target_fvl', 'samples_fvl');
end

%Source Sampling
idxs = sub2ind(size(source.luminance), samples.idxs(1,:), samples.idxs(2,:));
%Structs receive values
samples.fv = samples_fv(:,idxs);
target.fv = target_fv;
samples.fvl = samples_fvl;
target.fvl = target_fvl;

%Clear structured variables
clear target_fv samples_fv target_fvl samples_fvl;

%% Superpixel extraction
if (IP.SUPERPIXEL)
  disp('Superpixel extraction'); tic;

  try
    disp('Superpixel loading'); tic
    
    load(['./../temp/' dataName '_sps']);
    toc;
  catch
    disp('Loading failed. Recomputing superpixels.'); tic
    
    [src_sp, src_lin_sp, src_centrs, src_nSP] = ...
      SuperpixelExtraction(source.luminance, IP.nSuperpixels, 'turbo');
    [tgt_sp, tgt_lin_sp, tgt_centrs, tgt_nSP] = ...
      SuperpixelExtraction(target.image, IP.nSuperpixels, 'turbo');
    toc;
    
    save(['./../temp/' dataName '_sps'], 'src_sp', 'tgt_sp', 'src_lin_sp', 'tgt_lin_sp', ...
      'src_centrs', 'tgt_centrs', 'src_nSP', 'tgt_nSP');
  end
  source.sp = src_sp; target.sp = tgt_sp;
  source.lin_sp = src_lin_sp; target.lin_sp = tgt_lin_sp;
  source.sp_centroids = src_centrs; target.sp_centroids = tgt_centrs;
  source.nSuperpixels = src_nSP; target.nSuperpixels = tgt_nSP;
  clear src_sp tgt_sp src_lin_sp tgt_lin_sp src_centrs tgt_centrs src_nSP tgt_nSP
  
  if (OO.PLOT)
    %Show superpixels
    figure(figs.TargetSP); imshow(imoverlay(target.image, boundarymask(target.sp, 4), 'w')); 
    title('Target superpixels');
    figure(figs.SourceSP); imshow(imoverlay(source.image, boundarymask(source.sp, 4), 'w')); 
    title('Source superpixels');
  end

  %>Superpixel Labeling
  if (IP.COLOR_CLUSTERING)
    disp('Superpixel labeling'); tic;
    [source.sp_clusters, source.sp_chrom] = SuperpixelLabeling(source.lab, clusters.idxs, source.lin_sp, ...
      IP.LBL_MAJOR, OO.PLOT, source.sp, size(source.luminance));
    toc;
  end
     
  %> Superpixel Feature Aggregation
  disp('Superpixel feature aggregation'); tic;
  [target.fv_sp, source.fv_sp] = SuperpixelsFeatures(source, samples, target);
  
  target = rmfield(target, 'fv');
  samples = rmfield(samples, 'fv');
  toc;
  
  %> Saliency Feature Computation
  [ssi1, ssi2] = SaliencyFeature(source.luminance, source.sp, source.nSuperpixels);
  [tsi1, tsi2] = SaliencyFeature(target.luminance, target.sp, target.nSuperpixels);
  %Find unique superpixels indexes and concatenate their saliency values
  %onto the feature vector.
  [sp_idxs, src_idxs] = unique(source.lin_sp);
  [~, tgt_idxs] = unique(target.lin_sp);
  source.fv_sp = [source.fv_sp; ssi1(src_idxs)'; ssi2(src_idxs)'];
  target.fv_sp = [target.fv_sp; tsi1(tgt_idxs)'; tsi2(tgt_idxs)'];
  clear ssi1 ssi2 tti1 tti2;
end

%% Matching / Classification
disp('Feature matching / Classification in Feature Space'); tic;

%>Feature Space Distance Computation
PDs = CombinedPDists(source.fv_sp, target.fv_sp, FP.featsWeights);
neighbor_idxs = zeros(size(PDs));
neighbor_dists = zeros(size(PDs));
for i = 1:size(PDs,1)
  [neighbor_dists(i,:), neighbor_idxs(i,:)] = sort(PDs(i,:));
end
clear PDs

img_gen = {};
%>Matching:
[match_idxs, ~] = knnsearch(source.fv_sp', target.fv_sp');
img_gen{1, 1} = match_idxs;  img_gen{1,2} = 'match_idxs';

%>Classification:
%Classes assignment
neighbor_classes = source.sp_clusters(neighbor_idxs);

%kNN
labelsKNN = modeTies(neighbor_classes(:,1:IP.Kfs));
img_gen{2,1} = labelsKNN; img_gen{2,2} = 'labelsKNN';
%Predict Full
[labelsSPrF, labelsCPrF, doubtsPrF, scoresPrF, costsPrF] = ...
  PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, 0, IP.nClusters, clusters.mcCost);
img_gen{3,1} = labelsSPrF; img_gen{3,2} = 'labelsSPrF';
img_gen{4,1} = labelsCPrF; img_gen{4,2} = 'labelsCPrF';
%Predict Equality
[labelsSPrE, labelsCPrE, doubtsPrE, scoresPrE, costsPrE] = ...
  PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, IP.Kfs, IP.nClusters, clusters.mcCost);
img_gen{5,1} = labelsSPrE; img_gen{5,2} = 'labelsSPrE';
img_gen{6,1} = labelsCPrE; img_gen{6,2} = 'labelsCPrE';

if (OO.PLOT)
  imKNN = CreateLabeledImage(labelsKNN, target.sp, size(target.image));
  imSPF = CreateLabeledImage(labelsSPrF, target.sp, size(target.image));
  imDPF = CreateLabeledImage(labelsSPrF.*~doubtsPrF + -1*doubtsPrF, target.sp, size(target.image));
  imCPF = CreateLabeledImage(labelsCPrF, target.sp, size(target.image));
  imSPE = CreateLabeledImage(labelsSPrE, target.sp, size(target.image));
  imDPE = CreateLabeledImage(labelsSPrE.*~doubtsPrE + -1*doubtsPrE, target.sp, size(target.image));
  imCPE = CreateLabeledImage(labelsCPrE, target.sp, size(target.image));

  figure; imshow([imKNN zeros(size(imKNN)) zeros(size(imKNN));
                  imSPF imDPF imCPF;
                  imSPE imDPE imCPE], []); colormap jet;
  title('Locally assigned labels: [kNN - zeros - zeros] ; [Full > S D C] ; [Equality > S D C]');
  drawnow;
  
  clear imKNN imSPF imDPF imCPF imSPE imDPE imCPE;
end
 
toc;

%% Edge-Aware Labeling/Relabeling
disp('Edge-Aware Labeling/Relabeling'); tic;

try
    load (['./../temp/' dataName '_eaclusters']);
    figure(73); imshow(eaClustersImg, []); colormap 'jet'
catch
    [eaClusters, eaClustersImg] = EdgeAwareClustering(target);
    save (['./../temp/' dataName '_eaclusters'], 'eaClusters', 'eaClustersImg');
end

%Relabeling
relabelsKNN = EdgeAwareRelabeling(eaClusters, labelsKNN, []);
relabelsSPrE = EdgeAwareRelabeling(eaClusters, labelsSPrE, []);
relabelsCPrE = EdgeAwareRelabeling(eaClusters, labelsCPrE, []);
img_gen{7,1} = relabelsKNN; img_gen{7,2} = 'relabelsKNN';
img_gen{8,1} = relabelsSPrE; img_gen{8,2} = 'relabelsSPrE';
img_gen{9,1} = relabelsCPrE; img_gen{9,2} = 'relabelsCPrE';

%Costs Labeling
labelsEACPrE = EdgeAwareRelabeling(eaClusters, [], costsPrE);
img_gen{10,1} = labelsEACPrE; img_gen{10,2} = 'labelsEACPrE';
%Scores Labeling
labelsEASPrE = EdgeAwareRelabeling(eaClusters, [], scoresPrE);
img_gen{11,1} = labelsEASPrE; img_gen{11,2} = 'labelsEASPrE';


if (OO.PLOT)
  imKNN = CreateLabeledImage(relabelsKNN, target.sp, size(target.image));
  imSPE =  CreateLabeledImage(relabelsSPrE, target.sp, size(target.image));
  imCPE = CreateLabeledImage(relabelsCPrE, target.sp, size(target.image));  
  imEAPrEc = CreateLabeledImage(labelsEACPrE, target.sp, size(target.image));
  imEAPrEs = CreateLabeledImage(labelsEASPrE, target.sp, size(target.image));
  
  figure; imshow([imKNN imSPE imCPE;
                  zeros(size(imKNN)) imEAPrEs imEAPrEc], []); colormap jet;
  title('[Relabels> kNN - SE - CE] ; [EA Labels> - S C]');
  
  clear imKNN imSPE imCPE imEAPrEs imEAPrEc;
end

clear match_idxs labelsKNN labelsPrF labelsPrE;
clear relabelsKNN relabelsPrF relabelsPrE labelsEAPrFc labelsEAPrEc labelsEAPrFs labelsEAPrEs

%% Color transfer:
disp('Color transfer + Save'); tic

[tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelAvgScribble(source, target, img_gen{1,1});
    tgt_scribbled = lab2rgb(tgt_scribbled);
target.rgb = ColorPropagationLevin(tgt_scribbled, target.luminance, scribbles_mask);
imwrite(target.rgb, ['./../results/' batch_out dataName '_' img_gen{1,2} '.png'], 'png');

for i = 2:length(img_gen)
  [tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelFromClassScribble(source, target, ...
      neighbor_idxs, neighbor_classes, img_gen{i,1}, IP.Kfs);
  tgt_scribbled = lab2rgb(tgt_scribbled);
  target.rgb = ColorPropagationLevin(tgt_scribbled, target.luminance, scribbles_mask);
  
  imwrite(target.rgb, ['./../results/' batch_out dataName '_' img_gen{i,2} '.png'], 'png');
end

toc;