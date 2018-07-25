%% -----------------------------------------------------------------
% Exemplar-based colorization algorithm
% Author: Saulo Pereira
%-------------------------------------------------------------------

function initialize()
%   dir_el = dir('./../input');
%   
%   list = [];
%   for i = :
%     list = [list dir_el(i).name];
%   end
% 
%   batch_input = strsplit(list, '.in');
  
  batch_input = {'default'};

%TODO: remove function and use run.
  
  for i = 1:length(batch_input)
    close all;
    ColorizationPipeline(batch_input{i});
  end
end

function ColorizationPipeline(input_file)

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

  %% Superpixel extraction
  if (IP.SUPERPIXEL)
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
  if (IP.COLOR_CLUSTERING)
    disp('Source color clustering'); tic;
    clusters = ColorClustering( source.lab, IP.nClusters, IP.CL_CHANNELS, OO.PLOT);
   
    toc;
    
    if (IP.SUPERPIXEL)
      disp('Superpixel labeling'); tic;
      source.sp_clusters = SuperpixelLabeling(clusters.idxs, source.lin_sp, IP.LBL_MAJOR, ...
        OO.PLOT, source.sp, size(source.luminance));

      toc;
    end
    
  end
  
  %% Source sampling
	%TODO: REFACTORING
  disp('Source image sampling'); tic;

  if (IP.CLASSIFICATION && IP.SUPERPIXEL)
    disp('Class Rebalancing');

    %TODO: nao alterar o sp_clusters -> indexar utilizando o valid.
    [source.validSuperpixels, source.sp_clusters] = SuperpixelRebalSampling(source.sp_clusters, ...
      source.nSuperpixels, IP.nClusters);
  end
%   source.validSuperpixels = 1:source.nSuperpixels;
  
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
  try
    %TODO: criar mecanismo mais robusto para identificar necessidade de recalcular.
    %save no conjunto de parametros (FP). 
    load(['./../temp/' input_file(5:end)]);

  catch
    disp('Feature extraction'); tic;

    [target_fv, target_fvl] = FeatureExtraction(target.luminance, FP);
    [samples_fv, samples_fvl] = FeatureExtraction(source.luminance, FP);
    toc;
    
    save(['./../temp/' input_file(5:end)], 'target_fv', 'samples_fv', 'target_fvl', 'samples_fvl');
  end

  %Source Sampling
  idxs = sub2ind(size(source.luminance), samples.idxs(1,:), samples.idxs(2,:));
  samples.fv = samples_fv(:,idxs);
  target.fv = target_fv;
  samples.fvl = samples_fvl;
  target.fvl = target_fvl;
  
  %Clear structured variables
  clear target_fv samples_fv target_fvl samples_fvl;
  
  if (IP.SUPERPIXEL)
      disp('Superpixel feature averaging'); tic;
      [target.fv_sp, source.fv_sp] = SuperpixelFeatures(source, samples, target);

      target = rmfield(target, 'fv');
      samples = rmfield(samples, 'fv');

      toc;
  end

  % Principal components:
  if (IP.DIM_RED)
    disp('Dimensionality Reduction on Feature Space'); tic;

    if (~IP.SUPERPIXEL)
      [samples.fv, target.fv] = DimensionalityReduction(samples.fv, target.fv, IP.DIM_RED);
    else
      [source.fv_sp, target.fv_sp] = ...
        DimensionalityReduction(source.fv_sp, target.fv_sp, IP.DIM_RED);
    end
    toc;
  end

  %% Feature space analysis
  if (true)
    %>TEST: 180725----------------------------------------------
    %Cluster Distances
%     iiDists = clusterDist([source.fv_sp' source.sp_clusters']);
%     iiDistsN = ClusterDistNorm([source.fv_sp' source.sp_clusters']);

%     stem(iiDists(1,:)./iiDists(2,:), 'filled', 'MarkerSize', 3, 'LineStyle', 'none');
    [~,b] = FeatureCombinationSearchMedianColor(source, samples, target, clusters.mcCost, IP.nClusters);
    stem(b(1,:)./b(2,:), 'filled', 'MarkerSize', 3, 'LineStyle', 'none');
    pause;
    %-----------------------------------------------------------
    
    
    distsFSS = FeatureCombinationSearchMedianColor(source, samples, target, clusters.mcCost, IP.nClusters);
    [d_vals, d_idxs] = min(distsFSS);
    labels = -1*ones(1, source.nSuperpixels);
%     labels(source.validSuperpixels) = d_idxs;
    labels(d_vals > 1) = 0;
    
    [CI, C, ~, D] = kmeans(distsFSS', 4, 'Distance', 'sqEuclidean', ...
                          'Replicates', 5);
    
    tImg = CreateLabeledImage(labels, source.sp, size(source.luminance));
    
    error('Features Combination Test!');
  end
  
  if (IP.SUPERPIXEL && OO.ANALYSIS)
    %MOVE TO RELABELING
    [f_cluster, Cf, ~, Df] = kmeans(source.fv_sp', IP.nClusters, 'Distance', 'sqEuclidean', ...
                            'Replicates', 5);

    valid_f_cluster = zeros(source.nSuperpixels,1);
    valid_f_cluster(source.validSuperpixels) = f_cluster;
    valid_f_cluster(setdiff(1:source.nSuperpixels, source.validSuperpixels)) = -1;
    src_feat_labels = CreateLabeledImage(valid_f_cluster, source.sp, size(source.luminance));

    figure; imshow([src_col_labels src_feat_labels], []);
    title('Superpixel Labeling (left: colors, right: features)');
    colormap jet; drawnow;

    [f_cluster, Cf, ~, Df] = kmeans(target.fv_sp', IP.nClusters, 'Distance', 'sqEuclidean', ...
                            'Replicates', 5);

    tgt_feat_labels = CreateLabeledImage(f_cluster, target.sp, size(target.luminance));

    figure; imshow(tgt_feat_labels, []);
    title('Target Superpixel feature clustering');
    colormap jet; drawnow;

    target.fv_sp_labels = f_cluster;
  end

  if (OO.ANALYSIS && IP.COLOR_CLUSTERING)
    figure(figs.LabelsFS); title('Source: Labeled samples in feature space'); hold on; 
    figure(figs.LabelsImage); imshow(source.luminance); 
    title('Source: Labeled samples over image'); hold on;

    if (~IP.SUPERPIXEL)
      for i = 1:IP.nClusters
          instances = intersect(find(samples.clusters == i), source.validSuperpixels);
          figure(figs.LabelsFS); scatter(samples.fv(1,instances), samples.fv(2,instances),'.');
          figure(figs.LabelsImage); scatter(samples.idxs(2,instances), samples.idxs(1,instances),'.');
      end
    else
      for i = 1:IP.nClusters
        sp_instances = find(source.sp_clusters == i);
        figure(figs.LabelsFS); scatter(source.fv_sp(1,sp_instances), source.fv_sp(2,sp_instances), '.');
      end
    end
    figure(figs.LabelsFS); hold off;
    figure(figs.LabelsImage); hold off; 

    drawnow;
  end

  if (OO.ANALYSIS)
    figure; title('Target: Feature space distribution'); hold on;
    if (~IP.SUPERPIXEL) 
        scatter(target.fv(1,:), target.fv(2,:), '.k'); hold off;
    else
        scatter(target.fv_sp(1,:), target.fv_sp(2,:), '.k'); hold off;
    end
  end

  %% Matching / Classification
  disp('Feature matching / Classification in Feature Space'); tic;

  if (~IP.SUPERPIXEL && ~IP.CLASSIFICATION)
    [neighbor_idxs, neighbor_dists] = knnsearch(samples.fv', target.fv'); 
  elseif (~IP.SUPERPIXEL && IP.CLASSIFICATION)    
    [neighbor_idxs, neighbor_dists] = knnsearch(samples.fv', target.fv', 'K', IP.Kfs);
    neighbor_classes = samples.lin_idxs(neighbor_idxs,:);
    neighbor_classes = clusters.idxs(neighbor_classes);
    neighbor_classes = reshape(neighbor_classes, size(neighbor_idxs,1), size(neighbor_idxs,2));
    labels = mode(neighbor_classes,2);
  elseif (IP.SUPERPIXEL && ~IP.CLASSIFICATION)
    [neighbor_idxs, neighbor_dists] = knnsearch(source.fv_sp', target.fv_sp');
  elseif (IP.SUPERPIXEL && IP.CLASSIFICATION)    
    %TODO: Clean section.
    [neighbor_idxs, neighbor_dists] = knnsearch(source.fv_sp', target.fv_sp', ...
      'K', source.nSuperpixels); % Return all distances for further reference.
    neighbor_classes = source.sp_clusters(neighbor_idxs);
    
    [labels, scores] = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, IP.Kfs, IP.nClusters, ...
      clusters.mcCost);
    labels_m = modeTies(neighbor_classes(:,1:IP.Kfs));

    labeled_img = CreateLabeledImage(labels, target.sp, size(target.image));
    labeled_img_m = CreateLabeledImage(labels_m, target.sp, size(target.image));
    figure; imshow([labeled_img labeled_img_m], []); colormap jet;
    title('Predicted labels of each superpixel (left: NN posterior, right: NN mode)');
  end

  toc;

  %% Relabeling
  if (IP.SUPERPIXEL && IP.CLASSIFICATION && false)
    disp('Superpixel relabeling'); tic;

    %>TEST 180711:----------------------------
    %Linked to TEST 180710
    relabels = FeatureClustersRelabeling(target, labels);
    relabeled_img = CreateLabeledImage(relabels, target.sp, size(target.image));
    figure(100); imshow([labeled_img relabeled_img],[]); colormap jet;
    title('Before and after relabeling (Feat Cluster)');
    %-----------------------------------------
    
    %TODO: create purely spatial relabeling
    %Flow: FeatureClusterRelabeling -> SpatialRelabeling
%     relabels = ClassScoreSpatialRelabeling(target, IP.nClusters, IP.Kis, ...
%       neighbor_classes(:, 1:IP.Kfs));
    
  end

  %% Color transfer:
  disp('Color transfer'); tic

  switch IP.COL_METHOD
    case 0
      tgt_lab = CopyClosestFeatureColor(samples.ab, target, neighbor_idxs);
    case 1
      tgt_lab = CopyClosestFeatureInClassColor(samples.ab, target, neighbor_idxs, neighbor_classes, ...
        labels);
    case 2
      tgt_lab = CopyClosestSuperpixelAvgColor(source, target, neighbor_idxs);
    case 3
      tgt_lab = CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
        neighbor_classes, labels);
      
      %Relabeled
      tgt_lab_r = CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
        neighbor_classes, relabels);
    case 4
      [tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelAvgScribble(source, target, neighbor_idxs);
      tgt_scribbled = lab2rgb(tgt_scribbled);
      target.rgb = ColorPropagationLevin(tgt_scribbled, target.luminance, scribbles_mask);  
    case 5
      [tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelFromClassScribble(source, target, ...
        neighbor_idxs, neighbor_classes, labels);
      tgt_scribbled = lab2rgb(tgt_scribbled);
      target.rgb = ColorPropagationLevin(tgt_scribbled, target.luminance, scribbles_mask);
      
      %Relabeled
      if (exist('relabels'))
        [tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelFromClassScribble(source, target, ...
          neighbor_idxs, neighbor_classes, relabels);
        tgt_scribbled = lab2rgb(tgt_scribbled);
        target.rgb_r = ColorPropagationLevin(tgt_scribbled, target.luminance, scribbles_mask);
      end
    otherwise
      disp('Invalid COL_METHOD');
  end

  toc;

  % Color space reconversion
  if (~isfield(target, 'rgb'))
    target.rgb = lab2rgb(tgt_lab);
  end
  if (~isfield(target, 'rgb_r') && exist('tgt_lab_r'))
    target.rgb_r = lab2rgb(tgt_lab_r);
  end
  
  %% Show results
  if (OO.PLOT)
    figure; imshow(target.rgb);
    if(isfield(target, 'rgb_r'))
      figure; imshow(target.rgb_r); title('Relabeled');
    end
  end
  
  %% Save Output images
  if (OO.SAVE)
    disp('Saving output image');
    imwrite(target.rgb, ['./../results/' input_file '.png'], 'png');
    if(isfield(target, 'rgb_r'))
      imwrite(target.rgb_r, ['./../results/' input_file '_r.png'], 'png');
    end
  end

  %% Result Analysis
  
  if (OO.ANALYSIS)
    MatchingAnalysis(IP.COL_METHOD, figs, source, target, neighbor_idxs);
  end
  
end
