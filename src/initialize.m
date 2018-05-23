%% -----------------------------------------------------------------
% Colorization by example prototype.
% Author: Saulo Pereira
%
%TODO:
%-Estrutura de dados completa é cara para passar como argumento ?
%-Normalizar cada feature
%-Remover contagem de ties quando acabar a fase de debug.
%-Não repetir colorclustering
%-------------------------------------------------------------------

% clc
clear all; close all;

%% Setup
%TODO: arquivo
GRAPH =             false;
SAVE =              false;
SAMPLE_METHOD =     1;  %0 = brute-force, 1 = jittered-sampling, 2 = clustered-sampling
COL_METHOD =        0;  %0 = "regression", 1 = "classification"

%Parameters: 
nSamples =          2^10;
nClusters =         15;
features =          [false false false true true false false];

%% Input data (source and target images)
src_name = 'beach1_r.jpg';
tgt_name = 'beach2_r.jpg';

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
% tgt_lab = rgb2lab(cat(3, target.image, target.image, target.image));
% target.luminance = tgt_lab(:,:,1);

target.luminance = target.image;
source.luminance = luminance_remap(source.lab, target.luminance, src_name == tgt_name);

%% Source sampling
disp('Source image sampling'); tic;

switch SAMPLE_METHOD
    case 0
    [samples.idxs, samples.ab] = FullSampling(source.lab);

    case 1
    %Jittered sampling:
    [samples.idxs, samples_ab] = JitterSampleIndexes(source.lab, nSamples);
    %TODO: inverter na propria funcao.
    samples.idxs = [samples.idxs(2,:); samples.idxs(1,:)];
    samples.ab = samples_ab(2:3,:);
    
    case 2
    samples = ClusteredSampling(source.lab, nClusters, nSamples);
    
    otherwise
    disp('Invalid SAMPLE_METHOD');
end
samples.sourceSize = size(source.luminance);

if (GRAPH)
    figure(2); imshow(source.image); hold on;
    %Invert coordinates because it is a plot over an image.
    scatter(samples.idxs(2,:), samples.idxs(1,:), '.r');
    title('Samples from source');
    
    figure(1);
    scatter(samples.ab(1,:), samples.ab(2,:), 6, 'r');
    title('Lab chrominance distribution (total x sampled)');
    drawnow;
end

toc;
%% Feature extraction:
disp('Feature extraction'); tic;

[target.fv, target.fv_w] = FeatureExtraction(target.luminance, features);
[samples.fv, samples.fv_w] = FeatureExtraction(source.luminance, features, samples.idxs);

toc;

if (GRAPH)
    %Feature analysis
    d_ab = pdist(samples.ab');
%     D_ab = squareform(D_ab);
    d_fv = pdist(samples.fv');
    
    figure; scatter(d_fv, d_ab, '.');
    title('Relationship between distances in Feature and Color spaces');
end

%% Colorization:
disp('Color transfer'); tic

switch COL_METHOD
    case 0
    [tgt_lab, tiesIdx] = CopyClosestFeatureColor(samples, target);
    
    case 1
    clusters = ColorClustering(source.lab, nClusters, GRAPH);
    %TODO: melhorar forma de parametrizar esta funcao.
%     [tgt_lab, tiesIdx] = CopyClosestFeatureColor(samples, target, clusters);
    [tgt_lab, tiesIdx] = CopyClosestFeatureInClassColor(samples, target, clusters);
    
    otherwise
    disp('wtf');
end

toc;
%% Color space reconversion
tgt_rgb = lab2rgb(tgt_lab);

%% Show results
figure;
subplot(1,2,1); imshow(tgt_rgb); hold on;
if(~isempty(tiesIdx))
    scatter(tiesIdx(2,:), tiesIdx(1,:), '.k');
end
title('Colorized result (ties marked)');
subplot(1,2,2); imshow(tgt_rgb);

if (GRAPH)
    figure;
    abs = reshape(source.lab(:,:,2:3), size(source.lab,1)*size(source.lab,2), 2);
    scatter(abs(:,1), abs(:,2), '.'); hold on
    abs = reshape(tgt_lab(:,:,2:3), size(tgt_lab,1)*size(tgt_lab,2), 2);
    scatter(abs(:,1), abs(:,2), 6, 'g');
    title('Source x Target Lab chrominance distribution');
end

%% save images
% success = save_image(tgt_color, new_name, SAVE);
