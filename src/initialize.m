%% -----------------------------------------------------------------
% Basic colorization algorithm.
%
% Modified by: Saulo Pereira
%
%TODO:
%-Estrutura de dados completa é cara para passar como argumento ?
%-CADA FEATURE DEVE TER UMA DISTÂNCIA DIFERENTE.
% -os valores de distâncias devem estar entre 0-1 
% rotina colorclustering deve ser chamada nas porções de interesse
%-------------------------------------------------------------------

% clc
clear all; close all;
tic;

%% Setup
GRAPH = true;
SAVE = false;
SAMPLE_METHOD = 1; %0 = brute-force, 1 = jitter-sampling, 2 = clustered-sampling
COL_METHOD = 1; %0 = regression, 1 = classification

%TODO: include following parameters:
%-features weights
%-features activation
%-images paths/names

%TODO: Input from file
nSamples = 512;
nClusters = 15;

%% Input data (source and target images)
src_name = 'land1.jpg';
tgt_name = 'land1.jpg';

target = {};
source = {};

[source.image, target.image] = LoadImages(src_name, tgt_name, '../data/');

tic;
%% Color space conversion
source.lab = rgb2lab(source.image);
% target.lab = target.image;

if (GRAPH)
    figure(1); imshow([source.image source.lab]);
    title('RGB x Lab');
    
    abs = reshape(source.lab(:,:,2:3), size(source.lab,1)*size(source.lab,2), 2);
    figure(2); scatter(abs(:,1), abs(:,2), '.'); hold on
    title('Lab chrominance distribution');
    
    drawnow;
end

%% Map luminance to target luminance
%TODO: Corrigir (quebrando assumption dos valores de cor)
target.luminance = target.image;
source.luminance = luminance_remap(source.lab, target.luminance);

%% Source sampling
samples = {};

switch SAMPLE_METHOD
    case 0
    %Fully sampled
    sz = size(source.luminance);
    [i1, i2] = ind2sub(sz, 1:sz(1)*sz(2));
    samples.idxs = [i1; i2];
    samples_a = source.lab(:,:,2);
    samples_b = source.lab(:,:,3);
    samples.ab = [samples_a(:)'; samples_b(:)']; 

    case 1
    %Jittered sampling:
    [samples.idxs, samples_ab] = JitterSampleIndexes(source.lab, nSamples);
    samples.ab = samples_ab(2:3,:);
    
    case 2
%     clusters = ColorClustering(source.lab, nClusters, GRAPH);
%     test = ClusteredSampling(source.lab, clusters, nClusters, nSamples);
    
    otherwise
        disp('Invalid SAMPLE_METHOD');
end

if (GRAPH)
    figure(3); imshow(source.image); hold on;
    scatter(samples.idxs(1,:), samples.idxs(2,:), '.');
    title('Samples from source');
    
    figure(2);
    scatter(samples.ab(1,:), samples.ab(2,:), 6, 'r');
    title('Lab chrominance distribution (total x sampled)');
    drawnow;
end

%% Colorization

%Feature extraction:
target.fv = FeatureExtraction(target.luminance, true);
samples.fv = FeatureExtraction(source.luminance, true, samples.idxs);

%Colorization:
switch COL_METHOD
    case 0
    tgt_color = transfer_sample_fv(samples, target);
    
    case 1
    clusters = ColorClustering(source.lab, nClusters, GRAPH); 
    tgt_color = transfer_cluster_fv(samples, target, clusters);
    
    otherwise
    disp('wtf');
end
    
%% Color space reconversion
tgt_color = lab2rgb(tgt_color);

toc;
%% Show results
% figure;
% imshow(source.image)
% title('Source image');

figure;
imshow(tgt_color);
title('Colorized result');

%% save images
% success = save_image(tgt_color, new_name, SAVE);
toc;