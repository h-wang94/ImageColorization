%% -----------------------------------------------------------------
% Colorization by example.
% Modified by: Saulo Pereira
%
%TODO:
%-Estrutura de dados completa é cara para passar como argumento ?
%-CADA FEATURE DEVE TER UMA DISTÂNCIA DIFERENTE ?
% -os valores de distâncias devem estar entre 0-1 
% rotina colorclustering deve ser chamada nas porções de interesse
%-------------------------------------------------------------------

% clc
clear all; close all;
tic;

%% Setup
GRAPH = true;
SAVE = false;
SAMPLE_METHOD = 1;  %0 = brute-force, 1 = jittered-sampling, 2 = clustered-sampling
COL_METHOD = 1;     %0 = "regression", 1 = "classification"

%Parameters: 
nSamples = 2^10;
nClusters = 15;
features = [true true true];
% featuresWeights = 

%% Input data (source and target images)
src_name = 'beach2.jpg';
tgt_name = 'beach1.jpg';

target = {};
source = {};

[source.image, target.image] = LoadImages(src_name, tgt_name, '../data/');

tic;
%% Color space conversion
source.lab = rgb2lab(source.image);

if (GRAPH)
    figure(1); imshow([source.image source.lab]);
    title('RGB x Lab');
    
    abs = reshape(source.lab(:,:,2:3), size(source.lab,1)*size(source.lab,2), 2);
    figure(2); scatter(abs(:,1), abs(:,2), '.'); hold on
    title('Lab chrominance distribution');
    
    drawnow;
end

%% Map source luminance to target luminance
target.luminance = target.image;
source.luminance = luminance_remap(source.lab, target.luminance, src_name == tgt_name);

%% Source sampling
% TODO: colocar feature extraction aqui
samples = {};

switch SAMPLE_METHOD
    case 0
    [samples.idxs, samples.ab] = FullSampling(source.lab);
    samples.sourceSize = size(source.luminance);
    
    case 1
    %Jittered sampling:
    [samples.idxs, samples_ab] = JitterSampleIndexes(source.lab, nSamples);
    %TODO: inverter na propria funcao.
    samples.idxs = [samples.idxs(2,:); samples.idxs(1,:)];
    samples.ab = samples_ab(2:3,:);
    samples.sourceSize = size(source.luminance);
    
    case 2
    test = ClusteredSampling(source.lab, nClusters, nSamples);
    
    otherwise
        disp('Invalid SAMPLE_METHOD');
end

if (GRAPH)
    figure(3); imshow(source.image); hold on;
    %Invert coordinates because it is a plot over an image.
    scatter(samples.idxs(2,:), samples.idxs(1,:), '.');
    title('Samples from source');
    
    figure(2);
    scatter(samples.ab(1,:), samples.ab(2,:), 6, 'r');
    title('Lab chrominance distribution (total x sampled)');
    drawnow;
end

%% Colorization

%Feature extraction:
[target.fv, target.fv_w] = FeatureExtraction(target.luminance, features);
[samples.fv, samples.fv_w] = FeatureExtraction(source.luminance, features, samples.idxs);

%Colorization:
switch COL_METHOD
    case 0
    [tgt_color, tiesIdx] = CopyClosestFeatureColor(samples, target);
    
    case 1
    clusters = ColorClustering(source.lab, nClusters, GRAPH); 
    [tgt_color, tiesIdx] = CopyClosestFeatureColor(samples, target, clusters);
    
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
imshow(tgt_color); hold on;
if(length(tiesIdx))
    scatter(tiesIdx(2,:), tiesIdx(1,:), '.k');
end
title('Colorized result (ties marked)');

%% save images
% success = save_image(tgt_color, new_name, SAVE);
toc;