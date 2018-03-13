%% -----------------------------------------------------------------
% Basic colorization algorithm.
%
% Modified by: Saulo Pereira
%
%TODO:
%-Estrutura de dados completa é cara para passar como argumento ?
%-CADA FEATURE DEVE TER UMA DISTÂNCIA DIFERENTE.
% -os valores de distâncias devem estar entre 0-1 
%-------------------------------------------------------------------

% clc
clear all; close all;
tic;

%% Setup
GRAPH = true;
SAVE = false;
SAMPLE_METHOD = 1; % 0 for brute-force, 1, for jitter-sampling

%TODO: include following parameters:
%-features weights
%-features activation
%-images paths/names
%-

%TODO: Input from file

%% Input data (source and target images)
src_name = 'landscape1.jpg';
tgt_name = 'landscape2.jpg';

target = {};
source = {};

[source.image, target.image] = LoadImages(src_name, tgt_name, '../data/');

tic;
%% Color space conversion
source.lab = rgb2lab(source.image);
% target.lab = target.image;

if (GRAPH)
    figure; imshow([source.image source.lab]);
    title('RGB x Lab');
    
    abs = reshape(source.lab(:,:,2:3), size(source.lab,1)*size(source.lab,2), 2);
    figure; scatter(abs(:,1), abs(:,2), '.');
    title('Lab chrominance distribution');
end

%% Map luminance to target luminance
%TODO: Corrigir (quebrando assumption dos valores de cor)
target.luminance = target.image;
source.luminance = luminance_remap(source.lab, target.luminance);

%% Source sampling
%TODO:
% v2: em espaco de cor
samples = {};

[samples.idxs, samples_ab] = JitterSampleIndexes(source.lab, 256);
samples.ab = samples_ab(2:3,:);

if (GRAPH)
    figure;
    imshow(source.image); hold on;
    scatter(samples.idxs(1,:), samples.idxs(2,:), '.');
    title('Samples from source');
    drawnow;
end

%% Colorization

%Feature extraction:
target.fv = FeatureExtraction(target.luminance, true);
samples.fv = FeatureExtraction(source.luminance, true, samples.idxs);

%Colorization:
tgt_color = transfer_sample_fv(samples, target);

%TODO: separar feature extraction da otimização.

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