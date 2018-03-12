%% -----------------------------------------------------------------
% Basic colorization algorithm.
%
% Modified by: Saulo Pereira
%
%TODO:
%-Estrutura de dados completa é cara para passar como argumento ?
%-COLOCAR GUARDA DE TAMANHO EM TODAS AS FUNÇÕES.
%-------------------------------------------------------------------
% clc
clear all; close all;

%%
tic;
GRAPH = true;
SAVE = false;
SAMPLE_METHOD = 1; % 0 for brute-force, 1, for jitter-sampling

%TODO: Input from file

%% Input data (source and target images)
src_name = 'dog2.jpg'
tgt_name = 'dog1.jpg'

imgs = LoadImages(src_name, tgt_name, '../data/');

target = {};
source = {};
target.image = imgs.target_image;
source.image = imgs.source_image;

tic;
%% Color space conversion
source.lab = rgb2lab(source.image);
target.lab = target.image;

if (GRAPH)
    figure; imshow([source.image source.lab]);
    title('RGB x Lab');
    
    abs = reshape(source.lab(:,:,2:3), size(source.lab,1)*size(source.lab,2), 2);
    figure; scatter(abs(:,1), abs(:,2), '.');
    title('Lab chrominance distribution');
end
    
%% Map luminance to target luminance
%TODO: verificar necessidade e algoritmo.
source.luminance = luminance_remap(source.lab, target.lab);
target.luminance = target.lab;

%% Source sampling
%TODO:
% v2: em espaco de cor
%https://www.mathworks.com/matlabcentral/fileexchange/54340-image-sampling-algorithms
samples = {}

[samples.idxs, samples_ab] = JitterSampleIndexes(source.lab, 256);
samples.ab = samples_ab(1:2,:);

if (GRAPH)
    imshow(source.image); hold on;
    scatter(samples.idxs(1,:), samples.idxs(2,:), '.');
    title('Samples from source');
    drawnow;
end

% if SAMPLE_METHOD == 0
%     tgt_color = image_colorization_brute_force(target, gsource, source);
% %     new_name = strcat('bf_', img_name);
% elseif SAMPLE_METHOD == 1
%     tgt_color = image_colorization_jitter_sampling(target, source, GRAPH);
% %     new_name = strcat('jitter_', img_name);
% end

%% Colorization

%Feature extraction:
target.fv = FeatureExtraction(target.image, true);
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