function [inParams, ftsParams, outOptions] = InputAlgorithmParameters(filename)
%Input data from '.in' files on '../input' folder.
file_path = ['./../input/' filename '.in'];
fID = fopen(file_path, 'r');

text = fscanf(fID, '%s');
split = strsplit(text, ';');


for i = 1:(length(split) - 1)
    line = strsplit(split{i}, '=');
    val = line{end};
    
    switch line{1}
      case 'FOLDER'
        inParams.dataFolder = val;
      case 'SRC_NAME'
        inParams.sourceFile = val;
      case 'TGT_NAME'
        inParams.targetFile = val;

      case 'SAMPLE_METHOD'
        inParams.SAMPLE_METHOD = str2num(val);
      case 'COL_METHOD'
        inParams.COL_METHOD = str2num(val);
      case 'DIM_RED_METHOD'
        inParams.DIM_RED = str2num(val);
      case 'FEAT_SEL_METHOD'
        inParams.FEAT_SEL = str2num(val);
      case 'CL_CHANNELS'
        inParams.CL_CHANNELS = str2num(val);

      case 'N_SAMPLES'
        inParams.nSamples = str2num(val);
      case 'N_CLUSTERS'
        inParams.nClusters = str2num(val);
      case 'Kfs'
        inParams.Kfs = str2num(val);
      case 'Kis'
        inParams.Kis = str2num(val);
      case 'N_SUPERPIXELS'
        inParams.nSuperpixels = str2num(val);
      case 'LBL_MAJOR'
        inParams.LBL_MAJOR = str2num(val);

      case 'ON_FEATURES'
          feats_cells = strsplit(val, ',');
          feats = [];
          for f = 1:length(feats_cells)
              feats = [feats str2num(feats_cells{f})];
          end
          ftsParams.featsWeights = feats;
      case 'VEC_FEATURES'
          feats_cells = strsplit(val, ',');
          feats = [];
          for f = 1:length(feats_cells)
              feats = [feats str2num(feats_cells{f})];
          end
          ftsParams.vectorize = feats;
      case 'ON_STATS'
        feats_cells = strsplit(val, ',');
        feats = [];
        for f = 1:length(feats_cells)
              feats = [feats str2num(feats_cells{f})];
        end
        ftsParams.STATS = feats;  
      case 'FTS_STD_WINDOWSIZE'
          ftsParams.stdWS = str2num(val);
      case 'FTS_GABOR_WAVELENGHTS'
          ftsParams.gbWl = str2num(val);
      case 'FTS_GABOR_ORIENTATIONS' 
          ftsParams.gbOr = str2num(val);
      case 'FTS_DCT_WINDOWSIZE'
          ftsParams.dctWS = str2num(val);
      case 'FTS_DFT_WINDOWSIZE'
          ftsParams.dftWS = str2num(val);
      case 'FTS_SIFT_PATCHSIZE'
          ftsParams.siftPS = str2num(val);
      case 'FTS_SIFT_GRIDSPACING' 
          ftsParams.siftGs = str2num(val);
      case 'FTS_HARALICK_WINDOWSIZE'
          ftsParams.haraWS = str2num(val);

      case 'PLOT'
          outOptions.PLOT = str2num(val);
      case 'ANALYSIS'
          outOptions.ANALYSIS = str2num(val);
      case 'SAVE'
          outOptions.SAVE = str2num(val);

      otherwise
          warning([line{1} ' parameter not recognized']);
    end
end

fclose(fID);

%% Flow control flags
%COL_METHOD = 
% 0: Pixel + Matching
% 1: Pixel + KNN
% 2: Superpixel + Matching
% 3: Superpixel + KNN
% 4: Superpixel + Matching + scribbles
% 5: Superpixel + KNN + scribbles


inParams.SUPERPIXEL = inParams.COL_METHOD == 2 || inParams.COL_METHOD == 3 || ...
  inParams.COL_METHOD == 4 || inParams.COL_METHOD == 5;

inParams.COLOR_CLUSTERING = inParams.SAMPLE_METHOD == 2 || inParams.COL_METHOD == 1 || ...
  inParams.COL_METHOD == 3 || inParams.COL_METHOD == 5;

inParams.CLASSIFICATION = inParams.COL_METHOD == 1 || inParams.COL_METHOD == 3 || ...
  inParams.COL_METHOD == 5; 

%% Parameter consistency
  
if (inParams.SUPERPIXEL && inParams.SAMPLE_METHOD ~= 0);
    disp('Sampling method changed to full');
    inParams.SAMPLE_METHOD = 0;
end

end
