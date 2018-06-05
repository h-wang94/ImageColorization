function [inParams, outOptions] = InputAlgorithmParameters(filename)
%Input data from '.in' files on '../input' folder.
file_path = ['./../input/' filename '.in'];
fID = fopen(file_path, 'r');

text = fscanf(fID, '%s');
split = strsplit(text, ';');

for i = 1:length(split)
    line = strsplit(split{i}, '=');
    val = line(end);
    
    switch i
        case 1
            inParams.SAMPLE_METHOD = str2num(val{1});
        case 2
            inParams.COL_METHOD = str2num(val{1});
        case 3
            inParams.nSamples = str2num(val{1});
        case 4
            inParams.nClusters = str2num(val{1});
        case 5
            feats_cells = strsplit(val{1}, ',');
            feats = [];
            for f = 1:length(feats_cells)
                feats = [feats str2num(feats_cells{f})];
            end
            inParams.features = feats;
        case 6
            inParams.dataFolder = val{1};
        case 7
            inParams.sourceFile = val{1};
        case 8
            inParams.targetFile = val{1};
        case 9
            outOptions.PLOT = str2num(val{1});
        case 10
            outOptions.ANALYSIS = str2num(val{1});
        case 11
            outOptions.SAVE = str2num(val{1});
    end
end


%% Parameter consistency

if (inParams.COL_METHOD == 2 && inParams.SAMPLE_METHOD ~= 0);
    disp('Sampling method changed to full');
    inParams.SAMPLE_METHOD = 0;
end

end
