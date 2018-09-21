%Generate combination of input parameters
ParamList{1} =      2;
ParamList{2} =      5;
ParamList{3} =      0;
ParamList{4} =      0;
ParamList{5} =      2;	

ParamList{6} =			-1;
ParamList{7} =      {4, 5};
ParamList{8} =			{5, 7};
ParamList{9} =			-1;
ParamList{10} =     -1;
ParamList{11}	=     0.0;

ParamList{12} =     {'1,1,1,1,1,1,1,1,1', '1,1,1,1,1,1,1,0,1', '1,1,1,1,1,1,1,0,0', '1,1,0,1,1,1,1,0,0', '1,1,1,0,1,1,1,0,0'};
ParamList{13} =     'false, false, true, true, true, true, true, true';
ParamList{14} =     'true, true, true';

ParamList{15} =     -1;
ParamList{16} =     '2.^(1:0.2:2)';
ParamList{17} =     '0:-10:-170';
ParamList{18} =     -1;
ParamList{19} =     -1;
ParamList{20} =     8;
ParamList{21} =     1;
ParamList{22} =     -1;

ParamList{23} = 		'./../data/gupta/';
ParamList{24} = 		{'001_r.png', '002_r.png', '003_r.png', '004_r.png', '005_r.png', '006_r.png', '007_r.png', '008_r.png', '009_r.png', '010_r.png', 'I003_r.png'};
ParamList{25} = 		{'001_i.png', '002_i.png', '003_i.png', '004_i.png', '005_i.png', '006_i.png', '007_i.png', '008_i.png', '009_i.png', '010_i.png', 'I003_i.png'};

ParamList{26} =     false;
ParamList{27} =     false;
ParamList{28} =     true;

folder_name = 'DissertationExperiments/';


%% Name list
Names{1} = 'SAMPLE_METHOD'; 
Names{2} = 'COL_METHOD';
Names{3} = 'DIM_RED_METHOD';
Names{4} = 'FEAT_SEL_METHOD';
Names{5} = 'CL_CHANNELS';

Names{6} = 'N_SAMPLES';
Names{7} = 'N_CLUSTERS';
Names{8} = 'Kfs';	
Names{9} = 'Kis';
Names{10} = 'N_SUPERPIXELS';
Names{11} = 'LBL_MAJOR';

Names{12} = 'ON_FEATURES';
Names{13} = 'VEC_FEATURES';
Names{14} = 'ON_STATS';

Names{15} = 'FTS_STD_WINDOWSIZE';
Names{16} = 'FTS_GABOR_WAVELENGHTS';
Names{17} = 'FTS_GABOR_ORIENTATIONS';
Names{18} = 'FTS_DCT_WINDOWSIZE';
Names{19} = 'FTS_DFT_WINDOWSIZE';
Names{20} = 'FTS_SIFT_PATCHSIZE';
Names{21} = 'FTS_SIFT_GRIDSPACING';
Names{22} = 'FTS_HARALICK_WINDOWSIZE';

Names{23} = 'FOLDER';
Names{24} = 'SRC_NAME';
Names{25} = 'TGT_NAME';

Names{26} = 'PLOT';
Names{27} = 'ANALYSIS';
Names{28} = 'SAVE';

%% Generate combinations
numels = zeros(1, length(ParamList));
for nidx = 1:length(numels)
  if (iscell(ParamList{nidx }))
    numels(nidx) = length(ParamList{nidx}); 
  else
    numels(nidx) = 1;
  end
end

combs{1} = ones(1,length(numels));
combidx = 2;
over = 0;
while (~over)
  [combs{combidx}, over] = Increment(combs{combidx-1}, numels);
  combidx = combidx + 1;
end

%% Generate file for each combination
%Remove combinations of non-matching images
for c = 1:(length(combs) - 1)
  if (combs{c}(24) ~= combs{c}(25))
    combs{c} = [];
  end
end

mkdir(['./../input/' folder_name]);
for c = 1:(length(combs) - 1)
  if (isempty(combs{c}))
    continue;
  end;
  %Generate file name using the varying parameters
  fname = [];
  for i = length(numels):-1:1
    if (numels(i) > 1)
      %Write name of varying parameter and corresponding value.
      value = ParamList{i}{combs{c}(i)};
      if (isstring(value))
        fname = [fname Names{i}(1:3) value];
      else
        fname = [fname Names{i}(1:3) num2str(value)];
      end
    end
  end
  fname = [fname '.in'];
  fid = fopen(['./../input/' folder_name fname], 'w');
  
  %Write on file
  for i = 1:length(numels)
    value = ParamList{i};
    if (iscell(value))
      value = value{combs{c}(i)};
    end
    
    if (isstring(value))
      fprintf(fid, [Names{i} ' = ' value ';\n']);
    else
      fprintf(fid, [Names{i} ' = ' num2str(value) ';\n']);
    end
  end
  
  fclose(fid);
end