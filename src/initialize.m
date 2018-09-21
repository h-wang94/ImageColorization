%
clear all;
close all;

batch_folder = 'DissertationExperimentsSD/';

%%
% src_path = pwd;
input_folder = './../input/';
input_path = [input_folder batch_folder];

files = dir(input_path);

%Get filenames
for in_i = 3:length(files)
  list{in_i-2} = files(in_i).name;
end

%
for in_i = 61:length(list)
  batch_out = list{in_i};
  batch_out = [batch_out(1:end-3) '/'];
  mkdir(['./../results/' batch_out])
  copyfile([input_path list{in_i}], [input_folder 'default.in']);
  run('ColorizationSaulo');
  delete([input_folder 'default.in']);
end
