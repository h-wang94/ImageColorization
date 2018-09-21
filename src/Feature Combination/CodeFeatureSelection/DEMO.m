%This is a DEMO code for the algorithm proposed in [1].
% [1] Sujoy Paul and Swagatam Das. "Simultaneous feature selection and 
%     weighting–An evolutionary multi-objective optimization approach." 
%     Pattern Recognition Letters 65 (2015): 51-59.

%% Name of the dataset to execute
DatasetName = 'wbc_dataset';

%% Load the dataset
cd('Datasets');
load(DatasetName);
cd('..');

%% Call function to compute average accuracy and number of features selected
[Accuracy, NFeatures] = crossValidationAccuracy(dataset);
disp(['Average Accuracy: ' num2str(Accuracy*100)]);
disp(['Average no. of Features Selected: ' num2str(NFeatures)]);
