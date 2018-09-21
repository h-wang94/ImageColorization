function [Accs, Ws] = crossValidationAccuracy(dataset, K, mcCost)

% FUNCTION DESCRIPTION %
% This function is for computing the k-Fold Cross Validation (k-FCV) and
% its average over 'Niteration' number of times.
%
% INPUT %
% 'dataset': This is a matrix of dimension n x m+1, where n is the number
% of data points and m is the number of features per data point. The last
% column of the matrix contain the label of each data point.
%
% OUTPUTS %
% 'Accuracy': This is the k-FCV accuracy, averaged over 'Niteration'
% 'NFeatures': This is the average number of features selected by the
% feature selection algorithm.
%
% Written By: Sujoy Paul, Jadavpur University, India
% Email: paul.sujoy.ju@gmail.com


%% Initializations
Niteration = 2;
sz = size(dataset);
norm_dataset = zeros(sz(1),sz(2)-1);

%% Dataset Normalization

for ii=1:sz(1)
  for jj=1:sz(2)-1
    if max(dataset(:,jj))-min(dataset(:,jj))==0
      norm_dataset(ii,jj)=0;
    else
      norm_dataset(ii,jj)=(1+9*(dataset(ii,jj)-min(dataset(:,jj)))/(max(dataset(:,jj))-min(dataset(:,jj))));
    end
  end
end
dataset=[norm_dataset dataset(:,end)];

%% Feature Optimization
cvArgs.mcCost = mcCost;
cvArgs.Knn = K;
set(0,'userdata',cvArgs);

Accs = zeros(1,Niteration);
Ws = zeros(sz(2)-1,Niteration);
for iAvg=1:Niteration
  disp(['Avg. Iteration No. ' num2str(iAvg) '...  Running...']);
%   c = cvpartition(dataset(:,end),'KFold',kFCV);
  c = cvpartition(dataset(:,end),'holdout',0.1);
  outArgs = crossval(@findAccuracy, dataset(:,1:end-1),dataset(:,end), 'partition',c);
  
  Accs(iAvg) = outArgs.PredictAcc;
  Ws(:,iAvg) = outArgs.featsWeights;
end

end