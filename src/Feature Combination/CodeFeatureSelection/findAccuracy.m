%% Feature Selection and Classification Accuracy
function outArgs = findAccuracy(xtrain,ytrain,xtest,ytest)

% FUNCTION DESCRIPTION %
% This function is to compute the classification accuracy and the number 
% of selected features by the feature selection algorithm for the given
% training and testing dataset.
% INPUTS %%%%%%%%%%%%%%%
% 'xtrain': This is a matrix of size Tr x m, Tr being the no. of data
% points for training and m being the no. of features per data point.
% 'ytrain': This is a vector of size Tr x 1, containing the labels of each
% data point in the 'xtrain' matrix.
% 'xtest': This is a matrix of size Te x m, Te being the no. of data points
% for testing. 
% 'ytest': This is a vector of size Te x 1, containing the labels of each
% data point in the 'xtest' matrix.
% OUTPUT %%%%%%%%%%%%%%
% 'f' is a vector is size 1 x 2, s.t. f(1) contains the accuracy and f(2)
% contains the number of features selected by the feature selection method.
% Written By: Sujoy Paul, Jadavpur University, India
% Email: paul.sujoy.ju@gmail.com

%% Feature Selection
% Initializations of the parameters
Obj=2;                                 % Number of objective functions to be minimized
neighbours=20;                         % Number of nearest weights
idealpoint=inf*ones(1,Obj);            % Ideal Point
max_iteration=100;                     % Maximum number of iterations
max_val=10;                            % Maximum boundary of the search space
min_val=0;                             % Minimum boundary of the search space
F=0.7;                                 % Weighing Factor (DE)
Cr=0.95;                               % Crossover Rate (DE)

ldata = [xtrain ytrain];               % Learning Data
tdata = [xtest ytest];                 % Training Data
s=size(ldata);
D=s(2)-1;                              % Dimension of the search space
population=5*D;   %10                  % Number of subproblems

parent_fitness=zeros(population,Obj);
child_fitness=zeros(population,Obj);
% rand('state',40);
%Input arguments
udata = get(0, 'userdata');
kK = udata.Knn;

% Computing Inter and Intra Class distance vector (without weights)
Dist = clusterDist(ldata, udata.mcCost);
% Dist = ClusterDistNormMedian(ldata, udata.mcCost, 'pdist');

% Initialization of the weights of MOEA/D
weight=zeros(population,Obj);
for i=0:population-1
  weight(i+1,1)=i/(population-1);
  weight(i+1,2)=1-i/(population-1);
end

% Calculation of distance matrix and neighbouring weight vectors of MOEA/D
Neighbour_indices=zeros(population,neighbours);        % Neighbour indices to the each weight vector
Distance = squareform(pdist(weight));   % Inter weight distance matrix
for i=1:length(weight)
  [~, indices]=sort(Distance(i,:));
  Neighbour_indices(i,:)=indices(1:neighbours);
end

% Random initialization of the population
v=zeros(population,D);
u=zeros(population,D);
x=round(min_val+rand(population,D)*(max_val-min_val));     % Population of MOEA/D

% Function evaluation of the population and update the ideal point
for i=1:population
  parent_fitness(i,:)=weightedClusterDist(x(i,:),Dist);  % Compute objective function: Intra and Inter class distance
end
for i=1:Obj
  idealpoint(1,i)=min(min(parent_fitness(:,i)),idealpoint(1,i));
end

figure; hold on;
% Start Iterations
for iter=1:max_iteration
  for i=1:population
    % Reproduction
    % Mutation
    r=round(rand(1,3)*neighbours);
    while r(1)==r(2) || r(2)==r(3) || r(3)==r(1) || min(r)==0
      r=round(rand(1,3)*neighbours);
    end
    u(i,:)=x(Neighbour_indices(i,r(1)),:)+F*(x(Neighbour_indices(i,r(2)),:)-x(Neighbour_indices(i,r(3)),:));

    % Crossover
    for j=1:D
      if rand<=Cr || j==round(rand*D)
        v(i,j)=u(i,j);
      else
        v(i,j)=x(i,j);
      end
    end

    % Repair of the newly generated offsprings
    for j=1:D
      v(i,j)=min(v(i,j),max_val);
      v(i,j)=max(v(i,j),min_val);
    end
    v(i,:)=(v(i,:)>1).*v(i,:);

    %>FITNESS Compute objective function: Intra and Inter class distances
    child_fitness(i,:)=weightedClusterDist(v(i,:),Dist);  

    % Update ideal point
    for j=1:Obj
      idealpoint(1,j)=min(child_fitness(i,j),idealpoint(1,j));
    end

    % Update the neighbours
    for j=1:neighbours
      new_te=max(weight(Neighbour_indices(i,j),:).*abs(child_fitness(i,:)-idealpoint));
      old_te=max(weight(Neighbour_indices(i,j),:).*abs(parent_fitness(Neighbour_indices(i,j),:)-idealpoint));
      if new_te<=old_te
        x(Neighbour_indices(i,j),:)=v(i,:);
        parent_fitness(Neighbour_indices(i,j),:)=child_fitness(i,:);
      end
    end
  end
%   %TEST>
%   scatter(idealpoint(1,1), idealpoint(1,2),'.'); 
%   title(['Generation ' num2str(iter)]); drawnow;
end

% Best Compromise Solution selection
W=findWeights(x,parent_fitness);

% Select and weight the required features
t=repmat([W 1],size(tdata,1),1).*tdata;
t(:,W==0)=[];
l=repmat([W 1],size(ldata,1),1).*ldata;
l(:,W==0)=[];

%% Outputs
%Weighted features:
outArgs.featsWeights = W;

%Predict accuracy:
outArgs.PredictAccInit = PredictCrossValAccuracy(ldata, tdata, kK, udata.mcCost)
outArgs.PredictAcc = PredictCrossValAccuracy(l, t, kK, udata.mcCost)

%KNN accuracy:
outArgs.kNNAccInit = KNNAccuracy(ldata, tdata,kK)
outArgs.kNNAcc = KNNAccuracy(l,t,kK)

end