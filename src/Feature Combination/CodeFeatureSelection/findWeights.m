function W=findWeights(x,fitness)

% FUNCTION DESCRIPTION %
% This function is to compute the best compromise solution
%
% INPUT %
% 'x': It is a matrix of size P x m, where P is the total number of pareto
% optimal solutions and m is the total number of features per data point.
% each row of 'x' contains the trial weight vectors, which comprise the
% population of MOEA/D
% 'fitness': It is a matrix of size P x 2, Column 1 and 2 correspond to 
% objective values of function 1 (intra class distance with penalty) and 
% function 2 (inter class distance with penalty) respectively.
%
% OUTPUTS %
% Missing comment
%
% Written By: Sujoy Paul, Jadavpur University, India
% Email: paul.sujoy.ju@gmail.com

minFitness = min(fitness);
maxFitness = max(fitness);

Membership(:,1) = (maxFitness(1) - fitness(:,1)) ./ (maxFitness(1) - minFitness(1));
Membership(:,2) = (maxFitness(2) - fitness(:,2)) ./ (maxFitness(2) - minFitness(2));

[~,pos] = max(sum(Membership,2));
W = x(pos,:);



