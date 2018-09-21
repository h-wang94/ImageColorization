function f=weightedClusterDist(W,CDist)

% FUNCTION DESCRIPTION %
% This function is to compute the two objective functions
%
% INPUTS %
% 'W': It is the feature selection and weighting vector of size 1 x m, m
% being the total number of features
% 'CDist': It is a 2 x m matrix containg the intra and inter class distance
% in CDist(1,:) and CDist(2,:) respectively.
%
% OUTPUT %
% 'f' is a vector of size 1 x 2, s.t. f(1) and f(2) contain the two
% objective function's values
%
% Written By: Sujoy Paul, Jadavpur University, India
% Email: paul.sujoy.ju@gmail.com

Lambda1 = 500;
Lambda2 = 500;

% intra =  sum(W.*CDist(1,:)) + Lambda1*sum(W>0);
% inter = -sum(W.*CDist(2,:)) + Lambda2*sum(W>0);

intra =  sum(W.*CDist(1,:));
inter = -sum(W.*CDist(2,:));


f=[intra inter];