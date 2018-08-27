function [Fn, T] = IsotropicScaling(F)
%Performs isotropic scaling of a set of feature vectors.

ctroid = mean(F');
sig = std(F');

s = sqrt(size(F,1))/sqrt(sum(sig.^2));

%transformation matrix
    %multiply translation with scaling.
T = eye(size(F,1) + 1)*s;
T(end) = 1;
T(1:end-1,end) = -ctroid';

Fn = T*([F; ones(1,size(F,2))]);
Fn = Fn./repmat(Fn(end,:),size(Fn,1),1);

Fn = Fn(1:end-1,:);

end

