function mxc = MaxCyclicCrossCorr(x,y)

N = size(y,1);
% assert(~(N == 1))

mxc = zeros(1, N);

for i = 1:N
  [b] = xcorr(x,y(i,:));
  mxc(i) = max(b);
end

end

