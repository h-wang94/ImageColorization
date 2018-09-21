function y = modeTies(x)
%Mode with ties marked as -1
[y, ~, ties] = mode(x, 2);

for t = 1:length(ties)
  if(length(ties{t}) > 1)
    y(t) = -1;
  end
end

end