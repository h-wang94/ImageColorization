function [strInc, over] = Increment(str, strLimits)
  
stridx = 1;
while (stridx < length(str))
  aux = min(strLimits(stridx), str(stridx) + 1);
  if (aux ~= str(stridx))
    str(stridx) = aux;
    over = false;
    strInc = str;
    return
  else
    str(stridx) = 1;
    stridx = stridx + 1;
  end
end
  
  strInc = str;
  over = true;
end