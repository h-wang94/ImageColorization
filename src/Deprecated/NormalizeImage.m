function In = NormalizeImage(I)
In = (I - min(min(I)))/(max(max(I)) - min(min(I)));
end

