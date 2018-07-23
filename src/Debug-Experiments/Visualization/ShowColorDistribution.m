function ShowColorDistribution(rgb_img, lab_img)
%TODO

c_abs = reshape(lab_img, size(lab_img,1)*size(lab_img,2), 3);
c_rgb = reshape(rgb_img, size(rgb_img,1)*size(rgb_img,2), 3);

figure; hold on;
for i = 1:100:length(c_abs)
  scatter3(c_abs(i,2), c_abs(i,3), c_abs(i,1), '.', 'MarkerEdgeColor', c_rgb(i,:));
end
hold off;
xlabel('a'); ylabel('b'); zlabel('L');
title('Source Lab chrominance distribution (in colors)');

drawnow
end

