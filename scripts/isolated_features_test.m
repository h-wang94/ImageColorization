%Test the colorization from isolated KNN over isolated features

tgt_lab = CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
      neighbor_classes, labels);
figure; imshow(lab2rgb(tgt_lab)); title('Labeled by mode over combined features');

tgt_lab = CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
      neighbor_classes, labels_mm);
figure; imshow(lab2rgb(tgt_lab)); title('Labeled by mode over mode of each feature');

tgt_lab = CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
      neighbor_classes, labels_mt);
figure; imshow(lab2rgb(tgt_lab)); title('Labeled by mode over each feature distribution');


drawnow;

%%

for i = 1:length(target.fvl)
  tgt_lab = CopyClosestSuperpixelFromClassAvgColor(source, target, neighbor_idxs, ...
      neighbor_classes, labels_m(:,i));
  figure; imshow(lab2rgb(tgt_lab)); title(['Labeled by mode of feature ' num2str(i)]);

  drawnow;
end