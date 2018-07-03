%Test the colorization from isolated KNN over isolated features  
    l = 1;
    labels_m = [];
    labels_t = [];
    for i = 1:length(target.fvl)
      [nn_idx, nn_dist] = knnsearch(source.fv_sp(l:l+(target.fvl(i)-1),:)', ...
                                    target.fv_sp(l:l+(target.fvl(i)-1),:)', ...
                                    'K', IP.Kfs);
      classes = source.sp_clusters(nn_idx);

      labels_m_aux = modeTies(classes);
      labels_m = [labels_m labels_m_aux];
      labels_t = [labels_t classes];

      l = l + target.fvl(i);
    end
    %Mode of mode
    labels_mm = modeTies(labels_m);
    labels_mt = modeTies(labels_t);


%% 
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