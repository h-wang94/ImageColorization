function txt = SuperpixelClassificationVisualization(~, event_obj)
%TODO
pos = event_obj.Position;
udata = get(0, 'userdata');

%Argin
fSourceSP = udata.fSourceSP; 
t_sp = udata.targetSuperpixels;
s_sp = udata.sourceSuperpixels;
tgt_size = udata.targetSize;
neighbors_list = udata.neighborsList;

%Current cursor linear index (target)
lin_idx = sub2ind(tgt_size, pos(2), pos(1));

%Superpixel of current cursor position
idx_t_sp = t_sp(lin_idx);

%FS neighbors from source
neighbors = neighbors_list(idx_t_sp,:);

figure(fSourceSP); hold on;
hs = [];
for i = 1:length(neighbors)
    source_match = (s_sp == neighbors(i));    
    [rs, cs] = find(source_match);
    h = scatter(cs, rs, '.r'); 
    hs = [hs h];
end
hold off;
waitforbuttonpress;
for i = 1:length(neighbors)
    delete(hs(i));
end

end

