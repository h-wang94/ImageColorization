function txt = SuperpixelMatchVisualization(~, event_obj)
%TODO
pos = event_obj.Position;
udata = get(0, 'userdata');

fSourceSP = udata.fSourceSP; 
t_sp = udata.targetSuperpixels;
s_sp = udata.sourceSuperpixels;
tgt_size = udata.targetSize;
matches = udata.matchesList;

%Current cursor linear index (target)
lin_idx = sub2ind(tgt_size, pos(2), pos(1));

%Superpixel of current cursor position
idx_t_sp = t_sp(lin_idx);

%Matching superpixel from source
match = matches(idx_t_sp);

%Image show and Wait for click on Source image
figure(fSourceSP); hold on;
source_match = (s_sp == match);
[rs, cs] = find(source_match);
h = scatter(cs, rs, '.r'); hold off;
waitforbuttonpress; 
delete(h); 

%Return text tooltip.
txt = {['R = ' num2str(pos(2)) ', C = ' num2str(pos(1))]};

end

