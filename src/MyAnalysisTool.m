function txt = MyAnalysisTool(~, event_obj)
    pos = event_obj.Position;

    anal_argin = get(0, 'userdata');
    source = anal_argin.source;
    cd_anal = anal_argin.cddt_list;
    
    %Target linear index
    lin_idx = sub2ind(anal_argin.tgt_size, pos(2), pos(1));
    
    %Corresponding candidates for current pixel.
    candidates = cd_anal(:,lin_idx);
    src_size = size(source.luminance);
    [rs, cs] = ind2sub(src_size, candidates);
    
    figure(101);
    imshow(source.luminance); title('Source candidates'); hold on;
    scatter(cs(1), rs(1), '*g');
    scatter(cs(2:end), rs(2:end), '*r');
    hold off;
  
    txt = {['R = ' num2str(pos(2)) ', C = ' num2str(pos(1))]};
end