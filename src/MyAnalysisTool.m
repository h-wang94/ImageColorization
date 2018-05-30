function txt = MyAnalysisTool(~, event_obj)
    pos = event_obj.Position;

    udata = get(0, 'userdata');
%     source_img = udata.sourceImage;
    cd_anal = udata.cddt_list;
    tgt_size = udata.targetSize;
    src_size = udata.sourceSize;
    tgt_fv = udata.targetFS;
    
    fCandidatesImage = udata.fCandidatesImage;
    fCandidatesFS = udata.fCandidatesFS;
    
    %Current cursor linear index (target)
    lin_idx = sub2ind(tgt_size, pos(2), pos(1));
    
    %Corresponding candidates for current pixel.
    candidates = cd_anal(:,lin_idx);
    [rs, cs] = ind2sub(src_size, candidates);
    
    %Figure: Candidates from source
    figure(fCandidatesImage); hold on;
    h1 = scatter(cs(1), rs(1), '*g');
    h2 = scatter(cs(2:end), rs(2:end), '*r');
    hold off;
    
    %Feature: Candidates on feature space
    figure(fCandidatesFS); hold on;
    cursor_fv = tgt_fv(1:2, lin_idx);
    h = scatter(cursor_fv(1), cursor_fv(2), '*k');
    hold off;
    %Wait for click on Labels window
    waitforbuttonpress; 
    delete(h);  delete(h1); delete(h2);

    
    %Return text tooltip.
    txt = {['R = ' num2str(pos(2)) ', C = ' num2str(pos(1))]};
end