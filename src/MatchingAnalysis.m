function MatchingAnalysis(COL_METHOD, figs, source, target, neighbors_list)
%TODO:
%-argumento neighbors_list esta sendo usado para duas finalidades.

switch COL_METHOD
    case 1
    %Generate candidate source image
    figure(figs.CandidatesImage);
    imshow(source.image); title('Source candidates');

    %Generate cursor input image
    fig = figure(figs.AnalysisInput); 
    imshow(target.rgb); title('Colorized result (indexing)');
    datacursormode on;
    dcm_obj = datacursormode(fig);
    
    %Arguments for Analysis Tool
    AnalysisArguments.sourceSize = size(source.luminance);
    % AnalysisArguments.sourceImage = source.image;
    AnalysisArguments.cddt_list = neighbors_list;
    AnalysisArguments.targetSize = size(target.luminance);
    AnalysisArguments.targetFS = target.fv;
    AnalysisArguments.fCandidatesImage = figs.CandidatesImage;
    AnalysisArguments.fCandidatesFS = figs.LabelsFS;
    set(0,'userdata',AnalysisArguments);

    %Overwrite update function
    set(dcm_obj, 'UpdateFcn', @PixelClassificationVisualization)    
    
    case 2
    figure(figs.SourceSP);
    imshow(imoverlay(source.image, boundarymask(source.sp, 4), 'w')); 
    
    %Generate cursor input image
    fig = figure(figs.AnalysisInput); 
    imshow(imoverlay(target.rgb, boundarymask(target.sp, 4), 'w')); 
    title('Colorized superpixel result (indexing)');
    datacursormode on;
    dcm_obj = datacursormode(fig);
    
    AnalysisArguments.sourceSuperpixels = source.sp;
    AnalysisArguments.targetSuperpixels = target.sp;
    AnalysisArguments.matchesList = neighbors_list;
    AnalysisArguments.fSourceSP = figs.SourceSP;
    AnalysisArguments.targetSize = size(target.luminance);
    set(0,'userdata',AnalysisArguments);
    set(dcm_obj, 'UpdateFcn', @SuperpixelMatchVisualization)

    case 3
    figure(figs.SourceSP);
    imshow(imoverlay(source.image, boundarymask(source.sp, 4), 'w')); 
    
    %Generate cursor input image
    fig = figure(figs.AnalysisInput); 
    imshow(imoverlay(target.rgb, boundarymask(target.sp, 4), 'w')); 
    title('Colorized superpixel result (indexing)');
    datacursormode on;
    dcm_obj = datacursormode(fig);
   
    AnalysisArguments.sourceSuperpixels = source.sp;
    AnalysisArguments.targetSuperpixels = target.sp;
    AnalysisArguments.neighborsList = neighbors_list;
    AnalysisArguments.targetSize = size(target.luminance);
    AnalysisArguments.fSourceSP = figs.SourceSP;
    set(0,'userdata',AnalysisArguments);
    set(dcm_obj, 'UpdateFcn', @SuperpixelClassificationVisualization)
    
    otherwise
    error('Not recognized');
    
end

end

