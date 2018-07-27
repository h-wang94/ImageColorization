% Converts superpixel boundaries into segments
function segImage = get_segments_from_superpixel_boundaries(speed, boundaries)

    L = bwlabel(~boundaries,4);
    labels = zeros(size(boundaries,1),size(boundaries,2),8);
    speeds = zeros(size(boundaries,1),size(boundaries,2),8);
    mean_speed = zeros(size(boundaries,1),size(boundaries,2),8);
    
    label_ind = 1;
    
    for ofs_row = -1:1
        for ofs_col = -1:1
            if (ofs_row == 0 && ofs_col == 0)
                continue;
            end
            
            row_ind_org = 1+ofs_row:ofs_row+size(boundaries,1); 
            col_ind_org = 1+ofs_col:ofs_col+size(boundaries,2);
            row_ind = min(size(boundaries,1), max(1, row_ind_org));
            col_ind = min(size(boundaries,2), max(1, col_ind_org));
            labels(:,:,label_ind) = L(row_ind,col_ind);
            speeds(:,:,label_ind) = speed(row_ind,col_ind);
            labels(row_ind_org < 1, :,label_ind) = 0;
            labels(row_ind_org > size(boundaries,1), :,label_ind) = 0;
            labels(:, col_ind_org < 1,label_ind) = 0;
            labels(:, col_ind_org > size(boundaries,2),label_ind) = 0;
            label_ind = label_ind + 1;
        end
    end
    
    for i = 1:8
        ind = (labels == repmat(labels(:,:,i),[1,1,8]));
        mean_speed(:,:,i) = sum(speeds .* double(ind), 3) ./ ...
                            sum(double(ind), 3);
    end
    
    mean_speed(labels==0) = 0;
    [m,ind_segments] = max(mean_speed,[],3);
    [c,r] = meshgrid(1:size(boundaries,2), 1:size(boundaries,1));
    ind = sub2ind(size(labels),r(:),c(:),ind_segments(:));
    segImage = L;
    relabel = reshape(labels(ind(:)),[size(boundaries,1),size(boundaries,2)]);
    segImage(boundaries) = relabel(boundaries);
