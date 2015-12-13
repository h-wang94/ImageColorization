function fv = compute_fv(mat_struct)
    % currently use sds and luminance as FV;
    fv = [mat_struct.luminance(:), mat_struct.sds(:)];
end