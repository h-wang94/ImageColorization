function DCT_fvs = DCTWindowFeature(image, window_size)
%Compute DCT features for each pixel. For each pixel, returns the DCT
%coefficients of a square window of window_size dimensions centered on the
%pixel.
%TODO: Prototipar para usar mesma estrutura com outras transformadas.

if nargin < 2
    window_size = 7;
end

pad_frame = (window_size - 1)/2;
[kR, kC] = size(image);

padded_image = padarray(image, [pad_frame, pad_frame]);

%Computes the feature on windows centered on each pixel of the image.
DCT_fvs = zeros(window_size*window_size, kR*kC);
for i = (pad_frame+1):kR
    for j = (pad_frame+1):kC
        %index of the actual image (no pads)
        i_in = i - pad_frame;
        j_in = j - pad_frame;
        
        window = padded_image(i-pad_frame:i+pad_frame, ...
                              j-pad_frame:j+pad_frame);

        fv = dct2(window);
        DCT_fvs(:,(j_in-1)*kR + i_in) = fv(:);
    end
end

end