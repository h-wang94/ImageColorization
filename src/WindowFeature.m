function W_fvs = WindowFeature(image, type, window_size)
%Compute window-type features.
%Square window of window_size dimensions centered on the pixel.

if nargin < 3
    window_size = 7;
end

pad_frame = (window_size - 1)/2;
[kR, kC] = size(image);

padded_image = padarray(image, [pad_frame, pad_frame]);

%Computes the feature on windows centered on each pixel of the image.
W_fvs = zeros(window_size*window_size, kR*kC);
for i = (pad_frame+1):kR
    for j = (pad_frame+1):kC
        %index of the actual image (no pads)
        i_in = i - pad_frame;
        j_in = j - pad_frame;
        
        window = padded_image(i-pad_frame:i+pad_frame, ...
                              j-pad_frame:j+pad_frame);
        
        switch type
            case 'dct'
                fv = dct2(window);
            case 'dft'
                fv = abs(fft2(window));
            otherwise
                disp('Type not recognized');
        end
        
        W_fvs(:,(j_in-1)*kR + i_in) = fv(:);
    end
end

end