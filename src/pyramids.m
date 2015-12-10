function [gauss_stack, lap_stack] = pyramids(img, level)
   sigma = 2;
   prev_gauss = img;
   for i = 1:level
       %filter = fspecial('gaussian', 4, sigma);
        filter = fspecial('gaussian', 25, sigma);
        %filter = fspecial('gaussian', 100, sigma);
        %filter = fspecial('gaussian', 144, sigma);
        result = filter2(filter, img);
        if i == 1
            gauss_stack = result;
            %prev_gauss = img;
            lap_stack = prev_gauss - result;
        else
            gauss_stack = cat(ndims(img) + 1, gauss_stack, result);
            lap_stack = cat(ndims(img)+1, lap_stack, prev_gauss - result);
            prev_gauss = result;
        end
        
        sigma = sigma * 2;
        
   end
   % end
   lap_stack(:,:,level) = gauss_stack(:,:,level);
end