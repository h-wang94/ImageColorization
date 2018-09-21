function rgb_img = ColorPropagationLevin(scribbled_img, gray_img, scribble_mask)
%Colorization by Optimization.

%set solver=1 to use a multi-grid solver 
%and solver=2 to use an exact matlab "\" solver
solver=2; 

if (size(gray_img,3) > 1)
    gray_img = rgb2gray(gray_img);
end
gray_img(:,:,2) = gray_img;
gray_img(:,:,3) = gray_img(:,:,1);

colorIm = scribble_mask;

sgI=rgb2ntsc(gray_img);
scI=rgb2ntsc(scribbled_img);

ntscIm = zeros(size(scI));
ntscIm(:,:,1)=sgI(:,:,1);
ntscIm(:,:,2)=scI(:,:,2);
ntscIm(:,:,3)=scI(:,:,3);


max_d=floor(log(min(size(ntscIm,1),size(ntscIm,2)))/log(2)-2);
iu=floor(size(ntscIm,1)/(2^(max_d-1)))*(2^(max_d-1));
ju=floor(size(ntscIm,2)/(2^(max_d-1)))*(2^(max_d-1));
id=1; jd=1;
colorIm=colorIm(id:iu,jd:ju,:);
ntscIm=ntscIm(id:iu,jd:ju,:);

if (solver==1)
  rgb_img=getVolColor(colorIm,ntscIm,[],[],[],[],5,1);
  rgb_img=ntsc2rgb(rgb_img);
elseif (solver == 2)
  rgb_img=getColorExact(colorIm,ntscIm);
end

end

