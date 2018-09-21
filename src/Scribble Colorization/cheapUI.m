function []=cheapUI(im1);
strokeW=3;
gg=rgb2gray(im1);
figure(1);
hold off;
imshow(gg);
imMarked=im1;
imMarked(:,:,1)=gg;
imMarked(:,:,2)=gg;
imMarked(:,:,3)=gg;

[xx,yy]=meshgrid(1:256);xx=xx-128;yy=yy-128;

hh=atan2(yy,xx); hh=(hh+pi)/(2*pi); % hue of chooser
ss=sqrt(xx.^2+yy.^2); % saturation of chooser
ss=ss/128; ss=min(ss,1);
% quanitzed chooser looks better
ss=round(10*ss)/10;
hh=round(20*hh)/20;


colorChooser=zeros([256 256 3]);
colorChooser(:,:,1)=hh;
colorChooser(:,:,2)=ss;
colorChooser(:,:,3)=0.5;
colorRGB=hsv2rgb(colorChooser);

figure(3);
imshow(colorRGB);
figure(2);
clf;
figure(1);

g3=imMarked; % gray level in 3 channel format
disp('(d)rawline (a)pply solver (A)pply exact (c)hoose color (esc) quit');
hold on;
x=[];
y=[];
[sx,sy]=size(gg);
while(1)
    [xx,yy,button]=ginput(1);
    if (button==100); % 'd'
        mask=0;
        for i=2:length(x)
            x0=x(i-1);y0=y(i-1);
            x1=x(i);y1=y(i);
            mask=mask+drawLine([x0 y0],[x1 y1],strokeW,sx,sy);
        end
        mask=(mask>0);
        I=find(mask);
        medianC=zeros(3,1);
        for cc=1:3
            imC=im1(:,:,cc);
            medianC(cc)=median(imC(I));
            imMarked(:,:,cc)=mask*medianC(cc) ...
                +(1-mask).*imMarked(:,:,cc);
        end

        hold off;
        imshow(imMarked);
        x=[];y=[];
        hold on;
    elseif (button==99);% 'c'
        % first figure out luminances
        
        figure(3);
        imshow(colorRGB);
        [cx,cy]=ginput(1);cx=round(cx);cy=round(cy);
        medianC(1)=colorRGB(cy,cx,1);
        medianC(2)=colorRGB(cy,cx,2);
        medianC(3)=colorRGB(cy,cx,3);
         mask=0;
        for i=2:length(x)
            x0=x(i-1);y0=y(i-1);
            x1=x(i);y1=y(i);
            mask=mask+drawLine([x0 y0],[x1 y1],strokeW,sx,sy);
        end
        mask=(mask>0);
        for cc=1:3
            imC=im1(:,:,cc);
            imMarked(:,:,cc)=mask*medianC(cc) ...
                +(1-mask).*imMarked(:,:,cc);
        end

        figure(1);
        hold off;
        imshow(imMarked);
        x=[];y=[];
        hold on;


    elseif (button==97); %'a'
        disp('applying multigrid');
        nI=colorizeFun(g3,imMarked);
        figure(2);
        imshow(nI);
        figure(1);
        x=[];y=[];
     elseif (button==65); %'A'
        disp('running exact solver');
        nI=colorizeFun(g3,imMarked,2);
        figure(2);
        imshow(nI);
        figure(1);
        x=[];y=[];
    elseif (button==27); % esc
        break;
    else

    figure(1);
    plot(xx,yy,'x');
    x=[x;xx];
    y=[y;yy];
end
end


function [mask]=drawLine(x0,x1,strokeW,sx,sy)
[xG,yG]=meshgrid(1:sy,1:sx);
mask=zeros(sx,sy);
d=x1-x0;
step=strokeW/norm(d);
for t=0:step:1
    xn=x0+t*d;
    dImage=(xG-xn(1)).^2+(yG-xn(2)).^2;
    mask=mask+(dImage<strokeW^2);
end
mask=(mask>0);
