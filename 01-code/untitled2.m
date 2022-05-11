clc
clear all
% %% first demonstrate the imagesc function
% xx = [5 6 7 8];
% yy = [3 4 5];
% C  = [2  2  4  6;
%       8  10 12 14
%       16 18 20 22];
% figure (102)
% imagesc(xx,yy,C)


dx   =  1/100;
dy   =  dx;
xmin = -2;
xmax =  2;
ymin = -2;
ymax =  2;
x    = [xmin:dx:xmax];
y    = [ymin:dy:ymax];
[X,Y]= meshgrid(x,y);
Y    = flipud(Y);
z    = X+i*Y;
W1 = sqrt(1-z.^2);
W2 = -i*sqrt(i*(1-z)).*sqrt(i*(1+z));
W3 = -i*sqrt(-1+z.^2);
W4 = -i*sqrt(z-1).*sqrt(z-1);

mymap = jet;
figure(1);
    clf
    imagesc(x,y,angle(W1));
    axis square
    hold on
    plot([xmin,xmax],[0,0],'--k',[0,0],[ymin,ymax],'--k')
        xticklabels = xmin:1:xmax;
        xticks = linspace(xmin, xmax, numel(xticklabels));
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels(:));
        yticklabels = ymin:1:ymax;
        yticks = linspace(ymin,ymax, numel(yticklabels));
        set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)));
    colormap(mymap);

figure(2);
    clf
    imagesc(x,y,angle(W2));
    axis square
    hold on
    plot([xmin,xmax],[0,0],'--k',[0,0],[ymin,ymax],'--k')
        xticklabels = xmin:1:xmax;
        xticks = linspace(xmin, xmax, numel(xticklabels));
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels(:));
        yticklabels = ymin:1:ymax;
        yticks = linspace(ymin,ymax, numel(yticklabels));
        set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)));
    colormap(mymap);
figure(3);
    clf
    imagesc(x,y,angle(W3));
    axis square
    hold on
    plot([xmin,xmax],[0,0],'--k',[0,0],[ymin,ymax],'--k')
        xticklabels = xmin:1:xmax;
        xticks = linspace(xmin, xmax, numel(xticklabels));
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels(:));
        yticklabels = ymin:1:ymax;
        yticks = linspace(ymin,ymax, numel(yticklabels));
        set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)));
    colormap(mymap);
figure(4);
    clf
    imagesc(x,y,angle(W4));
    axis square
    hold on
    plot([xmin,xmax],[0,0],'--k',[0,0],[ymin,ymax],'--k')
        xticklabels = xmin:1:xmax;
        xticks = linspace(xmin, xmax, numel(xticklabels));
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels(:));
        yticklabels = ymin:1:ymax;
        yticks = linspace(ymin,ymax, numel(yticklabels));
        set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)));
    colormap(mymap);   
    
    
    
    
    