function s = gmcgse(gmc, gse, maxFlakeSize, filename, indFile)
%source filename;

data = dlmread (filename, " ");
%column: 1=Gmc, 2=Gse, 4=time(Secs), 5=assembly size, 6=mismatched tiles
radius = 0.5;
hold on;

for i=1:rows(data)

%y axis is Gmc
yCenter = (2*data(i,1))*radius;
%x axis is Gse
%xCenter = (2*data(i,2)-1)*radius;
xCenter = (3*(indFile+1))*radius;
%disp("gmc: "), disp(data(i,1)), disp(", gse: "), disp(data(i,2));
ratioSize = data(i,5)/maxFlakeSize;
ratioError = data(i,6)/data(i,5);
r = radius;
if (ratioSize < 0.05)
    plot(xCenter, yCenter, '*', "color", 'k');
elseif(ratioSize < 1)
    r = ratioSize*radius;
    drawcircle(xCenter, yCenter, r, ratioError, indFile);
else
    drawcircle(xCenter, yCenter, r, ratioError, indFile);
endif

endfor


endfunction

function ok = drawcircle(xc, yc, r, ratioError, indFile)
%disp("xcenter: "), disp(xc), disp(", ycenter: "), disp(yc);
if(mod(indFile,2) == 0)
theta = 0 : 0.25 : 2*pi;
x = r * cos(theta) + xc;
y = r * sin(theta) + yc;
else
    r -= 0.1;
    x = [xc-r xc+r xc+r xc-r];
    y = [yc-r yc-r yc+r yc+r];
endif

%color of circle depends on error. black 0,0,0 indicates no error. whitest - 1,1,1 indicates most error.
if(ratioError >= 0.5)
    ratioError = 1;
else
    ratioError = ratioError*2;
endif

color = [ratioError,ratioError,ratioError];
%c = repmat(red,rows(x),1);
fill(x,y,color);
endfunction

function drawMultiple(gmc, gse, maxFlakeSize, filenames)

%filenames = {"otm.dat", "sierpinski-compare.dat"}
%access filenames using filenames{1}, filenames{2} etc.


radius = 0.5;
axis image;
%xlim([0 (2*gse+1)*radius]);
%xlim([(2*gse-20)*radius (2*gse+20)*radius]);
xlim([-1 4*size(filenames,2)+1])
ylim([-2*radius (2*gmc+1)*radius]);
xlabel("Tile Assembly Model");
ylabel("Gmc");
grid on;
hold on;

   for indFile = 1:size(filenames,2)
        gmcgse(gmc, gse, maxFlakeSize, filenames{indFile}, indFile);
   endfor

hold off
endfunction
