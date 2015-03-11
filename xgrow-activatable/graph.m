function s = gmcgse(gmc, gse, maxFlakeSize, filename)
%source filename;

data = dlmread (filename, " ");
%column: 1=Gmc, 2=Gse, 4=time(Secs), 5=assembly size, 6=mismatched tiles
radius = 1;
hold on;

for i=1:rows(data)

%y axis is Gmc
yCenter = (2*data(i,1)-1)*radius;
%x axis is Gse
xCenter = (2*data(i,2)-1)*radius;
%disp("gmc: "), disp(data(i,1)), disp(", gse: "), disp(data(i,2));
ratioSize = data(i,5)/maxFlakeSize;
ratioError = data(i,6)/data(i,5);
r = radius;
if (ratioSize < 0.05)
    plot(xCenter, yCenter, '*', "color", 'k');
elseif(ratioSize < 1)
    r = ratioSize*radius;
    drawcircle(xCenter, yCenter, r, ratioError);
else
    drawcircle(xCenter, yCenter, r, ratioError);
endif

endfor

axis square;
xlim([0 (2*gse+1)*radius]);
ylim([0 (2*gmc+1)*radius]);
grid on;

endfunction

function ok = drawcircle(xc, yc, r, ratioError)
%disp("xcenter: "), disp(xc), disp(", ycenter: "), disp(yc);
theta = 0 : 0.25 : 2*pi;
x = r * cos(theta) + xc;
y = r * sin(theta) + yc;

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
