function s = errortime(gmc, gse, maxFlakeSize, filename, indFile, plotColour)
%source filename;

data = dlmread (filename, " ");
%column: 1=Gmc, 2=Gse, 4=time(Secs), 5=assembly size, 6=mismatched tiles
radius = 0.5;
hold on;

gmcs = [];
times = [];
prevY = -10;
xCenter = 0;
yCenter = 0;

for i=1:rows(data)
%Check the size of the assembly. If atleast 0.05% of the assembly has
%formed, then accept the data point, else not.
%accept=1;
%    accept=0;

if(data(i,5)/maxFlakeSize < 0.25)
    continue;
endif

%x axis is growth time
xCenter = log10(data(i,4));
%y axis is error rate (mismatches per tile)
if(data(i,6) != 0)
    prevY = yCenter;
    yCenter = log10(data(i,6)/data(i,5));
else
    yCenter = prevY;
endif

%if(accept)
gmcs = [gmcs xCenter];
times = [times yCenter];
%endif

endfor

if(mod(indFile, 2) == 1)
    plot(gmcs, times, '--*', 'color', plotColour(indFile, :), 'linewidth', 2, 'markersize', 6);
else 
    plot(gmcs, times, '-+', 'color', plotColour(indFile, :), 'linewidth',2,'markersize',6);
endif

endfunction

function drawMultipleErrors(gmc, gse, maxFlakeSize, filenames)

%filenames = {"otm.dat", "sierpinski-compare.dat"}
%access filenames using filenames{1}, filenames{2} etc.


radius = 0.5;
axis square;
xlim([-6 40]);
ylim([-3.2 0]);
xlabel("Growth Assembly Time (log_{10}(secs))", 'FontSize', 14);
ylabel("Error rate (log_{10}(mismatches/tile))", 'FontSize', 14);
set(gca,'FontName','Helvetica' ,'FontSize',12,'FontWeight','bold','linewidth',2);
grid on;
hold on;

plotColour = rand(20,3);

   for indFile = 1:size(filenames,2)
        errortime(gmc, gse, maxFlakeSize, filenames{indFile}, indFile, plotColour);
   endfor

   legend(filenames);
   legend location northeast;
   hl=findobj(gcf,'type','axes','tag','legend');
   set(hl,'position',[0.6 0.6 0.3 0.3]);
   legend boxoff;
   
   print ('-dsvg', 'sim-error.svg', '-S1200,1200');

hold off
endfunction
