function handle = assemblytime(gmc, gse, maxFlakeSize, filename, indFile, plotColour)
handle = 0;
%source filename;

data = dlmread (filename, " ");
%column: 1=Gmc, 2=Gse, 4=time(Secs), 5=assembly size, 6=mismatched tiles
radius = 0.5;
hold on;

gmcs = [];
times = [];

for i=1:rows(data)

%x axis is Gmc
xCenter = data(i,1);
%y axis is Assembly Time
yCenter = log10(data(i,4));

gmcs = [gmcs xCenter];
times = [times yCenter];

endfor

if(mod(indFile, 2) == 1)
    handle = plot(gmcs, times, '-.', 'color', plotColour(indFile,:), 'linewidth', 2, 'markersize',4);
else 
    handle = plot(gmcs, times, '-x', 'color', plotColour(indFile,:), 'linewidth', 2, 'markersize', 4);
endif

endfunction

function drawMultipleTimes(gmc, gse, maxFlakeSize, filenames)

%filenames = {"otm.dat", "sierpinski-compare.dat"}
%access filenames using filenames{1}, filenames{2} etc.

plotColour = rand(20,3);

radius = 0.5;
axis image;
xlim([0 101]);
ylim([-6 50]);
xlabel("Gmc", 'FontSize', 14);
ylabel("Assembly Growth Time(log_{10} secs)", 'FontSize', 14);
set(gca,'FontName','Helvetica' ,'FontSize',12,'FontWeight','bold','linewidth',2);
grid on;
hold on;

   for indFile = 1:size(filenames,2)
        handle = assemblytime(gmc, gse, maxFlakeSize, filenames{indFile}, indFile, plotColour);
   endfor

   legend(filenames);
   legend location northwest;
   hl=findobj(gcf,'type','axes','tag','legend');
   set(hl,'position',[0.15 0.5 0.3 0.3]);
   legend boxoff;
   
   print ('-dsvg', 'sim-growth.svg', '-S1600,1000');

   %{
   print('sim-growth.eps','-S1600,1000',
   '-deps',
   '-tight');
   %}

hold off
endfunction
