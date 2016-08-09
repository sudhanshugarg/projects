function [t,y] = processData(fname,col1,col2,left,right,options, yl, yn)
   [t,y] = readlhcfile(fname , col1, col2, left/60, right/60);
   t = t - t(1);

   %Can be optional
   if(options(1))
       lastNvalues = options(1);
   else
       lastNvalues = 10;
   endif

   %2 signifies whether to
   %read the value from yl or not.
   if(options(2))
        ait = min(yl);
    else
        ait = min(y);
   endif

   %3 signifies whether to
   %read the value from yn or not.
   if(options(3))
      aft = mean(yn(end-lastNvalues:end));
   else
      aft = mean(y(end-lastNvalues:end));
   endif
%  disp(aft);

   %4 signifies whether to normalize
   %according to conc.
   if(options(4))
       conc = options(4);
       y = mynormalize(y,ait,aft,conc);
   endif
   
   %5 signifies whether to smooth or not.
   if(options(5))
        y = smooth(y,21);
   endif
   %{
   %}

   %6 signifies whether to normalize
   %according to mintime.
   mintime = Inf;
   if(options(6))
       mintime = options(6);
   endif
   [t,y] = mycutrange(t,y,mintime);

   %7 signifies whether to zero or not.
   if(options(7))
        y = zeros(1,length(y));
   endif

endfunction


function s = solveode2(fn, init, shouldFit)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%Initial Concentrations - usually
    %%%%%% passed in through init.
    %%%%%% init has 6 values.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   init(1) = 0; %Init
   init(2) = 356*2/80; %H9
   init(3) = 0; %IH9
   init(4) = 9*100/80; %RC
   init(5) = 0; %TET
   init(6) = 0; %IH9FQ

   %conc = min([init(1), init(2)]);
   %conc = min([conc, init(4)]);
   conc = 10;

   mintime = Inf;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %DEFAULT OPTIONS
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %taking last values to mean
   options(1) = 10;
   %options 2 and 3: normalize from another curve
   %low and high respectively
   options(2) = 0;
   options(3) = 1;
   %concentration
   options(4) = conc;
   %smooth
   options(5) = 1;
   %mintime
   options(6) = mintime;
   %to zero or not.
   options(7) = 0;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %READ DATA IN
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   options(5) = 0;

   options(3) = 1;
   init(1) = 187*2/80; %L8R0
   [tl8r0h,yl8r0h] = readlhcfile('data/Aug8-RC10nM-H9-10nM-L8R0-Signal.csv' , 1, 2, 1520/60, 1540/60);
   [tl8r0,yl8r0] = processData('data/Aug8-RC10nM-H9-10nM-L8R0-Signal.csv' , 1, 2, 208, 1420, options, 0, yl8r0h);
   %tcalc = tl8r0;
   %ycalc = yl8r0;
   %k(1) = 2e-3

   init(1) = 54*2/80; %L7R1
   [tl7r1h,yl7r1h] = readlhcfile('data/Aug8-RC10nM-H9-10nM-L7R1-Signal.csv' , 1, 2, 6600/60, 7000/60);
   [tl7r1,yl7r1] = processData('data/Aug8-RC10nM-H9-10nM-L7R1-Signal.csv' , 1, 2, 120, 5940, options, 0, yl7r1h);
   %tcalc = tl7r1;
   %ycalc = yl7r1;
   %k(1) = 1.24e-3

   %Data collected Aug 9 2016
   init(1) = 232*2/80; %L6R2
   [tl6r2h,yl6r2h] = readlhcfile('data/Aug9-RC10nM-H9-10nM-L6R2-Signal.csv' , 1, 2, 2840/60, 2910/60);
   [tl6r2,yl6r2] = processData('data/Aug9-RC10nM-H9-10nM-L6R2-Signal.csv' , 1, 2, 190, 2750, options, 0, yl6r2h);
   tcalc = tl6r2;
   ycalc = yl6r2;
   %k(1) = 1.24e-3

   init(1) = 246*2/80; %L5R3
   [tl5r3h,yl5r3h] = readlhcfile('data/Aug9-RC10nM-H9-10nM-L5R3-Signal.csv' , 1, 2, 4000/60, 5200/60);
   [tl5r3,yl5r3] = processData('data/Aug9-RC10nM-H9-10nM-L5R3-Signal.csv' , 1, 2, 190, 3200, options, 0, yl5r3h);
   %tcalc = tl5r3;
   %ycalc = yl5r3;
   %k(1) = 1.24e-3

   init(1) = 202*2/80; %L7R0
   [tl7r0h,yl7r0h] = readlhcfile('data/Aug9-RC10nM-H9-10nM-L7R0-Signal.csv' , 1, 2, 1600/60, 2400/60);
   [tl7r0,yl7r0] = processData('data/Aug9-RC10nM-H9-10nM-L7R0-Signal.csv' , 1, 2, 240, 1480, options, 0, yl7r0h);
   %tcalc = tl7r0;
   %ycalc = yl7r0;
   %k(1) = 1.24e-3

   init(1) = 166*2/80; %L6R1
   [tl6r1h,yl6r1h] = readlhcfile('data/Aug9-RC10nM-H9-10nM-L6R1-Signal.csv' , 1, 2, 1580/60, 1620/60);
   [tl6r1,yl6r1] = processData('data/Aug9-RC10nM-H9-10nM-L6R1-Signal.csv' , 1, 2, 230, 1500, options, 0, yl6r1h);
   %tcalc = tl6r1;
   %ycalc = yl6r1;
   %k(1) = 1.24e-3

   init(1) = 82*2/80; %L5R2
   [tl5r2h,yl5r2h] = readlhcfile('data/Aug9-RC10nM-H9-10nM-L5R2-Signal.csv' , 1, 2, 2910/60, 3000/60);
   [tl5r2,yl5r2] = processData('data/Aug9-RC10nM-H9-10nM-L5R2-Signal.csv' , 1, 2, 260, 2820, options, 0, yl5r2h);
   %tcalc = tl5r2;
   %ycalc = yl5r2;
   %k(1) = 1.24e-3

   options(3) = 0;
   [trc,yrc] = processData('data/Aug8-RC-Test-Vincent-1to1pt5-10nM.csv' , 1, 2, 175, 760, options, 0, 0);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Calculate Fit.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   k = zeros(3,1);
   kmethod = zeros(3,1);
   kinterval = {};
   %kmethod(i) = 1 => using iterative, else default using doubling.
   %currently all unimolecular reaction use iterative methods

   %rate of hairpin + initiator strand
   k(1) = 2e-3;
   kmethod(1) = 0;
   kinterval{1} = [1e-4, 5e-3];

   %backward rate of hairpin:initiator (unimolecular)
   k(2) = 0;
   kmethod(2) = 0;
   kinterval{2} = [1e-2, 1e+1];

   % rate of rc being opened
   k(3) = 1.3e-3; %DONE
   kmethod(3) = 2; % don't calculate
   kinterval{3} = [1e-6, 1e-4];

   %5 is the TET curve, and 1000 is max number of iterations.
   %k is the initial set of values being passed.
   curve = 5;
   %64 chosen since assuming that the upper limit of any rate constant is
   %2^64 * 1e-16 ~ 2k /nMs, which is quite fast.

   if(ne(shouldFit,0))
   [kfit] = fit_rate_constants(tcalc,ycalc,fn,k,kmethod,kinterval,init,curve,64);
   else
       kfit = k;
   endif

   g = @(x,t) fn(x,t,kfit);
   try
   fitcalc = lsode(g, init, tcalc);
   catch
       printf("Error occured is lsode, continuing.\n");
       myflush();
   end_try_catch

   fit = fitcalc(:,curve);
   err = only_leastsquares(fit,ycalc);
   printf("error for this curve is %e\n",err);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %PLOT DATA
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %l = length(init);
   %i = 1;
   plotColor = rand(20,3);

   ops = {
            'tcalc', 'fit', 'FIT',
            'tl8r0', 'yl8r0', 'pL8R0',
            'tl7r1', 'yl7r1', 'pL7R1',
            'tl6r2', 'yl6r2', 'pL6R2',
            'tl5r3', 'yl5r3', 'pL5R3',
            'tl7r0', 'yl7r0', 'pL7R0',
            'tl6r1', 'yl6r1', 'pL6R1',
            'tl5r2', 'yl5r2', 'pL5R2',
            'trc', 'yrc', 'pRC'
        };

        hold on;
        len = size(ops,1);
        for i=1:len
            plot(eval([ops{i,1}]), eval([ops{i,2}]), 'color', plotColor(i,:));
            leg{i} = ops{i,3};
        endfor

        for i=1:9
            xc = 0.9*max(tl8r0);
            yc = (10-i)/10*max(yrc);
            %text(xc,yc,strcat("k(",mat2str(i),")=",mat2str(k(i))));
        endfor
        legend(leg);
        hold off;
        
       i++;
   s = 1;
        ax = gca();
        set(ax, 'fontsize', 15);
        filename = round(time() - 1400000000);
        filename = [num2str(filename) '.jpg']
        print (filename, '-djpg');
endfunction

function s = compare1_2Orders(k1,k2,conc,t)
B = conc*exp(-k1*(t-t(1)));
C = 1./(1/conc + k2*(t-t(1)));
plot(t,B,'r', t, C, 'g');
legend('First Order', 'Second Order');
endfunction

function idx = find_last_less_than(arr, value)
    i=1;
    while((i < length(arr)) && (arr(i) < value))
        i++;
    endwhile
    idx = i;
endfunction

function val = mymin (varargin)
  val = min ([varargin{:}]);
endfunction

function [smallest,avg] = myminmean(y,last)
   smallest = min(y);
   avg = mean(y(end-last:end));
endfunction

function normalized = mynormalize(y,ai,af,conc)
   normalized = (y-ai)*conc/(af-ai);
endfunction

function [tm,val] = mycutrange(t,y,mintime)
   i = find_last_less_than(t, mintime);
   tm = t(1:i);
   val = y(1:i);
endfunction
