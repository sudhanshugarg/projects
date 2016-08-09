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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   init(1) = 0; %A1
   init(2) = 0; %TB1
   init(4) = 0; %A2
   init(6) = 0; %A3
   init(14) = 0; %TA
   %{
   init(19) = 87.46; %B1
   init(21) = 92.95; %B2
   init(23) = 86.38; %B3
   %}
   init(19) = 19.20; %B1
   init(21) = 11.60; %B2
   init(23) = 9.80; %B3
   init(30) = 0; %TB
   init(35) = 15; %RC

   conc = min([init(19), init(21)]);
   conc = min([conc, init(23)]);

   mintime = Inf;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %DEFAULT OPTIONS
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %taking last values to mean
   options(1) = 10;
   %normalize from another curve
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
   %cuvette 3
   options(1) = 100;
   options(3) = 0;
   options(5) = 0;
   %{
   init(12) = 92.15; %TA1
   [t1x3,y1x3] = processData('system-b/PurifiedSysB-Apr18-200nM.csv' , 3, 4, 760, 7830,options,0,0);
   %}
   options(3) = 1;
   init(12) = 0; %TA1
   %max data.
   [t1x2,y1x2] = readlhcfile('M8-mar2016/M8-Set9-10per-20per-Mar2016/Mar9-rc-c-30nm-50ul-temp33pt1-b1-b2-b3m2-0xta1-laterTA1-20nM.csv' , 1, 2, 160000/60, 165000/60);
   [t0x2,y0x2] = processData('M8-mar2016/M8-Set9-10per-20per-Mar2016/Mar9-rc-c-30nm-50ul-temp33pt1-b1-b2-b3m2-0xta1-laterTA1-20nM.csv' , 1, 2, 250, 91800,options,0,y1x2);
   %{
   [t1x2,y1x2] = readlhcfile('system-b/PurifiedSysB-Apr18-200nM.csv' , 1, 2, 7900/60, 9900/60);
   [t0x2,y0x2] = processData('system-b/PurifiedSysB-Apr18-200nM.csv' , 1, 2, 767, 7836,options,0,y1x2);
   init(12) = 39.73; %TA1
   [tpt5x3,ypt5x3] = processData('system-b/PurifiedSysB-Apr18-100nM-pt5x-pt4x.csv' , 3, 4, 702, 7000, options, 0, 0);
   init(12) = 14.82; %TA1
   [tpt4x2,ypt4x2] = processData('system-b/PurifiedSysB-Apr18-100nM-pt5x-pt4x.csv' , 1, 2, 705, 7000, options, 0, 0);
   init(12) = 28.38; %TA1
   [tpt3x3,ypt3x3] = processData('system-b/PurifiedSysB-Apr18-100nM-pt3x-pt2x.csv' , 3, 4, 700, 7000, options, 0, 0);
   init(12) = 16.96; %TA1
   [tpt2x2,ypt2x2] = processData('system-b/PurifiedSysB-Apr18-100nM-pt3x-pt2x.csv' , 1, 2, 700, 7000, options, 0, 0);
   init(12) = 14.02; %TA1
   [tpt1x3,ypt1x3] = processData('system-b/PurifiedSysB-pt1x-ptoh5x.csv' , 3, 4, 928, 20000,options,0,0);
   init(12) = 4.74; %TA1
   [tptoh5x2,yptoh5x2] = processData('system-b/PurifiedSysB-pt1x-ptoh5x.csv' , 1, 2, 925, 80000,options,0,0);
   init(12) = 92.15; %TA1
    [tNA1x2,yNA1x2] = processData('sysb-may1/Check-Sys-B.csv' , 5, 6, 1900, 3300, options, 0, 0);
    [tNA1x2,yNA1x2] = processData('sysb-may1/Check-Sys-B-100nM-new.csv' , 5, 6, 1100, 2700, options, 0, 0);
    [tNA1x2,yNA1x2] = processData('sysb-may1/Check-Sys-B-200nM-new.csv' , 5, 6, 1750, 3100, options, 0, 0);
   init(1) = 95.77; init(4) = 88.02; init(6) = 93.12; init(14) = 100; init(30) = 100; init(2) = 95.47;
   %}
   %{
    [tCC0a1b2,yCC0a1b2] = processData('sysb-may1/Cross-catalytic-100nMTF.csv' , 1, 2, 1100, 8000, options, 0, 0);
    init(2) = 0;
    [tCC0a0b3,yCC0a0b3] = processData('sysb-may1/Cross-catalytic-100nMTF.csv' , 3, 4, 1100, 8000, options, 0, 0);
   init(2) = 6.78;
   [tCC0apt05b2,yCC0apt05b2] = processData('sysb-may1/Cross-catalytic-ptoh5x.csv' , 1, 2, 1400, 3700, options, 0, 0);
   init(2) = 0;
   init(12) = 4.74;
   [tCCpt05a0b3,yCCpt05a0b3] = processData('sysb-may1/Cross-catalytic-ptoh5x.csv' , 3, 4, 1400, 3700, options, 0, 0);
   tcalc = tCCpt05a0b3; ycalc = yCCpt05a0b3;
   %}
   tcalc = t0x2; ycalc = y0x2;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Calculate Fit.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   k = zeros(9,1);
   kmethod = zeros(9,1);
   kinterval = {};
   %kmethod(i) = 1 => using iterative, else default using doubling.
   %currently all unimolecular reaction use iterative methods
   %rate of hairpin + initiator strand
   k(1) = 1.62e-4; %DONE
   kmethod(1) = 1;
   kinterval{1} = [1e-5, 3.4e-4];

   % rate of displacement of initiator (unimolecular)
   k(2) = 5.69e-2; %DONE
   %k(2) = 50; %DONE
   kmethod(2) = 1;
   kinterval{2} = [1e-4, 1.3e-3];

   % rate of toehold exchange - reversibly
   k(3) = 7.6875e-5; %DONE
   kmethod(3) = 2; % don't calculate
   kinterval{3} = [1e-6, 1e-4];

   % rate of B2 displacing RC irreversibly - the toehold mediated strand displacement.
   k(4) = 1.6715e-1;
   kmethod(4) = 1;
   kinterval{4} = [1e-4, 1e-2];

   % rate of two hairpins opening up
   k(5) = 8.66e-8; %DONE, leak rate
   kmethod(5) = 1;
   kinterval{5} = [1e-10, 1e-3];

   % rate of just hairpin B2 opening up RC, via leak
   %k(6) = 5.76e-07; %DONE
   k(6) = 5.05e-14; %DONE
   kmethod(6) = 1;
   kinterval{6} = [1e-14, 1e-10];

   % rate of two single strands hybridizing.
   k(7) = 4e-3; %kta
   kmethod(7) = 2; % don't calculate
   kinterval{7} = [1e-3, 1e+3];

   % rate of duplex dehybridizing
   k(8) = 1.46615e-3; %krta
   kmethod(8) = 2;
   kinterval{8} = [1e-9, 1e-2];

   % rate of RC duplex dehybridizing
   %k(9) = 9.505e-7; %RC opening up spontaneously.
   k(9) = 0;
   kmethod(9) = 2;
   kinterval{9} = [0, 0];

   %33 is the TET curve, and 1000 is max number of iterations.
   %k is the initial set of values being passed.
   curve = 33;
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
   l = length(init);
   i = 1;
   plotColor = rand(20,3);

   ops = {
            'tcalc', 'ycalc', 'AB123 TA TB Expt',
            'tcalc', 'fit', 'fit for curve'
        };

        hold on;
        len = size(ops,1);
        for i=1:len
            plot(eval([ops{i,1}]), eval([ops{i,2}]), 'color', plotColor(i,:));
            leg{i} = ops{i,3};
        endfor

        for i=1:9
            xc = 0.9*max(tcalc);
            yc = (10-i)/10*max(ycalc);
            text(xc,yc,strcat("k(",mat2str(i),")=",mat2str(k(i))));
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
