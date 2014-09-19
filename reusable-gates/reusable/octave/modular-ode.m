function s = solveode(fn, init, tcalc, inIndex, inBits, outIndex, outBits, sequence, leg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%index contains the position in init
    %%%%where the input species that have to
    %%%%be added are located.
    %%%%index(1) and (2) are for bit A0 and A1.
    %%%%index(3) and (4) are for bit GA0 and GA1.
    %%%%index(5,6,7,8) B0,B1,GB0,GB1 resp.
    %%%% and so on.

    %%%%%%% inBits signifies number of input bits.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%Initial Concentrations - usually
    %%%%%%passed in through init.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iter = length(sequence)/inBits;

    plotColor = rand(20,3);
    marker = {'+','o','*','.','x','s','^','v'};
    kfit = 1;
    maxTime = tcalc(length(tcalc));
    yshift = 10;

    %sequence refers to order of inputs. B is first input,
    %while C is the second input.
    for j=1:iter
        bit = 1;
        if(j>1)
            for k=inBits-1:-1:0
                if(sequence(inBits*j-k) > sequence(inBits*(j-1)-k))% change to A1
                    init(inIndex(4*(bit-1)+3)) += 100; %GA0
                    init(inIndex(4*(bit-1)+2)) += 100; %A1
                elseif(sequence(inBits*j-k) < sequence(inBits*(j-1)-k))% change to A0
                    init(inIndex(4*(bit-1)+4)) += 100; %GA1
                    init(inIndex(4*(bit-1)+1)) += 100; %A0
                endif
                bit++;
            endfor
        endif

        g = @(x,t) fn(x,t,kfit);
        try
            fitcalc = lsode(g, init, tcalc);
        catch
            printf("Error occured in lsode, continuing.\n");
            myflush();
        end_try_catch

        ops = {};
        for i=1:outBits
            i0 = outIndex(4*(i-1)+1);
            i1 = outIndex(4*(i-1)+2);
            fit0 = fitcalc(:,i0);
            fit1 = fitcalc(:,i1);

            hold on;
            plot(tcalc+(j-1)*maxTime, fit0 + (i-1)*yshift, 'color', plotColor(2*i-1,:), 'marker', marker{2*i-1}, 'linewidth', 1);
            plot(tcalc+(j-1)*maxTime, fit1 + (i-1)*yshift, 'color', plotColor(2*i,:), 'marker', marker{2*i}, 'linewidth', 1);
            lg{2*i-1} = leg{i0};
            lg{2*i} = leg{i1};
            legend(lg, 'location', 'southeast');
            text(maxTime/2 + (j-1)*maxTime, 105, strcat(int2str(sequence(inBits*(j-1)+1: inBits*(j)))));
            hold off;
        endfor

        init = fitcalc(maxTime,:);

    endfor

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %PLOT DATA
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   i = 1;
   i++;
%   text(maxTime/2, 105, strcat(int2str(sequence(1:inBits))));
   xlabel("time(secs)","fontweight","bold");
   ylabel("conc(nM)","fontweight","bold");
   ax = gca();
   set(ax,"fontweight","bold","linewidth",2);
   set(ax, 'fontsize', 15);
   filename = round(time() - 1400000000);
   filename = [num2str(filename) '.jpg']
   print (filename, '-djpg', '-F:Helvetica:12');
   s=1;
endfunction

function [ret] = myflush()
    fflush(stdout);
    diary off;
    diary on;
endfunction

