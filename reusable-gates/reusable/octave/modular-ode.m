function s = solveode(fn, init, tcalc, inIndex, inBits, outIndex, outBits, sequence)

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
    kfit = 1;
    maxTime = tcalc(length(tcalc));

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
            ops = {
                'tcalc', 'fitcalc(:,14)', 'F0',
                'tcalc', 'fitcalc(:,17)', 'F1'
            };
            hold on;
            len = size(ops,1);
            for i=1:len
                plot(eval([ops{i,1}])+(j-1)*maxTime, eval([ops{i,2}]), 'color', plotColor(i,:));
                leg{i} = ops{i,3};
            endfor
            legend(leg);
            hold off;
            init = fitcalc(maxTime,:);
        catch
        printf("Error occured in lsode, continuing.\n");
        myflush();
        end_try_catch

    endfor

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %PLOT DATA
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   i = 1;
       i++;
        ax = gca();
        set(ax, 'fontsize', 15);
        filename = round(time() - 1400000000);
        filename = [num2str(filename) '.jpg']
        print (filename, '-djpg');
        s=1;
endfunction
