%% Not a single function file.

1;


function [A,err] = calculate_leastsquares(k, conc, x, y, ai, af)
    A = calculate_absorbance(k, conc, x, y, ai, af);
    diff = A .- y;
    err = sum(diff .^ 2);
endfunction

function [k,err,A] = calculate_rate_constant(conc, x, y, ai, af, iter)
    err1 = 1e21;
    err2 = 1e20;

    %assumption, the least value of k is incr+1.
    incr = 100;
    k = incr+1;
    threshold = 1e-6;
    low = k;
    high = 2*k;
    %keep doubling the value of k until the err is decreasing.
    %once the error starts to increase, we break.
    while(err1 > err2)
        [A,err1] = calculate_leastsquares(k-incr, conc, x, y, ai, af);
        [A,err2] = calculate_leastsquares(k+incr, conc, x, y, ai, af);
        low = k;
        high = 2*k;
        %disp(err1), disp(err2),disp("hullo:"), disp(k);
        
        if(abs(err2-err1) < threshold)
            break;
        elseif(err1 < err2) %increasing error - means min on the left
            low = k/2;
            high = k;
            break;
        elseif (err1 > err2) %decreasing error - means min on the right
            ;
        endif

        k = k*2;
    endwhile

    i = 0;
    %this is the error for low.
    while((i < iter) && (err1 > threshold))
        mid = (low+high)/2;
        [A,err1] = calculate_leastsquares(mid-incr, conc, x, y, ai, af);
        [A,err2] = calculate_leastsquares(mid+incr, conc, x, y, ai, af);
        %err2 should be less than err1. if so, then switch low to mid.
        %else switch low to mid. if err2 and err1 are pretty close, then break.
        if(abs(err2-err1) < threshold)
            disp(err2); disp(err1);
            disp("threshold achieved");
            break;
        elseif(err1 < err2) %increasing error - means min on the left
                high = mid;
        elseif (err1 > err2) %decreasing error - means min on the right
            low = mid;
        endif
        i++;
        %disp("low = "), disp(low), disp(" high = "), disp(high), disp(" err1 = "), disp(err1), disp(" err2 = "), disp(err2);
    endwhile
    k = mid;
    [A,err] = calculate_leastsquares(k, conc, x, y, ai, af);
    printf('%f',k);
    disp(i);
endfunction

function [R] = calculate_reactant_conc(k, conc, x)
%R is the concentration of the reactants at all times in x.
%prodconc is the concentration of the product.
    R = 1./(1/conc - k*(x(1)-x));
%    prodconc = a0 .- yconc;
endfunction

%% This function reads a file filename, chooses the columns xindex
%% for the x-axis, and yindex for the yaxis. It then converts 
%% the x-axis into seconds (from minutes), while it normalizes all
%% the y values to the value of the highest concentration, and returns
%% the same in [x,y].
function [x,y] = readfile(filename, xindex, yindex, tf, conc)
    csvf = csvread(filename);
    x = csvf(2:end,xindex);
    y = csvf(2:end,yindex);
    [v1,i] = max(y);
    x = x(i:end);
    y = y(i:end);
    [v2,i] = min(x);
    x = x(1:i-1);
    y = y(1:i-1);
    x = x-x(1);
    x = x * 60;
    y = y*conc/v1;

    i = 1;
    while((i < length(x)) && (x(i) < tf))
        i++;
    endwhile
    x = x(1:i);
    y = y(1:i);
    disp(i);
    return;
endfunction

function ret = whats(n)
disp(n);
endfunction

function [x,y] = readlhcfile(filename, xindex, yindex, t1, t2)
    csvf = csvread(filename);
    x = csvf(2:end,xindex);
    y = csvf(2:end,yindex);
    %eliminate terminal zeros.
    j = length(x);
    while(x(j) == 0)
        j--;
    endwhile
    x = x(1:j);
    y = y(1:j);

    i=1;
    while((i < length(x)) && (x(i) < t1))
        i++;
    endwhile

    j = i;
    while((j < length(x)) && (x(j) < t2))
        j++;
    endwhile

    if(x(j) > t2)
        j--;
    endif
    
    x = x * 60;
    x = x(i:j);
    y = y(i:j);
    return;

endfunction

function A = calculate_absorbance(k, conc, x, y, ai, af)
    R = calculate_reactant_conc(k, conc, x);
    %P is the product concentration, at all times in x.
    P = conc-R;
    % Now, we must convert concentration, to absorbance values,
    % to plot this correctly.
    % at prod. conc 0, the absorbance was ai, and at prod. concentration
    % conc, the absorbance was af. Hence..
    % conc - 0 / af - ai = y - 0 / a - ai
    A = ai + (af - ai)*P/conc;
    %disp(A);
endfunction

function ret = myplot(k, conc, x, y, ai, af)
%TODO hieu's values are inaccurate since he has calculated
% based on minutes instead of seconds. based on seconds, the
% rate constant changes. Lets figure that out as well.
    A = calculate_absorbance(k, conc, x, y, ai, af);
    plot(x,y,'-', x, A, '+');
endfunction


function least = errplot(k1, k2, inc, conc, x, y, ai, af)
   i = 1;
   while(k1 < k2)
       [A,err] = calculate_leastsquares(k1, conc, x, y, ai, af);
       ey(i) = err;
       ex(i) = k1;
       k1 = k1 + inc;
       i++;
   endwhile
   [v,j] = min(ey);
   least = ex(j);
   plot(ex,ey,'-');

endfunction

function [err] = only_leastsquares(y1, y2)
    diff = y1 - y2;
    err = sum(diff .^ 2);
endfunction

function [kfit] = fit_rate_constants(t, y, fn, kinit, kmethod, kinterval, init, curve, iter)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %BOOK KEEPING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    len = length(kinit);
    kfit = kinit;
    threshold = 1;

    %assumption, the least value of k is 1e-20, for any k.
    %assuming that 1e-20 ~ 0. Can change if required.
    %len = 1;
    numiter = 0;
    maxiter = 2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%Calculating initial error%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    g = @(x,t) fn(x,t,kinit);
    try
    fit = lsode(g, init, t);
    catch
    printf("Error occured in initial set of parameters in lsode\n");
    myflush();
    return;
    end_try_catch

    tet = fit(:,curve);
    init_err = only_leastsquares(tet,y);
    global_err = init_err;
    printf("Initial error based on kinit = %e\n",global_err);
    myflush();


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%Minimizing the global error%%%%%%
    %%%%%global_err initially corresponds
    %%%%%to kinit values. For each iteration
    %%%%% we get a new error. If this error
    %%%%% is lower than global_err, then
    %%%%% we set it as the new global_err
    %%%%% and the rate constants are set
    %%%%% into kfit.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error_decreased = 1;
    while((numiter<maxiter) || (error_decreased > 0))
    numiter++;
    %%flag to see if error was decreased even once in this iteration%%
    error_decreased = 0;

    for i = 1:len
        k = kfit;

        zeroerr = 1e-10;

        if(eq(kmethod(i),0))
            [prev_err, curr_err, low, high] = doublingMethod(t,y,fn,k,init,curve,kinterval{i}(1),i,iter);
        elseif(eq(kmethod(i),1))
            [prev_err, curr_err, low, high] = iterativeMethod(t,y,fn,k,init,curve,i,kinterval{i}(1),kinterval{i}(2));
        else
            printf("Skipping calculation for rate constant k(%d)\n",i);
            myflush();
            continue;
        endif

        %errl, errh, errm (err low, high, mid)
        %assume that the error low is decreasing, mid
        %is unknown, and high is increasing.
        %thus, errl should always be -ve initially, and
        %errh +ve.
        errl = prev_err;
        errm = Inf;
        errh = curr_err;
        printf("i = %d, low = %e, high = %e, errl = %e, errh = %e\n",i,low,high, prev_err, curr_err);
        myflush();

        %
        if(abs(abs(global_err)-abs(prev_err)) < threshold)
            printf("No change observed in error - reaction independent of change in rate constant k(%d)\n",i);
            continue;
        endif

        %this is the binary search error for low.
        for j =1:iter
            printf("binsrch iteration = %d\n",j);
            mid = (low+high)/2;
            printf("i = %d, low = %e, high = %e, mid = %e\n",i,low,high,mid);
            myflush();

            errm = isIncreasingError(t,y,fn,k,init,curve,mid,i);
            printf("errl = %f, errh = %f, errm = %f\n",errl, errh, errm);
            if((abs(abs(errm)-abs(errl)) < threshold) && (abs(abs(errm)-abs(errh)) < threshold))
                printf("number of iterations for rate constant %d=%d\n",i,j);
                myflush();
                break;
            endif
            if((errm > 0) && (errl < 0))
                high = mid;
                errh = errm;
            elseif((errm < 0) && (errh > 0))
                low = mid;
                errl = errm;
            else
                printf("error in calculations for rate constant %d\n",i);
                myflush();
                break;
            endif

            j++;
        endfor

        printf("for rate constant k(%d) = %e, global_err = %e\n",i,kfit(i),global_err);
        printf("for rate constant k(%d) = %e, new_err = %e\n",i,mid,errm);
        %%%%%%%%Check if the new error is lower, then set it.
        if(abs(errm) < abs(global_err))
            global_err = abs(errm);
            kfit(i) = mid;
            printf("global err decreased\n\n");
            error_decreased = 1;
        endif
        for iterator = 1:length(kfit)
            printf("k(%d) = %e\n",iterator, kfit(iterator));
        endfor
        myflush();

    endfor

    endwhile
endfunction

%returns positive value if error is decreasing
%else returns negative value if err is increasing
function [ret] = isIncreasingError(t,y,fn,k,init,curve,rate,i)

    err1 = calculateError(t,y,fn,k,init,curve,0.9*rate,i);
    err2 = calculateError(t,y,fn,k,init,curve,1.1*rate,i);

    %incase the values are too close, assume
    %that the function is still decreasing.
    ret = calculateError(t,y,fn,k,init,curve,rate,i);
    
    mythreshold = 1e2; %chosen arbitrarily
    if((err1 > err2) || abs(err1-err2) < mythreshold)
        ret*=-1;
    endif

    printf("isIncreasingError: rate = %e, err0.9 = %e, err = %e, err1.1 = %e, i=%d\n", rate, err1, ret, err2, i);
    myflush();
    %{
    if(ne(isdoubling,0) && (abs(err1-err2) < mythreshold))
        ret = -1;
    endif
    %}
endfunction

function [ret] = calculateError(t,y,fn,k,init,curve,rate,i)
    k(i) = rate;
    g = @(x,t) fn(x,t,k);
    
    try
    fit = lsode(g, init, t);
    tet = fit(:,curve);
    ret = only_leastsquares(tet,y);
    catch
    ret = Inf;
    printf("calculateError - error occurred lsode\n");
    myflush();
    end_try_catch

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keep doubling the value of k until the err is decreasing.
%once the error starts to increase, we break.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prev_err, curr_err, low, high] = doublingMethod(t,y,fn,k,init,curve,zeroerr,i,iter)
    printf("Using doubling method for rate constant k(%d)\n",i);
        prev_err = curr_err = Inf;
        low = zeroerr;
        high = 2*low;
        for j = 1:iter
            printf("doubling iteration = %d, low = %e\n",j,low);
            %increasing error - means passed min on the left
            % hence break
            prev_err = curr_err;
            curr_err = isIncreasingError(t,y,fn,k,init,curve,low,i);
            if(curr_err > 0)
                high = low;
                low = low/2;
                break;
            endif

            low = 2*low;
            j++;
        endfor
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Replacing doubling with for loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[prev_err, curr_err, low, high] = iterativeMethod(t,y,fn,k,init,curve,i,leftInterval, rightInterval)
printf("Using iterative method for rate constant k(%d)\n",i);
        min_error = Inf;
        prev_err = curr_err = Inf;
        bins = 30;
        diff = (rightInterval - leftInterval)/bins; %chosen arbitrarily
        for j=1:bins
            rate = (j-1)*diff + leftInterval;
            prev_err = curr_err;
            curr_err = calculateError(t,y,fn,k,init,curve,rate,i);
            if(curr_err < min_error)
                k(i) = rate;
                high = rate +diff;
                low = rate -diff;
                if(rate - diff < 0)
                    low = 0;
                endif
                min_error = curr_err;
            endif
            printf("iterative mode: k(%d) = %e, err=%e\n",i,rate,curr_err);
            myflush();
        endfor

        prev_err = isIncreasingError(t,y,fn,k,init,curve,low,i);
        curr_err = isIncreasingError(t,y,fn,k,init,curve,high,i);
endfunction

function [ret] = myflush()
    fflush(stdout);
    diary off;
    diary on;
endfunction
