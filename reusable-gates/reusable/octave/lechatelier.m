%test function, to solve the following set of reversible chemical reactions.
% X + Y + ZG <->{kc(1)} XYG + Z
% Z + MB <->{kc(2)} ZMB
%all concentrations are in nM.
%DONT SPLIT EQUATIONS BY NEWLINE

function F = lechatelier(r,c)
kc = [0.001589364572,31.91271327];

F(1) = kc(1)*r(1)^3 - r(1)^2*(kc(1)*c(1) + kc(1)*c(2) + kc(1)*c(3) - 1) + r(1)*(kc(1)*c(1)*c(2) + kc(1)*c(2)*c(3) + kc(1)*c(3)*c(1) + c(4) + c(5)) -r(1)*r(2) -c(4)*r(2) +c(4)*c(5) - kc(1)*c(1)*c(2)*c(3);


F(2) = kc(2)*r(2)^2 - r(2)*(kc(2)*c(5) + kc(2)*c(6) +1) + kc(2)*c(6)*r(1) - kc(2)*r(1)*r(2) + kc(2)*c(5)*c(6) - c(7);

endfunction

function F2 = quadraticTest(r)
F2(1) = r(1)^2 - r(2)^2 - 10;
F2(2) = r(1) - r(2) - 2;
endfunction

function F3 = rxn1(r)
c = [1000,1000,1000,0,0,1000,0];
kc = [0.001589364572,31.91271327];
F3(1) = kc(1)*r(1)^3 - r(1)^2*(kc(1)*c(1) + kc(1)*c(2) + kc(1)*c(3) - 1) + r(1)*(kc(1)*c(1)*c(2) + kc(1)*c(2)*c(3) + kc(1)*c(3)*c(1) + c(4) + c(5)) + c(4)*c(5) - kc(1)*c(1)*c(2)*c(3);
endfunction

function F4 = rxn2(r,c)
c = [1000,1000,1000,0,1000,1000,0];
kc = [0.001589364572,31.91271327];

F4(1) = kc(2)*r(1)^2 - r(1)*(kc(2)*c(5) + kc(2)*c(6) +1) + kc(2)*c(5)*c(6) - c(7);
endfunction

function [x,y] = twoval()
x = 2;
y = 3;
endfunction



%computes the minimum conc that isn't zero.
function [my1,my2] = getmin(c)
mymin = Inf;
my1 = 0;
my2 = 0;
for i = c
    if(i == 0)
        continue;
    elseif(i < mymin)
        mymin = i;
    endif
endfor
if(mymin == Inf)
    mymin = 0;
endif

my1 = 0.9*mymin;
my2 = 0.9*mymin;

if(c(1) < mymin)
    my1 = 0;
endif
if(c(2) < mymin)
    my2 = 0;
endif

endfunction

function s = singleconc(c,rdouble)
%[rdouble(1),rdouble(2)] = getmin(c);
fun = @(r) lechatelier(r,c);
r = fsolve(fun, rdouble);
disp(rdouble);
disp(r);
endfunction

function s = bestconc()
%variable is r
rdouble = [0,0];
rsingle = [0];
%fun = @quadraticTest;
c = [0,0,0,0,0,0,0];
%only varying the values c1, c2, c3, c6, to figure out the best combination of inputs.

conc = [10,50,100,500,1000,5000];
%conc = [1000,5000];
maxconc = [0,0,0,0,0,0,0];
maxval = -1;
allconc = [];
allvalues = [];

for c(1) = conc
    for c(2) = conc
        for c(3) = conc
            for c(6) = conc
%                disp("The current concentration is:");
%                disp(c);
                [rdouble(1),rdouble(2)] = getmin(c);
%                disp("this is rdouble");
%                disp(rdouble);
                fun = @(r) lechatelier(r,c);
                r = fsolve(fun, rdouble);
%                disp("Retrieved:");
%                disp(r);
                if (c(6) == 0)
                    continue;
                endif
                percent = (100*r(2))/c(6);
                if(percent > maxval)
                    maxval = percent;
                    maxconc = c;
                endif

                if(percent > 98)
                    %store this percentage.
                    allconc = [allconc; c];
                    allvalues = [allvalues; percent];
                endif
            endfor
        endfor
    endfor
endfor

final = [allconc, allvalues];
format short g;
disp(final);
endfunction
