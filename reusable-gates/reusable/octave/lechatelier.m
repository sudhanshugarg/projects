%test function, to solve the following set of reversible chemical reactions.
% X + Y + ZG <->{kc(1)} XYG + Z
% Z + MB <->{kc(2)} ZMB
%all concentrations are in nM.
%DONT SPLIT EQUATIONS BY NEWLINE

function F = lechatelier(r)
c = [1000,1000,1000,0,0,1000,0];
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

%variable is r
rdouble = [0,0];
rsingle = [0];
%fun = @quadraticTest;
fun = @lechatelier;
fr1 = @rxn1;
fr2 = @rxn2;
r = fsolve(fun, rdouble);
r1 = fsolve(fr1, rsingle);
r2 = fsolve(fr2, rsingle);
