function s = run_and()
source modular-ode.m; 
source my_and.m; 
init = zeros(22,1);
init(3) = 100;
init(6) = 100;
init(9) = 100;
init(12) = 100;
init(15) = 100;
init(18) = 100;
t = 0:1:500;
inBits = 2;
inIndex = [2,5,1,4,8,11,7,10];
outBits = 1;
outIndex = [14,17];
sequence = [0,0, 0,1, 1,1, 1,0, 0,0, 0,1, 1,0, 1,1];
%initially, 00
init(2) = init(8) = 100;
leg = {'GB0', 'B0', 'GB0-B0', 'GB1', 'B1', 'GB1-B1', 'GC0', 'C0', 'GC0-C0', 'GC1', 'C1', 'GC1-C1', 'GF0', 'F0', 'GF0-F0', 'GF1', 'F1', 'GF1-F1', 'B0-C0-GF0', 'B0-C1-GF0', 'B1-C0-GF0', 'B1-C1-GF1'};
s=solveode(@compute_and,init,t,inIndex,inBits,outBits,outIndex,sequence);
endfunction
