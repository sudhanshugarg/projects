function s = run_t_re()
source modular-ode.m; 
source t_re.m; 
init = zeros(56,1);

init(2) = 100;
init(8) = 200;
init(15) = 100;
init(18) = 100;
init(24) = 100;
init(31) = 100;
init(34) = 100;
init(39) = 100;
init(42) = 100;
init(49) = 100;
init(52) = 100;

t = 0:1:500;
inBits = 3;
inIndex = [1,2,4,5,7,8,10,11,23,24,26,27];
inMultiplier = [1,2,1];
outBits = 1;
outIndex = [47,48,50,51];
sequence = [0,0,0, 0,0,1, 0,1,1, 0,1,0, 1,1,0, 1,0,0, 1,0,1, 1,1,1, 0,1,1, 0,1,0, 0,0,0];
% sequence = [0,0,0, 0,0,1, 0,1,1, 0,1,0];

%initially, 00
conc = 100;
init(2) = inMultiplier(1)*conc;
init(8) = inMultiplier(2)*conc;
init(24) = inMultiplier(3)*conc;

leg = {
'GX0', 'X0', 'GX0-X0', 'GX1', 'X1', 'GX1-X1', 'GY0', 'Y0', 'GY0-Y0', 'GY1', 'Y1', 'GY1-Y1', 'GI0', 'I0', 'GI0-I0', 'GI1', 'I1', 'GI1-I1', 'X0-Y0-GI0', 'X0-Y1-GI0', 'X1-Y0-GI0', 'X1-Y1-GI1', 'GZ0', 'Z0', 'GZ0-Z0', 'GZ1', 'Z1', 'GZ1-Z1', 'GJ0', 'J0', 'GJ0-J0', 'GJ1', 'J1', 'GJ1-J1', 'Z0-GJ1', 'Z1-GJ0', 'GK0', 'K0', 'GK0-K0', 'GK1', 'K1', 'GK1-K1', 'J0-Y0-GK0', 'J0-Y1-GK0', 'J1-Y0-GK0', 'J1-Y1-GK1', 'GF0', 'F0', 'GF0-F0', 'GF1', 'F1', 'GF1-F1', 'I0-K0-GF0', 'I0-K1-GF1', 'I1-K0-GF1', 'I1-K1-GF1'
};
s=solveode(@compute_t_re,init,t,inIndex,inBits,inMultiplier,outIndex,outBits,sequence,leg);
endfunction
