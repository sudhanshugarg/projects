function s = run_t_re()
source modular-ode.m; 
source t_re.m; 
init = zeros(76,1);

init(2) = 100;
init(8) = 100;
init(13) = 100;
init(17) = 100;
init(23) = 100;
init(27) = 100;
init(35) = 100;
init(38) = 100;
init(44) = 100;
init(51) = 100;
init(54) = 100;
init(59) = 100;
init(62) = 100;
init(69) = 100;
init(72) = 100;

t = 0:1:500;
inBits = 3;
inIndex = [2,5,1,4,8,11,7,10,44,47,43,46];
outBits = 1;
outIndex = [68,71,67,70];
%sequence = [0,0,0, 0,1,0, 0,0,1, 1,1,0, 1,0,0, 1,1,1, 1,0,1, 0,1,0, 0,1,1];
sequence = [0,0,0, 0,0,1, 0,1,1, 0,1,0, 1,1,0, 1,0,0, 1,0,1, 1,1,1, 0,1,1, 0,1,0, 0,0,0];

%initially, 00
conc = 100;
init(2) = conc;
init(8) = conc;
init(44) = conc;

leg = {
'GX0', 'X0', 'GX0-X0', 'GX1', 'X1', 'GX1-X1', 'GY0', 'Y0', 'GY0-Y0', 'GY1', 'Y1', 'GY1-Y1', 'GY_1_0-Y_1_0-Y_1__0', 'Y0-GY_1_0', 'Y_1_0', 'Y_1__0', 'GY_2_0-Y_2_0', 'Y_1__0-GY_2_0', 'Y_2_0', 'GY_1_0', 'GY_1_0-Y_1_0', 'GY_2_0', 'GY_1_1-Y_1_1-Y_1__1', 'Y1-GY_1_1', 'Y_1_1', 'Y_1__1', 'GY_2_1-Y_2_1', 'Y_1__1-GY_2_1', 'Y_2_1', 'GY_1_1', 'GY_1_1-Y_1_1', 'GY_2_1', 'GI0', 'I0', 'GI0-I0', 'GI1', 'I1', 'GI1-I1', 'X0-Y_1_0-GI0', 'X0-Y_1_1-GI0', 'X1-Y_1_0-GI0', 'X1-Y_1_1-GI1', 'GZ0', 'Z0', 'GZ0-Z0', 'GZ1', 'Z1', 'GZ1-Z1', 'GJ0', 'J0', 'GJ0-J0', 'GJ1', 'J1', 'GJ1-J1', 'Z0-GJ1', 'Z1-GJ0', 'GK0', 'K0', 'GK0-K0', 'GK1', 'K1', 'GK1-K1', 'J0-Y_2_0-GK0', 'J0-Y_2_1-GK0', 'J1-Y_2_0-GK0', 'J1-Y_2_1-GK1', 'GF0', 'F0', 'GF0-F0', 'GF1', 'F1', 'GF1-F1', 'I0-K0-GF0', 'I0-K1-GF1', 'I1-K0-GF1', 'I1-K1-GF1'
};
s=solveode(@compute_t_re,init,t,inIndex,inBits,outIndex,outBits,sequence,leg);
endfunction
