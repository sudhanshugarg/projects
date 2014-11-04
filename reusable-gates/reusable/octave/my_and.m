function dxdt = compute_and(x,t)

% Assign derivatives
dx = zeros(22,1);


% Define reaction propensities
r_0 = 0.001*x(1)*x(2);
r_1 = 0.001*x(4)*x(5);
r_2 = 0.001*x(7)*x(8);
r_3 = 0.001*x(10)*x(11);
r_4 = 0.001*x(13)*x(14);
r_5 = 0.001*x(16)*x(17);
r_6 = 3E-05*x(2)*x(8)*x(15);
r_7 = 0.001*x(1)*x(19);
r_8 = 0.001*x(7)*x(19);
r_9 = 3E-05*x(2)*x(11)*x(15);
r_10 = 0.001*x(1)*x(20);
r_11 = 0.001*x(10)*x(20);
r_12 = 3E-05*x(5)*x(8)*x(15);
r_13 = 0.001*x(4)*x(21);
r_14 = 0.001*x(7)*x(21);
r_15 = 3E-05*x(5)*x(11)*x(18);
r_16 = 0.001*x(4)*x(22);
r_17 = 0.001*x(10)*x(22);

% Assign derivatives
dsp_1 = -r_0 - r_7 - r_10;
dsp_2 = -r_0 - r_6 + r_8 - r_9 + r_11;
dsp_3 = r_0 + r_7 + r_10;
dsp_4 = -r_1 - r_13 - r_16;
dsp_5 = -r_1 - r_12 + r_14 - r_15 + r_17;
dsp_6 = r_1 + r_13 + r_16;
dsp_7 = -r_2 - r_8 - r_14;
dsp_8 = -r_2 - r_6 + r_7 - r_12 + r_13;
dsp_9 = r_2 + r_8 + r_14;
dsp_10 = -r_3 - r_11 - r_17;
dsp_11 = -r_3 - r_9 + r_10 - r_15 + r_16;
dsp_12 = r_3 + r_11 + r_17;
dsp_13 = -r_4 + r_7 + r_8 + r_10 + r_11 + r_13 + r_14;
dsp_14 = -r_4 + r_6 + r_9 + r_12;
dsp_15 = r_4 - r_6 - r_9 - r_12;
dsp_16 = -r_5 + r_16 + r_17;
dsp_17 = -r_5 + r_15;
dsp_18 = r_5 - r_15;
dsp_19 = r_6 - r_7 - r_8;
dsp_20 = r_9 - r_10 - r_11;
dsp_21 = r_12 - r_13 - r_14;
dsp_22 = r_15 - r_16 - r_17;


dxdt = [dsp_1; dsp_2; dsp_3; dsp_4; dsp_5; dsp_6; dsp_7; dsp_8; dsp_9; dsp_10; dsp_11; dsp_12; dsp_13; dsp_14; dsp_15; dsp_16; dsp_17; dsp_18; dsp_19; dsp_20; dsp_21; dsp_22];

return
endfunction
