function sol = lbs_odes()

tfinal = 2000; % Final time for simulation
species = {'sp_1','sp_2','sp_3','sp_4','sp_5','sp_6','sp_7','sp_8','sp_9','sp_10','sp_11','sp_12','sp_13','sp_14','sp_15','sp_16','sp_17','sp_18','sp_19','sp_20','sp_21','sp_22','sp_23','sp_24','sp_25','sp_26','sp_27','sp_28','sp_29','sp_30','sp_31','sp_32','sp_33','sp_34','sp_35','sp_36','sp_37','sp_38','sp_39','sp_40','sp_41','sp_42','sp_43','sp_44','sp_45'}; % The list of all species
n = length(species);

% Assign initial conditions
x0 = zeros(n,1);
x0(1) = 95.77;
x0(4) = 88.02;
x0(6) = 93.12;
x0(14) = 100;
x0(19) = 87.46;
x0(21) = 92.95;
x0(23) = 86.38;
x0(30) = 100;
x0(35) = 100;

% Solve the ODEs
[t,x] = ode15s(@compute_odes,[0 tfinal],x0);

% Write out a solution structure to be returned by the function
sol.t = t;
for i = 1:n
  sol.(species{i}) = x(:,i);
end

% Produce a plot
figure;
plot(t, x)
legend('A1', 'TB1', 'A1TB1', 'A2', 'A12TB1', 'A3', 'A123TB1', 'A123', 'A12', 'A23', 'A31', 'TA1', 'TA2', 'TA', 'A123TA2', 'A123TB1TA2', 'A23TA2', 'A31TA2', 'B1', 'B1TA1', 'B2', 'B12TA1', 'B3', 'B123TA1', 'B123', 'B12', 'B23', 'B31', 'TB2', 'TB', 'B123TB2', 'B123TA1TB2', 'TET', 'FQ', 'RC', 'B2FQ', 'B12FQ', 'B23FQ', 'B12FQTA1', 'B123FQTA1', 'B123FQ', 'B123FQTA1TB2', 'B123FQTB2', 'Signal', 'SFQ')

return
endfunction

%%%

function dxdt = compute_odes(t,x)

% Assign derivatives
dx = zeros(45,1);


% Define reaction propensities
r_0 = 0.0004*x(1)*x(2);
r_1 = 0.0004*x(3)*x(4);
r_2 = 0.0004*x(5)*x(6);
r_3 = 5*x(7);
r_4 = 5.5E-08*x(1)*x(4);
r_5 = 5.5E-08*x(4)*x(6);
r_6 = 5.5E-08*x(6)*x(1);
r_7 = 0.0004*x(1)*x(10);
r_8 = 0.0004*x(4)*x(11);
r_9 = 0.0004*x(6)*x(9);
r_10 = 0.1*x(12)*x(13);
r_11 = 1E-07*x(14);
r_12 = 0.003*x(8)*x(14);
r_13 = 0.003*x(7)*x(14);
r_14 = 5*x(16);
r_15 = 0.003*x(10)*x(14);
r_16 = 0.003*x(11)*x(14);
r_17 = 0.0004*x(17)*x(1);
r_18 = 0.0004*x(18)*x(4);
r_19 = 0.0004*x(19)*x(12);
r_20 = 0.0004*x(20)*x(21);
r_21 = 0.0004*x(22)*x(23);
r_22 = 5*x(24);
r_23 = 5.5E-08*x(19)*x(21);
r_24 = 5.5E-08*x(21)*x(23);
r_25 = 5.5E-08*x(23)*x(19);
r_26 = 0.0004*x(19)*x(27);
r_27 = 0.0004*x(21)*x(28);
r_28 = 0.0004*x(23)*x(26);
r_29 = 0.1*x(2)*x(29);
r_30 = 1E-07*x(30);
r_31 = 0.003*x(25)*x(30);
r_32 = 0.003*x(24)*x(30);
r_33 = 5*x(32);
r_34 = 0.1*x(33)*x(34);
r_35 = 1E-07*x(35);
r_36 = 5.5E-08*x(21)*x(35);
r_37 = 0.0004*x(36)*x(19);
r_38 = 0.0004*x(36)*x(23);
r_39 = 0.007*x(26)*x(35);
r_40 = 0.007*x(22)*x(35);
r_41 = 0.0004*x(39)*x(23);
r_42 = 0.007*x(27)*x(35);
r_43 = 0.0004*x(37)*x(23);
r_44 = 0.0004*x(38)*x(19);
r_45 = 0.007*x(25)*x(35);
r_46 = 0.007*x(24)*x(35);
r_47 = 5*x(40);
r_48 = 0.003*x(40)*x(30);
r_49 = 0.003*x(41)*x(30);
r_50 = 0.007*x(31)*x(35);
r_51 = 0.007*x(32)*x(35);
r_52 = 5*x(42);
r_53 = 0.007*x(44)*x(35);

% Assign derivatives
dsp_1 = -r_0 - r_4 - r_6 - r_7 - r_17;
dsp_2 = -r_0 + r_3 + r_14 - r_29 + r_30 + r_31 + r_32 + r_48 + r_49;
dsp_3 = r_0 - r_1;
dsp_4 = -r_1 - r_4 - r_5 - r_8 - r_18;
dsp_5 = r_1 - r_2;
dsp_6 = -r_2 - r_5 - r_6 - r_9;
dsp_7 = r_2 - r_3 - r_13;
dsp_8 = r_3 + r_7 + r_8 + r_9 - r_12;
dsp_9 = r_4 - r_9;
dsp_10 = r_5 - r_7 - r_15;
dsp_11 = r_6 - r_8 - r_16;
dsp_12 = -r_10 + r_11 + r_12 + r_13 + r_15 + r_16 - r_19 + r_22 + r_33 + r_47 + r_52;
dsp_13 = -r_10 + r_11;
dsp_14 = r_10 - r_11 - r_12 - r_13 - r_15 - r_16;
dsp_15 = r_12 + r_14 + r_17 + r_18;
dsp_16 = r_13 - r_14;
dsp_17 = r_15 - r_17;
dsp_18 = r_16 - r_18;
dsp_19 = -r_19 - r_23 - r_25 - r_26 - r_37 - r_44;
dsp_20 = r_19 - r_20;
dsp_21 = -r_20 - r_23 - r_24 - r_27 - r_36;
dsp_22 = r_20 - r_21 - r_40;
dsp_23 = -r_21 - r_24 - r_25 - r_28 - r_38 - r_41 - r_43;
dsp_24 = r_21 - r_22 - r_32 - r_46;
dsp_25 = r_22 + r_26 + r_27 + r_28 - r_31 - r_45;
dsp_26 = r_23 - r_28 - r_39;
dsp_27 = r_24 - r_26 - r_42;
dsp_28 = r_25 - r_27;
dsp_29 = -r_29 + r_30;
dsp_30 = r_29 - r_30 - r_31 - r_32 - r_48 - r_49;
dsp_31 = r_31 + r_33 - r_50;
dsp_32 = r_32 - r_33 - r_51;
dsp_33 = -r_34 + r_35 + r_36 + r_39 + r_40 + r_42 + r_45 + r_46 + r_50 + r_51 + r_53;
dsp_34 = -r_34 + r_35;
dsp_35 = r_34 - r_35 - r_36 - r_39 - r_40 - r_42 - r_45 - r_46 - r_50 - r_51 - r_53;
dsp_36 = r_36 - r_37 - r_38;
dsp_37 = r_37 + r_39 - r_43;
dsp_38 = r_38 + r_42 - r_44;
dsp_39 = r_40 - r_41;
dsp_40 = r_41 + r_46 - r_47 - r_48;
dsp_41 = r_43 + r_44 + r_45 + r_47 - r_49;
dsp_42 = r_48 + r_51 - r_52;
dsp_43 = r_49 + r_50 + r_52;
dsp_44 = -r_53;
dsp_45 = r_53;


dxdt = [dsp_1; dsp_2; dsp_3; dsp_4; dsp_5; dsp_6; dsp_7; dsp_8; dsp_9; dsp_10; dsp_11; dsp_12; dsp_13; dsp_14; dsp_15; dsp_16; dsp_17; dsp_18; dsp_19; dsp_20; dsp_21; dsp_22; dsp_23; dsp_24; dsp_25; dsp_26; dsp_27; dsp_28; dsp_29; dsp_30; dsp_31; dsp_32; dsp_33; dsp_34; dsp_35; dsp_36; dsp_37; dsp_38; dsp_39; dsp_40; dsp_41; dsp_42; dsp_43; dsp_44; dsp_45];

return
endfunction
