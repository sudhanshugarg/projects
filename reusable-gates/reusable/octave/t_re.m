function sol = lbs_odes()

tfinal = 500; % Final time for simulation
species = {'sp_1','sp_2','sp_3','sp_4','sp_5','sp_6','sp_7','sp_8','sp_9','sp_10','sp_11','sp_12','sp_13','sp_14','sp_15','sp_16','sp_17','sp_18','sp_19','sp_20','sp_21','sp_22','sp_23','sp_24','sp_25','sp_26','sp_27','sp_28','sp_29','sp_30','sp_31','sp_32','sp_33','sp_34','sp_35','sp_36','sp_37','sp_38','sp_39','sp_40','sp_41','sp_42','sp_43','sp_44','sp_45','sp_46','sp_47','sp_48','sp_49','sp_50','sp_51','sp_52','sp_53','sp_54','sp_55','sp_56'}; % The list of all species
n = length(species);

% Assign initial conditions
x0 = zeros(n,1);
x0(2) = 100;
x0(8) = 200;
x0(15) = 100;
x0(18) = 100;
x0(24) = 100;
x0(31) = 100;
x0(34) = 100;
x0(39) = 100;
x0(42) = 100;
x0(49) = 100;
x0(52) = 100;

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
legend('GX0', 'X0', 'GX0-X0', 'GX1', 'X1', 'GX1-X1', 'GY0', 'Y0', 'GY0-Y0', 'GY1', 'Y1', 'GY1-Y1', 'GI0', 'I0', 'GI0-I0', 'GI1', 'I1', 'GI1-I1', 'X0-Y0-GI0', 'X0-Y1-GI0', 'X1-Y0-GI0', 'X1-Y1-GI1', 'GZ0', 'Z0', 'GZ0-Z0', 'GZ1', 'Z1', 'GZ1-Z1', 'GJ0', 'J0', 'GJ0-J0', 'GJ1', 'J1', 'GJ1-J1', 'Z0-GJ1', 'Z1-GJ0', 'GK0', 'K0', 'GK0-K0', 'GK1', 'K1', 'GK1-K1', 'J0-Y0-GK0', 'J0-Y1-GK0', 'J1-Y0-GK0', 'J1-Y1-GK1', 'GF0', 'F0', 'GF0-F0', 'GF1', 'F1', 'GF1-F1', 'I0-K0-GF0', 'I0-K1-GF1', 'I1-K0-GF1', 'I1-K1-GF1')

return
endfunction

%%%

function dxdt = compute_t_re(x,t,kfit)

% Assign derivatives
dx = zeros(56,1);


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
r_18 = 0.001*x(23)*x(24);
r_19 = 0.001*x(26)*x(27);
r_20 = 0.001*x(29)*x(30);
r_21 = 0.001*x(32)*x(33);
r_22 = 0.001*x(24)*x(34);
r_23 = 0.001*x(23)*x(35);
r_24 = 0.001*x(27)*x(31);
r_25 = 0.001*x(26)*x(36);
r_26 = 0.001*x(37)*x(38);
r_27 = 0.001*x(40)*x(41);
r_28 = 3E-05*x(30)*x(8)*x(39);
r_29 = 0.001*x(29)*x(43);
r_30 = 0.001*x(7)*x(43);
r_31 = 3E-05*x(30)*x(11)*x(39);
r_32 = 0.001*x(29)*x(44);
r_33 = 0.001*x(10)*x(44);
r_34 = 3E-05*x(33)*x(8)*x(39);
r_35 = 0.001*x(32)*x(45);
r_36 = 0.001*x(7)*x(45);
r_37 = 3E-05*x(33)*x(11)*x(42);
r_38 = 0.001*x(32)*x(46);
r_39 = 0.001*x(10)*x(46);
r_40 = 0.001*x(47)*x(48);
r_41 = 0.001*x(50)*x(51);
r_42 = 3E-05*x(14)*x(38)*x(49);
r_43 = 0.001*x(13)*x(53);
r_44 = 0.001*x(37)*x(53);
r_45 = 3E-05*x(14)*x(41)*x(52);
r_46 = 0.001*x(13)*x(54);
r_47 = 0.001*x(40)*x(54);
r_48 = 3E-05*x(17)*x(38)*x(52);
r_49 = 0.001*x(16)*x(55);
r_50 = 0.001*x(37)*x(55);
r_51 = 3E-05*x(17)*x(41)*x(52);
r_52 = 0.001*x(16)*x(56);
r_53 = 0.001*x(40)*x(56);

% Assign derivatives
dsp_1 = -r_0 - r_7 - r_10;
dsp_2 = -r_0 - r_6 + r_8 - r_9 + r_11;
dsp_3 = r_0 + r_7 + r_10;
dsp_4 = -r_1 - r_13 - r_16;
dsp_5 = -r_1 - r_12 + r_14 - r_15 + r_17;
dsp_6 = r_1 + r_13 + r_16;
dsp_7 = -r_2 - r_8 - r_14 - r_30 - r_36;
dsp_8 = -r_2 - r_6 + r_7 - r_12 + r_13 - r_28 + r_29 - r_34 + r_35;
dsp_9 = r_2 + r_8 + r_14 + r_30 + r_36;
dsp_10 = -r_3 - r_11 - r_17 - r_33 - r_39;
dsp_11 = -r_3 - r_9 + r_10 - r_15 + r_16 - r_31 + r_32 - r_37 + r_38;
dsp_12 = r_3 + r_11 + r_17 + r_33 + r_39;
dsp_13 = -r_4 + r_7 + r_8 + r_10 + r_11 + r_13 + r_14 - r_43 - r_46;
dsp_14 = -r_4 + r_6 + r_9 + r_12 - r_42 + r_44 - r_45 + r_47;
dsp_15 = r_4 - r_6 - r_9 - r_12 + r_43 + r_46;
dsp_16 = -r_5 + r_16 + r_17 - r_49 - r_52;
dsp_17 = -r_5 + r_15 - r_48 + r_50 - r_51 + r_53;
dsp_18 = r_5 - r_15 + r_49 + r_52;
dsp_19 = r_6 - r_7 - r_8;
dsp_20 = r_9 - r_10 - r_11;
dsp_21 = r_12 - r_13 - r_14;
dsp_22 = r_15 - r_16 - r_17;
dsp_23 = -r_18 - r_23;
dsp_24 = -r_18 - r_22;
dsp_25 = r_18 + r_23;
dsp_26 = -r_19 - r_25;
dsp_27 = -r_19 - r_24;
dsp_28 = r_19 + r_25;
dsp_29 = -r_20 + r_25 - r_29 - r_32;
dsp_30 = -r_20 + r_24 - r_28 + r_30 - r_31 + r_33;
dsp_31 = r_20 - r_24 + r_29 + r_32;
dsp_32 = -r_21 + r_23 - r_35 - r_38;
dsp_33 = -r_21 + r_22 - r_34 + r_36 - r_37 + r_39;
dsp_34 = r_21 - r_22 + r_35 + r_38;
dsp_35 = r_22 - r_23;
dsp_36 = r_24 - r_25;
dsp_37 = -r_26 + r_29 + r_30 + r_32 + r_33 + r_35 + r_36 - r_44 - r_50;
dsp_38 = -r_26 + r_28 + r_31 + r_34 - r_42 + r_43 - r_48 + r_49;
dsp_39 = r_26 - r_28 - r_31 - r_34 + r_44 + r_50;
dsp_40 = -r_27 + r_38 + r_39 - r_47 - r_53;
dsp_41 = -r_27 + r_37 - r_45 + r_46 - r_51 + r_52;
dsp_42 = r_27 - r_37 + r_47 + r_53;
dsp_43 = r_28 - r_29 - r_30;
dsp_44 = r_31 - r_32 - r_33;
dsp_45 = r_34 - r_35 - r_36;
dsp_46 = r_37 - r_38 - r_39;
dsp_47 = -r_40 + r_43 + r_44;
dsp_48 = -r_40 + r_42;
dsp_49 = r_40 - r_42;
dsp_50 = -r_41 + r_46 + r_47 + r_49 + r_50 + r_52 + r_53;
dsp_51 = -r_41 + r_45 + r_48 + r_51;
dsp_52 = r_41 - r_45 - r_48 - r_51;
dsp_53 = r_42 - r_43 - r_44;
dsp_54 = r_45 - r_46 - r_47;
dsp_55 = r_48 - r_49 - r_50;
dsp_56 = r_51 - r_52 - r_53;


dxdt = [dsp_1; dsp_2; dsp_3; dsp_4; dsp_5; dsp_6; dsp_7; dsp_8; dsp_9; dsp_10; dsp_11; dsp_12; dsp_13; dsp_14; dsp_15; dsp_16; dsp_17; dsp_18; dsp_19; dsp_20; dsp_21; dsp_22; dsp_23; dsp_24; dsp_25; dsp_26; dsp_27; dsp_28; dsp_29; dsp_30; dsp_31; dsp_32; dsp_33; dsp_34; dsp_35; dsp_36; dsp_37; dsp_38; dsp_39; dsp_40; dsp_41; dsp_42; dsp_43; dsp_44; dsp_45; dsp_46; dsp_47; dsp_48; dsp_49; dsp_50; dsp_51; dsp_52; dsp_53; dsp_54; dsp_55; dsp_56];

return
endfunction
