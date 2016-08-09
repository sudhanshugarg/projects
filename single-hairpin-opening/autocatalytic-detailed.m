function sol = lbs_odes()

tfinal = 10000; % Final time for simulation
species = {'sp_1','sp_2','sp_3','sp_4','sp_5','sp_6','sp_7','sp_8','sp_9','sp_10','sp_11','sp_12','sp_13','sp_14','sp_15','sp_16','sp_17','sp_18','sp_19','sp_20','sp_21','sp_22','sp_23','sp_24','sp_25','sp_26','sp_27','sp_28','sp_29','sp_30','sp_31','sp_32','sp_33','sp_34','sp_35','sp_36','sp_37','sp_38','sp_39','sp_40','sp_41','sp_42','sp_43','sp_44','sp_45','sp_46','sp_47'}; % The list of all species
n = length(species);

% Assign initial conditions
x0 = zeros(n,1);
x0(1) = 910;
x0(2) = 1000;
x0(4) = 920;
x0(6) = 930;
x0(12) = 950;
x0(14) = 940;
x0(19) = 960;
x0(21) = 970;
x0(23) = 980;
x0(30) = 990;
x0(37) = 900;

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
legend('A1', 'TB1', 'A1TB1', 'A2', 'A12TB1', 'A3', 'A123TB1', 'A123', 'A12', 'A23', 'A31', 'TA1', 'TA2', 'TA', 'A123TA2', 'A123TB1TA2', 'A23TA2', 'A31TA2', 'B1', 'B1TA1', 'B2', 'B12TA1', 'B3', 'B123TA1', 'B123', 'B12', 'B23', 'B31', 'TB2', 'TB', 'B123TB2', 'B123TA1TB2', 'B23TB2', 'B31TB2', 'TET', 'FQ', 'RC', 'B2FQ', 'B12FQ', 'B23FQ', 'B123FQ', 'B12FQTA1', 'B123FQTA1', 'B123FQTA1TB2', 'B123FQTB2', 'Signal', 'SFQ')

return
endfunction

%%%

function out = compute_odes(x,t,k)

% Assign derivatives

% Define reaction propensities
r_0 = k(1)*x(1)*x(2);
r_1 = k(2)*x(3)*x(4);
r_2 = k(3)*x(5)*x(6);
r_3 = k(4)*x(7);
r_4 = k(5)*x(1)*x(4);
r_5 = k(6)*x(4)*x(6);
r_6 = k(7)*x(6)*x(1);
r_7 = k(1)*x(1)*x(10);
r_8 = k(2)*x(4)*x(11);
r_9 = k(3)*x(6)*x(9);
r_10 = k(8)*x(12)*x(13);
r_11 = k(9)*x(14);
r_12 = k(10)*x(8)*x(14);
r_13 = k(11)*x(15)*x(12);
r_14 = k(10)*x(7)*x(14);
r_15 = k(11)*x(16)*x(12);
r_16 = k(4)*x(16);
r_17 = k(10)*x(10)*x(14);
r_18 = k(11)*x(17)*x(12);
r_19 = k(10)*x(11)*x(14);
r_20 = k(11)*x(18)*x(12);
r_21 = k(1)*x(17)*x(1);
r_22 = k(2)*x(18)*x(4);
r_23 = k(12)*x(19)*x(12);
r_24 = k(13)*x(20)*x(21);
r_25 = k(14)*x(22)*x(23);
r_26 = k(15)*x(24);
r_27 = k(16)*x(19)*x(21);
r_28 = k(17)*x(21)*x(23);
r_29 = k(18)*x(23)*x(19);
r_30 = k(12)*x(19)*x(27);
r_31 = k(13)*x(21)*x(28);
r_32 = k(14)*x(23)*x(26);
r_33 = k(19)*x(2)*x(29);
r_34 = k(20)*x(30);
r_35 = k(21)*x(25)*x(30);
r_36 = k(22)*x(31)*x(2);
r_37 = k(21)*x(24)*x(30);
r_38 = k(22)*x(32)*x(2);
r_39 = k(15)*x(32);
r_40 = k(21)*x(27)*x(30);
r_41 = k(22)*x(33)*x(2);
r_42 = k(21)*x(28)*x(30);
r_43 = k(22)*x(34)*x(2);
r_44 = k(12)*x(33)*x(19);
r_45 = k(13)*x(34)*x(21);
r_46 = k(23)*x(35)*x(36);
r_47 = k(24)*x(37);
r_48 = k(25)*x(21)*x(37);
r_49 = k(12)*x(38)*x(19);
r_50 = k(14)*x(38)*x(23);
r_51 = k(26)*x(26)*x(37);
r_52 = k(14)*x(39)*x(23);
r_53 = k(26)*x(22)*x(37);
r_54 = k(14)*x(42)*x(23);
r_55 = k(26)*x(27)*x(37);
r_56 = k(12)*x(40)*x(19);
r_57 = k(26)*x(25)*x(37);
r_58 = k(26)*x(24)*x(37);
r_59 = k(15)*x(43);
r_60 = k(21)*x(43)*x(30);
r_61 = k(22)*x(44)*x(2);
r_62 = k(21)*x(41)*x(30);
r_63 = k(22)*x(45)*x(2);
r_64 = k(26)*x(31)*x(37);
r_65 = k(26)*x(32)*x(37);
r_66 = k(15)*x(44);
r_67 = k(26)*x(46)*x(37);

% Assign derivatives
dsp_1 = -r_0 - r_4 - r_6 - r_7 - r_21;
dsp_2 = -r_0 + r_3 + r_16 - r_33 + r_34 + r_35 - r_36 + r_37 - r_38 + r_40 - r_41 + r_42 - r_43 + r_60 - r_61 + r_62 - r_63;
dsp_3 = r_0 - r_1;
dsp_4 = -r_1 - r_4 - r_5 - r_8 - r_22;
dsp_5 = r_1 - r_2;
dsp_6 = -r_2 - r_5 - r_6 - r_9;
dsp_7 = r_2 - r_3 - r_14 + r_15;
dsp_8 = r_3 + r_7 + r_8 + r_9 - r_12 + r_13;
dsp_9 = r_4 - r_9;
dsp_10 = r_5 - r_7 - r_17 + r_18;
dsp_11 = r_6 - r_8 - r_19 + r_20;
dsp_12 = -r_10 + r_11 + r_12 - r_13 + r_14 - r_15 + r_17 - r_18 + r_19 - r_20 - r_23 + r_26 + r_39 + r_59 + r_66;
dsp_13 = -r_10 + r_11;
dsp_14 = r_10 - r_11 - r_12 + r_13 - r_14 + r_15 - r_17 + r_18 - r_19 + r_20;
dsp_15 = r_12 - r_13 + r_16 + r_21 + r_22;
dsp_16 = r_14 - r_15 - r_16;
dsp_17 = r_17 - r_18 - r_21;
dsp_18 = r_19 - r_20 - r_22;
dsp_19 = -r_23 - r_27 - r_29 - r_30 - r_44 - r_49 - r_56;
dsp_20 = r_23 - r_24;
dsp_21 = -r_24 - r_27 - r_28 - r_31 - r_45 - r_48;
dsp_22 = r_24 - r_25 - r_53;
dsp_23 = -r_25 - r_28 - r_29 - r_32 - r_50 - r_52 - r_54;
dsp_24 = r_25 - r_26 - r_37 + r_38 - r_58;
dsp_25 = r_26 + r_30 + r_31 + r_32 - r_35 + r_36 - r_57;
dsp_26 = r_27 - r_32 - r_51;
dsp_27 = r_28 - r_30 - r_40 + r_41 - r_55;
dsp_28 = r_29 - r_31 - r_42 + r_43;
dsp_29 = -r_33 + r_34;
dsp_30 = r_33 - r_34 - r_35 + r_36 - r_37 + r_38 - r_40 + r_41 - r_42 + r_43 - r_60 + r_61 - r_62 + r_63;
dsp_31 = r_35 - r_36 + r_39 + r_44 + r_45 - r_64;
dsp_32 = r_37 - r_38 - r_39 - r_65;
dsp_33 = r_40 - r_41 - r_44;
dsp_34 = r_42 - r_43 - r_45;
dsp_35 = -r_46 + r_47 + r_48 + r_51 + r_53 + r_55 + r_57 + r_58 + r_64 + r_65 + r_67;
dsp_36 = -r_46 + r_47;
dsp_37 = r_46 - r_47 - r_48 - r_51 - r_53 - r_55 - r_57 - r_58 - r_64 - r_65 - r_67;
dsp_38 = r_48 - r_49 - r_50;
dsp_39 = r_49 + r_51 - r_52;
dsp_40 = r_50 + r_55 - r_56;
dsp_41 = r_52 + r_56 + r_57 + r_59 - r_62 + r_63;
dsp_42 = r_53 - r_54;
dsp_43 = r_54 + r_58 - r_59 - r_60 + r_61;
dsp_44 = r_60 - r_61 + r_65 - r_66;
dsp_45 = r_62 - r_63 + r_64 + r_66;
dsp_46 = -r_67;
dsp_47 = r_67;


out = [dsp_1; dsp_2; dsp_3; dsp_4; dsp_5; dsp_6; dsp_7; dsp_8; dsp_9; dsp_10; dsp_11; dsp_12; dsp_13; dsp_14; dsp_15; dsp_16; dsp_17; dsp_18; dsp_19; dsp_20; dsp_21; dsp_22; dsp_23; dsp_24; dsp_25; dsp_26; dsp_27; dsp_28; dsp_29; dsp_30; dsp_31; dsp_32; dsp_33; dsp_34; dsp_35; dsp_36; dsp_37; dsp_38; dsp_39; dsp_40; dsp_41; dsp_42; dsp_43; dsp_44; dsp_45; dsp_46; dsp_47];

endfunction
