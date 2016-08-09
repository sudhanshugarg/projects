function sol = lbs_odes()

tfinal = 10000; % Final time for simulation
species = {'sp_1','sp_2','sp_3','sp_4','sp_5','sp_6'}; % The list of all species
n = length(species);

% Assign initial conditions
x0 = zeros(n,1);
x0(1) = 10;
x0(2) = 10;
x0(4) = 10;

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
  legend('Init', 'H9', 'IH9', 'RC', 'TET', 'IH9FQ')

  return
  endfunction

  %%%

  function dxdt = compute_odes(x,t,k)

  % Assign derivatives
  dx = zeros(6,1);


  % Define reaction propensities
  r_0 = k(1)*x(1)*x(2);
  r_1 = k(2)*x(3);
  r_2 = k(3)*x(3)*x(4);

  % Assign derivatives
  dsp_1 = -r_0 + r_1;
  dsp_2 = -r_0 + r_1;
  dsp_3 = r_0 - r_1 - r_2;
  dsp_4 = -r_2;
  dsp_5 = r_2;
  dsp_6 = r_2;


  dxdt = [dsp_1; dsp_2; dsp_3; dsp_4; dsp_5; dsp_6];

  return
  endfunction
