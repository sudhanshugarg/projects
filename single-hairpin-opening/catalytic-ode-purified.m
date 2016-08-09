function xdot = arm3(x,t)
%This ode is to solve 2 differential equations x' = t^2 and y' = 2^t;
%In the initial conditions, we have to give the initial values of x and y, i.e.;
%call this function with x = lsode("f", [x0,y0], t);
%where t is the set of times over which the function has to be evalulated.
%k is a bimolecular rate constant here - 5e-5 nM-1 s-1, except k4 which is 1.
    k1 = k2 = k3 = 5e-4;
    k4 = 5e-1;
    %bimolecular leak rates;
    kl1 = kl2 = kl3 = 5e-8;
    xdot = zeros(11,1);
    xdot(1) = -k1*x(1)*x(2) -kl1*x(1)*x(4) -kl3*x(1)*x(6) -k1*x(10)*x(1);       %A1
    xdot(2) = -k1*x(1)*x(2)+k4*x(7);                                            %I
    xdot(3) = k1*x(1)*x(2) - k2*x(3)*x(4);                                      %A1I
    xdot(4) = -k2*x(3)*x(4) -kl1*x(1)*x(4) -kl2*x(4)*x(6) -k2*x(11)*x(4);       %A2
    xdot(5) = k2*x(3)*x(4) - k3*x(5)*x(6);                                      %A12I
    xdot(6) = -k3*x(5)*x(6) -kl2*x(4)*x(6) -kl3*x(1)*x(6) -k3*x(9)*x(6);        %A3
    xdot(7) = k3*x(5)*x(6) - k4*x(7);                                           %A123I
    xdot(8) = k4*x(7) + k3*x(9)*x(6) + k1*x(10)*x(1) + k2*x(11)*x(4);           %A123
    xdot(9) = kl1*x(1)*x(4) -k3*x(9)*x(6);                                      %A12
    xdot(10) = kl2*x(4)*x(6) -k1*x(10)*x(1);                                    %A23
    xdot(11) = kl3*x(1)*x(6) -k2*x(11)*x(4);                                    %A31
endfunction

function s = solveode(fn, init, t)
   %Leak
   [t1,y1] = readlhcfile('PurifiedSysB-Apr18-200nM.csv' , 1, 2, 450/60, 7800/60);
   t1 = t1 - t1(1);
   %Leak + Catalyst.
   [t2,y2] = readlhcfile('PurifiedSysB-Apr18-200nM.csv' , 1, 2, 7900/60, 9900/60);
   t2 = t2 - t2(1);
   t = t2;
   y1 = y1(1:length(y2));
   x = lsode(fn, init, t);
   %Remove initiator - now calculate only leak.
   init(2) = 0;
   x2 = lsode(fn, init, t);


   conc = init(1);
   ai = 6.5;
   af = 42;
   x = ai + (af - ai)*x/conc;
   x2 = ai + (af - ai)*x2/conc;

   l = length(init);
   i = 1;
   plotColor = 'brgkmcypw'; 
%   while( i <= l )
       c = sprintf('%s', plotColor(i));
%       plot (t,x(:,i),c, 'Linewidth', 2);
        h = figure(1);
    plot(
        t,x(:,1),plotColor(9),
        'Linewidth',3,
        t,x(:,2),plotColor(9),
        'Linewidth',3,
        t,x(:,3),plotColor(9),
        'Linewidth',3,
        t,x(:,4),plotColor(9),
        'Linewidth',3,
        t,x(:,5),plotColor(9),
        'Linewidth',3,
        t,x(:,6),plotColor(9),
        'Linewidth',3,
        t,x(:,7),plotColor(9),
        'Linewidth',3,
        t,x(:,8),plotColor(8),
        'Linewidth',1,
        t,x(:,9),plotColor(9),
        'Linewidth',3,
        t,x(:,10),plotColor(9),
        'Linewidth',3,
        t,x(:,11),plotColor(9),
        'Linewidth',3,
        t,x2(:,8),plotColor(3),
        'Linewidth',3,
        t,y1,plotColor(4),
        'Linewidth',2,
        t,y2,plotColor(5),
        'Linewidth',2
        );
        legend('A1','I','A1I','A2','A12I','A3','A123I','A123', 'A12', 'A23', 'A31', 'A123 Leak', 'Data Leak+Catalyst', 'Data Leak-Catalyst');
        
        FN = findall(h,'-property','FontName');
        set(FN,'FontName','/usr/share/fonts/truetype/ttf-lucida/LucidaSansDemiBold.ttf');
        FS = findall(h,'-property','FontSize');
        set(FS,'FontSize',100);

       i++;
%   endwhile
   s = 1;
endfunction

function s = compare1_2Orders(k1,k2,conc,t)
B = conc*exp(-k1*(t-t(1)));
C = 1./(1/conc + k2*(t-t(1)));
plot(t,B,'r', t, C, 'g');
legend('First Order', 'Second Order');
endfunction
