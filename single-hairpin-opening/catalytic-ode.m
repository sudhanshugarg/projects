function xdot = arm3(x,t)
%This ode is to solve 2 differential equations x' = t^2 and y' = 2^t;
%In the initial conditions, we have to give the initial values of x and y, i.e.;
%call this function with x = lsode("f", [x0,y0], t);
%where t is the set of times over which the function has to be evalulated.
%k is a bimolecular rate constant here - 5e-5 nM-1 s-1, except k4 which is 1.
    k1 = k2 = k3 = 4e-4;
    k4 = 5e-1;
    %bimolecular leak rates;
    kl1 = kl2 = kl3 = 5.5e-8;
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

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %READ DATA IN
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

   conc = 10;
   %cuvette 3
   [tl8r0,yl8r0] = readlhcfile('data/Aug8-RC10nM-H9-10nM-L8R0-Signal.csv' , 1, 2, 220/60, 1420/60);
   tl8r0 = tl8r0 - tl8r0(1);

   [tl8r0h,yl8r0h] = readlhcfile('data/Aug8-RC10nM-H9-10nM-L8R0-Signal.csv' , 1, 2, 1520/60, 1540/60);
   tl8r0h = tl8r0h - tl8r0h(1);

   [tl7r1,yl7r1] = readlhcfile('data/Aug8-RC10nM-H9-10nM-L7R1-Signal.csv' , 1, 2, 150/60, 5940/60);
   tl7r1 = tl7r1 - tl7r1(1);

   [tl7r1h,yl7r1h] = readlhcfile('data/Aug8-RC10nM-H9-10nM-L7R1-Signal.csv' , 1, 2, 6600/60, 7000/60);
   tl7r1h = tl7r1h - tl7r1h(1);

   [trc,yrc] = readlhcfile('data/Aug8-RC-Test-Vincent-1to1pt5-10nM.csv' , 1, 2, 175/60, 760/60);
   trc = trc - trc(1);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %Absorbance VALUES
   %pt0x2 means 0x initiator, cuvette 2.
   %IMPORTANT!!!!!!!!!! These values must be calculated 
   %before you shorten the data set.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

   ail8r0 = min(yl8r0);
   afl8r0 = mean(yl8r0h(1:end));
   ail7r1 = min(yl7r1);
   afl7r1 = mean(yl7r1h(1:end));
   airc = min(yrc);
   afrc = mean(yrc(end-100:end));

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Modify length to shortest data set.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   t = trc;
   yl8r0 = yl8r0(1:length(t));
   tl8r0 = tl8r0(1:length(t));
   yl7r1 = yl7r1(1:length(t));
   tl7r1 = tl7r1(1:length(t));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Calculate Fit.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %{
   fit1x = lsode(fn, init, t);
   %Vary Concentration of Initiator
   conc = min([init(1), init(4), init(6)]);
   
   % 40% initiator.
   init(2) = 39.73;
   fitpt5x = lsode(fn, init, t);
   % 15% initiator.
   init(2) = 14.82;
   fitpt4x = lsode(fn, init, t);
   % 28% initiator.
   init(2) = 28.38;
   fitpt3x = lsode(fn, init, t);
   % 17% initiator.
   init(2) = 16.96;
   fitpt2x = lsode(fn, init, t);
   % 14% initiator.
   init(2) = 14.02;
   fitpt1x = lsode(fn, init, t);
   % 5% initiator.
   init(2) = 4.74;
   fitptoh5x = lsode(fn, init, t);
   % 0% initiator.
   init(2) = 0;
   fit0x = lsode(fn, init, t);
   %}

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %NORMALIZE 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   yl8r0 = (yl8r0-ail8r0)*conc/(afl8r0-ail8r0);
   yl7r1 = (yl7r1-ail7r1)*conc/(afl7r1-ail7r1);
   yrc = (yrc-airc)*conc/(afrc-airc);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %SMOOTHEN DATA - For now.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   yl8r0 = supsmu(tl8r0, yl8r0);
   yl7r1 = supsmu(tl7r1, yl7r1);
   yrc = supsmu(trc, yrc);


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   %PLOT DATA
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   l = length(init);
   i = 1;
   plotColor = 'prgmcybk'; 
   c = sprintf('%s', plotColor(i));
   h = figure(1);

    plot(
    %{
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
%        t,fit0x(:,8),plotColor(1),
%        'Linewidth',3,
        t,y0x2,plotColor(2),
        'Linewidth',2,
%        t,y1x2,plotColor(2),
%        'Linewidth',2,
%        t,fit1x(:,8),plotColor(1),
%        'Linewidth',2,
        t,y1x3,plotColor(3),
        'Linewidth',2,
        t,ypt5x3,plotColor(4),
        'Linewidth',2,
%        t,fitpt5x(:,8),plotColor(1),
%        'Linewidth',2,
        t,ypt4x2,plotColor(5),
        'Linewidth',2,
%        t,fitpt4x(:,8),plotColor(1),
%        'Linewidth',2,
        t,ypt3x3,plotColor(6),
        'Linewidth',2,
%        t,fitpt3x(:,8),plotColor(3),
%        'Linewidth',2,
        t,ypt2x2,plotColor(7),
        'Linewidth',2,
%        t,fitpt2x(:,8),plotColor(3),
%        'Linewidth',2,
        t,ypt1x3,plotColor(2),
        'Linewidth',2,
%        t,fitpt1x(:,8),plotColor(3),
%        'Linewidth',2,
        t,yptoh5x2,plotColor(3),
        'Linewidth',2,
%        t,fitptoh5x(:,8),plotColor(3),
%        'Linewidth',2,
        t,yt1x1,plotColor(4),
        'Linewidth',2,
        t,yt1x2,plotColor(5),
        'Linewidth',2,
        t,yt1x3,plotColor(6),
        'Linewidth',2,
        t,ya1x1,plotColor(7),
        'Linewidth',2,
    %}
        tl8r0,yl8r0,plotColor(1),
        'Linewidth',2,
        tl7r1,yl7r1,plotColor(2),
        'Linewidth',2,
        trc,yrc,plotColor(3),
        'Linewidth',2,
        'linestyle',':');
        legend(
        'L8R0',
        'L7R1',
        'RC');
        
        FN = findall(h,'-property','FontName');
        set(FN,'FontName','/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf');
        FS = findall(h,'-property','FontSize');
        %set(FS,'FontSize',14);

       i++;
   s = 1;
endfunction

function s = compare1_2Orders(k1,k2,conc,t)
B = conc*exp(-k1*(t-t(1)));
C = 1./(1/conc + k2*(t-t(1)));
plot(t,B,'r', t, C, 'g');
legend('First Order', 'Second Order');
endfunction
