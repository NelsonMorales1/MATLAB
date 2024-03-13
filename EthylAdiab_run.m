clear;clc
close all;
global FTo vo To 

%Fc2h6=FT(1); Fc2h4=FT(2); Fh2=FT(3); Fch4=FT(4); T=FT(5) 

To=873.15;                    %Initial temperature in K
FTo=[100,0,0,100,To];         %Initial molar flowrates and temperature
vo=10;                        %Initial volumetric flowrate
sspan=[0,1000];              %Define residence time span
                           
%Run ODE solver
[s,FT]=ode45(@EthylAdiab,sspan,FTo);

%Plot
plot (s,FT(:,1),'+',s,FT(:,2),'*',s,FT(:,3),'o',s,FT(:,4),'d')
legend ('FC_2H_6','FC_2H_4','FH_2','FCH_4')
xlabel ('Residence Time(s)')
ylabel ('Flowrates')
title ('C_2H_4 Formation')

F1 = FT(:,1);
F2 = FT(:,2);
F3 = FT(:,3);
F4 = FT(:,4);

lasti = find(s==1000);

fprintf('The exit flow of ethane is %.2f\n', F1(lasti));
fprintf('The exit flow of ethylene is %.2f\n',F2(lasti));
fprintf('The exit flow of hydrogen is %.2f\n',F3(lasti));
fprintf('The exit flow of methane is %.2f\n',F4(lasti));

figure
X=(FTo(1)-FT(:,1))./(FTo(1));
plot(s,X)
xlabel('Residence Time (s)')
ylabel('C_2H_6 conversion')
title ('X of C_2H_6 in C_2H_4 Formation')
fprintf('The conversion of ethylene at 1,000s is %.2f\n',X(lasti));

figure
Y=FT(:,2)./(200-FT(:,1));
plot(s,Y)
xlabel('Residence Time (s)')
ylabel('Yield')
title('C_2H_4 Yield in C_2H_4 Formation')
fprintf('The yield of ethylene at 1,000s is %.2f\n',Y(lasti));

figure
L=FT(:,2)./(FT(:,2)+FT(:,3)+FT(:,4));                  %L is selectivity
plot(s,L)
xlabel('Residence Time (s)')
ylabel('Selectivity')
title('C_2H_4 Selectivity in C_2H_4 Formation')
fprintf('The selectivity of ethylene at 1,000s is %.2f\n',L(lasti));

figure
plot(s,FT(:,5))
xlabel('Residence Time (s)')
ylabel('Temperature')
title ('Temperature Change')
F5 = FT(:,5);
fprintf('The temperature at 1,000s is %.2f\n',F5(lasti));

