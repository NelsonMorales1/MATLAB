clear;clc
close all;
global Fio vo

%Fc2h6=F(1); Fc2h4=F(2); Fh2=F(3); Fch4=F(4)

Fio=[100,0,0,100];          %Initial molar flowrates
vo=10;                      %Initial volumetric flowrate
sspan=[0,700];              %Residence time span

%Run ODE solver
[s,F]=ode45(@C2H6_Hyd,sspan,Fio);

%Plot
figure
plot (s,F(:,1),'+',s,F(:,2),'*',s,F(:,3),'o',s,F(:,4),'d')
legend ('FC_2H_6','FC_2H_4','FH_2','FCH_4')
xlabel ('Residence Time')
ylabel ('Flowrates')
title ('C_2H_4 Formation')

F1 = F(:,1);
F2 = F(:,2);
F3 = F(:,3);
F4 = F(:,4);

lasti = find(s==700);

fprintf('The exit flow of ethane is %.2f\n', F1(lasti));
fprintf('The exit flow of ethylene is %.2f\n',F2(lasti));
fprintf('The exit flow of hydrogen is %.2f\n',F3(lasti));
fprintf('The exit flow of methane is %.2f\n',F4(lasti));

figure
X=(Fio(1)-F(:,1))./(Fio(1));
plot(s,X)
xlabel('Residence Time (s)')
ylabel('C_2H_6 conversion')
title ('X of C_2H_6 in C_2H_4 Formation')

fprintf('The conversion of ethylene at 700s is %.2f\n',X(lasti));

figure
Y=F(:,2)/(Fio(1)-F(:,1));
plot(s,Y)
xlabel('Residence Time (s)')
ylabel('Yield')
title('C_2H_4 Yield in C_2H_4 Formation')

fprintf('The yield of ethylene at 700s is %.2f\n',Y(lasti));

figure
L=F(:,2)./(F(:,1)+F(:,2)+F(:,3));                  %L is selectivity
plot(s,L)
xlabel('Residence Time (s)')
ylabel('Selectivity')
title('C_2H_4 Selectivity in C_2H_4 Formation')

fprintf('The selectivity of ethylene at 700s is %.2f\n',L(lasti));

