function dFTds=EthylAdiab (s,FT)
%(r1) C2H6 -> C2H4 + H2
%(r2) C2H6 + H2 -> 2CH4
%Fc2h6=FT(1); Fc2h4=FT(2); Fh2=FT(3); Fch4=FT(4); T=FT(5)<-Temperature

%Initial Conditions
global FTo vo To 

Fo=FTo(1)+FTo(2)+FTo(3)+FTo(4);        %Total initial total molar flow rate
Cto=Fo/vo;                             %Initial total concentration
Ft=FT(1)+FT(2)+FT(3)+FT(4);            %Total Flowrate

%Concentrations
C(1)=(FT(1)*(Cto/Ft))*(To/FT(5));
C(2)=(FT(2)*(Cto/Ft))*(To/FT(5));
C(3)=(FT(3)*(Cto/Ft))*(To/FT(5));
C(4)=(FT(4)*(Cto/Ft))*(To/FT(5));

%Energy Balance
Tr=298.15;                  %T is reference Tempt in K
DHr1=136360;                %Delta H of rx1 at Tr in J/mol
DHr2=65850;                 %Delta H of rx2 at Tr in J/mol
Cp1=113.70;                 %Cp in J/mol-K of C2H6
Cp2=88.48;                  %Cp of C2H4
Cp3=30.17;                  %Cp of H2
Cp4=66.35;                  %Cp of CH4
Ds1=120.7;                  %Delta S at Tr
Ds2=-12.9;

Cpint1=4.94*(FT(5)-Tr);
Cpint2=11.17*(FT(5)-Tr);
DH1=DHr1+Cpint1;                %Delta H at Trx of reaction 1
DH2=DHr2+Cpint2;                %Delta H at Trx of reaction 2
DS1=Ds1+Cpint1;                 %Delta S at Trx of reaction 1
DS2=Ds2+Cpint2;                 %Delta S at Trx of reaction 2
G1=DH1-FT(5)*DS1;               %Delta G at Trx of reaction 1
G2=DH2-FT(5)*DS2;               %Delta G at Trx of reaction 2

%Rate Constants
k1=46.52E13*exp((-272796.8/(8.314*FT(5))));
kf2=6.6E-20*(FT(5)^2.24)*exp(-3220/FT(5));
kc1=(exp(-G1/(8.314*FT(5))))/(8.314*FT(5));
kc2=(exp((-G2/(8.314*FT(5)))));
k2=kf2/kc2; 

%Rate laws
r1=k1*(C(1)-((C(2)*C(3))/kc1));
r2=0.5*k2*((C(4)^2)-((C(1)*C(3))/kc2));

%Net rates for each species
R(1)=-r1+r2;
R(2)=r1;
R(3)=r1+r2;
R(4)=-2*r2;

SumCpi=FT(1)*Cp1+FT(2)*Cp2+FT(3)*Cp3+FT(4)*Cp4;      %Summatory of FTi*Cpi               
FT(5)=(r1*-DH1+r2*-DH2)/SumCpi;                      %Adiabatic reactor Q=0

dFTds=[R(1);R(2);R(3);R(4);FT(5)].*vo;
end