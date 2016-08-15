%Algorithm for calculating the social utility of each state under ceratin
%action (a3)

function [ U ]  =  T_Utility( S4, a4, g4, Wm4, t4)
%Social Utility for a given state St
%   state = S4
%   action = a4, g4
%   time = t4
global theta1 theta2 theta3 L alpha Gcoeff deltarf Fex etha1
global thetaGE sai2temp sai2ocean sai2atmos nuG Geffective TP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Radiative Forcings change due to GE
RF0 = S4(9);
S4(9) = (deltarf * ((log(S4(4)) - log(596.4)) / log(2)) + Fex(t4)) * (1 - Geffective * g4);

%Atmospheric Temperature change due to GE
S4(2) = S4(2) + etha1 * (S4(9) - RF0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if TP == 7 && S4(10) == 1
    dm = 0.01;
else
    dm = 0;
end

%   Damage cost function
Damage = (1 - dm)/((1 + sai2temp * (Wm4 * S4(2)) ^ 2 + sai2ocean * (S4(5) - 1094) ^ 2 + sai2atmos * (S4(4) - 596.4) ^ 2) * (1 + nuG * g4 ^ 2));
Dam = 1 - Damage;

%   Abatement cost function
Abate = theta1(t4) * a4 ^ theta2;

%   Geoengineering cost function
%Geoengineer = Gcoeff * theta1(t4) * g4 ^ theta3;
Geoengineer = Gcoeff * thetaGE * g4 ^ theta3;

%   Net output after damage and abatement
Q = (1 - (Dam + Abate + Geoengineer)) * S4(7);

%   Consumption
C = (1 - 0.22) * Q;

%   Consumption per capita
c = C / L(t4) * 1000;

%   Utility per capita
u = 1 + c ^ (1 - alpha) / (1 - alpha);

%   Social utility
U = u * L(t4) * 10;
end

