%Algorithm for calculating the social utility of each state under ceratin
%action (a3)

function [ U ]  =  Rep_Utility( S4, a4, Wm4, t4)
%Social Utility for a given state St
%   state = S4
%   action = a4
%   time = t4
global theta1 theta2 L alpha
global sai1 sai2 sai3

%   Damage cost function
Damage = 1 / (1 + sai1 * S4(2) + sai2 * (Wm4 * S4(2)) ^ sai3);
Dam = 1 - Damage;

%   Abatement cost function
Abate = theta1(t4) * a4 ^ theta2;

%   Net output after damage and abatement
Q = (1 - (Dam + Abate)) * S4(7);

%   Consumption
C = (1 - 0.22) * Q;

%   Consumption per capita
c = C / L(t4) * 1000;

%   Utility per capita
u = 1 + c ^ (1 - alpha) / (1 - alpha);

%   Social utility
U = u * L(t4) * 10;
end

