function [S2, C2]  =  Rep_NextState( S1, a1, t1, W1, JJ )
%Finds the next state given the current state and actions
%   S1: Initial state (length 9) at time t
%   a1: GHG reduction rate
%   g1: Geoengineering rate
%   S2(1:9): Next state
%   C2(1:2): Damage and Abatement costs

global gama sai1 sai2 sai3 deltak theta1 theta2
global Sig Eland b11 b12 b21 b22 b23 b32 b33 etha4
global etha1 etha3 deltarf Fex L A deltaT

if JJ == 0
    bb11 = b11;
    bb12 = b12;
    cs = 1;
elseif JJ == 1
    bb11 = b11;
    bb12 = b12;
    cs = 3/4;
elseif JJ == 2
    bb11 = b11;
    bb12 = b12;
    cs = 3/5;    
elseif JJ == 3
    bb11 = b11;
    bb12 = b12;
    cs = 3/6;  
elseif JJ == 4
    bb12 = 0.75 * b12;
    bb11 = 100 - bb12;
    cs = 1;
elseif JJ == 5
    bb12 = 0.50 * b12;
    bb11 = 100 - bb12;
    cs = 1;   
elseif JJ == 6
    bb12 = 0.25 * b12;
    bb11 = 100 - bb12;
    cs = 1;    
end

S2 = zeros(9, 1);
C2 = zeros(2, 1);

%Industrial Emission
IE = Sig(t1 + 1) * (1 - a1) * S1(7);

%Total Emission
S2(8) = IE + Eland(t1 + 1);

%CO2 doubling coefficient
etha2 = deltarf / (deltaT/cs);

%Lower Ocean Concentration
S2(6) = b23 / 100 * S1(5) + b33 / 100 * S1(6);

%Upper Ocean Concentration
S2(5) = bb12 / 100 * S1(4) + b22 / 100 * S1(5)+ b32 / 100 * S1(6);

%Atmospheric Concentration
S2(4) = 10 * S2(8) + bb11 / 100 * S1(4) + b21 / 100 * S1(5);

%Lower Ocean Temperature
S2(3) = S1(3) + etha4 * (S1(2) - S1(3));

%Radiative Forcings
S2(9) = deltarf * ((log(S2(4)) - log(596.4)) / log(2)) + Fex(t1 + 1);

%Atmospheric Temperature
S2(2) = S1(2) + etha1 * (S2(9) - etha2 * S1(2) - etha3 * (S1(2) - S1(3)));

Damage = 1 / (1 + sai1 * S1(2) + sai2 * (W1 * S1(2)) ^ sai3);
Dam = 1 - Damage;
Abate = theta1(t1) * a1 ^ theta2;

C2(1) = Dam;
C2(2) = Abate;

%Economic output net of damage and abatement costs
Q = (1 - (Dam + Abate)) * S1(7);

%Capital at time t1+1
S2(1) = 10 * 0.22 * Q + (1 - deltak) ^ 10 * S1(1);

%Gross output at time t1+1
S2(7) = A(t1 + 1) * (S2(1) ^ gama) * L(t1 + 1)^(1 - gama);

end

