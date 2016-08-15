function [ value ] = T_fun( act )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global	T t coefm m S lambda deltaT aa0 ag0 ag1
global	TP TrigTemp
global	scenario nuG x Wm

%	Value of the Post-Decision state
[Sa10, cost1]= T_NextState(S(:, t), act(1), act(2), t, deltaT, nuG, Wm, 0);
[Sa11, cost1]= T_NextState(S(:, t), act(1), act(2), t, deltaT, nuG, Wm, TP);

if (x(m, t) == 0)
    if scenario == 1
        pt11 = 0;
    else
        pt11 = max(0, (min(Sa10(2), TrigTemp) - S(2, t))/(TrigTemp - S(2, t)));
    end
else
    pt11 = 1;
end

U1 = (1 - pt11) * T_Utility(Sa10, aa0, ag0, Wm, t + 1) + pt11 * T_Utility(Sa11, aa0, ag1, Wm, t + 1);

if (t < T-1)
    Sa200 = T_NextState(Sa10, aa0, ag0, t + 1, deltaT, nuG, Wm, 0);
    Sa201 = T_NextState(Sa10, aa0, ag0, t + 1, deltaT, nuG, Wm, TP);
    Sa211 = T_NextState(Sa11, aa0, ag1, t + 1, deltaT, nuG, Wm, TP);
    
    if scenario == 1
        pt22 = 0;
    elseif pt11 == 1
        pt22 = 1;
    else        
         pt22 = max(0, (min(Sa200(2), TrigTemp) - Sa10(2))/(TrigTemp - Sa10(2)));
    end
    U2 = (1 - pt11) * ((1 - pt22) * T_Utility(Sa200, aa0, ag0, Wm, t + 2) + pt22 * T_Utility(Sa201, aa0, ag1, Wm, t + 2)) + pt11 * T_Utility(Sa211, aa0, ag1, Wm, t + 2);
else
    U2 = U1;
end

Vb = coefm(m, t) * U2;

%	utility of the current state
W = T_Utility(S(:, t), act(1), act(2), Wm, t);

%   Approximate value of the current state
value = -(W + lambda * Vb);
end