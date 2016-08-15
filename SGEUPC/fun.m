function [ value ] = fun( act )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global	T t coefm m S lambda aa0 ag0 ag1
global	TP TrigTemp
global	scenario x Wm GE GEStrategy GETrigTemp

%	Value of the Post-Decision state
[Sa10, cost1]= NextState(S(:, t), act(1), act(2), t, Wm, 0);
[Sa11, cost1]= NextState(S(:, t), act(1), act(2), t, Wm, TP);

if (x(m, t) == 0)
    if scenario ~= 1
        pt11 = GE * GEStrategy * (S(2, t) > GETrigTemp);
    else
        pt11 = max(0, (min(Sa10(2), TrigTemp) - S(2, t))/(TrigTemp - S(2, t)));
    end
else
    pt11 = 1;
end

U1 = (1 - pt11) * Utility(Sa10, aa0, ag0, Wm, t + 1) + pt11 * Utility(Sa11, aa0, ag1, Wm, t + 1);

if (t < T-1)
    Sa200 = NextState(Sa10, aa0, ag0, t + 1, Wm, 0);
    Sa201 = NextState(Sa10, aa0, ag0, t + 1, Wm, TP);
    Sa211 = NextState(Sa11, aa0, ag1, t + 1, Wm, TP);
    
    if scenario ~= 1
        pt22 = GE * GEStrategy * (Sa10(2) > GETrigTemp);
    elseif pt11 == 1
        pt22 = 1;
    else        
         pt22 = max(0, (min(Sa200(2), TrigTemp) - Sa10(2))/(TrigTemp - Sa10(2)));
    end
    U2 = (1 - pt11) * ((1 - pt22) * Utility(Sa200, aa0, ag0, Wm, t + 2) + pt22 * Utility(Sa201, aa0, ag1, Wm, t + 2)) + pt11 * Utility(Sa211, aa0, ag1, Wm, t + 2);
else
    U2 = U1;
end

Vb = coefm(m, t) * U2;

%	utility of the current state
W = Utility(S(:, t), act(1), act(2), Wm, t);

%   Approximate value of the current state
value = -(W + lambda * Vb);
end