function [ prob, cum, ran ] = nuG_update(pr4, m4, t4, rand4)
global nuG_l nuG_u S action actionG x deltaT TP We TrigTemp

%Bayesian Update

y = zeros(1, 100);
ptt1 = zeros(1, 100);
ptt2 = zeros(1, 100);
prob = zeros(1, 100);
cum = zeros(1, 100);
sum = 0;
if t4>10 && x(m4, t4)==0
    ss=1;
end

for i = 1:1:100
    y(i) = nuG_l + (nuG_u - nuG_l) / 99 * (i - 1);
    if (t4 ~= 1)
        Sn = T_NextState(S(m4, :, t4 - 1), action(m4, t4 - 1), actionG(m4, t4 - 1), t4 - 1, deltaT, y(i), We(m4, t4), TP * x(m4, t4));
        ptt1(i) = max(0, (min(Sn(2), TrigTemp) - S(m4, 2, t4 - 1))/(TrigTemp - S(m4, 2, t4 - 1)));
        ptt2(i) = x(m4, t4) * x(m4, t4 + 1) + (1 - x(m4, t4)) * (ptt1(i) ^ x(m4, t4 + 1)) * ((1 - ptt1(i)) ^ (1 - x(m4, t4 + 1)));
        sum = sum + ptt2(i) * pr4(i);
    end
end


for i = 1:1:100
    if (sum ~= 0)
        prob(i) = ptt2(i) * pr4(i) / sum;
    else
        prob(i) = pr4(i);
    end    
    if (i > 1)
        cum(i) = cum(i - 1) + prob(i);
    else
        cum(i) = prob(i);
    end
    if (cum(i) > rand4)
        ran = y(i - 1);
        rand4 = 2;
    end       
end