
clc
clear;
initial;

global	K0 Tat0 Tlo0 Mat0 Mup0 Mlo0 a0 g0 aa0 ag0 Sig0 Eland
global	S L A T gama t m coefm deltaT lambda bprice Fex
global	deltarf sai2 sai3 theta1 theta2 theta3
global  pt1 GE TP x TrigTemp Wm We
global  scenario action actionG nuG Gcoeff Geffective nuG_l nuG_u
global	sai2temp sai2ocean sai2atmos GEStrategy

    GEStrategy = 0;
	prompt = '1:Deterministic, 2:Tipping Point, 3:Bayesian : ';
	scenario = input(prompt);
    
	if scenario == 3
        GE = 1;
    else
        prompt = '0:No Geoengineering, 1:Geoengineering: ';
        GE = input(prompt);
	end
    
    if GE == 1
        prompt = '0:Insurance, 1:Response: ';
        GEStrategy = input(prompt);
    end
    if scenario ~= 1
        prompt = 'Tipping point Type: 1:CS=4, 2:CS=5, 3:CS=6, 4:C=0.75, 5:C=0.5, 6:C=0.25, 7:dm = 10%:';
        TP = input(prompt);
    else
        TP = 0;
    end

%   Number of runs
    if scenario == 1
        mmax = 100;
    else
        mmax = 1000;
    end
   
%	Climate Sensitivity
    deltaT = 3;
   
%	Uncertainty in Geoengineering
    nuG = 0.03;
    
    % Lower bound
    nuG_l = 0.01;
    % Upper bound
    nuG_u = 0.06;
    % Mean
     nuG_m = log(nuG);
    % Standard Deviation
    nuG_sig = 1;
    
    nuG_path = nuG * ones(mmax, T);
    nuG_pdf = zeros(100, T);
    nuG_cdf = zeros(100, T); 
     
%   State space:
    S = zeros(mmax, 9, T);
    %S(m, 1, t): Capital "K" at time t
    %S(m, 2, t): Atmospheric Temperature "Tat" at time t
    %S(m, 3, t): Lower ocean Temperature "Tlo" at time t
    %S(m, 4, t): Atmospheric Carbon Concentration "Mat" at time t
    %S(m, 5, t): Upper Ocean Carbon Concentration "Mup" at time t
    %S(m, 6, t): Lower Ocean Carbon Concentration "Mlo" at time t
    %S(m, 7, t): Gross Economic Output "Y" at time t
    %S(m, 8, t): Emissions "E" at time t
    %S(m, 9, t): Radiative Forcings "F" at time t
    
    S_Temp = zeros(mmax, T);
    S_Con = zeros(mmax, T);
    S_RF = zeros(mmax, T);
    
%   Costs space:
    Costs = zeros(mmax, 3, T);
%   Costs(m, 1, t) = Damage Cost;
%   Costs(m, 2, t) = Abatement Cost;
%   Costs(m, 3, t) = Geoengineering Cost;

    Costs_Dm = zeros(mmax, T);
    Costs_Ab = zeros(mmax, T);
    Costs_Ge = zeros(mmax, T);

%   Optimal policy (GHG Reduction Rate %) "a" at state i at time t
    action = zeros(mmax, T);
    actionG = zeros(mmax, T);
    
%   Carbon Price
    Cprice = zeros(mmax, T);

%	Approximate value "V-" of the post-decision stateat time t (Vbar)
    Vbar = zeros(1, T);

%	True value "V^" at time t (Vhat)
    Vhat = zeros(1, T + 1);

%	Error in estimating the Value Function    
    error = zeros(mmax, T);
    
%	Coefficient of aaproximate value function    
    coefm = zeros(mmax + 1, T);
    
%	Compounents of aaproximate value function
    phi = zeros(1, T); 
    
%	Optimal Social Utiltiy
    OSU = zeros(mmax, T);
    
%	Total Social Welfare
    TSW = zeros(mmax, 1);

%	Weather Shocks
    Wm0 = 1;
    Wsd = 0.0068;
    We = ones(mmax, T);
    
%	Extreme Event
    % Binary variable detecting the Tipping Point
    x = zeros(mmax, T);
    % Probability of the Tipping Point
    pt = zeros(mmax, T);
    % Trigger Temperature
    TrigTemp = 4.27;
        
%   Initial State
    y0 = A(1) * (K0 ^ gama) * L(1) ^ (1 - gama);
    Em0 = Sig0 * (1 - a0) * y0 + Eland(1);
    F0 = (deltarf * ((log(Mat0) - log(596.4)) / log(2)) + Fex(1)) * (1 - Geffective * g0);
    S(:, :, 1) = ones(mmax, 1) * [K0, Tat0, Tlo0, Mat0, Mup0, Mlo0, y0, Em0, F0];
    Costs(:, 1, 1) = 1 - 1/((1 + sai2temp * (Wm0 * Tat0) ^ 2 + sai2ocean * (Mup0 - 1094) ^ 2 + sai2atmos * (Mat0 - 596.4) ^ 2) * (1 + nuG * g0 ^ 2));
    Costs(:, 2, 1) = theta1(1) * a0 ^ theta2;
    Costs(:, 3, 1) = Gcoeff * theta1(1) * g0 ^ theta3;
    
for m = 1:1:mmax
    if scenario == 2
        %	Generating a sample path based on truncated-lognormal distribution
        rand1 = randn(1, T);
        rand2 = rand(1, T);
        We(m, :) = Wm0 + Wsd * rand1;
        nuG_path(m, 1:T) = lognorm_trunc(rand2, nuG_l, nuG_u, nuG_m, nuG_sig );
    end     
    for t = 1:1:T
        Wm = We(m, t);
        nuG = nuG_path(m, t);
        if t == T
            action(m, t) = action(m, t - 1);
            actionG(m, t) = actionG(m, t - 1);
            Vhat(t) = T_Utility(S(m, :, t), action(m, t), actionG(m, t), Wm, t);
            Vbar(t) = coefm(m, t) * Vhat(t);
            phi(t) = Vhat(t);
        else
            if GE == 0	% No Geoengineering Case
                lb = [0; 0];
                ub = [1; 0];
            else        % Geoengineering Case
                lb = [0; 0];
                if GEStrategy == 1 && x(m, t) == 0
                    ub = [1, 0];
                else
                    ub = [1; 1.5];
                end
            end

            if t == 1   % Initial abatement rate
                lb(1) = a0;
                ub(1) = a0;
            end

            act0 = [a0; g0];
            myoptions = optimset('Display', 'off','FunValCheck','on','algorithm','sqp', 'MaxFunEvals', 1e5, 'TolX', 1e-9, 'TolCon', 1e-6, 'TolFun', 1e-4);
            [aopt, fval] = fmincon(@T_fun,act0,[],[],[],[],lb,ub,[], myoptions);
            Sa10 = T_NextState(S(m, :, t), aopt(1), aopt(2), t, deltaT, nuG, Wm, 0);
            Sa11 = T_NextState(S(m, :, t), aopt(1), aopt(2), t, deltaT, nuG, Wm, TP);
            if (x(m, t) == 0)
                if scenario == 1
                    pt1 = 0;
                else
                    pt1 = max(0, (min(Sa10(2), TrigTemp) - S(m, 2, t))/(TrigTemp - S(m, 2, t)));
                end
            else
                pt1 = 1;
            end
            U1 = (1 - pt1) * T_Utility(Sa10, aa0, ag0, Wm, t + 1) + pt1 * T_Utility(Sa11, aa0, ag0, Wm, t + 1);
              
            if (t < T-1)
                Sa200 = T_NextState(Sa10, aa0, ag0, t + 1, deltaT, nuG, Wm, 0);
                Sa201 = T_NextState(Sa10, aa0, ag0, t + 1, deltaT, nuG, Wm, TP);
                Sa211 = T_NextState(Sa11, aa0, ag0, t + 1, deltaT, nuG, Wm, TP);
                
                if scenario == 1
                    pt2 = 0;
                else
                     pt2 = max(0, (min(Sa200(2), TrigTemp) - Sa10(2))/(TrigTemp - Sa10(2)));
                end
                U2 = (1 - pt1) * ((1 - pt2) * T_Utility(Sa200, aa0, ag0, Wm, t + 2) + pt2 * T_Utility(Sa201, aa0, ag0, Wm, t + 2)) + pt1 * T_Utility(Sa211, aa0, ag0, Wm, t + 2);
            else
                U2 = U1;
            end
            
            Vhat(t) = -fval;
            phi(t) = U2;
            Vbar(t) = coefm(m, t) * U2;
            action(m, t) = aopt(1);
            actionG(m, t) = aopt(2);
            pt(m,t) = pt1;

            % Extreme event in the next state
            if (x(m, t) ~= 0)
                x(m, t + 1) = x(m, t);
                TPnxt = TP;
            elseif (pt(m, t) > rand)
                x(m, t + 1) = 1;
                TPnxt = TP;
            else
                x(m, t + 1) = 0;
                TPnxt = 0;
            end
            
            %   Bayesian update for nuG
            if scenario == 3
                if (t == 1)
                    for i = 1:1:100
                        yi = nuG_l + (nuG_u - nuG_l) / 99 * (i - 1);
                        nuG_cdf(i, t) = lognorm_trunc_inv(yi, nuG_l, nuG_u, nuG_m, nuG_sig);
                        if (i > 1)
                            nuG_pdf(i, t) = nuG_cdf(i, t) - nuG_cdf(i - 1, t);
                        else
                            nuG_pdf(i, t) = nuG_cdf(i, t);
                        end     
                    end
                    pri_pdf = nuG_pdf(:, t);
                else               
                    pri_pdf = nuG_pdf(:, t - 1);
                end
                [post_pdf, post_cdf, nuG_rand] = nuG_update(pri_pdf, m, t, rand);
                nuG_pdf(:, t) = post_pdf;
                nuG_cdf(:, t) = post_cdf;
                nuG_path(m, t + 1) = nuG_rand;
            end
                
            % Finding the next pre-decision state
            [S(m, :, t + 1), Costs(m, :, t + 1)] = T_NextState(S(m, :, t), action(m, t), actionG(m, t), t, deltaT, nuG_path(m, t + 1), We(m, t + 1), TPnxt); 
            Cprice(m, t) = bprice(t) * 1000 * action(m, t)^(theta2 - 1);
        end
        %	Updating coefficient of approximate value function
        if (t > 1)    
            %	Value approximation Error
            error(m, t - 1) = Vbar(t - 1) - Vhat(t);

            %   Bellman Update
            lm = phi(t - 1)^-2;
            coefm(m + 1, t - 1) = coefm(m, t - 1) - lm * error(m, t - 1) * phi(t - 1);
        end
        %   Social Utility
        OSU(m, t) = T_Utility(S(m, :, t), action(m, t), actionG(m, t), We(m, t), t);        
        %	Social Welfare
        TSW(m) = OSU(m, t) * lambda^(t - 1) + TSW(m);
        %   Costs
        Costs_Dm(m, t) = Costs(m, 1, t);
        Costs_Ab(m, t) = Costs(m, 2, t);
        Costs_Ge(m, t) = Costs(m, 3, t);
        %   Climate Parameters
        S_Temp(m, t) = S(m, 2, t);
        S_Con(m, t) = S(m, 4, t);
        S_RF(m, t) = S(m, 9, t);
        
    end
    if x(m, 60) == 0
        nuG_pdf_NT = nuG_pdf(:, :);
    end 
end

% Export the results into excel
if scenario == 1
    id1 = 3;
    id2 = 1;
    Temp = array2table(S_Temp(mmax, :));
    Con = array2table(S_Con(mmax, :));
    RF = array2table(S_RF(mmax, :));
    Act = array2table(action(mmax, :));
    ActG = array2table(actionG(mmax, :));
    CP = array2table(Cprice(mmax, :));
    Util = array2table(OSU(mmax, :));
    if GE == 1
        name = 'DGE';
    else
        name = 'DNG';
    end
else
    id1 = 3;
    id2 = mmax + 5;
    Temp = array2table(S_Temp);
    Con = array2table(S_Con);
    RF = array2table(S_RF);
    Act = array2table(action);
    ActG = array2table(actionG);
    CP = array2table(Cprice);
    Util = array2table(OSU);
    EE = array2table(x(:, 1:T));
    if scenario == 2
        if GE == 1
            name = 'TGE';
        else
            name = 'TNG';
        end
    else
        name = 'BGE';
    end
end
        
filename = 'blank.xlsx';
% Temperature
writetable(Temp, filename, 'Sheet', name, 'Range', 'B3', 'WriteVariableNames', false);
% Concentration
rng = strcat('B', num2str(id1 + 1 * id2));
writetable(Con, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% Radiative Forcing
rng = strcat('B', num2str(id1 + 2 * id2));
writetable(RF, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% Abatement
rng = strcat('B', num2str(id1 + 3 * id2));
writetable(Act, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% Geoengineering
rng = strcat('B', num2str(id1 + 4 * id2));
writetable(ActG, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% Carbon Price
rng = strcat('B', num2str(id1 + 5 * id2));
writetable(CP, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% Utility
rng = strcat('B', num2str(id1 + 6 * id2));
writetable(Util, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);

if scenario ~= 1
    % Extreme Event
    rng = strcat('B', num2str(id1 + 7 * id2));
    writetable(EE, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
end
