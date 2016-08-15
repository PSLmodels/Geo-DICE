
clc
clear;
initial;

global	K0 Tat0 Tlo0 Mat0 Mup0 Mlo0 a0 g0 aa0 ag0 ag1 Sig0 Eland
global	S L A T gama t m coefm deltaT lambda bprice Fex
global	deltarf theta1 theta2 theta3
global  pt1 GE TP x TrigTemp Wm We sai1 sai2 sai3
global  scenario action actionG nuG Gcoeff Geffective
global	sai2temp sai2ocean sai2atmos GEStrategy GETrigTemp

    scenario = 0;
    GE = 0;
    GEStrategy = 0;
    SenAn = 0;
    SenPar = 0;
    GETrigTemp = 0;
    k3 = 1;
    WShock = 0;
    StochVar = 0;
    
%   Number of runs
    mmax = 100;
    
	prompt = '0:Deterministic, 1:Tipping Point, 2:Stochastic: ';
	scenario = input(prompt);
    if scenario ~= 0
        mmax = 200;
        prompt = 'Weather Shocks 0:off, 1:on: ';
        WShock = input(prompt); 
        prompt = '1:nuG, 2:deltaT: 3:nuG+deltaT: ';
        StochVar = input(prompt); 
    end
	prompt = '0:No Geoengineering, 1:Geoengineering: ';
    GE = input(prompt);
    if GE == 1
        prompt = '0:Unconstrained, 1:Reparation: ';
        GEStrategy = input(prompt);
    end
    if GEStrategy == 0
        if GE == 1 && scenario == 0
            prompt = '0:No SenAn, 1:SenAn: ';
            SenAn = input(prompt);
            if SenAn == 1
                k3 = 3;
                prompt = '1:DiscRate, 2:nuG, 3:Geffective, 4:Gcoeff, 5:Damcomp, 6:AbateCost, 7:Climate Sensitivity: ';
                SenPar = input(prompt);
                if SenPar == 7
                    k3 = 9;
                    mmax = 200;
                    prob = zeros(k3, mmax);
                end
            end
        end
    else
        if scenario ~= 1
            prompt = 'Trigger Temperature: ';
            GETrigTemp = input(prompt);
        end
    end
    if scenario == 1
        prompt = 'Tipping point Type: 1:CS=4, 2:CS=5, 3:CS=6, 4:C=0.75, 5:C=0.5, 6:C=0.25, 7:dm = 1%:';
        TP = input(prompt);
    else
        TP = 0;
        TPnxt = 0;
    end


    
    %	Climate Sensitivity
    deltaT = 3;
    
    % Lower bound
    deltaT_l = 0.5;
    % Upper bound
    deltaT_u = 20;
    % Mean
	deltaT_m = log(deltaT);
    % Standard Deviation
	deltaT_sig = 0.7;
	deltaT_path = deltaT * ones(mmax, T);
    
    deltaTk = [1.1, 0.7; 1.1, 0.7; 1.1, 0.7; 1.1, 0.7; 1.1, 0.7; 1.1, 0.7; 1.1, 0.7; 1.1, 0.7; 1.1, 0.7]; 
    
    rok = 1.5 * ones(1, 9);
    nuGk = 0.03 * ones(1, 9);
    Geffectivek = 1.0 * ones(1, 9);
    Gcoeffk = 1.0 * ones(1, 9);
    Damk = [0.8, 0.1, 0.1; 0.8, 0.1, 0.1; 0.8, 0.1, 0.1; 0.8, 0.1, 0.1; 0.8, 0.1, 0.1; 0.8, 0.1, 0.1; 0.8, 0.1, 0.1; 0.8, 0.1, 0.1; 0.8, 0.1, 0.1];
    AbateCost = 2.8 * ones(1, 9);
    
    if SenPar == 1
    %   Discount Rate
        rok = [1.5, 0.75, 3.0];
    elseif SenPar == 2
    %   SGE Damage
        nuGk = [0.03, 0.015, 0.06];
    elseif SenPar == 3
    %   Effectiveness of SGE
        Geffectivek = [1.0, 0.5, 2.0];
    elseif SenPar == 4
    %   Coefficient of SGE cost function
        Gcoeffk = [1.0, 0.5, 2.0];
    elseif SenPar == 5
    %   Damage cost function
        Damk = [0.8, 0.1, 0.1; 0.5, 0.4, 0.1; 0.5, 0.1, 0.4];
    elseif SenPar == 6
    %   Damage cost function
        AbateCost = [2.8, 1.4, 5.6];
    elseif SenPar == 7
    %   Climate Sensitivity
        deltaTk = [0.9, 0.3; 0.9, 0.5; 0.9, 0.7; 1.1, 0.3; 1.1, 0.5; 1.1, 0.7; 1.3 0.3; 1.3, 0.5; 1.3, 0.7];
        mmax = 200;
    end
%%Calibration of new damage coefficients
%   This is DICE's damage specification in initial period, based just on
%   temperature
%   Initial temp = 0.7307 degrees above preindustrial
    origdamage = sai1 * .7307 + sai2 * .7307 ^ sai3;

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
     
%   State space:
    S = zeros(10, T);
    %S(1, t): Capital "K" at time t
    %S(2, t): Atmospheric Temperature "Tat" at time t
    %S(3, t): Lower ocean Temperature "Tlo" at time t
    %S(4, t): Atmospheric Carbon Concentration "Mat" at time t
    %S(5, t): Upper Ocean Carbon Concentration "Mup" at time t
    %S(6, t): Lower Ocean Carbon Concentration "Mlo" at time t
    %S(7, t): Gross Economic Output "Y" at time t
    %S(8, t): Emissions "E" at time t
    %S(9, t): Radiative Forcings "F" at time t
    %S(10, t): tipping point reached at time t
    
    S_Temp = zeros(mmax, T);
    K_Temp = zeros(3, T);
    
    S_Con = zeros(mmax, T);
    K_Con = zeros(3, T);
    
    S_RF = zeros(mmax, T);
    K_RF = zeros(3, T);
    
%   Costs space:
    Costs = zeros(3, T);
%   Costs(1, t) = Damage Cost;
%   Costs(2, t) = Abatement Cost;
%   Costs(3, t) = Geoengineering Cost;

    Costs_Dm = zeros(mmax, T);
    Costs_Ab = zeros(mmax, T);
    Costs_Ge = zeros(mmax, T);

%   Optimal policy (GHG Reduction Rate %) "a" at state i at time t
    action = zeros(mmax, T);
    K_act =zeros(3, T);
    
    actionG = zeros(mmax, T);
    K_actg= zeros(3, T);
    
%   Carbon Price
    Cprice = zeros(mmax, T);
    K_Cprice= zeros(3, T);

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
    Wm = 1;
    Wm0 = 1;
    Wsd = 0.0068;
    We = Wm0 * ones(mmax, T);
    
%	Extreme Event
    % Binary variable detecting the Tipping Point
    x = zeros(mmax, T);
    % Probability of the Tipping Point
    pt = zeros(mmax, T);
    % Trigger Temperature
    TrigTemp = 4.27;
        

    
for k = 1:1:k3
        
%   Social time preference factor
    ro = rok(k) * ones(1,T);
    R = zeros(1,T);
    R(1) = 1;
    for i = 2:T
        R(i) = R(i - 1)/(1 + ro(i - 1) / 100)^10;
    end
    alpha = 2;
    lambda = 1/(1+ro(1)/100)^10;
    
%   In new specification, 80% of this should be directly from temperature
    sai2temp = Damk(k, 1) * origdamage / 0.7307 ^ 2;
%   Initial concentration = 1255; preindustrial = 1094; difference = 161
    sai2ocean = Damk(k, 2) * origdamage / 161 ^ 2;
%   Initial concentration = 808.9; preindustrial = 596.4; difference = 212.5
    sai2atmos = Damk(k, 3) * origdamage / 212.5 ^ 2;
    
    %   SGE Damage
    nuG = nuGk(k);
%   Effectiveness of SGE
    Geffective = Geffectivek(k);
%   Coefficient of SGE cost function
    Gcoeff = Gcoeffk(k);
    
%   Coefficient of Abatement cost function
    theta2 = AbateCost(k);

%   Initial State
    y0 = A(1) * (K0 ^ gama) * L(1) ^ (1 - gama);
    Em0 = Sig0 * (1 - a0) * y0 + Eland(1);
    F0 = (deltarf * ((log(Mat0) - log(596.4)) / log(2)) + Fex(1)) * (1 - Geffective * g0);
    S(:, 1) = [K0, Tat0, Tlo0, Mat0, Mup0, Mlo0, y0, Em0, F0, 0];
    Costs(1, 1) = 1 - 1/((1 + sai2temp * (Wm0 * Tat0) ^ 2 + sai2ocean * (Mup0 - 1094) ^ 2 + sai2atmos * (Mat0 - 596.4) ^ 2) * (1 + nuG * g0 ^ 2));
    Costs(2, 1) = theta1(1) * a0 ^ theta2;
    Costs(3, 1) = Gcoeff * theta1(1) * g0 ^ theta3;
    
    for m = 1:1:mmax
        %	Generating a sample path for stochastic parameter
        if WShock == 1
            randW = randn(1, T);
            We(m, :) = Wm0 + Wsd * rand1;
        end
        %	Generating a sample path for uncertain parameters
        rand1 = rand(1, T);
        rand2 = rand(1)*ones(1, T);
        if StochVar == 1
            nuG_path(m, 1:T) = lognorm_trunc(rand1, nuG_l, nuG_u, nuG_m, nuG_sig );
        elseif StochVar == 2 || SenPar == 7
            deltaT_m = deltaTk(k, 1);
            deltaT_sig = deltaTk(k, 2);
            deltaT_path(m, 1:T) = lognorm_trunc(rand2, deltaT_l, deltaT_u, deltaT_m, deltaT_sig );
        elseif StochVar == 3
            nuG_path(m, 1:T) = lognorm_trunc(rand1, nuG_l, nuG_u, nuG_m, nuG_sig );
            deltaT_path(m, 1:T) = lognorm_trunc(rand2, deltaT_l, deltaT_u, deltaT_m, deltaT_sig );
        end
             
        for t = 1:1:T 
            Wm = We(m, t);
            nuG = nuG_path(m, t);
            deltaT = deltaT_path(m, t);
            if t == T
                action(m, t) = action(m, t - 1);
                actionG(m, t) = actionG(m, t - 1);
                Vhat(t) = Utility(S(:, t), action(m, t), actionG(m, t), Wm, t);
                Vbar(t) = coefm(m, t) * Vhat(t);
                phi(t) = Vhat(t);
            else
                if GE == 0      % No Geoengineering Case
                    lb = [0; 0];
                    ub = [1; 0];
                    ag1 = 0;
                else            % Geoengineering Case
                    lb = [0; 0];
                    if GEStrategy == 1 && x(m, t) == 0
                        ub = [1; 0];
                        ag1 = 0.1;
                    else
                        ub = [1; 1.5];
                        ag1 = 0;
                    end
                end

                if t == 1   % Initial abatement rate
                    lb(1) = a0;
                    ub(1) = a0;              
                end

                act0 = [a0; g0];
                myoptions = optimset('Display', 'off','FunValCheck','on','algorithm','sqp', 'MaxFunEvals', 1e5, 'TolX', 1e-9, 'TolCon', 1e-6, 'TolFun', 1e-4);
                [aopt, fval] = fmincon(@fun, act0,[],[],[],[],lb,ub,[], myoptions);
                Sa10 = NextState(S(:, t), aopt(1), aopt(2), t, Wm, 0);
                Sa11 = NextState(S(:, t), aopt(1), aopt(2), t, Wm, TP);
                if (x(m, t) == 0)
                    if scenario ~= 1
                        pt1 = GE * GEStrategy * (S(2, t) > GETrigTemp);
                    else
                        pt1 = max(0, (min(Sa10(2), TrigTemp) - S(2, t))/(TrigTemp - S(2, t)));
                    end
                else
                    pt1 = 1;
                end
                U1 = (1 - pt1) * Utility(Sa10, aa0, ag0, Wm, t + 1) + pt1 * Utility(Sa11, aa0, ag1, Wm, t + 1);

                if (t < T-1)
                    Sa200 = NextState(Sa10, aa0, ag0, t + 1, Wm, 0);
                    Sa201 = NextState(Sa10, aa0, ag0, t + 1, Wm, TP);
                    Sa211 = NextState(Sa11, aa0, ag1, t + 1, Wm, TP);

                    if scenario ~= 1
                        pt2 = GE * GEStrategy * (Sa10(2) > GETrigTemp);
                    else
                         pt2 = max(0, (min(Sa200(2), TrigTemp) - Sa10(2))/(TrigTemp - Sa10(2)));
                    end
                    U2 = (1 - pt1) * ((1 - pt2) * Utility(Sa200, aa0, ag0, Wm, t + 2) + pt2 * Utility(Sa201, aa0, ag1, Wm, t + 2)) + pt1 * Utility(Sa211, aa0, ag1, Wm, t + 2);
                else
                    U2 = U1;
                end

                Vhat(t) = -fval;
                phi(t) = U2;
                Vbar(t) = coefm(m, t) * U2;
                action(m, t) = aopt(1);
                actionG(m, t) = aopt(2);
                pt(m,t) = pt1;

                % Tipping point in the next state
                if scenario == 1
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
                else
                    if GEStrategy == 1
                        if (x(m, t) ~= 0)
                            x(m, t + 1) = x(m, t);
                        elseif (S(2, t) > GETrigTemp)
                            x(m, t + 1) = 1;
                        else
                            x(m, t + 1) = 0;
                        end
                    end
                end

                % Finding the next pre-decision state
                [S(:, t + 1), Costs(:, t + 1)] = NextState(S(:, t), action(m, t), actionG(m, t), t, We(m, t + 1), TPnxt); 
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
            OSU(m, t) = Utility(S(:, t), action(m, t), actionG(m, t), We(m, t), t);        
            %	Social Welfare
            TSW(m) = OSU(m, t) * lambda^(t - 1) + TSW(m);
            %   Costs
            Costs_Dm(m, t) = Costs(1, t);
            Costs_Ab(m, t) = Costs(2, t);
            Costs_Ge(m, t) = Costs(3, t);
            %   Climate Parameters
            S_Temp(m, t) = S(2, t);
            S_Con(m, t) = S(4, t);
            S_RF(m, t) = S(9, t);
        end
    if SenPar == 7
        prob(k, m) = deltaT;
    end
    end
    if SenPar == 7
        K_Temp(k, :) = median(S_Temp);
        K_Con(k, :) = median(S_Con);
        K_RF(k, :) = median(S_RF);
        K_act(k, :) = median(action);
        K_actg(k, :) = median(actionG);
        K_Cprice(k, :) = median(Cprice);
        
    else
        for t = 1:1:T
            K_Temp(k, t) = S_Temp(m, t);
            K_Con(k, t) = S_Con(m, t);
            K_RF(k, t) = S_RF(m, t);
            K_act(k, t) = action(m, t);
            K_actg(k, t) = actionG(m, t);
            K_Cprice(k, t) = Cprice(m, t);
        end
    end
end