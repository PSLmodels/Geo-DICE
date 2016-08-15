
clc
clear;
initial;

global	K0 Tat0 Tlo0 Mat0 Mup0 Mlo0 a0 aa0 Sig0 Eland
global	S L A T gama t m coefm deltaT lambda bprice Fex
global	deltarf sai2 theta1 theta2
global  pt1 TP x TrigTemp Wm We
global  scenario action

	prompt = '1:Deterministic, 2:Tipping Point: ';
	scenario = input(prompt);

    if scenario ~= 1
        prompt = 'Tipping point Type: 1:CS=4, 2:CS=5, 3:CS=6, 4:C=0.75, 5:C=0.5, 6:C=0.25:';
        TP = input(prompt);
    else
        TP = 0;
    end

%   Number of runs
    if scenario == 1
        mmax = 100;
    else
        mmax = 200;
    end
   
%	Climate Sensitivity
    deltaT = 3;  
     
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
    Costs = zeros(mmax, 2, T);
%   Costs(m, 1, t) = Damage Cost;
%   Costs(m, 2, t) = Abatement Cost;

    Costs_Dm = zeros(mmax, T);
    Costs_Ab = zeros(mmax, T);

%   Optimal policy (GHG Reduction Rate %) "a" at state i at time t
    action = zeros(mmax, T);
    
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
    F0 = (deltarf * ((log(Mat0) - log(596.4)) / log(2)) + Fex(1));
    S(:, :, 1) = ones(mmax, 1) * [K0, Tat0, Tlo0, Mat0, Mup0, Mlo0, y0, Em0, F0];
    Costs(:, 1, 1) = 1 - 1/(1 + sai1 * Tat0 + sai2 * (Wm0 * Tat0) ^ 2);
    Costs(:, 2, 1) = theta1(1) * a0 ^ theta2;
    
for m = 1:1:mmax
    if scenario == 2
        %	Generating a sample path based on truncated-lognormal distribution
        rand1 = randn(1, T);
        We(m, :) = Wm0 + Wsd * rand1;
    end     
    for t = 1:1:T
%         Wm = We(m, t);
        Wm = 1;
        
        if t == T
            action(m, t) = action(m, t - 1);
            Vhat(t) = Rep_Utility(S(m, :, t), action(m, t), Wm, t);
            Vbar(t) = coefm(m, t) * Vhat(t);
            phi(t) = Vhat(t);
        else
            lb = 0;
            ub = 1;

            if t == 1   % Initial abatement rate
                lb = a0;
                ub = a0;
            end

            act0 = a0;
            myoptions = optimset('Display', 'off','FunValCheck','on','algorithm','sqp', 'MaxFunEvals', 1e5, 'TolX', 1e-9, 'TolCon', 1e-6, 'TolFun', 1e-4);
            [aopt, fval] = fmincon(@Rep_fun,act0,[],[],[],[],lb,ub,[], myoptions);
            Sa10 = Rep_NextState(S(m, :, t), aopt, t, Wm, 0);
            Sa11 = Rep_NextState(S(m, :, t), aopt, t, Wm, TP);
            if (x(m, t) == 0)
                if scenario == 1
                    pt1 = 0;
                else
                    pt1 = max(0, (min(Sa10(2), TrigTemp) - S(m, 2, t))/(TrigTemp - S(m, 2, t)));
                end
            else
                pt1 = 1;
            end
            U1 = (1 - pt1) * Rep_Utility(Sa10, aa0, Wm, t + 1) + pt1 * Rep_Utility(Sa11, aa0, Wm, t + 1);
              
            if (t < T-1)
                Sa200 = Rep_NextState(Sa10, aa0, t + 1, Wm, 0);
                Sa201 = Rep_NextState(Sa10, aa0, t + 1, Wm, TP);
                Sa211 = Rep_NextState(Sa11, aa0, t + 1, Wm, TP);
                
                if scenario == 1
                    pt2 = 0;
                else
                     pt2 = max(0, (min(Sa200(2), TrigTemp) - Sa10(2))/(TrigTemp - Sa10(2)));
                end
                U2 = (1 - pt1) * ((1 - pt2) * Rep_Utility(Sa200, aa0, Wm, t + 2) + pt2 * Rep_Utility(Sa201, aa0, Wm, t + 2)) + pt1 * Rep_Utility(Sa211, aa0, Wm, t + 2);
            else
                U2 = U1;
            end
            
            Vhat(t) = -fval;
            phi(t) = U2;
            Vbar(t) = coefm(m, t) * U2;
            action(m, t) = aopt;
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
        
            % Finding the next pre-decision state
            [S(m, :, t + 1), Costs(m, :, t)] = Rep_NextState(S(m, :, t), action(m, t), t, We(m, t), TPnxt); 
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
        OSU(m, t) = Rep_Utility(S(m, :, t), action(m, t), We(m, t), t);        
        %	Social Welfare
        TSW(m) = OSU(m, t) * lambda^(t - 1) + TSW(m);
        %   Costs
        Costs_Dm(m, t) = Costs(m, 1, t);
        Costs_Ab(m, t) = Costs(m, 2, t);
        %   Climate Parameters
        S_Temp(m, t) = S(m, 2, t);
        S_Con(m, t) = S(m, 4, t);
        S_RF(m, t) = S(m, 9, t);        
    end
end