
clc
clear;
initial;

global	K0 Tat0 Tlo0 Mat0 Mup0 Mlo0 a0 g0 aa0 ag0 ag1 Sig0 Eland
global	S L A T gama t m coefm deltaT lambda bprice Fex
global	deltarf sai2 sai3 theta1 theta2 theta3
global  pt1 GE TP x TrigTemp Wm We
global  scenario action actionG nuG Gcoeff Geffective
global	sai2temp sai2ocean sai2atmos GEStrategy

    GEStrategy = 0;
	prompt = '1:Deterministic, 2:Tipping Point: ';
	scenario = input(prompt);
	prompt = '0:No Geoengineering, 1:Geoengineering: ';
    GE = input(prompt);
    if GE == 1 && scenario ==2
        prompt = '0:Unconstrained, 1:Reparation: ';
        GEStrategy = input(prompt);
    end
    if scenario ~= 1
        prompt = 'Tipping point Type: 1:CS=4, 2:CS=5, 3:CS=6, 4:C=0.75, 5:C=0.5, 6:C=0.25, 7:dm = 1%:';
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
    %S(10, t): Tipping point stateus at time t
    
    S_Temp = zeros(mmax, T);
    S_Con = zeros(mmax, T);
    S_RF = zeros(mmax, T);
    
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
    S(:, 1) = [K0, Tat0, Tlo0, Mat0, Mup0, Mlo0, y0, Em0, F0, 0];
    Costs(1, 1) = 1 - 1/((1 + sai2temp * (Wm0 * Tat0) ^ 2 + sai2ocean * (Mup0 - 1094) ^ 2 + sai2atmos * (Mat0 - 596.4) ^ 2) * (1 + nuG * g0 ^ 2));
    Costs(2, 1) = theta1(1) * a0 ^ theta2;
    Costs(3, 1) = Gcoeff * theta1(1) * g0 ^ theta3;
    
%      Geffective = 0;
for m = 1:1:mmax
    if scenario == 2
        %	Generating a sample path based on truncated-lognormal distribution
        rand1 = randn(1, T);
        rand2 = rand(1, T);
%         We(m, :) = Wm0 + Wsd * rand1;
%         nuG_path(m, 1:T) = lognorm_trunc(rand2, nuG_l, nuG_u, nuG_m, nuG_sig );
    end     
    for t = 1:1:T        
        Wm = We(m, t);
        nuG = nuG_path(m, t);
        if t == T
            action(m, t) = action(m, t - 1);
            actionG(m, t) = actionG(m, t - 1);
            Vhat(t) = T_Utility(S(:, t), action(m, t), actionG(m, t), Wm, t);
            Vbar(t) = coefm(m, t) * Vhat(t);
            phi(t) = Vhat(t);
        else
            if GE == 0	% No Geoengineering Case
                lb = [0; 0];
                ub = [1; 0];
                ag1 = 0;
            else        % Geoengineering Case
                lb = [0; 0];
                if GEStrategy == 1 && x(m, t) == 0
                    if TP == 2
                        GEx = [0,0.107434584490395,0.140499937319662,0.179012611699957,0.222135175180202,0.269043198281727,0.318675828070613,0.341933180113276,0.368814673275469,0.410183630394222,0.449466528830037,0.482362631024878,0.511112711055709,0.534741048219690,0.551992472094033,0.561320109314122,0.561158519948835,0.549648994041184,0.546080427369739,0.531835991952055,0.504018806909382,0.482523423921699,0.463942740932249,0.447898697058971,0.433954011542929,0.421695942431020,0.410807077455643,0.401001712114981,0.392168022377709,0.384030211171632,0.376336220624992,0.369282277704972,0.362292085340308,0.356048596004393,0.349821706282337,0.343542339584937,0.337855797642749,0.332341344518111,0.327378775711458,0.320743637031588,0.315960445720148,0.310558263065344,0.305443717190098,0.299647988353383,0.294161116935071,0.288493611675500,0.282512470027056,0.275728226999515,0.269004054484226,0.260712907288535,0.252879581151832,0.242580191349021,0.230954214752703,0.217984915226522,0.201243457521829,0.181698275600982,0.156563826963907,0,0,0];
                    elseif TP == 5
                        GEx = [0,0.107874433869796,0.139691217389115,0.176362309432007,0.217101139653028,0.261085418578513,0.307412687048827,0.343302271519424,0.369140017162491,0.411705014280397,0.452550745074851,0.487657647492930,0.519459992737786,0.547285629623799,0.570318453879630,0.587430700269060,0.597460384521840,0.598822804594195,0.589685299628902,0.578284689722598,0.553122372898349,0.527137246342213,0.504200050487731,0.484267807129318,0.466745217832522,0.451373583051741,0.437730129311616,0.425334203958543,0.414147953336748,0.403794818546084,0.394183359993001,0.385060025710366,0.376332315054719,0.368370459998334,0.360505330262459,0.352982998724717,0.345766845561754,0.338996899298213,0.331995427147212,0.325044191964469,0.318847558749195,0.312237946397734,0.306263962880220,0.299614834559796,0.292609889742286,0.286360108582102,0.280412233527236,0.272641225523983,0.266291173534864,0.256533490736244,0.247643989445073,0.238149657187816,0.226375621170583,0.213175196362304,0.197040909463361,0.178913416275185,0.156002586747024,0,0,0];
                    elseif TP == 7
                        GEx = [0,0.216137005957665,0.259597401613045,0.309709807611373,0.365503385585065,0.398610434427858,0.366676612822494,0.337940756153797,0.312127336771672,0.317002680981912,0.346397975349267,0.369322416456496,0.387427618456471,0.399468815965731,0.403816382789668,0.398699993923124,0.395136848763170,0.395990509100569,0.396319469452275,0.387997036931702,0.370564171014800,0.357030149155434,0.345731143274946,0.336172921161989,0.328057532252538,0.321044939318711,0.315060760425527,0.309611966577275,0.304915627523532,0.300463257536970,0.296573088379362,0.292787643346475,0.289219226154214,0.285995519729736,0.282423722885242,0.279754969063726,0.276722313750419,0.273630310621410,0.270866266890900,0.267247102396828,0.264643511551467,0.261437143907404,0.258880325301679,0.254557851499634,0.251764332563209,0.247674782635341,0.243501349004839,0.238681615175159,0.234119700008785,0.228569154291758,0.221678321678322,0.213780628275803,0.204249333175160,0.193499347035369,0.179452508559149,0.162336799947444,0.139131155303030,0,0,0];
                    end
                    ub = [1, 0];
                     %ag1 = GEx(t+1);
                    ag1 = 0.1;
                else
                    ub = [1; 1.5];
                    ag1 = 0;
                end
            end

            if t == 1   % Initial abatement rate
                lb(1) = a0;
                ub(1) = a0;
%                 lb(1) = g0;
%                 ub(1) = g0;                
            end

            act0 = [a0; g0];
            myoptions = optimset('Display', 'off','FunValCheck','on','algorithm','sqp', 'MaxFunEvals', 1e5, 'TolX', 1e-9, 'TolCon', 1e-6, 'TolFun', 1e-4);
            [aopt, fval] = fmincon(@T_fun,act0,[],[],[],[],lb,ub,[], myoptions);
            Sa10 = T_NextState(S(:, t), aopt(1), aopt(2), t, deltaT, nuG, Wm, 0);
            Sa11 = T_NextState(S(:, t), aopt(1), aopt(2), t, deltaT, nuG, Wm, TP);
            if (x(m, t) == 0)
                if scenario == 1
                    pt1 = 0;
                else
                    pt1 = max(0, (min(Sa10(2), TrigTemp) - S(2, t))/(TrigTemp - S(2, t)));
                end
            else
                pt1 = 1;
            end
            U1 = (1 - pt1) * T_Utility(Sa10, aa0, ag0, Wm, t + 1) + pt1 * T_Utility(Sa11, aa0, ag1, Wm, t + 1);
              
            if (t < T-1)
                Sa200 = T_NextState(Sa10, aa0, ag0, t + 1, deltaT, nuG, Wm, 0);
                Sa201 = T_NextState(Sa10, aa0, ag0, t + 1, deltaT, nuG, Wm, TP);
                Sa211 = T_NextState(Sa11, aa0, ag1, t + 1, deltaT, nuG, Wm, TP);
                
                if scenario == 1
                    pt2 = 0;
                else
                     pt2 = max(0, (min(Sa200(2), TrigTemp) - Sa10(2))/(TrigTemp - Sa10(2)));
                end
                U2 = (1 - pt1) * ((1 - pt2) * T_Utility(Sa200, aa0, ag0, Wm, t + 2) + pt2 * T_Utility(Sa201, aa0, ag1, Wm, t + 2)) + pt1 * T_Utility(Sa211, aa0, ag1, Wm, t + 2);
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
        
            % Finding the next pre-decision state
            [S(:, t + 1), Costs(:, t + 1)] = T_NextState(S(:, t), action(m, t), actionG(m, t), t, deltaT, nuG_path(m, t + 1), We(m, t + 1), TPnxt); 
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
        OSU(m, t) = T_Utility(S(:, t), action(m, t), actionG(m, t), We(m, t), t);        
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
    if x(m, 50) == 0
        mtst = m;
    end
end

% %coefx = 
% 
% % Export the results into excel
% if scenario == 1
%     id1 = 3;
%     id2 = 1;
%     Temp = array2table(S_Temp(mmax, :));
%     Con = array2table(S_Con(mmax, :));
%     RF = array2table(S_RF(mmax, :));
%     Act = array2table(action(mmax, :));
%     ActG = array2table(actionG(mmax, :));
%     CP = array2table(Cprice(mmax, :));
%     Util = array2table(OSU(mmax, :));
%     if GE == 1
%         name = 'DGE';
%     else
%         name = 'DNG';
%     end
% else
%     id1 = 3;
%     id2 = mmax + 5;
%     Temp = array2table(S_Temp);
%     Con = array2table(S_Con);
%     RF = array2table(S_RF);
%     Act = array2table(action);
%     ActG = array2table(actionG);
%     CP = array2table(Cprice);
%     Util = array2table(OSU);
%     EE = array2table(x(:, 1:T));
% 	if GE == 1
% 		name = 'TGE';
% 	else
% 		name = 'TNG';
%     end
% end
%         
% filename = 'blank.xlsx';
% % Temperature
% writetable(Temp, filename, 'Sheet', name, 'Range', 'B3', 'WriteVariableNames', false);
% % Concentration
% rng = strcat('B', num2str(id1 + 1 * id2));
% writetable(Con, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% % Radiative Forcing
% rng = strcat('B', num2str(id1 + 2 * id2));
% writetable(RF, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% % Abatement
% rng = strcat('B', num2str(id1 + 3 * id2));
% writetable(Act, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% % Geoengineering
% rng = strcat('B', num2str(id1 + 4 * id2));
% writetable(ActG, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% % Carbon Price
% rng = strcat('B', num2str(id1 + 5 * id2));
% writetable(CP, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% % Utility
% rng = strcat('B', num2str(id1 + 6 * id2));
% writetable(Util, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% 
% if scenario == 2
%     % Extreme Event
%     rng = strcat('B', num2str(id1 + 7 * id2));
%     writetable(EE, filename, 'Sheet', name, 'Range', rng, 'WriteVariableNames', false);
% end
