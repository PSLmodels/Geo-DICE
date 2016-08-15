global T gama sai1 sai2 sai3 nuG deltak k0 pb pr pd theta1 theta2 theta3 bprice
global Sig0 Sigg0 Siggd Sigga Sigg Sig E0 Eland Mat0 Mlo0 Mup0
global b11 b12 b21 b22 b23 b32 b33 F2000 F2100 Tat0 Tlo0 
global etha1 etha3 etha4 deltarf Fex L0 Lg0 Lg LA L Ag Ag0 Agd A A0 pai ro
global R alpha K0 lambda a0 g0 Gcoeff Geffective aa0 ag0 thetaGE
global sai2temp sai2ocean sai2atmos

% OUTPUT AND CAPITAL ACCUMULATION
%   T = Total Time horizon
T = 60;
time = 1:T;

% OUTPUT AND CAPITAL ACCUMULATION

%   capital share
gama = 0.3;
%   damage coefficient on temperature
sai1 = 0.0;
%   damage coefficient on temperature squared
sai2 = 0.0028388;
%   Exponent on damages
sai3 = 2.0;
%   Geoengineering damages multiplier
nuG = 0.03;
%   rate of depreciation (percent per year)
deltak = 0.1;
%   Initial capital stock ($ trillion)
k0 = 137;
%   cost of backstop 2005
pb = 1.17;


backrat = 2.0;
gback = 0.05;
limmiu = 1.0;
bprice = zeros(1,T);
for  i=1:T
bprice(i) = pb*((backrat-1+exp(-gback*(i-1)))/backrat);
end

%   ratio initial to final backstop cost
pr = 2;
%   initial cost decline backstop
pd = 0.05;
%   Exponent of control cost function
theta2 = 2.8;
%   Exponent of geoengineering cost function
theta3 = 2;
%   Coefficient on GE cost function - fixed
thetaGE = 0.0027;
%   Coefficient of geoengineering cost function
Gcoeff = 1;
%   Effectiveness of geoengineering
Geffective = 1;

%%Calibration of new damage coefficients
%   This is DICE's damage specification in initial period, based just on
%   temperature
%   Initial temp = 0.7307 degrees above preindustrial
origdamage = sai1 * .7307 + sai2 * .7307 ^ sai3;
%   In new specification, 80% of this should be directly from temperature
sai2temp = 0.80 * origdamage / 0.7307 ^ 2;
%   10% should be directly from upper ocean concentration
%   Initial concentration = 1255; preindustrial = 1094; difference = 161
sai2ocean = 0.10 * origdamage / 161 ^ 2;
%   10% should be directly from atmospheric concentration
%   Initial concentration = 808.9; preindustrial = 596.4; difference =
%   212.5
sai2atmos = 0.10 * origdamage / 212.5 ^ 2;



% EMISSIONS

%    Initial Sig
Sig0 = 0.1342;
%    Initial growth rate of Sig (percent per decade)
Sigg0 = -7.3;
%    Rate of decrease in the growth rate of Sig (percent per year)
Siggd = 0.3;
%    Accleration parameter of growth rate of Sig
Sigga = 0.0;
%    Growth rate of Sig (percent per decade)
Sigg = (Sigg0 / 100) * exp(-10 * (Siggd / 100) * (time - 1) - 10 * Sigga * (time - 1).^2);
Sigg(1) = Sigg0;
%    Sig (industrial CO2 emissions/output -- MTC/$1000)
Sig = zeros(1, T);
Sig(1) = Sig0;
for i = 2:T
    Sig(i) = Sig(i - 1) / (1 - Sigg(i));
end
%    Initial carbon emissions from land use change (GTC per year)
E0 = 1.1;
%    Carbon emissions from land use change (GTC per year)
Eland = exp(log(E0) + log(0.9) * (time - 1));
%   Abatement cost function coefficient
theta1 = (pb * Sig / theta2).*((pr - 1 + exp(-pd * (time - 1))) / pr);

% CONCENTRATIONS

%    Initial atmospheric concentration of CO2 (GTC)
Mat0 = 808.9;
%    Initial concentration of CO2 in biosphere/shallow oceans (GTC)
Mup0 = 1255;
%    Initial concentration of CO2 in deep oceans (GTC)
Mlo0 = 18365;
%    Carbon cycle transition coefficients (percent per decade)
%       atmosphere to atmosphere (b11)
b11 = 81.0712;
%       biosphere/shallow oceans to atmosphere (b21)
b21 = 9.7213;
%        atmosphere to biosphere/shallow oceans (b12)
b12 = 18.9288;
%       biosphere/shallow oceans to biosphere/shallow oceans (b22)
b22 = 85.2787;
%        deep oceans to biosphere/shallow oceans (b32)
b32 = 0.3119;
%       biosphere/shallow oceans to deep oceans (b23)
b23 = 5;
%       deep oceans to deep oceans (b33)
b33 = 99.6881;

% TEMPERATURE

%   2000 forcings other ghg
F2000 = -0.06;
%   2100 forcings other ghg
F2100 = 0.3;
%   Initial atmospheric temperature (deg. C above 1900)
Tat0 = 0.7307;
%   Initial temperature of deep oceans (deg. C above 1900)
Tlo0 = 0.0068;
%   Speed of adjustment parameter for atmospheric temperature
etha1 = 0.22;
%   FCO22x
deltarf = 3.8;
%   Coefficient of heat loss from atmosphere to oceans
etha3 = 0.3;
%   Coefficient of heat gain by deep oceans
etha4 = 0.05;
%   Exogenous forcing (Watts per square meter)
Fex = F2000 + 0.1 *(F2100 - F2000) * (time - 1);
Fex(12:T) = F2100;

% POPULATION

%   Initial population (millions)
L0 = 6514;
%   Init rate of pop growth (per decade)
Lg0 = 0.35;
%   Asymptotic population
LA = 8600;
%   Growth factor population
Lg = ((exp(Lg0 * (time - 1)) - 1))./exp(Lg0 * (time - 1));
%   Population (millions)
L = L0 * (1 - Lg) + LA * Lg;

% PRODUCTIVITY

%   Initial level of total factor productivity
A0 = 0.02722;
%    Init rate of prod growth (per decade)
Ag0 = 0.092;
%   Rate of decline in productivity growth rate (percent per year)
Agd = 0.001;
%   Rate of growth of productivity (percent per decade)
Ag = Ag0 * exp(-Agd * 10 * (time - 1));
%   Total factor productivity
A = zeros(1,T);
A(1) = A0;
for i = 2:T
    A(i) = A(i - 1)/(1 - Ag(i - 1));
end

% PARTICIPATION

%   PARTFRACT (participation rate)
pai = ones(1, T);

% WELFARE

%   Social rate of time preference (percent per year)
ro = 1.5 * ones(1,T);
%   Social time preference factor
R = zeros(1,T);
R(1) = 1;
for i = 2:T
    R(i) = R(i - 1)/(1 + ro(i - 1) / 100)^10;
end
alpha = 2;
lambda = 1/(1+ro(1)/100)^10;

% Initialization

%   Initial Capital ($trill)
K0 = 137;

% Initial Actions
a0 = 0.005;
g0 = 0.005;

% Approximate Actions
aa0 = 0;
ag0 = 0;