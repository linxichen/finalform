% quarterly
bbeta = 0.99; % if you don't know beta... good luck
ttau = 0.1; % search cost
aalpha = 0.25; % y = a*k^aalpha*l^v
v = 0.5; % labor share
aalpha0 = 2; % search elasticity
ddelta = .1/4;
pphi = 0.000000; % price of disinvestment relative to investment
MC = 1; % How many units of consumption goods is needed for 1 inv good
rrhox = 0.95; % persistence of idio TFP
ppsi = 0.00; % quadratic cost of investment adjustment
rrhoz = rrhox; % persistence of agg TFP
ssigmaz = 0.01; % std of z innov
ssigmax_low = 0.04; % low std of x innov
ssigmax_high= 0.04*3; % high std of x innov
Pssigmax = [1-.05 .05; .08 1-.08]; % Transition prob of ssigmax

% Notational difference
nu = v;
ppsi_n = 1;
xi_n = ppsi_n;
C_ss = 1;

