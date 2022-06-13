
% Load empirical parameters
P = load('./data/CA-counties-P.csv');   % Matrix of county visitation probs
s0 = load('./data/CA-counties-s0.csv'); % Vector of county populations
scale = 2.3654e-07;                     % Ensure that R0 = 2.5
A = scale * (P * P');                   % Approx. matrix of inter-county contact rates
nNodes = size(A, 1);

% Allocation settings
budget = 0.1;
delOffset = 2;

% Simulation settings  
tf = 10000;
nPts = 20000; 
tt = linspace(0, tf, nPts);

% Infection initial settings
i0 = s0 * 0.001; 
x0 = [s0; zeros(nNodes, 1); i0; zeros(nNodes, 1)];
