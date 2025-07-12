function S = get_WC_Dynamics(G, D, stimulated_oscillator, c5, CS)
% =====================================================================
%
% This function initializes the Wilson-Cowan dynamics, runs it, then
%  calculates the pair-wise cognitive system synchrony matrix, S.
%
%  Parameters
%  ----------
%                     G: The streamline density, matrix of size 128x128.
%                     D: The time-delay matrix, matrix of size 128x128.
% stimulated_oscillator: The brain region being stimulated, integer.
%                    c5: The "critical" point, the value just before the
%                         active state.
%                    CS: Cognitive system assignments, container map.
%
%
% Returns
% -------
%         S: The pair-wise cognitive system synchrony matrix, matrix 9x9.
% ======================================================================

%% Initialize Wilson-Cowan

% Time parameters.
simulation_time = 1500; %Total simulation time, ms
dt = 1e-3; %Time step, ms

% Stimulation
stimulation_strength = 1.15; %Strength of stimulation
stimulation_start_time = 300; %Start of the stimulation, ms
Stim_P = [stimulated_oscillator,...
          stimulation_strength,...
          stimulation_start_time,...
          stimulation_time]; % Stimulation on the excitatory population.
Stim_Q = []; % No stimulation on the inhibitory population.

% Coupling parameters
c6 = c5/4; % Roughly 1 inhibitory neuron for every 4 excitatory neurons.

%% Run Wilson-Cowan model.
Dynamics = wc_coupled_stochastic(G, D, time, dt, c5, c6, Stim_P, Stim_Q);

% Grab time and excitatory/inhibitory activity.
t = Dynamics.t(1000/dt:end);
E = Dynamics.e(1000/dt:end,:);
I = Dynamics.i(1000/dt:end,:);

% Cognitive system container map keys (change if necessary). 
key = {'Att', 'Aud', 'CO', 'FP', 'DM', 'MS', 'SC', 'VT', 'V'};

% Get sizes
N = size(E,2); % Number of brain regions (n=128).
M = length(key); % Number of cognitive systems (n=9).

%% get phase of the limit cycle.
phi = atan2(I,E); % Geometric phase.

%% kuramoto order parameter & Pair-wise cognitive system synchrony matrix.

S = zeros(M,M);
for i = 1:M
    for j = 1:M
        S(i,j) = mean(abs(mean(exp(1j.*phi(:,[CS(key{i}),CS(key{j})])),'all')));
    end
end