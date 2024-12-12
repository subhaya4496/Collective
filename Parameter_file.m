function [N, nu, D_R, D_T, dt, nT, sig, sizedif, r_cut, r_cutE, r_list, pre_F, pre_T, delta, Nsave] = Parameter_file(alpha)

%  Parameters

N = 100;                    % Number of cells in simulation
nu = 0.1;                   % Poisson's Ratio
D_T= 1;                   % Sqrt of Translational diffusion constant
D_R= 1;                   % Sqrt of Rotation diffusion constant
dt = 0.0001;                  % Time-step duration
nT = 1000000;                % Number of time steps
sig = 1.0;                  % Size of cell
sizedif = 0.0;              % Difference in size for two cell types
delta = 0.0;                % Buffer length                                              
Nsave = 1000;                % Number of steps after which a figure is saved

% rsigA = 0.5*(sig - sizedif);      % Size of cell type A
% rsigB = 0.5*(sig + sizedif);      % Size of cell type B
% v0A = v0;                   % Velocity of cell type A
% v0B = v0;                   % Velocity of cell type B

% EHSDAB = (rsigA + rsigB);
% EHSDAA = (rsigA + rsigA);
% EHSDBB = (rsigB + rsigB);

MAXHSD = sig + sizedif;     % This is equal to EHSDBB

%% Cut off lengths
r_cut =  1.0*MAXHSD;        % Maximum cut-off LJ and LS  
r_cutE = 7.0*MAXHSD;        % Cut off - Elasticity                                              
r_list = 7.5*MAXHSD;        % Cut off - List

%% Prefactors
pre_F = alpha*(1.+nu)*(3./(16.*pi));     % force prefactor term
pre_T = -alpha*(1.+nu)/(8.*pi);          % torque prefactor term

