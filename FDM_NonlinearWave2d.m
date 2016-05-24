%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           A FDM solver for the Nonlinear 2D Wave Equation
%                  by Manuel Diaz, NHRI, 2016.05.22
%
%           p_tt - c^2*(p_xx + p_rr) + a*p_ttt + b*p_tt = 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Hallaj, Ibrahim M., and Robin O. Cleveland. "FDTD simulation of
%     finite-amplitude pressure and temperature fields for biomedical
%     ultrasound." The J.Acoust.Soc.Am. 105.5 (1999): L7-L12.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% Parameters
     c = 1.0;   % Advection Coef.
    nx = 200;	% number of space gridpoints
    ny = 200;	% number of space gridpoints
   cfl = 0.5;	% CFL stability condition
tFinal = 8.0;	% final time
plotfigs = true;

% Build mesh
ax=-1; bx=1; x=linspace(ax,bx,nx); dx=x(2)-x(1);
ay=-1; by=1; y=linspace(ay,by,ny); dy=y(2)-y(1);

% Initial time step
dt0=cfl*min(dx/c,dy/c);

% Solution Arrays
up=zeros(nx,ny); % previous state of wave
uo=zeros(nx,ny); % present state of wave
un=zeros(nx,ny); % next state of wave

% plot region
region = [ax,bx,ay,by,-1.2,1.2];