%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            A basic MUSCL solver for Hyperbolic system equations
%                      by Manuel Diaz, NTU, 29.04.2015
%
%                       [p]    [ 0 , r c^2 ][p]   
%                       [u]t + [ 1 / r , 0 ][u]x = 0,
%
%	MUSCL based numerical schemes extend the idea of using a linear
%	piecewise approximation to each cell by using slope limited left and
%	right extrapolated states. This results in the following high
%	resolution, TVD discretisation scheme.   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Refs:
%   [1] LeVeque, Randall J., and Randall J. Le Veque. Numerical methods for
%   conservation laws. Vol. 132. Basel: Birkhäuser, 1992. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %clc; close all;
global c r

%% Parameters
     dx = 0.02; % cell size
    cfl = 0.50;	% CFL number
   tEnd = 10.0; % Final time
limiter ='MC';  % MC, MM, VA.
fluxMth ='LF';	% LF.
plot_fig= true; % plot figures

% Discretize spatial domain
a=-1; b=11; xc=a+dx/2:dx:b; nx=length(xc);

% Constants
c=1; r=1;

% Set IC
%u0 = zeros(size(xc));
%u0 = sin(2*pi*xc); 
u0 = zeros(size(xc)) + (xc>-0.5 & xc<0.5);
%u0 = exp(-5*xc.^2);
%u0 = max(1-20*xc.^2,0).^2;
%u0 = cos(5*pi*xc)+ 3*(xc>-0.3 & xc<0.3);
p0 = zeros(size(xc));
%p0 = u0;

% Set q-array & adjust grid for ghost cells
nx=nx+2; q0=[p0;u0]; zero=[0;0]; q0=[zero,q0,zero];

% Boundary Conditions in ghost cells
q0(:,1)=q0(:,2); q0(:,nx)=q0(:,nx-1);   % Natural BCs

% Plot region
region1 = [a,b,-1.2,1.2];
region2 = [a,b,-1.2,1.2];

% Initial time step
dt0=cfl*dx/c;

% print to terminal dx and dt
fprintf('using dx: %1.2e and dt: %1.2e\n',dx,dt0);

%% Solver Loop

% Load IC
q=q0; t=dt0; it=0; dt=dt0;

tic
while t < tEnd
    
    % RK2 1st step
    qs = q - dt*MUSCL_NonlinearWaveRes1d(q,c,dx,nx,limiter,fluxMth,t);
  
    % RK2 2nd step  / update q
    q = (q + qs - dt*MUSCL_NonlinearWaveRes1d(qs,c,dx,nx,limiter,fluxMth,t+dt/2))/2;
        
    % compute conserved properties
    p=q(1,:); u=q(2,:);
    
    % Update dt and time
    dt=cfl*dx/c; if t+dt>tEnd; dt=tEnd-t; end; t=t+dt; it=it+1;
    
    % Plot figure
    if plot_fig ~= false
        if rem(it,10) == 0
            subplot(2,1,1); plot(xc,p(2:nx-1),'.r'); axis(region1); grid minor;
            subplot(2,1,2); plot(xc,u(2:nx-1),'.r'); axis(region2); grid minor;
            drawnow
        end
    end
end
cputime = toc; fprintf('CPUtime: %1.2f\n',cputime);

%% Post Processing

% Remove ghost cells
q=q(:,2:nx-1); nx=nx-2; 

% compute flow properties
p=q(1,:); u=q(2,:);

% Plots results
figure(1);
subplot(2,1,1); plot(xc,p,'.r',xc,p0,':b'); xlabel('x'); ylabel('p'); 
axis(region1); legend('MUSCL','IC','location','southwest'); grid on; grid minor;
title('SSP-RK2 MUSCL for Wave System Eqns.')
subplot(2,1,2); plot(xc,u,'.r',xc,u0,':b'); xlabel('x'); ylabel('u'); 
axis(region2); legend('MUSCL','IC','location','southwest'); grid on; grid minor;