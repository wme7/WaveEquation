%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           A WENO solver for a System of Hyperbolic PDE's
%                  by Manuel Diaz, NHRI, 2016.05.20
%
%   A Numerical solver based on weno scheme for solving the second order
%   wave equation q_tt + c^2 q_xx = 0, with c > 0, as a 1st order system
%
%                    [p]     [ 0 , r c^2 ][p]   
%                    [u]_t + [ 1 / r , 0 ][u]_x = 0,
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] LeVeque, Randall J., and Randall J. Le Veque. Numerical methods for
% conservation laws. Vol. 132. Basel: Birkhäuser, 1992. 
% [2] Shu, Chi-Wang. "High order weighted essentially nonoscillatory
% schemes for convection dominated problems." SIAM review 51.1(2009):82-126. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %clc; close all;
global c r vis eps s t0

%% Parameters
     dx = 0.02; % cell size
    cfl = 0.90;	% CFL number
   tEnd = 10.0; % Final time
plot_fig= true;

% Discretize spatial domain
a=-1; b=11; xc=a+dx/2:dx:b; nx=length(xc);

% Coefs
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
nx=nx+4; q0=[p0;u0]; zero=[0;0]; q0=[zero,zero,q0,zero,zero];

% Boundary Conditions in ghost cells
q0(:,1 )=q0(:, 3  ); q0(:, 2  )=q0(:, 3  ); % Natural BCs
q0(:,nx)=q0(:,nx-2); q0(:,nx-1)=q0(:,nx-2); 

% Plot region
region1 = [a,b,-1.2,1.2];
region2 = [a,b,-1.2,1.2];

% Initial time step
dt0=cfl*dx/c;

% print to terminal dx and dt
fprintf('using dx: %1.2f and dt: %1.2f\n',dx,dt0);

%% Solver Loop

% Load IC
q=q0; t=dt0; it=0; dt=dt0;

tic
while t < tEnd
    % RK Initial step
    qo = q;
    
    % 1st stage
    dF=WENO5_NonlinearWaveRes1d(q,c,nx,dx,'LF',t); q = qo-dt*dF;
    
    % 2nd Stage
    dF=WENO5_NonlinearWaveRes1d(q,c,nx,dx,'LF',t+dt/2); q = 0.75*qo+0.25*(q-dt*dF);
    
    % 3rd stage
    dF=WENO5_NonlinearWaveRes1d(q,c,nx,dx,'LF',t+dt); q = (qo+2*(q-dt*dF))/3;
        
    % compute conserved properties
    p=q(1,:); u=q(2,:);
    
    % Update dt and time
    dt=cfl*dx/c; if t+dt>tEnd; dt=tEnd-t; end; t=t+dt; it=it+1;
    
    % Plot figure
    if plot_fig ~= false;
        if rem(it,10) == 0
            subplot(2,1,1); plot(xc,p(3:nx-2),'.-r'); axis(region1); grid minor;
            subplot(2,1,2); plot(xc,u(3:nx-2),'.-r'); axis(region2); grid minor;
            drawnow
        end
    end
end
cputime = toc; fprintf('CPUtime: %1.2f\n',cputime);

%% Post-process

% Remove ghost cells
q=q(:,3:nx-2); nx=nx-4; 

% compute flow properties
p=q(1,:); u=q(2,:);

% Plots results
figure(1);
subplot(2,1,1); plot(xc,p,'.-r',xc,p0,':b'); xlabel('x'); ylabel('p'); 
axis(region1); legend('WENO5','IC','location','southwest'); grid on; grid minor;
title('SSP-RK3 WENO5 for Wave System Eqns.')
subplot(2,1,2); plot(xc,u,'.-r',xc,u0,':b'); xlabel('x'); ylabel('u'); 
axis(region2); legend('WENO5','IC','location','southwest'); grid on; grid minor;