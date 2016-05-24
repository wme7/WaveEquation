%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           A WENO solver for a System of Hyperbolic PDE's
%                  by Manuel Diaz, NHRI, 2016.05.22
%
%   A Numerical solver based on weno scheme for solving the second order
%   wave equation q_tt + c^2 q_xx = 0, with c > 0, as a 1st order system
%
%                    [u]     [ 0 , c ][u]   
%                    [v]_t + [ c , 0 ][v]_x = 0,
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

%% Parameters
	 nx = 4000;	% Number of cells/Elements
    cfl = 0.90;	% CFL number
   tEnd = 50.0;	% Final time
plot_fig= true; % plot figures

% Discretize spatial domain
a=-50; b=50; dx=(b-a)/nx; xc=a+dx/2:dx:b;

% Coefs
c=1;

% Set IC
%u0 = zeros(size(xc));
%u0 = sin(2*pi*xc); 
%u0 = zeros(size(xc)) + (xc>-0.5 & xc<0.5);
u0 = exp(-10*xc.^2);
%u0 = max(1-20*xc.^2,0).^2;
%u0 = cos(5*pi*xc)+ 3*(xc>-0.3 & xc<0.3);
%v0 = zeros(size(xc));
v0 = zeros(size(xc));
%v0 = u0;

% Set q-array & adjust grid for ghost cells
nx=nx+4; q0=[u0;v0]; zero=[0;0]; q0=[zero,zero,q0,zero,zero];

% Boundary Conditions in ghost cells
q0(:,1 )=q0(:, 3  ); q0(:, 2  )=q0(:, 3  ); % Natural BCs
q0(:,nx)=q0(:,nx-2); q0(:,nx-1)=q0(:,nx-2); 

% Plot region
region = [a,b,-1.2,1.2];

% Initial time step
dt0=cfl*dx/c;

%% Solver Loop

% Load IC
q=q0; t=dt0; it=0; dt=dt0; 

% print to terminal dx and dt
fprintf('using dx: %1.2f and dt: %1.2f\n',dx,dt0);

tic
while t < tEnd
    % RK Initial step
    qo = q;
    
    % 1st stage
    dF=WENO5_WaveRes1d(q,c,nx,dx,'LF');     q = qo-dt*dF;
    q(:,1)=qo(:,3); q(:, nx )=qo(:,nx-2); % Neumann BCs
    q(:,2)=qo(:,3); q(:,nx-1)=qo(:,nx-2); % Neumann BCs
    
    % 2nd Stage
    dF=WENO5_WaveRes1d(q,c,nx,dx,'LF');     q = 0.75*qo+0.25*(q-dt*dF);
    q(:,1)=qo(:,3); q(:, nx )=qo(:,nx-2); % Neumann BCs
    q(:,2)=qo(:,3); q(:,nx-1)=qo(:,nx-2); % Neumann BCs
    
    % 3rd stage
    dF=WENO5_WaveRes1d(q,c,nx,dx,'LF');     q = (qo+2*(q-dt*dF))/3;
    q(:,1)=qo(:,3); q(:, nx )=qo(:,nx-2); % Neumann BCs
    q(:,2)=qo(:,3); q(:,nx-1)=qo(:,nx-2); % Neumann BCs
        
    % compute conserved properties
    u=q(1,:); v=q(2,:);
    
    % Update dt and time
    dt=cfl*dx/c; if t+dt>tEnd; dt=tEnd-t; end; t=t+dt; it=it+1;
    
    % Plot figure
	if plot_fig ~= false
        if rem(it,10) == 0
            subplot(2,1,1); plot(xc,u(3:nx-2),'.-r'); axis(region); grid on; grid minor;
            subplot(2,1,2); plot(xc,v(3:nx-2),'.-r'); axis(region); grid on; grid minor;
            drawnow
        end
    end
end
cputime = toc; disp(['CPUtime: ',num2str(cputime)])

%% Post-process

% Remove ghost cells
q=q(:,3:nx-2); nx=nx-4; 

% compute flow properties
u=q(1,:); v=q(2,:);

% Plots results
figure(1);
subplot(2,1,1); plot(xc,u,'.-r',xc,u0,':b'); xlabel('x'); ylabel('u'); 
axis(region); legend('MUSCL','IC','location','southeast'); grid on; grid minor;
title('SSP-RK3 WENO5 for Wave System Eqns.')
subplot(2,1,2); plot(xc,v,'.-r',xc,v0,':b'); xlabel('x'); ylabel('v'); 
axis(region); legend('MUSCL','IC','location','southeast'); grid on; grid minor;