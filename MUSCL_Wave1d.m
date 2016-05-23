%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          A basic MUSCL solver for Hyperbolic system equations
%                   by Manuel Diaz, NHRI, 2016.05.20
%
%                       [u]    [ 0 , c ][u]   
%                       [v]t + [ c , 0 ][v]x = 0,
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

%% Parameters
cfl     = 0.8;	% CFL number
tEnd    = 20.0;	% Final time
nx      = 1000;  % Number of cells/Elements
limiter ='MM';  % MC, MM, VA.
fluxMth ='LF';	% LF.
plot_fig= true; % plot figures

% Coefs
c = 1;

% Discretize spatial domain
a=-50; b=50; dx=(b-a)/nx; xc=a+dx/2:dx:b;

% Set IC
%u0 = zeros(size(xc));
%u0 = sin(2*pi*xc); 
%u0 = zeros(size(xc)) + (xc>-0.2 & xc<0.2);
u0 = exp(-5*xc.^2);
%u0 = max(1-20*xc.^2,0).^2;
%u0 = cos(5*pi*xc)+ 3*(xc>-0.3 & xc<0.3);
%v0 = zeros(size(xc));
v0 = zeros(size(xc));
%v0 = u0;

% Set q-array & adjust grid for ghost cells
nx=nx+2; q0=[u0;v0]; zero=[0;0]; q0=[zero,q0,zero];

% Boundary Conditions in ghost cells
q0(:,1)=q0(:,2); q0(:,nx)=q0(:,nx-1);   % Natural BCs

% Plot region
region = [a,b,-1.2,1.2];

% Initial time step
dt0=cfl*dx/c;

% print to terminal dx and dt
fprintf('using dx: %1.2f and dt: %1.2f\n',dx,dt0);

%% Solver Loop

% Load IC
q=q0; t=0; it=0; dt=dt0; 

tic
while t < tEnd
    
    % RK2 1st step
    qs = q - dt*MUSCL_WaveRes1d(q,c,dx,nx,limiter,fluxMth);
    
    qs(:,1)=qs(:,2); 
    qs(:,nx)=qs(:,nx-1);   % Natural BCs
    
    % RK2 2nd step  / update q
    q = (q + qs - dt*MUSCL_WaveRes1d(qs,c,dx,nx,limiter,fluxMth))/2;
    
    q(:,1)=q(:,2); 
    q(:,nx)=q(:,nx-1);   % Natural BCs
        
    % compute conserved properties
    u=q(1,:); v=q(2,:);
    
    % Update dt and time
    dt=cfl*dx/c; if t+dt>tEnd; dt=tEnd-t; end; t=t+dt; it=it+1;
    
    % Plot figure
    if plot_fig ~= false
        if rem(it,10) == 0
            subplot(2,1,1); plot(xc,u(2:nx-1),'.r'); axis(region);
            subplot(2,1,2); plot(xc,v(2:nx-1),'.r'); axis(region);
            drawnow
        end
    end
end
cputime = toc; fprintf('CPUtime: %1.2f\n',cputime);

%% Post Processing

% Remove ghost cells
q=q(:,2:nx-1); nx=nx-2; 

% compute flow properties
u=q(1,:); v=q(2,:);

% Plots results
figure(1);
subplot(2,1,1); plot(xc,u,'.r',xc,u0,':b'); xlabel('x'); ylabel('u'); 
axis(region); legend('MUSCL','IC'); title('SSP-RK2 MUSCL for Wave System Eqns.')
subplot(2,1,2); plot(xc,v,'.r',xc,v0,':b'); xlabel('x'); ylabel('v'); 
axis(region); legend('MUSCL','IC');