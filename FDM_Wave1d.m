%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           A FDM solver for the Second-order Wave Equation
%                  by Manuel Diaz, NHRI, 2016.05.20
%
%                   q_tt + c^2 q_xx = 0, with c > 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES: Example uses homogeneous Dirichlet b.c. (fixed ends) 
%      or homogeneous Neumann b.c. (loose ends) or any combination of them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% Parameters
     c = 1.0;   % Advection Coef.
    nx = 200;	% number of space gridpoints without boundaries
   cfl = 0.9;	% CFL stability condition
tFinal = 8.0;	% final time
plotfigs = true;

% Build mesh
a=-1; b=1; x=linspace(a,b,nx)'; dx=x(2)-x(1); 

% Initial Displacement / IC: u(x,0) = f(x)
%f = @(x) max(1-20*x.^2,0).^2;
%f = @(x) sin(2*pi*x);
f = @(x) exp(-10*x.^2);
%f= @(x) zeros(size(x)) + (x>-0.2 & x<0.2);

% Initial Velocity / IC: u_t(x,0) = g(x)
g = @(x) -[diff(f(x));0]/dx/2;

% Initial time step
dt0=cfl*dx/c; steps=ceil(tFinal/dt0);

% Build Array U(x,t) = U(space,time)
U = zeros(nx,steps); % we want to save all data

% plot region
region = [a,b,-1.2,1.2];

%% Solver Loop
t=3*dt0; dt=dt0; it=0; u0=f(x); U(:,1)=u0; U(:,2)=u0+dt*g(x);

%while t < tFinal
for n=3:steps
    
    % Compute and update [u_old,u_new]
    for i = 2:nx-1
        U(i,n) = 2*U(i,n-1)-U(i,n-2)+(c*dt/dx)^2*(U(i-1,n-1)-2*U(i,n-1)+U(i+1,n-1));
    end
    
    % Update dt and iteration counter
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt; it=it+1;
    
    % plot figure
    if plotfigs~=false
    plot(x,u0,'b:',x,U(:,n),'r.-'); axis(region);
    title(sprintf('time t=%0.2f',t)); drawnow
    end
end

%% Post processing
figure(2);
[X,T] = meshgrid(x,linspace(0,tFinal,steps));
surf(U,'EdgeColor','none'); %colormap hot;
axis tight; axis equal;

