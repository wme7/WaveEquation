%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           A Leapfrog solver for the Second-order Wave Equation
%                  by Manuel Diaz, NHRI, 2016.05.20
%
%                   q_tt + c^2 q_xx = 0, with c > 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES: Example uses homogeneous Dirichlet b.c. (fixed ends) 
%      or homogeneous Neumann b.c. (loose ends) or any combination of them. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by Benjamin Seibold 02/2007 -- http://www-math.mit.edu/~seibold/
% modifs by Manuel Diaz 05/2016 -- http://www-github.com/wme7/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feel free to modify for teaching and learning.
clc; clear; close all;

%% Parameters
     c = 1.0;   % Advection Coef.
    nx = 200;	% number of space gridpoints without boundaries
   cfl = 0.9;	% CFL stability condition
tFinal = 8.0;	% final time

% Build mesh
a=-1; b=1; x=linspace(a,b,nx+2)'; x=x(2:end-1); dx=x(2)-x(1); 

% Initial Displacement / IC: u(x,0) = f(x)
%f = @(x) max(1-20*x.^2,0).^2;
%f = @(x) sin(2*pi*x);
f = @(x) exp(-10*x.^2);
%f= @(x) zeros(size(x)) + (x>-0.2 & x<0.2);

% Initial Velocity / IC: u_t(x,0) = g(x)
g = @(x) -[diff(f(x));0]/dx/2;

% Initial time step
dt0=cfl*dx/c; 

% Build Coefs Matrix
I=eye(nx); R=diag(ones(1,nx-1),1); r=dt0/dx; A=r^2*R+2*(1-r^2)*I+r^2*R';

% Impose Neumann BCs
%A( 1 , 1 ) = 2-r^2;
%A(end,end) = 2-r^2;

% plot region
region = [a,b,-1.2,1.2];

%% Solver Loop
t=3*dt0; dt=dt0; it=0; u0=f(x); U=[u0,u0+dt*g(x)];

while t < tFinal

   % Compute and update [u_old,u_new]
   U = [U(:,2),A*U(:,2)-U(:,1)];
   
   % Update dt and iteration counter
   if t+dt>tFinal; dt=tFinal-t; r=dt/dx; A=r^2*R+2*(1-r^2)*I+r^2*R'; A(end,end) = 2-r^2; end; 
   t=t+dt; it=it+1;
   
   % plot figure
   plot(x,u0,'b:',x,U(:,2),'r.-'); axis(region); 
   title(sprintf('time t=%0.2f',t)); drawnow
end
