%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           A FDM solver for the Second-order Wave Equation
%                  by Manuel Diaz, NHRI, 2016.05.20
%
%               q_tt + c^2 (q_xx + q_yy) = 0, with c > 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES: Example uses homogeneous Dirichlet b.c. (fixed ends) 
%      or homogeneous Neumann b.c. (loose ends) or any combination of them.
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

%% Solver Loop
t=3*dt0; dt=dt0; it=0;

%for k=1:n
while t < tFinal
    
    for j=2:ny-1
        for i=2:nx-1
            un(i,j)= 2*un(i,j)-up(i,j) + ...
                (c*dt/dx)^2*(uo(i+1,j)-2*uo(i,j)+uo(i-1,j)) + ...
                (c*dt/dy)^2*(uo(i,j+1)-2*uo(i,j)+uo(i,j-1));
        end
    end
    
    % update info
    up=uo; uo=un;

    % Update dt and iteration counter
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt; it=it+1;
    
    % Perturbation in time
    %uo(nx/2,ny/2) = exp(-((t+8)/tFinal)^2);  
    uo(nx/2,ny/2) = sin(0.25*pi*t);
    
    % plot figure
    if plotfigs~=false
        imagesc(x,y,uo); drawnow;
    end
end