%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         A FDM solver for the Second-order Westervelt Equation
%                 by Manuel A. Diaz, NHRI, 2016.05.22
%
%    p_xx-1/(c^2)*p_tt+delta/(c0^4)*p_ttt+beta/(r0*c0^4)*p^2_tt = 0,  (1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Hallaj, Ibrahim M., and Robin O. Cleveland. "FDTD simulation of
%     finite-amplitude pressure and temperature fields for biomedical
%     ultrasound." The J.Acoust.Soc.Am. 105.5 (1999): L7-L12.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %clc; close all;

%% Parameters
     cfl = 0.80;	% CFL stability condition
  tFinal = 32E-6;	% final time
plotfigs = true;	% plot figures?

% Thermofluid Physical parameters
dx=0.1E-3; % m
c0=1600; % m/s
r0=1100; % kg/m^3 
alpha=4.5; % Np/m, absortion
beta=5.5; % 1+B/(2A), where B/A nonlinarity parameter of fluid
f=1E6; % 1 MHz 
omega=2*pi*f; % angular freq
delta=2*c0^2*alpha/omega^2; % sound viscosity

% Bioheat model parameters
kt=0.6; % W/(m.K)
Ct=3800; % J/(kg.K)
Cb=3800; % J/(kg.K)
wb=0.5; % kg/(m^3 s)

% Build mesh
ax=-2.56E-2; bx=2.56E-2; x=ax+dx/2:dx:bx; nx=length(x);

% Initial time step
dt0=cfl*dx/c0;

% Solution Arrays
pmmmmm=zeros(1,nx); % previous state of wave
pmmmm=zeros(1,nx); % previous state of wave
pmmm=zeros(1,nx); % previous state of wave
pmm=zeros(1,nx); % previous state of wave
pm=zeros(1,nx); % previous state of wave
po=zeros(1,nx); % present state of wave
pn=zeros(1,nx); % next state of wave

% plot region
region = [ax,bx,-2.2,2.2];

%% Solver Loop
t=6*dt0; dt=dt0; it=0; 

% coefs
a=(c0*dt0/dx)^2/12; b=(delta/(8*dt0*c0^2)); c=beta/(r0*c0^2);

method = 2;
while t < tFinal
    
    switch method
        case 1
            for i=3:nx-2
                pn(i)= 2*pn(i)-pm(i) + ...
                    a*(-po(i+2)+16*po(i+1)-30*po(i)+16*po(i-1)-po(i-2)) + ...
                    b*(6*po(i)-23*pm(i)+34*pmm(i)-24*pmmm(i)+8*pmmmm(i)-pmmmmm(i)) + ...
                    c*(2*po(i)^2-5*pm(i)^2+4*pmm(i)^2-pmmm(i)^2);                
            end
        case 2
            for i=3:nx-2
                pn(i)= 2*pn(i)-pm(i) + ...
                    a*(-po(i+2)+16*po(i+1)-30*po(i)+16*po(i-1)-po(i-2)) + ...
                    b*(6*po(i)-23*pm(i)+34*pmm(i)-24*pmmm(i)+8*pmmmm(i)-pmmmmm(i)) + ...
                    beta*dt^2/(12*r0*dx^2)*(-po(i+2)^2+16*po(i+1)^2-30*po(i)^2+16*po(i-1)^2-po(i-2)^2);
            end
    end
    
    % update info
    pmmmmm=pmmmm; pmmmm=pmmm; pmmm=pmm; pmm=pm; pm=po; po=pn;
    
    % Update dt and iteration counter
    if t+dt>tFinal; dt=tFinal-t; a=(c0*dt/dx)^2/12; b=(delta/(8*dt*c0^2)); end; t=t+dt; it=it+1; 
    
    % Perturbation in time
    po(nx/2) = sin(omega*t)*exp(-((t-4E-6)/2E-6)^2); 
    
    % plot figure
    if plotfigs~=false
    plot(x,po,'r.-'); axis(region); grid minor; 
    title(sprintf('time t=%0.2f',t)); drawnow
    end
end
