%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         A FDM solver for the Second-order Westervelt Equation
%                 by Manuel A. Diaz, NHRI, 2016.05.20
%
% (q_xx+q_rr)-1/(c^2)*p_tt+delta/(c0^4)*p_ttt+beta/(r0*c0^4)*p^2_tt=0, (1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Hallaj, Ibrahim M., and Robin O. Cleveland. "FDTD simulation of
%     finite-amplitude pressure and temperature fields for biomedical
%     ultrasound." The J.Acoust.Soc.Am. 105.5 (1999): L7-L12.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: 
%  The numerical solution of Eq.(1) was calculated on a polar cylindrical
%  grid. The source was modeled as a spherical cap (bowl) having azimuthal
%  symmetry about the axis of the source. The acoustic and temperature
%  fields need only be computed using a two dimensional spatial grid, x, in
%  the axial direction and r in the radial direction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %clc; close all;

%% Parameters
     cfl = 0.50;	% CFL stability condition
  tFinal = 32E-6;	% final time
plotfigs = true;	% plot figures?

% Thermofluid Physical parameters
dx=0.1E-3; % m
dr=0.1E-3; % m
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
ar=-1.00E-2; br=4.12E-2; r=ar+dr/2:dr:br; nr=length(r);

% Build device
[X,R]=meshgrid(x,r);
faceA = (X.^2+R.^2-3.00E-2^2) > 0;
faceB = (X.^2+R.^2-3.01E-2^2) < 0;
faceC = abs(X)<0.02;
device= find(faceA.*faceB.*faceC);
clear X R faceA faceB faceC;

% Initial time step
dt0=cfl*min(dx/c0,dr/c0);

% Solution Arrays
pmmmmm=zeros(nr,nx); % previous state of wave
pmmmm=zeros(nr,nx); % previous state of wave
pmmm=zeros(nr,nx); % previous state of wave
pmm=zeros(nr,nx); % previous state of wave
pm=zeros(nr,nx); % previous state of wave
po=zeros(nr,nx); % present state of wave
pn=zeros(nr,nx); % next state of wave

% plot region
region = [ax,bx,ar,br,-1.2,1.2];

%% Solver Loop
t=3*dt0; dt=dt0; it=0; 

% coefs
a=(c0*dt0/dx)^2/12; b=(c0*dt0/dr)^2/12; c=(delta/(8*dt0*c0^2)); d=beta/(r0*c0^2);

while t < tFinal
        
    for j=3:nr-2
        for i=3:nx-2
            pn(i,j)= 2*pn(i,j)-pm(i,j) + ...
                a*(-po(i,j+2)+16*po(i,j+1)-30*po(i,j)+16*po(i,j-1)-po(i,j-2)) + ...
                b*(-po(i+2,j)+16*po(i+1,j)-30*po(i,j)+16*po(i-1,j)-po(i-2,j)) + ...
                c*(6*po(i,j)-23*pm(i,j)+34*pmm(i,j)-24*pmmm(i,j)+8*pmmmm(i,j)-pmmmmm(i,j)) + ...
                d*(2*po(i,j)^2-5*pm(i,j)^2+4*pmm(i,j)^2-pmmm(i,j)^2); 
        end
    end   
    
    % update info
    pmmmmm=pmmmm; pmmmm=pmmm; pmmm=pmm; pmm=pm; pm=po; po=pn;

    % Update dt and iteration counter
    if t+dt>tFinal; dt=tFinal-t; % update coefs 
    a=(c0*dt/dx)^2/12; b=(c0*dt/dr)^2/12; c=(delta/(8*dt*c0^2)); end; 
    t=t+dt; it=it+1;
    
    % Perturbation in time
    po(device) = sin(1E6*2*pi*t)*exp(-((t-4E-6)/1.5E-6)^2); 
    
    % plot figure
    if plotfigs~=false
        imagesc(x,r,po,[-1,1]); ylabel('r'); xlabel('x'); axis('square'); 
        title(sprintf('time t=%1.2e',t)); drawnow
    end
end