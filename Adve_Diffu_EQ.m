% PDE Equation for solving
%Ct=-u*Cx+kCxx +qs   Advection equation with source term
% Ct=time derivative
% Cx= spatial derivative
% u= Velocity
% qs= Source term

%% Domain 
% Space
% Define total length of the domain (Lx) in meters
Lx=10;
% Divide the total length of domain into small parts (dx)
dx=0.5;
% Calculate the number of nodes
nx=fix(Lx/dx); % fix command outputs interger
% store the individual node points in a variable x
x=linspace(0,Lx,nx);

% Time
% setup the time in domain

% define the total time to be run
T=5; % seconds
% the interval of time that should be used in Calculation is given by CFL
% CFL Number
CFL=0.05;     % ( CFL=v*Dt/Dx)
v=1; % Velocity is 1 (given)
% define the time interval now
dt=CFL*dx/v;
% Calculate the total number of time steps
nt=fix(T/dt);

%% Field Arrays
%Variables
C=zeros(nx,1);

% parameters ( here parameters are speed)
u=zeros(nx,1);
k=zeros(nx,1);

%% Initail conditions 
 t=0;
 C(:)=0;
 u(:)=1;
 k(:)=2;

%% Time stepping loop
 % for starting the simulation this part is the one that performs the
 % analysis at every time step 
 % thus loop is give for every time step

for n=1:nt

    % Boundary conditions
 % for periodic boundary ( where the convection occurs periodically)
 C(end)=C(end-1)
 C(1)=C(end)

    % Source 
     if n==1
         C(1)=1; % You can set it up for whatever the source value is 
    % here the source value is 1
    % else
         %C(1)=0
     end

    

    % Solution
    t=t+dt;
    Cn=C
    for i= 2:nx-1
        % Advection term
         A=u(i)*((Cn(i)-Cn(i-1))/dx);
         % Diffusion term
         D=k(i)*(Cn(i+1)-2*Cn(i)+Cn(i-1))/dx^2
        % Eulers method
        C(i)=Cn(i)+dt*((-A)+(D));
    end
    % Check convergence

    % Visulaize at selected steps
    %clf;
    plot(x,C,'-o')
    title( sprintf('t= %.2f',t)); % .2f denotes the 2 decimal places
    axis([0 Lx 0 0.2]);
    shg; pause(0.01);

    % shg makes the current figure visible and places it in front of all other figures on the screen

end





