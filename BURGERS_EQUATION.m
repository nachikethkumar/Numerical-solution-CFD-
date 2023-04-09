%% Burger PDE Equation for solving %%
%Ct=-u*Cx+kCxx   Advection and diffusion equation ( Burger Equation)
%Along with Initail condition
% Ct=time derivative
% Cx= spatial derivative
% u= Velocity

clear all; clc;
%% Domain Creation

% Space (Spatial Condition)
% Define total length of the domain (Lx) in meters
Lx=6;
% Divide the total length of domain into small parts (dx)
dx=0.3;
% Calculate the number of nodes
nx=fix(Lx/dx); % fix command outputs interger
% store the individual node points in a variable x
x=linspace(0,Lx,nx);

% Time (Temporal condition)
% setup the time in domain

% define the total time to be run
T= 0.06%0.52 %0.06 ; %2 0.06seconds
% the interval of time that should be used in Calculation is given by CFL
% CFL Number
CFL=0.08 %0.1 0.005;     % ( CFL=v*Dt/Dx)
v=1; % Velocity is 1 (given)
% define the time interval now
dt=0.002 %CFL*dx/v;
% Calculate the total number of time steps
nt=fix(T/dt);

%% Field Arrays
%Variables
u=zeros(nx,1);

% parameters ( here parameters are speed and Viscosity)
Vel=zeros(nx,1);% Initialize speed term
k=zeros(nx,1); % Initialize viscosity term

%% Initail conditions 
 t=0;
 Vel(:)=1
 Phi=zeros(nx,1)
 dPhi=zeros(nx,1)
 k(:)=0.02;
 for i=1:length(x)
 Phi(i,1)=(exp(-(x(i)^2)/(4*k(i))))+(exp(-((x(i)-(2*pi))^2)/(4*k(i))))
 dPhi(i,1)=((-0.5*(1/k(i))*x(i))*exp(-(x(i)^2)/(4*k(i))))+((-0.5*((x(i)-(2*pi))/k(i)))*exp((-(x(i)-(2*pi))^2)/(4*k(i))))
 u(i,1)=(-(2*k(i)*dPhi(i))/Phi(i))+4
 end
plot(x,u,'-k',LineWidth=2)
hold on

%% Analytical Solution 

uany=u
for n=1:nt
    t=(n-1)*dt
   
    for i = 2:nx-1
      
  phi_any(i,1) = exp((-(x(i)-4*t).^2)/(4*k(i))/(t+1))...
                             + exp(-(x(i)-4*t-2*pi).^2/(4*k(i))/(t+1));
    
  dphi_any(i,1) = (-(x(i)-4*t)*0.5/k(i)/(t+1))*...
                     exp(-(x(i)-4*t).^2/(4*k(i))/(t+1))...
                  +(-(x(i)-4*t-2*pi)/k(i)/(t+1))*0.5*...
                     exp(-(x(i)-4*t-2*pi).^2/(4*k(i))/(t+1));
            
    uany(i,1) = (-(2*k(i))*(dphi_any(i,1)/phi_any(i,1))) + 4;
    uany(end)=uany(end-1)
    uany(1)=u(1)
    end
end
 
%% Time stepping loop
 % for starting the simulation this part is the one that performs the
 % analysis at every time step 
 % thus loop is give for every time step
 t=0
Cn=u
%C=zeros(nx,1)
C=u
for n=1:nt

 % Boundary conditions
 % for periodic boundary ( where the convection occurs periodically)
 Cn(end)=Cn(end-1)
Cn(1)=Cn(end)
C(1)=Cn(end)
C(end)=Cn(end-1)
 
    % Solution
    t=t+dt;
    %Cn=u
    %u(i)=C(i)
    %Cn(i)=u(i)
    for i= 2:nx-1
        % Advection term
         A=Vel(i)*((Cn(i)-Cn(i-1))/dx);
         % Diffusion term
         %D=k(i)*(Cn(i+1)-2*Cn(i)+Cn(i-1))/dx^2
         D=k(i)*(Cn(i+1)-2*Cn(i)+Cn(i-1))/dx^2
        % Eulers method
        %C(i)=Cn(i)+dt*((-A)+(D));

         C(i)=Cn(i)-(Cn(i)*dt*(Cn(i)-Cn(i-1))/dx)+(k(i)*dt*(Cn(i+1)-2*Cn(i)+Cn(i-1))/(dx*dx));
         %C(i)=Cn(i)-(Vel(i)*dt*(Cn(i)-Cn(i-1))/dx)+(k(i)*dt*(Cn(i+1)-2*Cn(i)+Cn(i-1))/(dx*dx));

         Cn(i)=C(i)
    end

    

    % Visulaize at selected steps
    %clf;
  
    plot(x,C,'-o');
    grid on
   %refreshdata(h);
    %hold on
    %title( sprintf('t= %.2f',t)); % .2f denotes the 2 decimal places
    %hold on
    %plot(x,u,'-c')
    title({['Conserved parameter vs grid point (Viscocity co-efficient = ',num2str(k(1)),')'];[sprintf('t= %.2f',t)]})
    xlabel('Grid points')
    ylabel('Conserved parameter')
    legend('Initial Condition','Numerical Solution')
     axis([0 Lx 0.5 7.5]);
    shg; pause(0.1);
    
    %drawnow;
    % shg makes the current figure visible and places it in front of all
    % other figures on the screen.
  
end

figure(2)
plot(x,u,'-k',LineWidth=2)
hold on
plot(x,uany,'g--',LineWidth=2)
hold on
plot(x,C,'-o');
grid on
title({['Conserved parameter vs grid point (Viscocity co-efficient = ',num2str(k(1)),')'];[sprintf('t= %.2f',t)]})
xlabel('Grid points')
ylabel('Conserved parameter')
legend('Initial Condition','Analytical Solution','Numerical Solution')
axis([0 Lx 0.5 7.5]);



