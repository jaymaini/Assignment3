clear all
clearvars
clearvars -GLOBAL
clc
close all
format shorte

set(0, 'DefaultFigureWindowStyle', 'docked')
global C

%Physics Constants
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s²
C.am = 1.66053892e-27;


%Time interval and simulation time 
dt = 0.01e-12;
nt = 100;
%Number of particles - can be run with multiple
np = 2000;
maxXBound = 200e-9;
maxYBound = 100e-9;
x = rand(np,1)*maxXBound;
y = rand(np,1)*maxYBound;

T = 300;
m = 0.26*C.m_0;
v_th = sqrt(2*C.kb*T/m);

%New stuff for assignment #3

%Hand calculations
Voltage = 0.1;
EField = Voltage/maxXBound
E_force = EField * C.q_0

Voltage = 0.5;
EField = Voltage/maxXBound;
E_force = EField * C.q_0;
E_accel = E_force / m;


%Gives -1 or 1
%2*randi([0 1], np, 1)-1

vx = randn(np, 1)*v_th/sqrt(2);
vy = randn(np, 1)*v_th/sqrt(2);
%vy = (2*randi([0 1], np, 1)-1)*v_th/sqrt(2);
%Other Variables
t = 0;

mtc = 0.3e-12; %Mean time between collisions
%Calculate scattering probability
P_scat = 1 - exp(-dt/mtc);
%Time before scattering
T_scat = zeros(1,np);
Scatter_Time_Array = 0;
    

%Visualize a certain number of the electroncs
visuals = 7;
watcher = randi(np,visuals);
 
for i=2:nt
    %Advance the electron in time updating its position and velocity
    t(i) = t(i-1)+dt;
    %Increment Time before scattering
    T_scat = T_scat + dt;
    
    x(:,i) = x(:,i-1) + vx(:,1)*dt;
    y(:,i) = y(:,i-1) + vy(:,1)*dt;
    
    vx(:,1) = vx(:,1) + (1/2)*E_accel*dt;
    
    Jx(i) = C.q_0 * (1e15) * (1e14) * mean(vx(:,1));
    
    hold on;
    for j=1:visuals
        %Create line segments to plot
        
        plot([x(watcher(j),i-1) x(watcher(j),i)],[y(watcher(j),i-1) y(watcher(j),i)])
        hold on
    end
    xlim([0 maxXBound])
    ylim([0 maxYBound])
   
    %Plot Temperature over simulation time
    magVel = sqrt(vx.^2+vy.^2);
    avgV(i,:) = mean(magVel.^2);
    %Calculate Temperature based off average velocity
    t_thermal = (m*(avgV))./(2*C.kb);
    
    %Scattering --------------------------------------------------
    scatter = P_scat > rand(np,1);
    
    vx(scatter) = randn(length(vx(scatter)), 1)*v_th/sqrt(2);
    vy(scatter) = randn(length(vy(scatter)), 1)*v_th/sqrt(2);
    
    %Now affect the Time before scattering:
    Scatter_Time_Array = [Scatter_Time_Array T_scat(scatter)];
    T_scat(scatter) = 0;
    avgScatter = mean(Scatter_Time_Array);
    avgMFP = avgScatter*v_th;
   
  %Reflections for y axis and transmissions for the x asis
  vy((y(:,i)>maxYBound)) = -vy((y(:,i)>maxYBound));
  vy((y(:,i)<0)) = -vy((y(:,i)<0));
  %Perform vertical shift here
 
  x(x>maxXBound) = x(x>maxXBound)-maxXBound;
  x(x<0) = x(x<0)+maxXBound;
 
    
  pause(0.01);
  
end

figure

subplot(2,2,1)
plot(Jx)
title('Current Density over Time')
xlabel('Time')
ylabel('J (A/nm^2')

%Temperature map

%Temperature Calculation
magVel = sqrt(vx.^2+vy.^2);
temperatures = (m*(magVel.^2))./(2*C.kb);
subplot(2,2,2)
%Create Matrix for mapping temperatures to positions
xv = linspace(min(x(:,i)), max(x(:,i)), 100);
yv = linspace(min(y(:,i)), max(y(:,i)), 50);
[Xgrid,Ygrid] = meshgrid(xv, yv);
Z = griddata(x(:,i),y(:,i),temperatures(:,1),Xgrid,Ygrid);
imagesc(xv,yv,Z),colorbar,title('Temperature Map');
axis([0,200e-9,0,100e-9]); 

%Electron Density map
ptsX = linspace(0, 200e-9, 100);
ptsY = linspace(0, 100e-9, 50);
N = histcounts2(y(:,i), x(:,i), ptsY, ptsX);
subplot(2, 2, 3)
imagesc(ptsY,ptsX,N),colorbar,title('Electron Density Map');






