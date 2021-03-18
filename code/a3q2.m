% clear all
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
nt = 500;
%Number of particles - can be run with multiple
np = 10000;
maxXBound = 200e-9;
maxYBound = 100e-9;

%Create the bottleneck - top
bbox_x = [80e-9 120e-9];
bbox_y = [60e-9 100e-9];
%Bottom (uses same X coordinates)
bbox_y2 = [0 40e-9];

%Diffusive option - enable to switch from specular barriers to diffusive
diffusive = 0;

%Initialize positions 
x = rand(np,1)*maxXBound;
y = rand(np,1)*maxYBound;

%For all particles in the region, re-initialize position 
viols = (x(:,1) > bbox_x(1)) & (x(:,1) < bbox_x(2)) & (y(:,1) > bbox_y(1)) & (y(:,1) < bbox_y(2));
while(sum(viols) > 0)
    x(viols) = rand(sum(viols),1)*maxXBound;
    y(viols) = rand(sum(viols),1)*maxYBound;
    
    viols = (x(:,1) > bbox_x(1)) & (x(:,1) < bbox_x(2)) & (y(:,1) > bbox_y(1)) & (y(:,1) < bbox_y(2));
end

%Do the same for the second box
viols = (x(:,1) > bbox_x(1)) & (x(:,1) < bbox_x(2)) & (y(:,1) > bbox_y2(1)) & (y(:,1) < bbox_y2(2));
while(sum(viols) > 0)
    x(viols) = rand(sum(viols),1)*maxXBound;
    y(viols) = rand(sum(viols),1)*maxYBound;
    
    viols = (x(:,1) > bbox_x(1)) & (x(:,1) < bbox_x(2)) & (y(:,1) > bbox_y2(1)) & (y(:,1) < bbox_y2(2));
end

T = 300;
m = 0.26*C.m_0;
Voltage = 0.8;

%New stuff for assignment #3
%Obtain E fields and find acceleration in each direction
[Ex Ey] = EFieldFinder(Voltage,0.33);
E_force_x = (1e9)*Ex' * C.q_0;
E_force_y = (1e9)*Ey' * C.q_0;
E_accel_x = E_force_x ./ m;
E_accel_y = E_force_y ./ m;

%Velociies
v_th = sqrt(2*C.kb*T/m);
vx = randn(np, 1)*v_th/sqrt(2);
vy = randn(np, 1)*v_th/sqrt(2);

t = 0;
mtc = 0.2e-12;               %Mean time between collisions
P_scat = 1 - exp(-dt/mtc); % scattering probability
T_scat = zeros(1,np);       %Time before scattering
Scatter_Time_Array = 0;
    
%Visualize a certain number of the electroncs
visuals = 7;
watcher = randi(np,visuals);
figure 
for i=2:nt
    %Advance the electron in time updating its position and velocity
    t(i) = t(i-1)+dt;
    T_scat = T_scat + dt; %Increment Time before scattering
    x(:,i) = x(:,i-1) + vx(:,1)*dt;
    y(:,i) = y(:,i-1) + vy(:,1)*dt;
    
    [x_bin, edge_x] = discretize(x(:,i),200);
    [y_bin, edge_y] = discretize(y(:,i),100);
    
    %APPLY ACCELERATION DUE TO E FIELD
    vx(:,1) = vx(:,1) + (1/2)*E_accel_x(sub2ind(size(E_accel_x),x_bin,y_bin))*dt;
    vy(:,1) = vy(:,1) + (1/2)*E_accel_y(sub2ind(size(E_accel_y),x_bin,y_bin))*dt;
    
    hold on;
    for j=1:visuals
        %Create line segments to plot      
        plot([x(watcher(j),i-1) x(watcher(j),i)],[y(watcher(j),i-1) y(watcher(j),i)])
        hold on
    end
    xlim([0 maxXBound])
    ylim([0 maxYBound])

    %Scattering --------------------------------------------------
    scatter = P_scat > rand(np,1);
    
    vx(scatter) = randn(length(vx(scatter)), 1)*v_th/sqrt(2);
    vy(scatter) = randn(length(vy(scatter)), 1)*v_th/sqrt(2);
    %After scattering, average temperature fluctuates but remains averaged
    %around 300K.
    
    %Now affect the Time before scattering:
    Scatter_Time_Array = [Scatter_Time_Array T_scat(scatter)];
    T_scat(scatter) = 0;
    avgScatter = mean(Scatter_Time_Array);
        
  %Reflections for y axis and transmissions for the x asis
  vy((y(:,i)>maxYBound)) = -vy((y(:,i)>maxYBound));
  vy((y(:,i)<0)) = -vy((y(:,i)<0));

  x(x>maxXBound) = x(x>maxXBound)-maxXBound;
  x(x<0) = x(x<0)+maxXBound;
    
  
  %Draw the bottleneck
  if i == 2
      x1 = bbox_x(1);
      x2 = bbox_x(2);
      y1 = bbox_y(1);
      y2 = bbox_y(2);
      boxX = [x1, x1, x2, x2, x1];
        boxY = [y1, y2, y2, y1, y1];
        fill(boxX, boxY, 'w');
        
      x1 = bbox_x(1);
      x2 = bbox_x(2);
      y1 = bbox_y2(1);
      y2 = bbox_y2(2);
      boxX = [x1, x1, x2, x2, x1];
        boxY = [y1, y2, y2, y1, y1];
        fill(boxX, boxY, 'w');  
  end
  
  %Reflections for bottlenecks - shift the particle to avoid getting them
  %stuck!
   
   %Top left edge
   viol_y = (y(:,i) > bbox_y(1));
   viol_x = ((x(:, i-1) < bbox_x(1)) & (x(:, i-1) > 0) & (x(:,i) > bbox_x(1)) & (x(:,i) < maxXBound));
   if (diffusive)
       vx(viol_y & viol_x) = randn(sum(viol_y & viol_x), 1)*v_th/sqrt(2);
   else
       vx(viol_y & viol_x) = -vx(viol_x & viol_y);
       %Shift
       x((viol_y & viol_x), i) = (x((viol_y & viol_x), i-1));
   end
   
   %Top right edge 
   viol_x = ((x(:, i-1) > bbox_x(2)) & (x(:,i) < bbox_x(2)) & (x(:,i) > 0.5*maxXBound) & (x(:, i-1) < maxXBound));
   if (diffusive)
       vx(viol_y & viol_x) = randn(sum(viol_y & viol_x), 1)*v_th/sqrt(2);
   else
     vx(viol_y & viol_x) = -vx(viol_x & viol_y);
     x((viol_y & viol_x), i) = (x((viol_y & viol_x), i-1));
   end
   
   %Horizontal bottleneck - top edge
   viol_y = y(:,i-1) < bbox_y(1) & (y(:,i) > bbox_y(1));
   viol_x = (x(:, i) > bbox_x(1)) & (x(:,i) < bbox_x(2));
   if (diffusive)
       vy(viol_y & viol_x) = randn(sum(viol_y & viol_x), 1)*v_th/sqrt(2);
   else
     vy(viol_y & viol_x) = -vy(viol_x & viol_y);
     y((viol_y & viol_x), i) = (y((viol_y & viol_x), i-1));
   end
   
   %Bottom left edge
   viol_y = (y(:,i) < bbox_y2(2));
   viol_x = ((x(:, i-1) < bbox_x(1)) & (x(:, i-1) > 0) & (x(:,i) > bbox_x(1)) & (x(:,i) < maxXBound));
   if (diffusive)
       vx(viol_y & viol_x) = randn(sum(viol_y & viol_x), 1)*v_th/sqrt(2);
   else
     vx(viol_y & viol_x) = -vx(viol_x & viol_y);
     x((viol_y & viol_x), i) = (x((viol_y & viol_x), i-1));
   end
   
   %Bottom right edge 
   viol_x = ((x(:, i-1) > bbox_x(2)) & (x(:,i) < bbox_x(2)) & (x(:,i) > 0.5*maxXBound)  & (x(:, i-1) < maxXBound));
   if (diffusive)
       vx(viol_y & viol_x) = randn(sum(viol_y & viol_x), 1)*v_th/sqrt(2);
   else
     vx(viol_y & viol_x) = -vx(viol_x & viol_y);
      x((viol_y & viol_x), i) = (x((viol_y & viol_x), i-1));
   end
   
   %Horizontal bottleneck - bottom edge
   viol_y = y(:,i-1) > bbox_y2(2) & (y(:,i) < bbox_y2(2));
   viol_x = (x(:, i) > bbox_x(1)) & (x(:,i) < bbox_x(2));
   if (diffusive)
       vy(viol_y & viol_x) = randn(sum(viol_y & viol_x), 1)*v_th/sqrt(2);
   else
      vy(viol_y & viol_x) = -vy(viol_x & viol_y);
      y((viol_y & viol_x), i) = (y((viol_y & viol_x), i-1));
   end
   
  pause(0.01);
  
end
    figure
    %Electron Density map
    ptsX = linspace(0, 200e-9, 100);
    ptsY = linspace(0, 100e-9, 50);
    N = histcounts2(y(:,i), x(:,i), ptsY, ptsX);
    imagesc(ptsY,ptsX,N),colorbar,title('Electron Density Map');







