%Main code

clear all;
clc;
close all;

%% Initilize
% Parameters
tic
D_o = 0.2; %Size of square cylinder
vel_o = 0.05; %Velocity (u)
nu = 10^-4.; %Kinematic Viscosity 
Lx_o = 5.; %Length in x
Ly_o = 1; %Length in y
xc_o = 0.7; %Distance of square cylinder from x
yc_o = 0.5; %Distance of square cylinder from y
dx_o = 0.01; %Increment in x
dy_o = 0.01; %Increment in y
N = (Lx_o/dx_o)+1; %Number of grid points in x
M = (Ly_o/dy_o)+1; %Number of grid points in y
tmax_o = 100.; %Maximum time needed
Re = (vel_o*D_o)/nu; %Reynolds number based on smallest dimension?
cfl = 0.3; %CFL condition
sig = 0.6; %Sigma condition

% Non-Dimensionalize parameters
vel_scale = vel_o; %new velocity scale
l_scale = D_o; %new length scale

D = D_o/l_scale;   % Dimensionless diameter of square
Lx = Lx_o/l_scale;   % Dimensionless length of domain 
Ly = Ly_o/l_scale;   % Dimensionless height of domain
xc = xc_o/l_scale;   % Dimensionless location of square (x)
yc = yc_o/l_scale;   % Dimensionless location of square (x)
dx = dx_o/l_scale; % Dimensionless delta x
dy = dy_o/l_scale; % Dimensionless delta y
vel = vel_o/vel_scale;   % Dimensionless u velocity

tmax = tmax_o * (vel_scale/l_scale); %Dimensionless maximum time [0,5,10,20,40,100];

%% Time plot variables

T_o = [5,10,20,40,100]; %Time plots needed
T = T_o *(vel_scale/l_scale); %Dimensionless time plots
ti = 2; %Plots start at 2, plot 1 taken up by t=0

%% Setup of variables and initilize

u = zeros(N,M); %initilize u
v = zeros(N,M); %initilize u
F1C = zeros(N,M); %initilize F1C
F2C = zeros(N,M); %initilize F2C
F1V = zeros(N,M); %initilize F1V
F2V = zeros(N,M); %initilize F2V
P = zeros(N,M); %initilize P

u = u+vel; %inital condition, u starts with vel everywhere
PO = 0; %inlet pressure condition
a = [0, -5/9,-153/128]; %a coefficient - not needed in loop

%% Creating VOF

%Square cylinder occurs at
xi = (round((xc -(D/2))/dx)) : round(((xc +(D/2))/dx)); %x values of square cylinder
yj = (round((yc -(D/2))/dy)) : round(((yc +(D/2))/dy)); %y values of square cylinder

VOF = ones(N,M); %Initlize with all ones
for i = xi(1):xi(end)
    for j = yj(1):yj(end)   
VOF (i,j) = 0; %Add zeros to parts that the square cyclinder occurs
    end
end

%% Creating grid of x and y points

[x,y] = meshgrid(0:dx:Lx,0:dy:Ly); %meshgrid of x and y values at each point
starty = 0:dy:Ly; %strealines start from all y points
startx = zeros(1,M); %x starts from 0 for streamline function

%Creating t = 0 plot
figure(1)
u = u .* VOF; %removing velocity terms from inside square cylinder
v = v .* VOF; %removing velocity terms from inside square cylinder
streamline(x,y,u',v',startx,starty); %ploting steamline function
xlabel('Dimensionless Length, x')
ylabel('Dimensionless Height, y')

%% Calculate A

A1 = A(dx,dy,N,M); %Calculate A matrix

%% Loop in time

t=0; %Time starting at 0
nt = 1; %Number of time points start
while t < tmax %start of time loop
dt = min([cfl*min([min(dx./max((abs(u)))),min(dy./max((abs(v))))]),sig*min([dx^2/nu,dy^2/nu])]); %calcuation of dt with CFL and sigma conditions
fprintf('t = %f\n',t) %print total time

b = [(1/3)*dt, (15/16)*dt, (8/15)*dt]; %coefficient b
c = [3/dt, 12/(5*dt), 4/dt]; %coefficient c

G1 = zeros(N,M); %Reset after every timestep
G2 = zeros(N,M); %Reset after every timestep

for k = 1:3 %Loop in coefficient values / Runge-Kutta 3rd order
    
[F1C,F2C,F1V,F2V] = F(F1C,F2C,F1V,F2V,Re,u,v,dx,dy,N,M); %Convective and Viscous flux

    G1 = (a(k)*G1) + F1C + F1V; %Sum of Convective and Viscous flux with amendment
    G2 = (a(k)*G2) + F2C + F2V; %Sum of Convective and Viscous flux with amendment
    W1 = u + (b(k)*G1); %Initial correction
    W2 = v + (b(k)*G2); %Initial correction

%i=1 Inlet Boundary condition (velocity)
    W1(1,:) = 1;
    W2(1,:) = 0;  

%j=1 Bottom Boundary condition (velocity)
    W2(:,1) = 0; 

%j=M Top Boundary condition (velocity)
    W2(:,M) = 0;    
    
    W1 = W1 .* VOF; %removing velocity terms from inside square cylinder
    W2 = W2 .* VOF; %removing velocity terms from inside square cylinder

B1 = B(PO,c,k,W1,W2,dx,dy,N,M); %calculating c(k)*dWi/dXi
B2 = reshape(B1,((N-1)*M),1); %reshape to convert to column matrix

    P1 = A1\B2; %Solve for pressure
    P1 = reshape(P1,N-1,M); %Reshape to add inlet pressure
    P = [zeros(1,M);P1]; %Add inlet pressure

    u = U(W1,c,k,P,dx,dy,N,M,VOF); %new u velocity
    v = V(W2,c,k,P,dx,dy,N,M,VOF); %new v velocity

end

t = t+dt; %total time

%% Plotting each graph and time needed
if t >= T(ti) %Finding plotted time needed against total time

figure(ti) %new figure for every plot
streamline(x,y,u',v',startx,starty); %streamline plot
xlabel('Dimensionless Length, x')
ylabel('Dimensionless Height, y')
pause(1) %Pause to show figures
ti=ti+1; %next figure and time for new plot
end

vel_u_st(nt) = u(240,38); %plotting u velocity at i=240 j=30
vel_v_st(nt) = v(240,38); %plotting v velocity at i=240 j=30
dt_st(nt) = dt; %saving each dt
nt = nt+1; %total time steps
end
toc

%% Strouhal number

DT_st = mean(dt_st(1:end)); %average dt of all time
nt_st = 19/DT_st; %number of timestep vortex shedding is stable at
VEL_u_st = vel_u_st(nt_st:end); %u velocity from when vortex shedding is stable
fs=1/DT_st; %sampling frequency
X=abs(fft(VEL_u_st))/length(VEL_u_st); %fourier transform 

F = 0:(fs/(length(VEL_u_st)-1)):fs; %frequency domain

FF = F(2:(end)); %remove first 0 value from frequency domain
XX = X(2:(end)); %remove first value from fourier transform domain

[Y,N] = max(XX); %finding maximum vortex shedding value and position
St = FF(N); %use vortex shedding position to find frequency and thus Strouhal number

%% Change of u in time

figure(ti) %New figure 1 after the last
plot(linspace(0,(nt)*(DT_st/(vel_scale/l_scale)),nt-1),vel_u_st) %Create timescale within plot to match velocity vector
xlabel('Time [s]')
ylabel('Dimensionless Velocity, u')
