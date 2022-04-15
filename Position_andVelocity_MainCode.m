clc; clear;
%---------------------------------------------------------------%
% This script runs the function propagateOnCircle_Visic_Zorana  %
% to find a table of size N x 7 where each column 1 represents  %
% a vector with N components of equally spaced time values from %
% t0 to tf, columns 2 through 4 represent the (x,y,z) position  %
% of the satelite at each time value, and columns 5 through 7   %
% represent the (v_x,v_y,v_z) inertial velocity of the satelite %
% at eacch time value. Both position and inertial velocity are  %
% represented in the Earth-Centered Inertial (ECI) Cartesian    %
% coordinates.                                                  %
%                                                               %
% Once this is done, the code will run the function ode113 with %
% the same time grid as the propagateOnCircle function. This    %
% function will have a relative error tolerance (epsi). Then,   %
% the two functions will be plotted on the same graph in order  %
% to compare the two.                                           %
%---------------------------------------------------------------%
%---------------------------------------------------------------%
% Given Conditions %
%---------------------------------------------------------------%
rv0 = [-5613.97603835865; -2446.44383433555; 2600.48533877841]; %initial position (ECI) [km]
vv0 = [2.12764777374332; -7.13421216656605; -2.1184067703542]; %initial inerital velocity (ECI) [km * s^-1]
t0 = 33.2; %initial time [min]
tf = 100.2; %final time [min]
mu = 398600; %gravitational parameter of Earth [km^3 * s^-2]
N = 45; %number of time intervals
epsi = (10)^-8; %error tolerance of ode113
opts = odeset('RelTol',epsi);
azimuth = 45; %first plot view azimuth
elevation = 45; %first plot view elevation
x = 1; %the position of the x-dimension in each matrix row
y = 2; %the position of the y-dimension in each matrix row
z = 3; %the position of the z-dimension in each matrix row

%---------------------------------------------------------------%
% Transforming time to seconds, and creating a time span between%
% the two that is evenly spaced and length N                    %
%---------------------------------------------------------------%
t0 = t0*60; %inital time [sec]
tf = tf*60; %final time [sec]
tspan = linspace(t0,tf,N); %time interval for ode113

%---------------------------------------------------------------%
% Creating a vector p that contains the initial posiiton and    %
% inertial velocity of the satelite                             %
%---------------------------------------------------------------%
p0(x:z) = rv0;
p0(x+3:z+3) = vv0;

%---------------------------------------------------------------%
% Calling the function propagateOnCircle_Visic_Zorana using     %
% given and derived information.                                %
%---------------------------------------------------------------%
[tv,rv1,vv1] = propagateOnCircle_Visic_Zorana(rv0,vv0,t0,tf,mu,N);

%---------------------------------------------------------------%
% Calling the function ode113 using given and derived           %
% information.                                                  %
%---------------------------------------------------------------%
[t, p] = ode113(@(t,p) twoBodyOde_Visic_Zorana(t,p,mu), tspan, p0, opts);
rv2 = p(:,1:3);
vv2 = p(:,4:6);

%---------------------------------------------------------------%
% Creating a table of the propagateOnCircle_Visic_Zorana values %
%---------------------------------------------------------------%
T = table(tv,rv1(:,x),rv1(:,y),rv1(:,z),vv1(:,x),vv1(:,y),vv1(:,z),'VariableNames',{'Time [sec]', 'x', 'y', 'z', 'v_x', 'v_y', 'v_z'});
T(:,:)

%---------------------------------------------------------------%
% Plotting the two sets of position data obtained from the two  %
% methods used above. Set the view to (azimuth, elevation).     %
%---------------------------------------------------------------%
figure(1)
plot3(rv1(:,x),rv1(:,y), rv1(:,z),'-o',rv2(:,x),rv2(:,y),rv2(:,z),'-s')
view(azimuth,elevation);
xl = xlabel('$x$ (km)');
yl = ylabel('$y$ (km)');
zl = zlabel('$z$ (km)');
ll = legend('Propagate on Circle','ODE113','Location','none');
tit = title('Propogation and Integration of a Given Orbit over Time');
set(xl,'FontSize',16,'Interpreter','LaTeX');
set(yl,'FontSize',16,'Interpreter','LaTeX');
set(zl,'FontSize',16,'Interpreter','LaTeX');
set(ll,'FontSize',16,'Interpreter','LaTeX');
set(tit,'FontSize',16,'Interpreter','LaTeX');
set(gca,'FontSize',16,'TickLabelInterpreter','LaTeX');
grid on;

%---------------------------------------------------------------%
% Plotting the interial velocity data obtained from the         %
% propagation method with the corresponding velocity vector at  %
% each value of tv.                                             %
%---------------------------------------------------------------%
quiverScale = 2;
figure(2)
pp = quiver3(rv1(:,x),rv1(:,y),rv1(:,z),vv1(:,x),vv1(:,y),vv1(:,z),quiverScale,'-o');
view(azimuth,elevation);
xl = xlabel('$x$ (km)');
yl = ylabel('$y$ (km)');
zl = zlabel('$z$ (km)');
tit = title('Position and Velocity of the Given Orbit over Time');
set(xl,'FontSize',16,'Interpreter','LaTeX');
set(yl,'FontSize',16,'Interpreter','LaTeX');
set(zl,'FontSize',16,'Interpreter','LaTeX');
set(tit,'FontSize',16,'Interpreter','LaTeX');
set(gca,'FontSize',16,'TickLabelInterpreter','LaTeX');
grid on;