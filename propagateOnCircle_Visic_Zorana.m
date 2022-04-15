function [tv,rv0tof,vv0tof] = propagateOnCircle_Visic_Zorana(rv0,vv0,t0,tf,mu,N)
%-----------------------------------------------------------%
% This function propogates the position and inertial        %
% velocity of a satelite orbiting Earth in a CIRCULAR orbit %
% in Earth-Ceneterd Ineritial (ECI) coordinates.            %
%-----------------------------------------------------------%
% Inputs:                                                   %
% -The initial position, rv0 (ECI)                          %
% -The intiial ineritial vecloity, vv0 (ECI)                %
% -The initial time, t0                                     %
% -The final time, tf                                       %
% -The gravitational parameter, mu                          %
% -The number of time intervals, N                          %
%                                                           %
% Outputs:                                                  %
% -a column vector of length N that contains equally        %
% spaced values of time from t0 to tf as a function of      %
% the time interval N [tv]                                  %
% -a matrix of size N x 3 whose rows contains the           %
% position expressed in ECI for the N time values in tv     %
% [rv0tof]                                                  %
% -a matrix of size N x 3 whose rows contain the inertial   %
% velocity expressed in ECI for the N time values in tv     %
% [vv0tof]                                                  %
%-----------------------------------------------------------%

%-----------------------------------------------------------%
% Creating the column vector tv where each value is spaced  %
% deltat time units apart from t0 to tf                     %
%-----------------------------------------------------------%
tv = linspace(t0,tf,N); %row vector of time values
tv = tv.'; %column vector of time values

%-----------------------------------------------------------%
% Using the initial position and velocity to extract all    %
% orbital elements using the code rv2oe_Visic_Zorana.m      %
% In this case because the orbit is CIRCULAR a = r which is %
% constant where a represents the semi-major axis and r     %
% represents the magnitude of the position vector at any    %
% point in time t. Then using the oe vector to solve for    %
% the transformation matrices to convert our solutions from %
% the basis {w_1,w_2,w_3} into ECI cartesian coordinates.   %
%-----------------------------------------------------------%
oe = rv2oe_Visic_Zorana(rv0,vv0,mu); %orbital elements vector
a = oe(1); %semi-latus rectum = radius
bOmega = oe(3); %longitude of the ascending node
inc = oe(4); %inclination
lOmega = oe(5); %arguement of the periapsis
nu = oe(6); %true anomaly
Tin = [ cos(bOmega) -sin(bOmega) 0; sin(bOmega) cos(bOmega) 0; 0 0 1]; 
Tnq = [ 1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
Tqp = [ cos(lOmega) -sin(lOmega) 0; sin(lOmega) cos(lOmega) 0; 0 0 1];
Tpw = [cos(nu) -sin(nu) 0; sin(nu) cos(nu) 0; 0 0 1];

%-----------------------------------------------------------%
% Using the given values for a, mu, and the tv vector to    %
% determine the angle theta of the spacecraft from the      %
% initial position and the rate of change of the angle      %
% theta. These values will be tabulated for each time value %
% in the vector tv and used to solve for rv0tof and vv0tof  %
%-----------------------------------------------------------%
thetadot = sqrt(mu/(a^3));
rv0tof = zeros(N,3);
vv0tof = zeros(N,3);
rv0tof(1,:) = rv0';
vv0tof(1,:) = vv0';
for ii = 2:N
 theta = (tv(ii) - t0)*thetadot;
 pos = [(a*cos(theta)); (a*sin(theta)); 0]; %in the basis {w_1,w_2,w_3}
 vel = [-(a*thetadot*sin(theta)); (a*thetadot*cos(theta)); 0]; %in the basis {w_1,w_2,w_3}
 rvECI = Tin*Tnq*Tqp*Tpw*pos; %in the ECI basis
 vvECI = Tin*Tnq*Tqp*Tpw*vel; %in the ECI basis
 rv0tof(ii,:) = rvECI';
 vv0tof(ii,:) = vvECI';
end
end
