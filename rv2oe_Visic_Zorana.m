function oe = rv2oe_Visic_Zorana(rv,vv,mu)
%-------------------------------------------------------------------------%
%Inputs:                                                                 
% rv: Cartesian planet-centered inertial (PCI) position (3 by 1)        
% vv: Cartesian planet-centered inertial (PCI) velocity (3 by 1)        
% mu: gravitational parameter of centrally acting body                  
%Outputs:                                                               
% oe(1): semi-major axis (a)                                            
% oe(2): eccentricity (e)                                               
% oe(3): longitude of the ascending node (big omega) [rad]              
% oe(4): inclination (inc) [rad]                                        
% oe(5): arguement of the periapsis (little omega) [rad]                
% oe(6): true anomaly (nu) [rad]                                        
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Part (1): Preparation of orbital elements calculations                
%                                                                       
% Calculating the angular momentum vector (hv), angular momentum        
% magnitude (h), semi-latus rectum (p), magnitude of the position vector 
% (r), and the eccentricity vector (ev). Also creating                  
% a matrix to represent the planet-centered inertial (PCI) basis
% where the matrix (Ib) represents the basis {Ix, Iy, Ix}. Also, finding 
% the vector (nv) which is used to calculate the longitude of the ascending
% node (bOmega), the orbital inclination (inc), and the argument of the
% periapsis (lOmega). This section is used to neatly set up the rest of 
% the orbital elements calculations in Part (2).
%-------------------------------------------------------------------------%
hv = cross(rv,vv); %angular momentum vector
h = norm(hv); %magnitude of the angular momentum
p = (h^2)/mu; %semi-latus rectum
r = norm(rv); %magnitude of the position
ev = ((cross(vv,hv))/mu) - (rv/r); %eccentricity vector
Ix = [1; 0; 0]; %vector which represents the PCI basis z-direction
Iy = [0; 1; 0]; %vector which represents the PCI basis z-direction
Iz = [0; 0; 1]; %vector which represents the PCI basis z-direction
nv = cross(Iz,hv); %calculating the line of nodes n
n = norm(nv); %magnitude of the line of nodes n

%-------------------------------------------------------------------------%
% Part (2): Calculating the orbital elements
% 
% Using the values derived in Part (1) to solve for all of the orbital
% elements oe = [a; e; bOmega; inc; lOmega; nu].
%-------------------------------------------------------------------------%
e = norm(ev); %eccentricity
a = p/(1-(e^2)); %semi-major axis
bOmega = atan2((dot(nv,Iy)),(dot(nv,Ix))); %longitude of the ascending node [rad]
if bOmega < 0 %correcting if bOmega is negative
 bOmega = bOmega + (2*pi);
end
inc = atan2((dot(-hv,(cross(Iz,nv)))),(n*(dot(hv,Iz)))); %inclination [rad]
lOmega = atan2((dot(ev,(cross(hv,nv)))),(h*(dot(ev,nv)))); %argument of the periapsis [rad]
if lOmega < 0 %correction if lOmega is negative
 lOmega = lOmega + (2*pi);
end
nu = atan2((dot(rv,cross(hv,ev))),(h*(dot(rv,ev)))); %true anomaly [rad]
if nu < 0 %correcting if nu is negative
 nu = nu + (2*pi);
end

%-------------------------------------------------------------------------%
% Part (3): Construction the orbital elements vector to be output from the
% function
%-------------------------------------------------------------------------%
oe = [a; e; bOmega; inc; lOmega; nu]; %orbital elements vector
end
