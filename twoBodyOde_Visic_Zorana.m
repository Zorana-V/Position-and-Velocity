function pdot = twoBodyOde_Visic_Zorana(t,p,mu)
%-%-------------------------------------------------%-%
%-% Computes the right-hand side dynamics           %-%
%-% for the two-body differential equation.         %-%
%-%-------------------------------------------------%-%
%-% Inputs:                                         %-%
% t: current time                                   %-%
% p: 6 by 1 column matrix containing the            %-%
% state of the spacecraft at time t                 %-%
% p(1:3): Planet-centered position                  %-%
% p(4:6): Planet-centered inertial position         %-%
%-% Output:                                         %-%
%-% pdot: 6 by 1 column matrix containing the       %-% 
%-% rate of change of the state at time t           %-% 
%-%-------------------------------------------------%-%
%-% ASSUMPTION: CANONICAL UNITS ARE USED (MU=1)     %-%
%-%-------------------------------------------------%-%
%-%-------------------------------------------------%-%
%-% Gravitational Parameter = mu = 1                %-%
%-%-------------------------------------------------%-%
%-%-------------------------------------------------%-%
%-% Extracting the position vector from the p       %-%
%-% vector in order to determine the magnitude of   %-%
%-% the radius vector                               %-%
%-%-------------------------------------------------%-%
rv = p(1:3);
vv = p(4:6);
r = norm(rv);

%-%-------------------------------------------------%-%
%-% Creating a new vector pdot which contains the   %-%
%-% time derivates for p                            %-%
%-%-------------------------------------------------%-%
rvdot = vv;
vvdot = -(mu/(r^3))*rv;
pdot = [rvdot; vvdot];
end
