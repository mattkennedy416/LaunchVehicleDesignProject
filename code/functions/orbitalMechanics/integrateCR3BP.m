function stateSolution = integrateCR3BP(X, tf)
%INTEGRATECR3BP Summary of this function goes here
%   Detailed explanation goes here


sunMass = 1.989E30;
earthMass = 5.972E24; 
moonMass = 7.34767309E22;

lunarDist = 384402000; % distance from moon to earth
au2m = 1.496e+11;

mu = moonMass / (earthMass + moonMass);


options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10); % set the integration options
[t, stateSolution] = ode45(@stateDerivative, [0, tf], X, options, mu); % and perform the integration

% lets just add t as a 7th column
stateSolution(:,7) = t;

end





%% needed functions
function Xdot = stateDerivative(dt, X, mu)
% function for the ode45 integrator to use
% returns the derivative of our state variable X
%global mu

p = X(1:3);
v = X(4:6);

[dudx,dudy,dudz] = GradU(p(1),p(2),p(3), mu);

dv = [2*v(2) + dudx;
    -2*v(1) + dudy;
    dudz];

dp = X(4:6);

Xdot = [ dp; dv];

end


function [dudx,dudy,dudz] = GradU(x,y,z,mu)

% solved via mathematica

r13 = sqrt( (mu+x).^2 + y.^2 + z.^2 );
r23 = sqrt( (mu+x-1).^2 + y.^2 + z.^2 );

dudx = x - mu*(mu+x-1) / r23^3 - (1-mu).*(mu+x) / r13^3;
dudy = y - y*mu / r23^3 - y*(1-mu) / r13^3;
dudz = -z*mu / r23^3 - z*(1-mu) / r13^3;

% elliptical restricted three body problem equations of motion:
% per: http://www.cds.caltech.edu/~marsden/wiki/uploads/projects/surf/Gawlik_Surf07.pdf
% where f(t) is the true anomaly and e is the eccentricity
%   - this parallels the equation of orbit r=P/(1+e*cos(f))
% global e f
% ER3BP_factor = 1/(1 + e*cos(f));
% dudx = dudx * ER3BP_factor;
% dudy = dudy * ER3BP_factor;
% dudz = dudz * ER3BP_factor;

end


