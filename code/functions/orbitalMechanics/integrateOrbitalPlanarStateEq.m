function [t, X] = integrateOrbitalPlanarStateEq(orbit, startTime, tSpan)
%INTEGRATEORBITALPLANARSTATEEQ FOR EARTH - Integrate the base equations of motion,
%without atmospheric considerations
%   startTime will take the values calculated by keplerian equations at
%   this time to use as initial conditions

[~, ind] = min(abs(orbit.time - startTime));
initialNu = orbit.nu(ind);

X0(1,1) = atan( e*sin(initialNu)/(1+e*cos(initialNu)) );
X0(2,1) = norm(orbit.v(ind,:));
X0(3,1) = norm(orbit.r(ind,:));
X0(4,1) = initialNu;


X0(1,1) = X0(1,1)*180/pi; % for some stupid reason MATLAB gives the wrong answer if we use radians instead of degrees
X0(4,1) = X0(4,1)*180/pi;

options = odeset('RelTol',10e-10,'AbsTol',10e-10);
[t,X] = ode45(@stateDeriv, tSpan, X0, options);

X(:,1) = X(:,1) * pi/180;
X(:,4) = X(:,4) * pi/180;

end




function Xdot = stateDeriv(t, X)

%global beta rEarth A cD m liftToDragRatio muEarth

rEarth = 6371;
muEarth = 398600.4415;

gamma = X(1) * pi / 180;
v = X(2); % km/s
r = X(3); % km
theta = X(4);

rho = atmosphericDensity( (r - rEarth) ); % kg/km^3

g = (9.81/1000) * ( rEarth / r)^2; % km/s^2
D = 0.5*rho*v^2 * A * cD;
L = liftToDragRatio * D;

Xdot = double(zeros(4,1));

Xdot(1,1) = (1/v)*( L/m + (v^2/r - g)*cos(gamma) );
Xdot(2,1) = (v^2/r - g)*sin(gamma) - D/m;
Xdot(3,1) = v*sin(gamma);
Xdot(4,1) = (v/r)*cos(gamma) * 180 / pi;

end


function rho = atmosphericDensity(h)
    
if h > 150
    rho = 3.875E-9 * exp( -h / 59.06 );
else
    rho = 1.226*exp(-h/7.524); % units of kg/m^3
end

rho = rho * 1E9; % units of kg / km^3
end

