
global mu e f

sunMass = 1.989E30;
earthMass = 5.972E24; 
moonMass = 7.34767309E22;

lunarDist = 384402000; % distance from moon to earth
au2m = 1.496e+11;

mu = moonMass / (earthMass + moonMass);

R_SOI = (earthMass/sunMass)^(2/5)*(au2m/lunarDist); % need earth-sun orbit normalized to earth-moon orbit


% mu = 0.012277471; % earth moon
% mu = 1E-6;
% mu = 0.04;
% e = 1; % eccentricity for ER3BP 
% f = 0; % true anomaly for ER3BP


[L1,L2,L3,L4,L5] = findLagrangePoints(mu);

%% initial state specified in homework problem
% p = [1.061692; 0; 0];
% v = [0; 0.403877; 0];

% % semi-stable L2 orbit
% p = [0.825; 0; 0];
% v = [0; 0.103561313; 0];

% % can we try l5?
% L5 = [0.5 - mu;
%     -sqrt(3)/2;
%     0];
% p = L5 + [+0.05; 0; 0];
% [dudx,dudy,dudz] = GradU(p(1),p(2),p(3),mu);
% 
% v = [0.05; 0.1; 0];
% 
%  p = [0.5;0;0];
%  v = [0.0;0.5;0];
%  
%  p = [0.233724738610180; -0.241905389892218; 0];
%  v = -[0.247849857342514; -0.013605124443645; 0];


% test backwards compatibility (ie if we take some solution point and
% reverse the velocity, we should go back to where we came from right?) lol nope
% p = transpose([0.568431071839008,-0.941343232283432,0]);
% v = -transpose([-0.00322362599497506,-0.205405503517633,0]);


% compare 3 body dynamics to the linearized solutions
% p = L5 + [-0.1; -0.1; 0];
% v = [0.05; 0.; 0];


% play around with some manifold tubes
% p = L1 + [-0.0664268271886969; 0.0362220992775682; 0];
% v = [
% p = L1 + [-0.1; 0.02; 0];
% v = -[-0.1276;0.1349; 0];



% p = L1 + [0.002; -0.007172998214837; 0];
% v = [0.004668771548911; -0.016744544992909; 0];

% p = L2 + [-4E-4; -0.010399774027823; 0];
% v = [-0.004830093696154; 5.396474160970063E-4; 0];

% p = [1.0922; 0.3510; 0];
% v = [0.0884; -0.4724; 0];
% % v = [0.0909; -0.4729; 0];
% 
% 
% 

% p = [-1; 0; 0];
% v = [0; 0.1; 0];

nu = 0; % orbital angle of secondary body
A_r2i = [cos(nu), -sin(nu); sin(nu), cos(nu)];

a = 0.1;
n = a^(-1/2);

p_inertial = [0; -a];
v_inertial = [n*1.33; 0];

p = inv(A_r2i) * p_inertial;
v = inv(A_r2i) * v_inertial - [-p(2); p(1)];

p(3) = 0;
v(3) = 0;

% p = [0; 0.99*R_SOI; 0];
% v = [2.95; 0; 0];

% p = [0; 0.1; 0];
% v = [3; 0; 0];




X = [p; v]; % stack together for integration

tf = 10; % final integration time

%% calculate integration
options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10); % set the integration options
[~, stateSolution] = ode45(@stateDerivative, [0, tf], X, options); % and perform the integration

% check the conservation of energy using the Jacobi Integral
x0 = stateSolution(1,1);
y0 = stateSolution(1,2);
z0 = stateSolution(1,3);
c0 = norm(stateSolution(1,4:6))^2 - 2*potential(x0,y0,z0);

xf = stateSolution(end,1);
yf = stateSolution(end,2);
zf = stateSolution(end,3);
cf = norm(stateSolution(end,4:6))^2 - 2*potential(xf,yf,zf);

% assert(abs(cf-c0) < 1E-10)


%% calculate a map of the pseudo-potential
outerVal = 4;
resolution = 0.01;
[pot_x,pot_y] = meshgrid(-outerVal:resolution:outerVal,-outerVal:resolution:outerVal);
pot_z = 0; % just show the z=0 plane

potentialMap = potential(pot_x,pot_y,pot_z);
potentialMap(potentialMap>3) = nan; % at the bodies themselves this goes to infinity;

%% and create the plots
x = stateSolution(:,1);
y = stateSolution(:,2);

figure(1)
hold on
plot(0,0, 'k*')
plot(1,0, 'c*')
contour(pot_x,pot_y,potentialMap,100)
plot(x, y, 'r')
legend('Earth','Moon','Potential Map', 'Orbit Path')
xlabel('Distance from Earth (DU)')
ylabel('Distance from Earth (DU)')


plot(L1(1),L1(2), 'r*')
plot(L2(1),L2(2), 'r*')
plot(L3(1),L3(2), 'r*')
plot(L4(1),L4(2), 'r*')
plot(L5(1),L5(2), 'r*')

rMoon = 1.7371E6 / lunarDist;
rectangle('Position',[1-rMoon, -rMoon, 2*rMoon, 2*rMoon], 'Curvature',[1,1])



%% needed functions
function Xdot = stateDerivative(dt, X)
% function for the ode45 integrator to use
% returns the derivative of our state variable X
global mu

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


function U = potential(x, y, z)

global mu

w=1;
U = 0.5*(x.^2 + y.^2)*w^2 + (1-mu)./sqrt( (x + mu).^2 + y.^2 + z.^2 ) + mu ./ sqrt( (x-(1-mu)).^2 + y.^2 + z.^2);

% elliptical restricted three body problem equations of motion:
% per: http://www.cds.caltech.edu/~marsden/wiki/uploads/projects/surf/Gawlik_Surf07.pdf
% where f(t) is the true anomaly and e is the eccentricity
%   - this parallels the equation of orbit r=P/(1+e*cos(f))
% global e f
% ER3BP_factor = 1/(1 + e*cos(f));
% U = U * ER3BP_factor;

end
