
global mu e f

sunMass = 1.989E30;
earthMass = 5.972E24; 
moonMass = 7.34767309E22;

lunarDist = 384402000; % distance from moon to earth
au2m = 1.496e+11;

mu = moonMass / (earthMass + moonMass);

R_SOI = (earthMass/sunMass)^(2/5)*(au2m/lunarDist); % need earth-sun orbit normalized to earth-moon orbit

rAtmosphere = 1 / (lunarDist/1000/(6371 + 100));

[L1,L2,L3,L4,L5] = findLagrangePoints(mu);



nu = 0; % orbital angle of secondary body
A_r2i = [cos(nu), -sin(nu); sin(nu), cos(nu)];

a = rAtmosphere; % starting altitude


% start with some approximate solution and then optimize
n = a^(-1/2);
p_inertial = [0; -a];
v_inertial = [n; 0];

p0 = inv(A_r2i) * p_inertial;
v0 = inv(A_r2i) * v_inertial - [-p0(2); p0(1)];

p0(3) = 0;
v0(3) = 1;

v0Norm = v0 / norm(v0);

rDist = 1 / (lunarDist/1000/(1737.1 + 100)); % final orbit around moon

multiplier = 1.2;
multiplierIter = 0.1;
breakError = rDist / 100;
tf = 1; % final integration time

% v0(1) = v0(1)*1.2545;
prevMinVal = inf;
for iter = 1:50
    v = v0;
    p = p0;
    v = v0Norm * norm(v0) * multiplier; % multiply the entire velocity vector not just the x component

    X = [p; v]; % stack together for integration

    

    %% calculate integration
    options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10); % set the integration options
    [~, stateSolution] = ode45(@stateDerivative, [0, tf], X, options); % and perform the integration

    target = [1+rDist,0,0];

    dist = abs(stateSolution(:,1:3) - target);
    dist = sqrt(dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2);
   
    [minVal, minInd] = min(dist);
    
    if minVal < rDist && stateSolution(minInd,1) > 1
        break
    end

    % now want to check if it turned back or went through the moon
    % so compare the distance from origin to moon closest approach and the origin to last point?
    % - well we really just want to look at the max x right?
    %   - no this can get stuck in cases where we go right below the moon, closest point is <1 but max is >1
    %   - should we specify the closest point be >1 or <1?
    
    % so I think the REAL metric is if it's above or below the axis at x=1
    % just want the first crossing I think
    for n = 2:size(stateSolution,1)
        if stateSolution(n-1,1) < target(1) && stateSolution(n,1) >= target(1)
            
            if stateSolution(n,2) > 0 && multiplierIter > 0
                multiplierIter = -multiplierIter / 2;
            elseif stateSolution(n,2) < 0 && multiplierIter < 0
                multiplierIter = -multiplierIter / 2;
            end
            
            break;
        end
        
        if n == size(stateSolution,1) && multiplierIter < 0 % didn't find a crossing (reach the moon) make sure we're increasing
            multiplierIter = -multiplierIter / 2;            
        end
        
    end
    

    %if max(stateSolution(:,1)) > 1 && multiplierIter > 0 && (stateSolution(minInd,1) > 1 || stateSolution(minInd,2) > 0) % decrease velocity
%     if (max(stateSolution(:,1)) > 1 || stateSolution(minInd,1) > 1)  && multiplierIter > 0
%         multiplierIter = -multiplierIter / 2;
%     elseif (max(stateSolution(:,1)) <= 1 || stateSolution(minInd,1) < 1)  && multiplierIter < 0% increase velocity
%         multiplierIter = -multiplierIter / 2;
%     end
%     if stateSolution(minInd,1) > 1 && multiplierIter > 0 % decrease velocity
%         multiplierIter = -multiplierIter / 5;
%     elseif stateSolution(minInd,1) <= 1  && multiplierIter < 0% increase velocity
%         multiplierIter = -multiplierIter / 5;
%     end

    multiplier = multiplier + multiplierIter;

end



%% and then can we get an approximate delta V to capture into at the moon
% X0 = stateSolution(minInd,:);
% % and we're going to get the most bang for cost just firing in the +y direction at the periapse right?
% 
% multiplier = 0.9;
% multiplierIter = -0.05;
% tf = 1;
% 
% X = X0;
% X(5) = X(5) * multiplier;
% 
%     %% calculate integration
%     options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10); % set the integration options
%     [~, stateSolution] = ode45(@stateDerivative, [0, tf], X, options); % and perform the integration
%     
%     % and this time for our goal, we want to maximize the number of close
%     % passes to the moon, indicating several orbits
%     target = [1,0,0];
%     
%     dist = abs(stateSolution(:,1:3) - target);
%     dist = sqrt(dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2);
%     
%     T = islocalmin(dist);
%     numMinimums = sum(T);
    
    
    







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


figure(2)
plot3(stateSolution(:,1), stateSolution(:,2), stateSolution(:,3))



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
