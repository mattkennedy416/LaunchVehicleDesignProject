

perigeeLG = 3000; % km
apogeeLG = 70000; % km

moonMass = 7.34767309E22; % kg
earthMass = 5.972E24; 
sunMass = 1.989E30;

muMoon = 4667.9; % km^3/s^2

muMoonEarth = moonMass / (earthMass + moonMass); % mass ratios
muEarthSun = earthMass / (earthMass + sunMass);


orientationLG = defineOrientation(pi/2, 0, 0);

orbitLunarGateway = keplerianOrbits_inPlane(perigeeLG, apogeeLG, muMoon);
orbitLunarGateway = keplerianOrbits_inertial(orbitLunarGateway, orientationLG); % this isn't working ...

% figure(1)
% hold on
% plot3(orbitLuarGateway.r(:,1), orbitLuarGateway.r(:,2), orbitLuarGateway.r(:,3))
% xlabel('x')
% ylabel('y')
% zlabel('z')
% 
% plot3([0],[0],[0], 'o')
% plot3([-384400], [0], [0], 'o')



% %% need to convert from inertial values to relative values
lunarDist = 384402; % km - distance from moon to earth
earthDist = 1.496E8; % km - 1au

soiMoonEarth = (moonMass)^(2/5)*lunarDist;

% ** need the 3d transform between CR3BP and inertial

nu = 0; % orbital angle of secondary body
A_r2i = [cos(nu), -sin(nu); sin(nu), cos(nu)];

p = inv(A_r2i) * orbitLunarGateway.r(1,:)';
v = inv(A_r2i) * orbitLunarGateway.v(1,:)' - [p(3); -p(2); p(1)];


% p = [0; 1; 0];
% v = [2.95; 0; 0.4];

X = [p; v]; % stack together for integration

tf = 2; % final integration time

stateSolution = integrateCR3BP(X, tf);

plot3(stateSolution(:,1), stateSolution(:,2), stateSolution(:,3))
