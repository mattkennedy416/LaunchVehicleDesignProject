function output = keplerianOrbits_inPlane(perigee, apogee, mu)
%KEPLERIANORBITS_INPLANE Calculate the in-plane orbit of a body
%   Detailed explanation goes here


a = 0.5*(perigee + apogee);
e = 1 - perigee/a;
P = a*(1-e^2);

n = sqrt(mu/a^3); % mean motion
T = 2*pi/n; % period


time = 0:T/1000:T;
time = time';


nu(1) = 0;
for i = 2:length(time)
    
    E0 = 2*atan( tan(nu(1)/2) * sqrt( (1-e)/(1+e) ) );
    M0 = E0 - e*sin(E0);
    M1 = M0 + n*(time(i) - time(1));
    
    %E1 = M1 + e*sin(time(n-1));
    E1 = 0;
    for ii = 1:6
        E1 = M1 + e*sin(E1);
    end
    
    nu(i,1) = 2*atan( sqrt( (1+e)/(1-e) ) * tan(E1/2) );  
    
end

r_pqr = [P*cos(nu)./(1 + e*cos(nu)), P*sin(nu)./(1 + e*cos(nu)),zeros(length(nu),1)];

v_pqr = [-sqrt(mu/P)*sin(nu),sqrt(mu/P)*(e + cos(nu)),zeros(length(nu),1)];


output = struct();
output.perigee = perigee;
output.apogee = apogee;
output.mu = mu;
output.a = a;
output.e = e;
output.P = P;
output.T = T;
output.time = time;
output.nu = nu;
output.r_pqr = r_pqr;
output.v_pqr = v_pqr;

%plot(r_pqr(:,1), r_pqr(:,2))

% output.plot = @plotOrbit;


end



function plotOrbit()

plot(r_pqr(:,1), r_pqr(:,2))

end


