
function [L1,L2,L3,L4,L5] = findLagrangePoints(mu)

if mu > 1E-5
    [L1,L2,L3,L4,L5] = findLagrangePoints_CR3BP(mu);
else
    [L1,L2,L3,L4,L5] = findLagrangePoints_HillProblem(mu);
end

end



function [L1,L2,L3,L4,L5] = findLagrangePoints_HillProblem(mu)

% so as mu goes to zero, we need higher and higher tolerances to find the
% exact points using the zero crossing method
% but fortunately, they also converge to known locations

L3 = [-1; 0; 0];

% xL1 = -(1/3)^(1/3) - (1/3)*(mu/3)^(1/3) + (1/9)*(mu/3)^(3/2);
% xL2 = (1/3)^(1/3) + (1/3)*(mu/3)^(1/3) + (1/9)*(mu/3)^(3/2);


L1 = [1 - (mu/3)^(1/3); 0; 0];
L2 = [1 + (mu/3)^(1/3); 0; 0];

L4 = [0.5 - mu; +sqrt(3)/2; 0];
L5 = [0.5 - mu; -sqrt(3)/2; 0];

end



function [L1,L2,L3,L4,L5] = findLagrangePoints_CR3BP(mu)



% find zero crossings of the x-gradiant of the pseudo-potential function
% this is using the knowledge that L1, L2, and L3 will lie somewhere along y=0; z=0
% LOGIC:
%   - start at specified minimum x
%   - calculate gradiant of the pseudo-potential for (x,y,z)=(x,0,0)
%   - if the current step's gradient has a different sign from previous
%   step's gradient, we found a zero crossing
%   - reduce step size and turn back to zoom in on the zero crossing 
%   - once resolution threshold is met, print value and continue on

% define a search range in distance units
xMin = -1.5; 
xMax = 1.5;

y = 0; % known property of l1,l2,l3
z = 0; % known property of l1,l2,l3

dx = 0.01;
xResolution = dx / 1E8;

x = xMin;
dx_mod = 1;
dudx_prev = nan;
iter = 1; % track iterations needed to complete
zeroCrossings = [];
while (1)
    
    [dudx, ~, ~] = GradU(x,y,z);
    
%     fprintf("%f -- %f\n", ~isnan(dudx_prev), sign(dudx_prev) ~= sign(dudx))
    if ~isnan(dudx_prev) && sign(dudx_prev) ~= sign(dudx)
       % found a zero crossing! 
       dx_mod = -dx_mod/10;% go back and decrease resolution
    end
    
    if abs(dx*dx_mod) < xResolution
        % be satisfied with the value of this root - move on
%         fprintf("Found root at x=%f\n", x);
        zeroCrossings(length(zeroCrossings)+1) = x;
        dx_mod = 1; % reset
        dudx_prev = nan; % reset
    else
        dudx_prev = dudx; % continue search 
    end
    
    x = x + dx*dx_mod; % update the next x position
    iter = iter + 1; % track iterations if needed
    
    if x > xMax
        break
    end
end

L3 = [zeroCrossings(1); 0; 0];
L1 = [zeroCrossings(3); 0; 0];
L2 = [zeroCrossings(5); 0; 0];

L4 = [0.5 - mu; +sqrt(3)/2; 0];
L5 = [0.5 - mu; -sqrt(3)/2; 0];


end


function [dudx,dudy,dudz] = GradU(x,y,z)
global mu

% solved via mathematica

r13 = sqrt( (mu+x).^2 + y.^2 + z.^2 );
r23 = sqrt( (mu+x-1).^2 + y.^2 + z.^2 );

dudx = x - mu*(mu+x-1) / r23^3 - (1-mu).*(mu+x) / r13^3;
dudy = y - y*mu / r23^3 - y*(1-mu) / r13^3;
dudz = -z*mu / r23^3 - z*(1-mu) / r13^3;

end


