function output = keplerianOrbits_inertial(orbit, orientation)
% KEPLERIANORBITS_INERTIAL Calculate an orbit's orientation in space.

r = zeros(length(orbit.time),3);
v = zeros(length(orbit.time),3);

for n = 1:length(r)
    r(n,:) = orientation.rotation_pqr2I * orbit.r_pqr(n,:)';
    v(n,:) = orientation.rotation_pqr2I * orbit.v_pqr(n,:)';
end

output = orbit;
output.r = r;
output.v = v;

output.orientation = orientation;



end