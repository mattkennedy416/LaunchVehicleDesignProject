function output = defineOrientation(inclination, rightAscension, argumentOfPeriapsis)
% DEFINEORIENTATION Calculate and hold orbit orientation information.

output = struct();
output.i = inclination;
output.OMEGA = rightAscension;
output.omega = argumentOfPeriapsis;

output.rotation_pqr2I = R3(-rightAscension) * R1(-inclination) * R3(-argumentOfPeriapsis);

end




function rot = R1(val)

rot = [1, 0, 0;
    0, cos(val), sin(val);
    0, -sin(val), cos(val)];

end

function rot = R3(val)

rot = [cos(val), sin(val), 0;
    -sin(val), cos(val), 0;
    0, 0, 1];

end