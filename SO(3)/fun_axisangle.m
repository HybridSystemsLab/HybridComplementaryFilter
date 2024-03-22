function output = fun_axisangle(angle, axis)
% output is the rotation corresponding to the input axis angle
% angle in radians, axis should be a unit vector

if norm(axis) > 0.9999 && norm(axis) < 1.0001
    axiscross = fun_cross(axis);
    output = eye(3) + axiscross * sin(angle) + axiscross^2. * (1. - cos(angle));
else
    error("Axis is not a unit vector")
end
end



