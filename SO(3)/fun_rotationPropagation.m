function output = fun_rotationPropagation(R, omega, h)
% given rotatation matrix, angular velocity and time step,
% this function calculates the rotation matrix at the next time step.

if norm(omega) < 0.0001
    output = R;
else
    angle = h*norm(omega);
    axis = omega/norm(omega);
    incrementalRotation = fun_axisangle(angle, axis);
    output = R * incrementalRotation;
end
end