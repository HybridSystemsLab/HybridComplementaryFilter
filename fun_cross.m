function output = fun_cross(v)
% returns the cross product matrix for v \in \R3
output = zeros(3,3);
output(1,2) = -v(3);
output(1,3) = v(2);
output(2,1) = v(3);
output(2,3) = -v(1);
output(3,1) = -v(2);
output(3,2) = v(1);
end