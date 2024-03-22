function output = fun_potential(R)
% potential function
    output = trace(eye(3) - R);
    output = output/4.;
end