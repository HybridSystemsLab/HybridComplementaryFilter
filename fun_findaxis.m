function axis = fun_findaxis(R)
% finds the axis of rotation for a given rotation matrix

potential = fun_potential(R);

if potential <= 0.00001
    axis = [1.; 0.; 0.;]; % setting arbitrary axis since all axis can be considered
    
elseif potential > 0.99999
    [V, D] = eig(R);
    if (real(D(1,1))) >= 0.99999999
        axis = V(:,1);
    elseif (real(D(2,2))) >= 0.99999999
        axis = V(:,2);
    elseif (real(D(3,3))) >= 0.99999999
        axis = V(:,3);
    else
        error("No eigenvalue 1 for 180 deg rotation")
    end
    
else
    skew = fun_skewSymmetricMatrix(R); % sin theta * cross(axis)
    axis = fun_vex(skew);
    axis = axis/norm(axis);
    
end
    