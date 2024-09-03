function output = fun_passiveComplementaryFilter(Rhat, Ry, Omega_y, kp, h)
% This function propogates the estimate of R_hat
% output = Rhat(t+1) \in SO(3)

Rtilde_y = Rhat'*Ry;

% if fun_potential(Rtilde_y) > 0.9999
%     error("Stuck at the 180 deg point")
% end

w = fun_vex(fun_skewSymmetricMatrix(Rtilde_y));
Omega = Omega_y + kp*w;

output = fun_rotationPropagation(Rhat, Omega, h);
end