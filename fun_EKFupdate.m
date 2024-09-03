function [Rhat_Kp1givenKp1, P_Kp1givenKp1] = fun_EKFupdate(Rhat_Kp1givenK, P_Kp1givenK, Ry, Rk)
% Rk is the covariance of v, where v is such that Ry = exp(v_x)

Hk = eye(3);

% yk = Ry - Rhat_Kp1givenK;
yk = eye(3) - Ry' * Rhat_Kp1givenK;
Sk = Hk * P_Kp1givenK * Hk' + Rk;
Kk = P_Kp1givenK * Hk' * inv(Sk);

Rhat_Kp1givenKp1 = Rhat_Kp1givenK + Kk * yk;

[U, S, V] = svd(Rhat_Kp1givenKp1);

Rhat_Kp1givenKp1 = U * diag([1, 1, det(U*V')]) * V'; % Projecting onto SO(3)

P_Kp1givenKp1 = (eye(3) - Kk*Hk)*P_Kp1givenK;
end