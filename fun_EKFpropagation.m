function [Rhat_Kp1givenK, P_Kp1givenK] = fun_EKFpropagation(Rhat_KgivenK, P_KgivenK, omega_y, deltaT)

Rhat_Kp1givenK = fun_rotationPropagation(Rhat_KgivenK, omega_y, deltaT);
Fk = eye(3);
P_Kp1givenK = Fk * P_KgivenK * Fk';
end