function [Rhat_Kp1givenKp1, P_Kp1givenKp1] = fun_EKF(Rhat_KgivenK, P_KgivenK, Ry, omega_y, Rk, deltaT)

[Rhat_Kp1givenK, P_Kp1givenK] = fun_EKFpropagation(Rhat_KgivenK, P_KgivenK, omega_y, deltaT);
[Rhat_Kp1givenKp1, P_Kp1givenKp1] = fun_EKFupdate(Rhat_Kp1givenK, P_Kp1givenK, Ry, Rk);
end