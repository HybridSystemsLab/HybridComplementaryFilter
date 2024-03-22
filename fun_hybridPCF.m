function [output_R, output_q, output_jumps] = fun_hybridPCF(q, Rhat, Ry, norm_Rstar, Omega_y, kp, c0, c1, h, jumps)
% returns (Rhat(+), q(+), jumps(+)) for the hybrid filter on SO(2)

axis = [1;0;0];
r = axis;

if inFlowC0(Rhat, Ry, q, c0) == 1
    output_R = fun_passiveComplementaryFilter(Rhat, Ry, Omega_y, kp, h);
    output_q = q;
    output_jumps = jumps;

elseif inFlowC1(Rhat, Ry, q, c1) == 1
    theta_star = acos(1. - 2.*norm_Rstar);
    Rstar = fun_axisangle(theta_star, r);
    R1_y = Rstar*Ry;
    output_R = fun_passiveComplementaryFilter(Rhat, R1_y, Omega_y, kp, h);
    output_q = q;
    output_jumps = jumps;
   

elseif inJumpD0(Rhat, Ry, q, c0) == 1
    output_R = Rhat;
    output_q = 1 - q;
    output_jumps = jumps+1;
    
elseif inJumpD1(Rhat, Ry, q, c1) == 1
    output_R = Rhat;
    output_q = 1 - q;
    output_jumps = jumps+1;

else
    error("(Rtilde, q) does not lie in any of the flow or jump sets.")
end
end