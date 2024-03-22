function [output_R, output_jumps] = fun_hybridPCF(Rhat, Ry, Omega_y, kp, c, h, jumps)
% returns (Rhat(+), q(+), r(+)) for the hybrid observer on SO(2)

axis = [1;0;0];
r = axis;

if inFlow(Rhat, Ry, c) == 1
    output_R = fun_passiveComplementaryFilter(Rhat, Ry, Omega_y, kp, h);
    output_jumps = jumps;
   
elseif inJump(Rhat, Ry, c) == 1
    fprintf('in Jump D0 \n')
    output_R = Ry;
    output_jumps = jumps+1;

else
    error("Rtilde does not lie in the flow or jump sets.")
end
end