function output = inFlowC1(Rhat, Ry, q, c1)
% output = 1 if (Rtilde_y, q, r) \in FlowC1
% output = 0 otherwise

Rtilde_y = Rhat' * Ry;

if q == 1 && fun_potential(Rtilde_y) >= c1 % && fun_evalJumpCondition(r, Rtilde_y, Rhat, norm_Rstar) >= delta
    output = 1;
else
    output = 0;
end
end