function output = inJumpD0(Rhat, Ry, q, c0)
% output = 1 if (Rtilde_y, q) \in D0
% output = 0 otherwise

Rtilde_y = Rhat' * Ry;

if q == 0 && fun_potential(Rtilde_y) >= c0
    output = 1;
else
    output = 0;
end
end