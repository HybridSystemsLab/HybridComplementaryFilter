function output = inJumpD1(Rhat, Ry, q, c1)
% output = 1 if (Rtilde_y, q) \in D1
% output = 0 otherwise

Rtilde_y = Rhat' * Ry;

if q == 1 && fun_potential(Rtilde_y) <= c1
    output = 1;
else
    output = 0;
end
end