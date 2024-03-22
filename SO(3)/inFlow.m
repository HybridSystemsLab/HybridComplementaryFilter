function output = inFlow(Rhat, Ry, c)
% output = 1 if (Rtilde_y, q) \in C0
% output = 0 otherwise

Rtilde_y = Rhat' * Ry;

if fun_potential(Rtilde_y) <= c
    output = 1;
else
    output = 0;
end
end