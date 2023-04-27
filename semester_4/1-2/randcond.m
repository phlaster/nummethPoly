function Mout = randcond(side, cnum)
% Random square matrix with fixed condition number
    [Q, R] = qr(rand(side));
    D = diag(linspace(1, cnum, side));
    Mout = Q'*D*Q;
end

