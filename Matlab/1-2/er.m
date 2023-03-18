function out = er(A, B, X)
% Погрешность решения СЛАУ
    out = norm(B-A*X)/norm(B);
end

