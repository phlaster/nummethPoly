function A = fix_cond(n, k)
    % Задание случайной матрицы с фиксированным числом обусловленности
    R = randn(n, n);
    [U, ~, V] = svd(R);
    S_new = diag(linspace(1, k, n));
    A = U * S_new * V';
end




