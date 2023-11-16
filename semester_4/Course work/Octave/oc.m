clc;
F = b5d(5)
spy(F)
figure
spy(A)
cond(full(A))

droptol = 1e-5; %1e-2
options = struct('type', 'ict', 'droptol', droptol);
iL = ichol(A, options);
L = chol(A, 'lower');
spy(L)

P = iL*iL';
%norm(A - P, 1)
cond(full(inv(P)*A))
subplot(1,3,2)
spy(iL)
subplot(1,3,3)
spy(L)

% ПРЕДОБУСЛОВЛЕННЫЙ МЕТОД СОПРЯЖЕННЫХ ГРАДИЕНТОВ
% 5-и диагональная матрица (МКР для уравнения Пуассона)
% конструирование матрицы для 5-ти точечного разностного оператора Лапласа
n = 30; I = ones(n,1); D = spdiags(-I, 0, n, n); T = spdiags([I -2*I I], [-1 0 1], n, n);
A = kron(D, T) + kron(T, D);
nnz(A)
% создание правой части системы для точного решения [1, 2, 3, ...]
x_ex = (1:n*n)'; b = A * x_ex;
% в цикле изменяем порог droptol в неполном разложении Холецкого R'R = AN
% и используем факторизованный предобусловливатель R'R
DROPTOL = logspace(-6, 0, 50); ITER = []; TIME = [];
for droptol = DROPTOL
    droptol
    options = struct('type', 'ict', 'droptol', droptol);
    t = cputime;
    iL = ichol(A, options);
    nnz(iL)
    [x, flag, relres, iter] = pcg(A, b, 1e-5, 1000, iL', iL);
    t = cputime - t;
    TIME = [TIME t];
    ITER = [ITER iter];
end
norm(x - x_ex) / norm(x_ex);
figure('Color', 'w')
loglog(DROPTOL, ITER)
xlabel('Порог'), ylabel('Количество итераций')
% сопряженные градиенты без предобусловливания
t = cputime;
[x, flag, relres, iter] = pcg(A, b, 1e-5, 1000);
t = cputime - t;
hold all
loglog([DROPTOL(1) DROPTOL(end)], [iter iter])
grid on
legend('c предобусловливателем', 'без предобусловливателя')
