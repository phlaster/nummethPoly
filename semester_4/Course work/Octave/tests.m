clc;
F = b5d(4);
spy(F)
full(F)

C = chol(F, "lower");
spy(C)
full(C)

iC = ichol(F);
spy(iC)
full(iC)

options = struct('type', 'ict', 'droptol', 0.3);
iC_t = ichol(F, options);
spy(iC_t)
full(iC_t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [
    7.7   4.36  7.79   5.77  0.40
    4.36  7.26  8.17   3.29  0
    7.79  8.17  12.87  7.23  1.34
    5.77  3.29  7.23   7.87  1.36
    0.40  0     1.34   1.36  1.20
]
C = chol(A, 'lower');
full(C)
iC = ichol(sparse(A));
full(iC)
iCt = ichol(sparse(A), struct('type', 'ict', 'droptol', 0.05));
full(iCt)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ПРЕДОБУСЛОВЛЕННЫЙ МЕТОД СОПРЯЖЕННЫХ ГРАДИЕНТОВ
% 5-и диагональная матрица (МКР для уравнения Пуассона)
% конструирование матрицы для 5-ти точечного разностного оператора Лапласа
n=30;
A = b5d(30);
% создание правой части системы для точного решения [1, 2, 3, ...]
x_ex = (1:n*n)';
b = A * x_ex;
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
