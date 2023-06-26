clc; clear;
% Решение с помощью fzero
% Полиномиальная функция
f1 = @(x) polyval([1,-1,-2,3,-3], x);
X1_0 = [0.5, 2.73205];
X1 = fzero(f1, X1_0);
sprintf('Точный корень f1: %.15f', X1)

% Трансцендентная функция
f2 = @(x) 3*exp(x) + polyval([-5, -3], x);
X2_0 = [0.5, 1.5];
X2 = fzero(f2, X2_0);
sprintf('Точный корень f2: %.15f', X2)


% Зависимости числа итераций от заданной точности
% Для fzero
n_iters_1 = zeros(1, 15);
n_iters_2 = zeros(1, 15);
eps = zeros(1, 15);
for i = 1:15
    eps(i) = 10^(-i);
    options = optimset(tolX=eps(i));
    [~, fval, exitflag, output_1] = fzero(f1, X1_0, options);
    [~, fval, exitflag, output_2] = fzero(f2, X2_0, options);
    n_iters_1(i) = output_1.iterations;
    n_iters_2(i) = output_2.iterations;
end



% Для МПД (посчитано на C++)
bisect1 = readtable('bisect1.csv');
bisect2 = readtable('bisect2.csv');

% Для метода Ньютона (посчитано на C++)
newton1 = readtable('newton1.csv');
newton2 = readtable('newton2.csv');

% Сабплот для МПД
subplot(1,2,1);
semilogy(bisect1.steps, bisect1.eps)
grid on
hold on
semilogy(bisect2.steps, bisect2.eps)
xlabel('Количество итераций'); ylabel('Точность');
legend("bisect f1","bisect f2");

% Сабплот для быстрых методов
subplot(1,2,2);
semilogy(n_iters_1, eps);
hold on
semilogy(n_iters_2, eps);
semilogy(newton1.steps, newton1.eps)
semilogy(newton2.steps, newton2.eps)
xlabel('Количество итераций'); ylabel('Точность');
legend("fzero f1", "fzero f2", "newton f1", "newton f2");
grid on

