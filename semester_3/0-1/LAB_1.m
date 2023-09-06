clc; clear;
% Решение с помощью fzero
% Полиномиальная функция
f1 = @(x) polyval([1,-1,-2,3,-3], x);
X1_0 = [-2.414, -1.693];
X1 = fzero(f1, X1_0);
sprintf('Exact root on chosen interval: %.15f', X1)

% Трансцендентная функция
f2 = @(x) 3*exp(x) + polyval([-5, -3], x);
X2_0 = [0.6, 1.5];
X2 = fzero(f2, X2_0);
sprintf('Exact root on chosen interval: %.15f', X2)


X1_ERR = [0.5, 3];
X3 = fzero(f1, X1_0);
sprintf('Exact root on chosen (unappropriate) interval: %.15f', X3)

% Зависимости числа итераций от заданной точности
% Для fzero
n_iters_1 = zeros(1, 15);
n_iters_2 = zeros(1, 15);
eps = zeros(1, 15);
for i = 1:15
    eps(i) = 10^(-i);
    options = optimset(tolX=eps(i));
    [~, ~, ~, output_1] = fzero(f1, X1_0, options);
    [~, ~, ~, output_2] = fzero(f2, X2_0, options);
    n_iters_1(i) = output_1.iterations;
    n_iters_2(i) = output_2.iterations;
end

% ГРАФИКИ
semilogx(eps, n_iters_1, 'red');
hold on; grid on;
semilogx(eps, n_iters_2, 'blue');
xlabel('eps'); ylabel('N iters'); title('fzero'); legend("f1", "f2");
saveas(gcf, 'fzero.png'); hold off;


% Для МПД
bisect1 = readtable('bisect1.csv');
bisect2 = readtable('bisect2.csv');
semilogx(bisect1.eps, bisect1.steps, 'red');
grid on; hold on;
semilogx(bisect2.eps, bisect2.steps, 'blue');
title('bisection method'); xlabel('eps'); ylabel('N iters'); legend("f1","f2");
saveas(gcf, 'mpd.png'); hold off;


% Для метода Ньютона
newton1 = readtable('newton1.csv');
newton2 = readtable('newton2.csv');
semilogx(newton1.eps, newton1.steps, 'red');
grid on; hold on;
semilogx(newton2.eps, newton2.steps, 'blue');
title('Newton method'); xlabel('eps'); ylabel('N iters'); legend("f1","f2");
saveas(gcf, 'newton.png'); hold off;

% Сломанный Ньютон
broken_limits = readtable('broken_newton(limits).csv');
broken_limits_x0 = readtable('broken_newton(limits+x0).csv');
subplot(2,1,1);
semilogx(newton1.eps, newton1.steps);
grid on; hold on;
semilogx(broken_limits.eps, broken_limits.steps, '-.');
semilogx(broken_limits_x0.eps, broken_limits_x0.steps, 'x');
title('broken Newton conversion'); xlabel('eps'); ylabel('N iters');

subplot(2,1,2);
loglog(newton1.eps, newton1.err);
grid on; hold on;
loglog(broken_limits.eps, broken_limits.err, '-.');
loglog(broken_limits_x0.eps, broken_limits_x0.err, 'x');
title('broken Newton errors'); xlabel('eps'); ylabel('err');
legend("No violations","0 ∈ f'", "0 ∈ f'; f(x_0)*f''(x_0)<0", 'Location', 'west');
saveas(gcf, 'broken_newton.png'); hold off;
