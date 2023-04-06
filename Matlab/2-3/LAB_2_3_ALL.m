clc;
% Лабораторная работа №2. (часть 3)
% Устойчивость проблемы поиска собственных значений
% Задача. Исследовать устойчивость поиска собственных значений на
% матрице
% БАЗА (0)
%1. Создать матрицу 3х3. В качестве ε взять 10-3
N = 3;
eps = 1e-3;
A = zeros(N);

for i=1:N-1
    A(i, i+1) = 1.0;
end
A(N, 1) = eps;

% 2. С помощью функции eig вычислить собственные числа
eig(A);

%3. Вычислить точные собственные значения по формуле 𝜆 = 3√𝜀. Учесть,
% что результат – комплексное число с модулем 0.1
r = nthroot(eps, N);
eigens = [];
for i=1:N
    [x,y] = pol2cart(i*2*pi/N,r);
    eigens = [eigens; x+1i*y];
end
eigens;

% 4. Вычислить фактическую ошибку для всех собственных значений
errors(1:N) = abs(eigens - eig(A));

% МИНИМУМ (+1)
% 5. Изменить размер матрицы на 4х4, 10х10,100х100. За ε брать 10-n,
% где n – размер матрицы. Проследить за изменением ошибки в вычислении
% собственных значений
% ДОСТАТОЧНО (+1)
% 6. В цикле по размеру матрицы от 3 до 100 создать массив ошибок
% 7. Вывести на график зависимость ошибки от размера матрицы
clear;
err_progression = [];
for size = 3:100
    eps = 10^-size;
    A = zeros(size);
    for i=1:size-1
        A(i, i+1) = 1.0;
    end
    A(size, 1) = eps;
    r = nthroot(eps, size);
    eigens = [];
    for i=1:size
        [x,y] = pol2cart(i*2*pi/size,r);
        eigens = [eigens; x+1i*y];
    end
    % Комплексная сортировка
    errors(1:size) = abs(cplxpair(eigens) - cplxpair(eig(A)));
    err_progression = [err_progression max(errors)];
end
err_progression;

semilogy(3:100, err_progression)
grid on
hold on
xlabel('Размер матрицы')
ylabel('Ошибка')
title("Сравнение")


%%%%%%%%%
%%%%%%%%%
% Q-R алгоритм
size = 5
eps = 10^-size;
A = zeros(size);
for i=1:size-1
    A(i, i+1) = 1.0;
end
A(size, 1) = eps;
A

for i=1:5 % Несколько применений к одной и той же матрице
    [Q, R] = qr(A);
    A = R * Q   % Посмотреть, как эпсилон бегает по ненулевым элементам
end
% Как видим, QR алгоритм не работает для такой матрицы

