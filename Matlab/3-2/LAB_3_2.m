clc;
% Лабораторная работа №3. (часть 2)
% Решение интегралов в МatLab. Методы Симпсона, Гаусса, Чебышева
% Задача. Использование встроенных функции MatLab и методов Симпсона,
%Гаусса и Чебышева. Вычислить интеграл от функции
% g(x) = f(x)*sin(2*(x-a)*4^)
% где f(x) – полином из части 1, (вариант 23)
% f(x) = x^5 - 3.2*x^3 + 9.5*x^2 - 7*x - 7.5
% a – левая граница интервала ВАРИАНТ 23
a = -0.7;
b = 1.6;
% fplot(@(x) g(x,a),[a,b])
N = 2000; % точки разбиения
% БАЗА (0) Решение интеграла с заданным числом разбиений методом Симпсона. 
% Создать функцию для метода Симпсона:
% Баловство, т.к. тут отрезок не разбивается, так что ошибка очень большая 
simp3 = s3_f(a, b); % для разбиения на 1 кусок :)
% Ньютон-Лейбниц, т.к. первообразная взята ручками:
prescise = pervoobr(b) - pervoobr(a);
% Соответственно, ошибка > 1
e = abs(simp3 - prescise);

% То же самое для полинома f, но уже с разбиением отрезка на N-1 кусков
otrezok = linspace(a,b,N);
simpsonFull_f = 0;
error_f = []
for i = 1:N-1
    simpsonFull_f = simpsonFull_f + s3_f(otrezok(i), otrezok(i+1));
    error_f = [error_f abs(simpsonFull_f - prescise)];
end
error_f(end) % Уже вменяемая ошибка
% loglog(error_f); % Просто поглядеть, как приближаемся к точному значению


% МИНИМУМ (+1) Применение функции integral. Построение графика
clc;
error_g = [];
for N = 1:2000 % Цикл по измельчению разбиения
   simpsonFull_g = 0;
   otrezok = linspace(a,b,N); % Измельчаем
   for i = 1:N-1
        simpsonFull_g = simpsonFull_g + s3_g(a, otrezok(i), otrezok(i+1));
    end
    
    % Итоговая ошибка вычисления для каждого разбиения
    error_g = [error_g abs(simpsonFull_g - integral(@(x) g(x, a), a, b))];
end

error_g(end); % Пришли сюда
loglog(error_g); % Смотрим


% ДОСТАТОЧНО (+1) Решение интеграла с помощью формул Гаусса с заданными
% узлами и коэффициентами (для 4-х узлов)
c0 = sqrt(30)/36;
c1 = 2*sqrt(30)/35;
p1 = 1/2;
p2 = 3/7; % константы метода для 4х слагаемых
param = (b+a)/2 + (b-a)/2;
s4 = @(f, x) (b-a)/2 * ( ...
        (p1 - c0)*f(-sqrt(p2 + c1) * param) +...
        (p1 + c0)*f(-sqrt(p2 - c1) * param) +...
        (p1 + c0)*f( sqrt(p2 - c1) * param) +...
        (p1 - c0)*f( sqrt(p2 + c1) * param)...
     );
s4(@(x) g(x, a), 0); % Вроде работает...


% Далее to do то же, что и в МИНИМУМе...
% error_g = [];
% for N = 1:2000
%    simpsonFull_g = 0; 
%  otrezok = linspace(a,b,N);
%    for i = 1:N-1
%         simpsonFull_g = simpsonFull_g + s3_g(a, otrezok(i), otrezok(i+1));
%     end
%     error_g = [error_g abs(simpsonFull_g - integral(@(x) g(x, a), a, b))];
% end
% 
% error_g;
% loglog(error_g);
