% Функция f1 из задания:
f1 = @(x) cot(x) + x.^2;
nNodes = 5;
x = linspace(0.5, 2.75, nNodes);
y = f1(x);

% Вызов функции csape с различными граничными условиями:
% По умолчанию функция применяет конечные условия Лагранжа к каждому концу
% данных и сопоставляет производные концов сплайна с производной кубического
% полинома, который соответствует последним четырем точкам данных на каждом конце.
spline_default = csape(x,y);

% Полный кубический сплайн: проходит через все заданные точки данных.
spline_complete = csape(x,y,'complete');

% Условие не-узла требует, чтобы первая и последняя производные были
% одинаковы на конечных границах сплайна.
spline_not_a_knot = csape(x,y,'not-a-knot');

% Периодический кубический сплайн, то есть значение и производная на концах
% сплайна должны быть равны.
spline_periodic = csape(x,y,'periodic');

% Кубический сплайн с условием второй производной которое требует,
% чтобы первая и последняя вторые производные были равны нулю
% на конечных границах сплайна.
spline_second = csape(x, y, 'second');

% Кубический сплайн с использованием вариационного подхода, который
% минимизирует интеграл от квадрата первой производной по всей
% области интерполяции.
spline_variational = csape(x, y, 'variational');



% Оценка значений сплайнов на новом наборе точек для отображения результатов
xx = linspace(0.5, 2.75, nNodes*100);
yy = f1(xx);
yy_default = fnval(spline_default, xx);
yy_complete = fnval(spline_complete,xx);
yy_not_a_knot = fnval(spline_not_a_knot,xx);
yy_periodic = fnval(spline_periodic,xx);
yy_second = fnval(spline_second,xx);
yy_variational = fnval(spline_variational,xx);

errs_default = abs(yy-yy_default);
errs_complete = abs(yy-yy_complete);
errs_not_a_knot = abs(yy-yy_not_a_knot);
errs_periodic = abs(yy-yy_periodic);
errs_second = abs(yy-yy_second);
errs_variational = abs(yy-yy_variational);

% Визуализация результатов
figure;
subplot(2,1,1);
hold on;
plot(x,y,'o');
plot(xx,yy, '--');
plot(xx,yy_default);
plot(xx,yy_complete);
plot(xx,yy_not_a_knot);
plot(xx,yy_periodic);
plot(xx,yy_second);
plot(xx,yy_variational);
grid on;

legend('Узлы','Точно','По умолчанию','Полный','Не-узел', 'Периодический', '2-я производная', 'Вариационный');
xlabel('x');
ylabel('y');
title('Кусочно-кубические сплайны с различными граничными условиями');

subplot(2,1, 2);
hold on;
semilogy(xx,errs_default);
semilogy(xx,errs_complete);
semilogy(xx,errs_not_a_knot);
semilogy(xx,errs_periodic);
semilogy(xx,errs_second);
semilogy(xx,errs_variational);
set(gca,'YScale','log')

grid on;

legend('По умолчанию','Полный','Не-узел', 'Периодический', '2-я производная', 'Вариационный');
xlabel('x');
ylabel('eps');
title('Профили ошибок различных граничных условий');
