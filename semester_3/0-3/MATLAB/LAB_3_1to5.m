% Лабораторная работа №3. Приближение данных по методу наименьших квадратов
% A) Приближение данных полиномиальным МНК
% 1. Приближение данных полиномиальным МНК при помощи функции `polyfit`
%1.1
a=0.5; b=4.5; x=a:0.05:b;
p = [1, -8.1, 19.04, -10.56];
p_legend = 'y = x**3 - 8.1x**2 + 19.04x - 10.56'
y = polyval(p, x);
plot(x, y, 'b')
legend(p_legend)
title('График полинома')

%1.2
x = linspace(a,b,150);

%1.3
y=polyval(p,x)+randn(size(x));

%1.4
p1=polyfit(x,y,1);

%1.5
t=linspace(a,b,1001); y1=polyval(p1,t);
subplot(2,1,1);

plot(x,y,'b.', 'DisplayName', 'Данные');
hold on;
plot(t,y1,'r', 'DisplayName', 'Прибл полином ст.1');
xlabel('x');
ylabel('y');
title('Приближение полиномом 1-й степени');
legend;

err1=y-polyval(p1,x);
subplot(2, 1, 2); 
plot(x,err1,'r.', 'DisplayName', 'Ошибка');
xlabel('x');
ylabel('Err');
title('Ошибка приближения');
legend;



%1.6
p2=polyfit(x,y,2); y2=polyval(p2,t);
subplot(2,1,1);

plot(x,y,'b.', 'DisplayName', 'Данные');
hold on;
plot(t,y2,'g', 'DisplayName', 'Прибл. полином ст.2');
xlabel('x');
ylabel('y');
title('Приближение полиномом 2-й степени');
legend;

err2=y-polyval(p2,x);

subplot(2, 1, 2); 
plot(x,err2,'g.', 'DisplayName', 'Ошибка');
xlabel('x');
ylabel('Err');
title('Ошибка приближения');
legend;



%1.7
p3=polyfit(x,y,3); y3=polyval(p3,t);
subplot(2,1,1);

plot(x,y,'b.', 'DisplayName', 'Данные');
hold on;
plot(t,y3,'g', 'DisplayName', 'Полином 3 степени', color='green');
xlabel('x');
ylabel('y');
% title('Приближение полиномом 3-й степени');
legend;

err3=y-polyval(p3,x);

subplot(2, 1, 2); 
plot(x,err3,'g.', 'DisplayName', 'Ошибка 3 степени', color='green');
xlabel('x');
ylabel('Err');
% title('Ошибка приближения');
legend;

hold on;

p4=polyfit(x,y,4); y4=polyval(p4,t);
subplot(2,1,1);

plot(t,y4,'g', 'DisplayName', 'Полином 4 степени', color='blue');
xlabel('x');
ylabel('y');
title('Приближение полиномами 3 и 4 степени');
legend;

err4=y-polyval(p4,x);

subplot(2, 1, 2); 
plot(x,err4,'g.', 'DisplayName', 'Ошибка 4 степени', color='blue');
xlabel('x');
ylabel('Err');
title('Ошибки приближения');
legend;


%2 cftool МНК:
% Степень: I       II       III      IV
% p1 =  -0.3316  -0.5113   0.9649   0.0659
% p2 =   2.128    2.225   -7.748    0.3059
% p3 =           -0.3762  17.97    -5.505
% p4 =                    -9.586   15
% p5 =                             -8.35
% Наиболее близкие значения при приближении 3-й степени.
% Если степень полинома на 1 меньше количества точек, то ошибка, в идеале,
% Равна нулю, т.к. такой полином проходит через каждую точку данных.
% В действительности, остаётся ошибка округления машинной арифметики.


%3 cftool c выбросами
y1 = y; y1(1:10:end) = 5;
%3.3
% bisquare I       II      III      IV
% p1 =  -0.5086  -0.544    1      -1.623e-15
% p2 =   2.933    2.282   -8.1     1
% p3 =           -0.0615  19.04   -8.1
% p4 =                   -10.56   19.04
% p5 =                           -10.56

%3.4
% LAR     I       II      III      IV
% p1 =  -0.8056  -0.4      1      -7.968e-14
% p2 =   3.747    1.318   -8.1     1
% p3 =            1.347   19.04   -8.1
% p4 =                   -10.56   19.04
% p5 =                           -10.56

%4 cftool сглаживающие сплайны

%5.1
y=polyval(p,x)+randn(size(x))*0.5;
plot(x,y,'b.'); hold on

%5.2
pp = csaps(x, y, 1);

%5.3
yt = fnval(pp, t); plot(t, yt,'r')















