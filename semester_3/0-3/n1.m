% Лабораторная работа №3. Приближение данных по методу наименьших квадратов
% ВАРИАНТ 70
% A) Приближение данных полиномиальным МНК
% 1. Приближение данных полиномиальным МНК при помощи функции `polyfit`
%1.1
a=0.5; b=4.5; x=a:0.05:b;
p = [1, -8.1, 19.04, -10.56];
p_legend = 'y = x^3 - 8.1x^2 + 19.04x - 10.56';
y = polyval(p, x);

plot(x, y, 'b'); grid on;
legend(p_legend); title('Polynomial');
saveas(gcf, '1.png'); hold off;

%1.2
x = linspace(a,b,150); t = linspace(a,b,1001); rng(1);
y = polyval(p,x) + randn(size(x));

% P 1
p1 = polyfit(x,y,1); y1 = polyval(p1,t);
err1 = y - polyval(p1,x);

subplot(2,1,1);
plot(x, y, 'b.');
hold on; grid on;
plot(t, y1, 'r');
xlabel('x'); ylabel('y');
title('Approximation with polynomial with degree 1'); hold off;

subplot(2, 1, 2);
plot(x, err1,'r.'); grid on;
xlabel('x'); ylabel('Err'); title('max error ', max(abs(err1)));
saveas(gcf, '2.png'); 


% P 2
p2 = polyfit(x,y,2); y2 = polyval(p2,t);
err2 = y - polyval(p2, x);

subplot(2,1,1);
plot(x, y, 'b.');
hold on; grid on;
plot(t, y2, 'r');
xlabel('x'); ylabel('y');
title('Approximation with polynomial with degree 2'); hold off;

subplot(2, 1, 2);
plot(x, err2,'r.'); grid on;
xlabel('x'); ylabel('Err'); title('max error ', max(abs(err2)));
saveas(gcf, '3.png');



% P 3
p3 = polyfit(x,y,3); y3 = polyval(p3,t);
err3 = y - polyval(p3, x);

subplot(2,1,1);
plot(x, y, 'b.');
hold on; grid on;
plot(t, y3, 'r');
xlabel('x'); ylabel('y');
title('Approximation with polynomial with degree 3'); hold off;

subplot(2, 1, 2);
plot(x, err3,'r.'); grid on;
xlabel('x'); ylabel('Err'); title('max error ', max(abs(err3)));
saveas(gcf, '4.png');



% P 4
p4 = polyfit(x,y,4); y4 = polyval(p4,t);
err4 = y - polyval(p4, x);

subplot(2,1,1);
plot(x, y, 'b.');
hold on; grid on;
plot(t, y4, 'r');
xlabel('x'); ylabel('y');
title('Approximation with polynomial with degree 4'); hold off;

subplot(2, 1, 2);
plot(x, err4,'r.'); grid on;
xlabel('x'); ylabel('Err'); title('max error ', max(abs(err4)));
saveas(gcf, '5.png');