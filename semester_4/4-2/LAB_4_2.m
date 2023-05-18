clc;
% Вариант 11:
% xy^2y' = x^2 + y^3;           x in [1.1,3];   y^3 = 3x^2*(x-1)
% xy'/y  = x^2/(3x^2*(x-1)) + 1                 y^2 = cbrt(sqr(3x^2*(x-1)))
% y'/y   = 1/x   + x/(3x^2*(x-1))               y   = (3x^2*(x-1)).^1/3
% y'     = y/x + (3x^2*(x-1)).^(1/3) / (3x^2*(x-1))
a = 1.1;
b = 3.0;
ns = [10 1000];
for i = 1:2
    n = ns(i);
    h = (b-a)/n;

    x = linspace(a, b, n);


    % Явный метод

    y_presc = [];
    y_appr  = [];
    y_i = y(a);
    for x_i = x(1:end-1)
        y_presc = [y_presc y(x_i)];
        y_appr  = [y_appr  y_i];
        y_i = y_i + h * dif_y(x_i, y_i);
    end

    subplot(2, 2, i)
    plot(x(1:n-1),y_presc);
    xlabel('X')
    ylabel('Y')
    title(n)
    hold on;
    grid on;
    plot(x(1:n-1),y_appr);
    legend("Точное","Численное")
end

errors = [];
hs = []
for n = 10:1000
    h = (b-a)/n;
    hs = [hs h];
    x = linspace(a, b, n);
    y_i = y(a);
    maxerr = 0;
    for x_i = x(1:end-1)
        y_i = y_i + h * dif_y(x_i, y_i);
        if abs(y_i-y(x_i)) > maxerr
            maxerr = abs(y_i-y(x_i));
        end
    end
    errors = [errors maxerr];
        
end

subplot(2, 2, 3)
loglog((b-a)./(10:1000), errors);
xlabel('n')
ylabel('err')
title("По количеству узлов")
hold on;
grid on;
% semilogy(10:1000, hs);
% legend("err", "h")
% Похоже, ошибка в исходном преобразовании формул. Графики неверные!
    
    
    
% Неявный метод
% ...

