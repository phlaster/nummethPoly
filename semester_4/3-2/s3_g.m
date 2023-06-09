function result = s3_g(a0, a, b)
% a0 - левая граница отрезка интегрирования
% a - левая граница отрезка разбиения
% g - объявлена глобально, чтобы не передавать функцию в функцию
    result = (b-a)/6 *(...
                  g(a, a0) +...
                4*g((a+b)/2, a0) +...
                  g(b, a0)...
              );
end

