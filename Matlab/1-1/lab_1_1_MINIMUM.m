clc;
% МИНИМУМ (+1)
% 4. Обратный алгоритм Катхилла-Макки.
% 4.1. Создать случайную разреженную матрицу размерности 150
% с плотностью 0,01 (лучше по формуле). Сколько у нее ненулевых элементов?
S = sprandsym(150, 0.06, 0.8, 1); 
nonzero = nnz(S); % Ненулевые элементы
subplot(3,2,1), spy(S), title("S = sprand(150)");

% 4.2. Применить разложение Холецкого:
S_chol = chol(S); % для этого была создана ПО матрица в стр. 89
subplot(3,2,2), spy(S_chol), title("chol(S)");

% 4.3. Найти матрицу перестановок по обратному алгоритму Катхилла-Макки:
rkm = symrcm(S);
S_rkm = S(rkm, rkm);
subplot(3,2,3), spy(S_rkm), title("symrcm(S)");

% 4.4. Применить разложение Холецкого:
chol_rkm = chol(S_rkm);
subplot(3,2,4), spy(chol_rkm), title("chol(symrcm(S))");


% 5. Алгоритм минимальной степени
% Проделать аналогичные действия, используя для создания матрицы
% перестановок функцию symamd
sap = symamd(S);
S_sap = S(sap, sap);
subplot(3,2,5), spy(S_sap), title("symamd(S)");

chol_sap = chol(S_sap);
subplot(3,2,6), spy(chol_sap), title("chol(symamd(S))");