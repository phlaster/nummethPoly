x = load('3-1.txt');
x_r = load('3-1_runge.txt');
x_R(:,2) % выбрать строки с максимальным n

loglog(x)
hold on;
loglog(x_R(:,1).*2, x_R(:,2))