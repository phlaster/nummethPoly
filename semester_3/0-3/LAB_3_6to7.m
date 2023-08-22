%6 cftool нелинейные модели МНК
%6.1
a0=-3; b0=3;
x=a0:0.05:b0;

a=2.0; b=0.8; c=2.8; d=0.3;
y = exp(a-b*x.^2) .* cos(sqrt(d+c * x.^2));

plot(x,y,'b');

%6.2
x=linspace(a0,b0,150);
y=exp(a-b*x.^2) .* cos(sqrt(d+c * x.^2)) + randn(size(x));


%7 fit: приближение нелинейными моделями МНК
NlModel = 'exp(a-b*x*x).*cos(sqrt(d+c*x*x))';
StartPoint=[2,1,3,0];
f1 = fit(x',y',NlModel,'Start', StartPoint)
plot(f1,x,y)






