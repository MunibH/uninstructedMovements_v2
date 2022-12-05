
close all;
clear x y z out
x = 5*sin(1:0.01:100)';
y = x + randn(size(x));

% ----------------------


N = 51;

out = mySmooth(y,N,'reflect');

% ----------------------

figure;
subplot(3,1,1); plot(x);
xlim([-100 200])
subplot(3,1,2); plot(y)
xlim([-100 200])
subplot(3,1,3); plot(out); %hold on; plot(out,'.','MarkerSize',6)
xlim([-100 200])






