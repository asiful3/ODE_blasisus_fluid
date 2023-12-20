
clc
close all

[eta,f] = ode45(@blas, [0, 10], [0, 0, 0.33]);

figure(1)
plot(eta, f(:,1),'DisplayName','F(eta)')
set(gca,'FontSize', 14)
xlabel('eta')
ylabel('u')
title('Boundary Layer')
grid
axis([0 10 0 5])
lgd = legend
hold on;
plot(eta, f(:,2),'DisplayName','G(eta)')
plot(eta, f(:,3),'DisplayName','H(eta)')
hold off;

g = eta - f(:,1)

figure(2)
plot (eta, g)
lgd = legend
