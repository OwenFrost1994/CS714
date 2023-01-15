clear
clc
close all

% %% 3
% a = -3:0.01:1;
% b = -2:0.01:2;
% 
% for i = 1:1:size(a,2)
%     for j = 1:1:size(b,2)
%         z = a(i) + b(j)*1i;
%         R(j,i) = abs(1 + z + z^2/2);
%     end
% end
% 
% figure(1)
% contourf(a,b,R,[0 1])
% hold on
% grid on
% colorbar
% xlabel({'Re(z)'},'FontSize',20,'Interpreter','latex');
% xlim([-3,1])
% ylabel({'Im(z)'},'FontSize',20,'Interpreter','latex');
% ylim([-2,2])
% axis square
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[100 100 800 800])
% 
% %% 4b
% a = -100:0.1:100;
% b = -100:0.1:100;
% 
% for i = 1:1:size(a,2)
%     for j = 1:1:size(b,2)
%         z = a(i) + b(j)*1i;
%         R(j,i) = abs((5*z+12)/((4-z)*(3-z)));
%     end
% end
% 
% figure(2)
% contourf(a,b,R,[0 1])
% hold on
% grid on
% colorbar
% xlabel({'Re(z)'},'FontSize',20,'Interpreter','latex');
% xlim([-100,100])
% ylabel({'Im(z)'},'FontSize',20,'Interpreter','latex');
% ylim([-100,100])
% axis square
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[100 100 800 800])

%% 4c
t0 = 0;
t1 = 3;
u0 = 1.5;
lambda = -1E6;
%TR-BDF method
k0 = 0.01;
n = (t1-t0)/k0+1;
T(1) = t0;
U(1) = u0;
for i = 2:1:n
    T(i) = T(i-1)+k0;
    U(i)=0.5*exp(lambda*T(i))+cos(T(i));
end

ks = [0.3 0.2 0.1];
figure(3)
plot(T,U,'r-','LineWidth',3)
hold on;
legendlist = {'exact solution'};
for j = 1:1:3
    k = ks(j);
    N = floor((t1-t0)/k)+1;
    u(1) = u0;
    t(1) = t0;
    for i = 2:1:N
        t(i) = t(i-1)+k;
        fun = lambda*(u(i-1)-cos(t(i-1)))-sin(t(i-1));
        fun12 = -lambda*cos(t(i-1)+0.5*k)-sin(t(i-1)+0.5*k);
        Us = (u(i-1)+k*(fun+fun12)/4)/(1-k*lambda/4);
        fun2 = -lambda*cos(t(i))-sin(t(i-1));
        u(i) = (4*Us-u(i-1)+k*fun2)/(3-lambda*k);
    end
    plot(t,u,'--','LineWidth',1)
    hold on;
    legendlist = [legendlist strcat('k=',num2str(k))];
    Error(j) = abs(u(end)-U(end));
end
xlabel({'t'},'FontSize',20,'Interpreter','latex');
ylabel({'error'},'FontSize',20,'Interpreter','latex');
axis square
title('Trapezoidal');
legend([legendlist],'Location','northeast')
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[1000 200 800 800])

%% 4b
a = -5:0.1:5;
b = -5:0.1:5;

for i = 1:1:size(a,2)
    for j = 1:1:size(b,2)
        z = a(i) + b(j)*1i;
        R(j,i) = abs((1+z/2)/(1-z/2));
    end
end

figure(4)
contourf(a,b,R,[0 1])
hold on
grid on
colorbar
xlabel({'Re(z)'},'FontSize',20,'Interpreter','latex');
xlim([-5,5])
ylabel({'Im(z)'},'FontSize',20,'Interpreter','latex');
ylim([-5,5])
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[100 100 800 800])
