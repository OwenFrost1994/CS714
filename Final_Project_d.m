clear
clc
close all

L = 1;%length
H = 1;%hight
T = 2;%time
deltat = 0.02;%time step
n = 90;
Re = 10000;
deltax = 1/n;%space step x
deltay = 1/n;%space step y
N = n+1;
K = T/deltat+1;

Ax = sparse(N*N,N*N);
Ay = sparse(N*N,N*N);
Bx = sparse(N*N,N*N);
By = sparse(N*N,N*N);
I = eye(N*N);
W0 = sparse(N*N,1);
W1 = sparse(N*N,1);
P0 = sparse(N*N,1);
U0 = sparse(N*N,1);
V0 = sparse(N*N,1);
for j = 1:1:N
    for i = 1:1:N
        x(i,j) = deltax*(i-1);
        y(i,j) = deltay*(j-1);
        if i == 1 || j == 1 || i == N || j == N
            Ax(i+(j-1)*N,i+(j-1)*N) = 1;
        else
            Ax(i+(j-1)*N,i+(j-1)*N+1) = 1;
            Ax(i+(j-1)*N,i+(j-1)*N-1) = 1;
            Ax(i+(j-1)*N,i+(j-1)*N) = -2;
        end
        if i == 1
            Bx(i+(j-1)*N,i+(j-1)*N+1) = 2;
            Bx(i+(j-1)*N,i+(j-1)*N) = -2;
        else if i == N
                Bx(i+(j-1)*N,i+(j-1)*N-1) = -2;
                Bx(i+(j-1)*N,i+(j-1)*N) = 2;
            else
                Bx(i+(j-1)*N,i+(j-1)*N+1) = 1;
                Bx(i+(j-1)*N,i+(j-1)*N-1) = -1;
            end
        end
        if i == 1 || j == 1 || i == N || j == N
            Ay(i+(j-1)*N,i+(j-1)*N) = 1;
        else
            Ay(i+(j-1)*N,i+j*N) = 1;
            Ay(i+(j-1)*N,i+(j-2)*N) = 1;
            Ay(i+(j-1)*N,i+(j-1)*N) = -2;
        end
        if j == 1
            By(i+(j-1)*N,i+(j)*N) = 2;
            By(i+(j-1)*N,i+(j-1)*N) = -2;
        else if j == N
                By(i+(j-1)*N,i+(j-2)*N) = -2;
                By(i+(j-1)*N,i+(j-1)*N) = 2;
            else
                By(i+(j-1)*N,i+(j)*N) = 1;
                By(i+(j-1)*N,i+(j-2)*N) = -1;
            end
        end
        W0(i+(j-1)*N,1) = 32*(y(i,j)^2-y(i,j)^3)*(1-6*x(i,j)+6*x(i,j)^2)+32*(1-3*y(i,j))*x(i,j)^2*(1-x(i,j))^2;
        ub(i) = 16*x(i,j)^2*(1-x(i,j))^2;
    end
end
Ax = Ax/deltax^2;
Ay = Ay/deltay^2;
A = Ax+Ay;
clearvars Ax
clearvars Ay
Bx = Bx/deltax/2;
By = By/deltay/2;

P0 = A\W0;
P0 = Thomsnz(N,P0,W0,deltax,deltay,ub);
U0 = -By*P0;
V0 = Bx*P0;
Ws = W0;
Ps = P0;
Us = U0;
Vs = V0;

for k = 2:1:K
    f1 = -U0.*(Bx*W0)-V0.*(By*W0)+A*W0/Re;
    
    w2 = W0+deltat*f1/2;
    p2 = A\w2;
    f2 = -(-By*p2).*(Bx*w2)-(Bx*p2).*(By*w2)+A*w2/Re;
    
    w3 = W0+deltat*f2/2;
    p3 = A\w3;
    f3 = -(-By*p3).*(Bx*w3)-(Bx*p3).*(By*w3)+A*w3/Re;
    
    w4 = W0+deltat*f3;
    p4 = A\w4;
    f4 = -(-By*p4).*(Bx*w4)-(Bx*p4).*(By*w4)+A*w4/Re;
    
    W1 = W0 + deltat*(f1+2*f2+2*f3+f4)/6;
    P1 = A\W1;
    P1 = Thomsnz(N,P1,W1,deltax,deltay,ub);
    U1 = -By*P1;
    V1 = Bx*P1;
    W0 = W1;
    P0 = P1;
    U0 = U1;
    V0 = V1;
end

figure(1)
subplot(1,2,1)
contourf(x,y,reshape(Ws,N,N))
hold on
grid on
colorbar
xlabel({'x'},'FontSize',20,'Interpreter','latex');
xlim([0,1])
ylabel({'y'},'FontSize',20,'Interpreter','latex');
ylim([0,1])
title(strcat('\omega, T=',num2str(0),'s'))
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
subplot(1,2,2)
contourf(x,y,reshape(Ps,N,N))
hold on
quiver(x,y,reshape(Us,N,N),reshape(Vs,N,N),'-r')
grid on
colorbar
xlabel({'x'},'FontSize',20,'Interpreter','latex');
xlim([0,1])
ylabel({'y'},'FontSize',20,'Interpreter','latex');
ylim([0,1])
title(strcat('\psi, T=',num2str(0),'s'))
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[100 100 1200 600])

figure(2)
subplot(1,2,1)
contourf(x,y,reshape(Ps,N,N))
hold on
quiver(x,y,reshape(Us,N,N),reshape(Vs,N,N),'-r')
grid on
colorbar
xlabel({'x'},'FontSize',20,'Interpreter','latex');
xlim([0,1])
ylabel({'y'},'FontSize',20,'Interpreter','latex');
ylim([0,1])
title(strcat('\psi, T=',num2str(0),'s'))
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[100 100 1200 1200])

subplot(1,2,2)
contourf(x,y,reshape(P1,N,N))
hold on
quiver(x,y,reshape(U1,N,N),reshape(V1,N,N),'-r')
grid on
colorbar
xlabel({'x'},'FontSize',20,'Interpreter','latex');
xlim([0,1])
ylabel({'y'},'FontSize',20,'Interpreter','latex');
ylim([0,1])
title(strcat('\psi, T=',num2str(T),'s'))
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[100 100 1200 600])

a = min(P1);
b = max(P1);
level1 = -exp(linspace(log(-a),log(0.00001),10));
level2 = exp(linspace(log(0.00001),log(b)-log(10),10));
level3 = linspace(0.1*b,b,10);
level = [level1,level2,level3];

figure(3)
subplot(1,2,1)
contourf(x,y,reshape(W1,N,N))
hold on
grid on
colorbar
xlabel({'x'},'FontSize',20,'Interpreter','latex');
xlim([0,1])
ylabel({'y'},'FontSize',20,'Interpreter','latex');
ylim([0,1])
title(strcat('\omega, T=',num2str(T),'s'))
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
subplot(1,2,2)
contourf(x,y,reshape(P1,N,N),level)
hold on
grid on
colorbar
xlabel({'x'},'FontSize',20,'Interpreter','latex');
xlim([0,1])
ylabel({'y'},'FontSize',20,'Interpreter','latex');
ylim([0,1])
title(strcat('\psi, T=',num2str(T),'s'))
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[100 100 1200 600])