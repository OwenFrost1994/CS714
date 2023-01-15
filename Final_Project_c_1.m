clear
clc
close all

L = 1;%length
H = 1;%hight
T = 0.2;%time
deltat = 0.02;%time step
ns = [8,16,32,64];
for s = 1:1:size(ns,2)
    n = ns(s);
    Re = 1000;
    deltax = 1/n;%space step x
    deltay = 1/n;%space step y
    N = n+1;
    K = T/deltat+1;

    Ax = zeros(N*N,N*N);
    Ay = zeros(N*N,N*N);
    Bx = zeros(N*N,N*N);
    By = zeros(N*N,N*N);
    I = eye(N*N);
    w = zeros(N*N,1);
    p = zeros(N*N,1);
    u = zeros(N*N,1);
    v = zeros(N*N,1);
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
            w(i+(j-1)*N,1) = -2*pi^2*exp(-2*pi^2*0/Re)*sin(pi*x(i,j))*sin(pi*y(i,j));
        end
    end
    Ax = Ax/deltax^2;
    Ay = Ay/deltay^2;
    A = Ax+Ay;
    Bx = Bx/deltax/2;
    By = By/deltay/2;
    
    p(:,1) = A\w(:,1);
    
    p(:,1) = Thoms(N,p,w,deltax,deltay);
    u(:,1) = -By*p(:,1);
    v(:,1) = Bx*p(:,1);
    for k = 2:1:K
        f1 = -u(:,k-1).*(Bx*w(:,k-1))-v(:,k-1).*(By*w(:,k-1))+A*w(:,k-1)/Re;

        w2 = w(:,k-1)+deltat*f1/2;
        p2 = A\w2;
        f2 = -(-By*p2).*(Bx*w2)-(Bx*p2).*(By*w2)+A*w2/Re;
        
        w3 = w(:,k-1)+deltat*f2/2;
        p3 = A\w3;
        f3 = -(-By*p3).*(Bx*w3)-(Bx*p3).*(By*w3)+A*w3/Re;

        w4 = w(:,k-1)+deltat*f3;
        p4 = A\w4;
        f4 = -(-By*p4).*(Bx*w4)-(Bx*p4).*(By*w4)+A*w4/Re;

        w(:,k) = w(:,k-1) + deltat*(f1+2*f2+2*f3+f4)/6;
        p(:,k) = A\w(:,k);
        p(:,k) = Thoms(N,p(:,k),w(:,k),deltax,deltay);
        u(:,k) = -By*p(:,k);
        v(:,k) = Bx*p(:,k);
    end
    for j = 1:1:N
        for i = 1:1:N
            w1(i+(j-1)*N,1) = -2*pi^2*exp(-2*pi^2*T/Re)*sin(pi*x(i,j))*sin(pi*y(i,j));
            p1(i+(j-1)*N,1) = exp(-2*pi^2*T/Re)*sin(pi*x(i,j))*sin(pi*y(i,j));
        end
    end
    err(s)=norm(abs(p(:,end)-p1),'inf');
end

figure(1)
subplot(1,2,1)
z = reshape(w(:,1),N,N);
contourf(x,y,z)
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
z = reshape(p(:,1),N,N);
contourf(x,y,z)
hold on
quiver(x,y,reshape(u(:,1),N,N),reshape(v(:,1),N,N),'-r')
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
subplot(2,2,1)
z = reshape(w(:,end),N,N);
contourf(x,y,z)
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
subplot(2,2,2)
z = reshape(p(:,end),N,N);
contourf(x,y,z)
hold on
quiver(x,y,reshape(u(:,end),N,N),reshape(v(:,end),N,N),'-r')
grid on
colorbar
xlabel({'x'},'FontSize',20,'Interpreter','latex');
xlim([0,1])
ylabel({'y'},'FontSize',20,'Interpreter','latex');
ylim([0,1])
title(strcat('\psi, T=',num2str(T),'s'))
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[100 100 1200 1200])
subplot(2,2,3)
z = reshape(w1,N,N);
contourf(x,y,z)
hold on
grid on
colorbar
xlabel({'x'},'FontSize',20,'Interpreter','latex');
xlim([0,1])
ylabel({'y'},'FontSize',20,'Interpreter','latex');
ylim([0,1])
title('\omega, Exact')
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
subplot(2,2,4)
z = reshape(p1,N,N);
contourf(x,y,z)
hold on
quiver(x,y,reshape(u(:,end),N,N),reshape(v(:,end),N,N),'-r')
grid on
colorbar
xlabel({'x'},'FontSize',20,'Interpreter','latex');
xlim([0,1])
ylabel({'y'},'FontSize',20,'Interpreter','latex');
ylim([0,1])
title('\psi, Exact')
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[100 100 1200 1200])

figure(6)
subplot(1,2,1)
plot(ns,err,'-o','LineWidth',2)
hold on;
xlabel({'n'},'FontSize',20,'Interpreter','latex');
ylabel({'error'},'FontSize',20,'Interpreter','latex');
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
subplot(1,2,2)
plot(1./ns,err,'-o','LineWidth',2)
hold on;
xlabel({'h'},'FontSize',20,'Interpreter','latex');
ylabel({'error'},'FontSize',20,'Interpreter','latex');
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[1000 100 1500 600])


X = zeros(size(ns,2),2);
X(:,1) = 1;
X(:,2) = log(1./ns)';

b = log(err)';

k = X\b;

K = k(1);
p = k(2);

Errt = [ns',1./ns',err'];