clear
clc
close all

% %% 2(b)
% L = 1;%length
% T = 0.2;%time
% deltax = 0.05;%space step
% deltat = 0.001;%time step
% N = L/deltax+1;
% M = T/deltat+1;
% 
% 
% e = ones(N,1);
% A = deltat*spdiags([e -2*e e],-1:1,N,N)/deltax^2/2;
% I = eye(N);
% 
% for i = 1:1:N
%     x(i) = (i-1)*deltax;
%     if x(i)<0.6 && x(i)>0.4
%         U(i,1) = 2;
%     else
%         U(i,1) = 1;
%     end
% end
% t(1) = 0;
% for j = 2:1:M
%     t(j) = (j-1)*deltat;
%     U(:,j) = (I-A)\(I+A)*U(:,j-1);
%     U(1,j) = 0;
%     U(N,j) = 0;
% end
% 
% figure(1)
% [T1,X] = meshgrid(t,x);
% contourf(T1,X,U,10)
% hold on
% grid on
% colorbar
% xlabel({'t'},'FontSize',20,'Interpreter','latex');
% ylabel({'x'},'FontSize',20,'Interpreter','latex');
% axis square
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[100 0 800 800])
% 
% figure(2)
% J = [1,2,3,4,5,6,10,15,20,25,30,35,40,45,50,55];
% for j = 1:1:16
%     subplot(4,4,j)
%     plot(x,U(:,J(j)),'--','LineWidth',1)
%     hold on;
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
%     ylim([0,3])
%     axis square
%     title(strcat('t=',num2str(t(J(j)))));
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
% end
% set(gcf,'position',[1000 200 1000 1000])

% %% 2(c)
% L = 1;%length
% T = 0.5;%time
% Deltax  = [0.1,0.05,0.01,0.005,0.002];
% figure(3)
% for h = 1:1:5
%     deltax = Deltax(h);%space step h
%     deltat = 0.001;%time step k
%     N = L/deltax+1;
%     M = T/deltat+1;
%     
%     e = ones(N,1);
%     A = deltat*spdiags([e -2*e e],-1:1,N,N)/deltax^2/2;
%     I = eye(N);
%     for i = 1:1:N
%         x(i) = (i-1)*deltax;
%         if x(i)<0.6 && x(i)>0.4
%             U(i,1) = 2;
%         else
%             U(i,1) = 1;
%         end
%     end
%     
%     t(1) = 0;
%     for j = 2:1:M
%         t(j) = (j-1)*deltat;
%         U(:,j) = (I-A)\(I+A)*U(:,j-1);
%         U(1,j) = 0;
%         U(N,j) = 0;
%     end
%     
%     J = [1,3,5,7];
%     for j = 1:1:size(J,2)
%         subplot(size(Deltax,2),size(J,2),j+(h-1)*size(J,2))
%         plot(x,U(:,J(j)),'--','LineWidth',1)
%         hold on;
%         xlabel({'x'},'FontSize',20,'Interpreter','latex');
%         ylabel({'U'},'FontSize',20,'Interpreter','latex');
%         ylim([0,3])
%         axis square
%         title(strcat('h=',num2str(Deltax(h)),',','t=',num2str(t(J(j)))));
%         set(gca, 'FontName','Times New Roman','FontSize', 20);
%     end
% end
% set(gcf,'position',[1000 100 1000 1200])

%% 2(d)
% L = 1;%length
% T = 0.2;%time
% 
% deltax = 0.01;
% deltat = 0.002*deltax;
% N = floor(L/deltax)+1;
% M = floor(T/deltat)+1;
% Tr = (M-1)*deltat;
% for j = 1:1:M
%         t0(j) = (j-1)*deltat;
%         for i = 1:1:N
%             x0(i) = (i-1)*deltax;
%             u0(i,j) = 0;
%             for n = 1:1:10
%                 an = (-2*cos(n*pi*0.4)+2*cos(n*pi*0)-2*cos(n*pi*1)+2*cos(n*pi*0.6)-4*cos(n*pi*0.6)+4*cos(n*pi*0.4))/(n*pi);
%                 u0(i,j) = u0(i,j)+an*exp(-n^2*pi^2*t0(j))*sin(n*pi*x0(i));
%             end
%         end
% end
% figure(4)
% [T0,X0] = meshgrid(t0,x0);
% surf(T0,X0,u0,'EdgeColor','none')
% hold on
% grid on
% xlabel({'t'},'FontSize',20,'Interpreter','latex');
% ylabel({'x'},'FontSize',20,'Interpreter','latex');
% axis square
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[100 100 800 800])

L = 1;%length
T = 0.05;%time
r = 0.2;
Deltax  = [0.1,0.05,0.04,0.02,0.01];
figure(5)
for h = 1:1:5
    deltax = Deltax(h);%space step h
    deltat = r*deltax^2;%time step k
    N = floor(L/deltax)+1;
    M = floor(T/deltat)+1;
%     deltat/(deltax^2)
    e = ones(N,1);
    A = deltat*spdiags([e -2*e e],-1:1,N,N)/deltax^2/2;
    I = eye(N);
    for i = 1:1:N
        x(i) = (i-1)*deltax;
        if x(i)<0.6 && x(i)>0.4
            U(i,1) = 2;
        else
            U(i,1) = 1;
        end
    end
    Tr = (M-1)*deltat;
    for i = 1:1:N
        u(i,1) = 0;
        for n = 1:1:100
            an = (-2*cos(n*pi*0.4)+2*cos(n*pi*0)-2*cos(n*pi*1)+2*cos(n*pi*0.6)-4*cos(n*pi*0.6)+4*cos(n*pi*0.4))/(n*pi);
            u(i,1) = u(i,1)+an*exp(-n^2*pi^2*Tr)*sin(n*pi*x(i));
        end
    end
    t(1) = 0;
    for j = 2:1:M
        t(j) = (j-1)*deltat;
        U(:,j) = (I-A)\(I+A)*U(:,j-1);
        U(1,j) = 0;
        U(N,j) = 0;
    end
    E = abs(U(:,end)-u);
    err(h) = norm(E,'inf');
%     figure(5+h)
    subplot(2,size(Deltax,2),h)
%     subplot(1,2,1)
    [T1,X] = meshgrid(t,x);
    contourf(T1,X,U,10)
%     surf(T0,X0,u0,'EdgeColor','none');
%     hold on
%     surf(T1,X,U,'EdgeColor','none');
%     hold on
    colorbar
    grid on
    xlabel({'t'},'FontSize',20,'Interpreter','latex');
    ylabel({'x'},'FontSize',20,'Interpreter','latex');
    axis square
    set(gca, 'FontName','Times New Roman','FontSize', 20);
    
    subplot(2,size(Deltax,2),h+size(Deltax,2))
%     subplot(1,2,2)
    plot(x,U(:,end),'--b','LineWidth',1)
    hold on;
    plot(x,u,'-r','LineWidth',3)
    hold on;
    xlabel({'x'},'FontSize',20,'Interpreter','latex');
    ylabel({'U'},'FontSize',20,'Interpreter','latex');
    axis square
    title(strcat('t=',num2str(Tr)));
    legend(strcat('\Deltax=',num2str(deltax),',\Deltat=',num2str(deltat)),'exact soln')
    set(gca, 'FontName','Times New Roman','FontSize', 20);
%     set(gcf,'position',[1000 100 1000 500])
end
set(gcf,'position',[1000 100 2000 1000])

figure(6)
plot(Deltax,err,'-o','LineWidth',1)
hold on;
xlabel({'h'},'FontSize',20,'Interpreter','latex');
set(gca, 'XScale', 'log')
ylabel({'error'},'FontSize',20,'Interpreter','latex');
set(gca, 'YScale', 'log')
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[1000 100 800 800])

X = zeros(5,2);
X(:,1) = 1;
X(:,2) = log(Deltax)';

b = log(err)';

k = X\b;

K = k(1);
p = k(2);

