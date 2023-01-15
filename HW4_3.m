clear
clc
close all
% %% 3(b)
% L = 1;%length
% H = 1;%hight
% T = 0.1;%time
% deltax = 0.05;%space step x
% deltay = 0.05;%space step y
% deltat = 0.001;%time step
% N = L/deltax+1;
% M = H/deltay+1;
% K = T/deltat+1;
% 
% A = zeros(N*M,N*M);
% B = zeros(N*M,N*M);
% I = eye(N*M);
% for j = 1:1:M
%     y(j) = deltay*(j-1);
%     for i = 1:1:N
%         x(i) = deltax*(i-1);
%         A(i+(j-1)*N,i+(j-1)*N) = -2;
%         B(i+(j-1)*N,i+(j-1)*N) = -2;
%         if i == 1
%             A(i+(j-1)*N,i+(j-1)*N+1) = 1;
%         else if i == N
%                 A(i+(j-1)*N,i+(j-1)*N-1) = 1;
%             else
%                 A(i+(j-1)*N,i+(j-1)*N+1) = 1;
%                 A(i+(j-1)*N,i+(j-1)*N-1) = 1;
%             end
%         end
%         if j == 1
%             B(i+(j-1)*N,i+j*N) = 1;
%         else if j == N
%                 B(i+(j-1)*N,i+(j-2)*N) = 1;
%             else
%                 B(i+(j-1)*N,i+j*N) = 1;
%                 B(i+(j-1)*N,i+(j-2)*N) = 1;
%             end
%         end
%         %%initial condition
%         U(i+(j-1)*N,1) = sin(2*pi*x(i))*sin(2*pi*y(j));
%     end
% end
% for i = 1:1:N
%     A(i,:) = 0;
%     A(i+(M-1)*N,:) = 0;
%     B(i,:) = 0;
%     B(i+(M-1)*N,:) = 0;
% end
% for j = 1:1:M
%     A(1+(j-1)*N,:) = 0;
%     A(N+(j-1)*N,:) = 0;
%     B(1+(j-1)*N,:) = 0;
%     B(N+(j-1)*N,:) = 0;
% end
% A = A*deltat/deltax^2/2;
% B = B*deltat/deltay^2/2;
% 
% t(1) = 0;
% for j = 2:1:K
%     t(j) = (j-1)*deltat;
%     Us = (I-A)\(I+B)*U(:,j-1);
%     U(:,j) = (I-B)\(Us+A*U(:,j-1));
% end
% 
% figure(7)
% [X,Y] = meshgrid(x,y);
% Ts = [1,2,4,7,11,16,21,26,31];
% for i = 1:1:9
%     subplot(3,3,i)
%     contourf(X,Y,reshape(U(:,Ts(i)),N,M)',20,'LineWidth',0)
%     colorbar
%     caxis([-1 1])
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     ylabel({'y'},'FontSize',20,'Interpreter','latex');
%     axis square
%     title(strcat('t=',num2str(t(Ts(i))),'s'));
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
% end
% set(gcf,'position',[500 0 1400 1200])

% %% 3(c)
% L = 1;%length
% H = 1;%hight
% T = 0.005;%time
% figure(8)
% hs = [0.2,0.1,0.05,0.02,0.01];
% for h = 1:1:size(hs,2)
%     deltax = hs(h);%space step x
%     deltay = hs(h);%space step y
%     deltat = 0.001;%time step
%     N = L/deltax+1;
%     M = H/deltay+1;
%     K = T/deltat+1;
%     
%     A = zeros(N*M,N*M);
%     B = zeros(N*M,N*M);
%     I = eye(N*M);
%     for j = 1:1:M
%         y(j) = deltay*(j-1);
%         for i = 1:1:N
%             x(i) = deltax*(i-1);
%             A(i+(j-1)*N,i+(j-1)*N) = -2;
%             B(i+(j-1)*N,i+(j-1)*N) = -2;
%             if i == 1
%                 A(i+(j-1)*N,i+(j-1)*N+1) = 1;
%             else if i == N
%                     A(i+(j-1)*N,i+(j-1)*N-1) = 1;
%                 else
%                     A(i+(j-1)*N,i+(j-1)*N+1) = 1;
%                     A(i+(j-1)*N,i+(j-1)*N-1) = 1;
%                 end
%             end
%             if j == 1
%                 B(i+(j-1)*N,i+j*N) = 1;
%             else if j == N
%                     B(i+(j-1)*N,i+(j-2)*N) = 1;
%                 else
%                     B(i+(j-1)*N,i+j*N) = 1;
%                     B(i+(j-1)*N,i+(j-2)*N) = 1;
%                 end
%             end
%             %%initial condition
%             U(i+(j-1)*N,1) = sin(2*pi*x(i))*sin(2*pi*y(j));
%         end
%     end
%     for i = 1:1:N
%         A(i,:) = 0;
%         A(i+(M-1)*N,:) = 0;
%         B(i,:) = 0;
%         B(i+(M-1)*N,:) = 0;
%     end
%     for j = 1:1:M
%         A(1+(j-1)*N,:) = 0;
%         A(N+(j-1)*N,:) = 0;
%         B(1+(j-1)*N,:) = 0;
%         B(N+(j-1)*N,:) = 0;
%     end
%     A = A*deltat/deltax^2/2;
%     B = B*deltat/deltay^2/2;
%     
%     t(1) = 0;
%     for j = 2:1:K
%         t(j) = (j-1)*deltat;
%         Us = (I-A)\(I+B)*U(:,j-1);
%         U(:,j) = (I-B)\(Us+A*U(:,j-1));
%     end
%     
%     [X,Y] = meshgrid(x,y);
%     subplot(2,3,h)
%     contourf(X,Y,reshape(U(:,end),N,M)',20,'LineWidth',0)
%     colorbar
%     caxis([-1 1])
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     ylabel({'y'},'FontSize',20,'Interpreter','latex');
%     axis square
%     title(strcat('h=',num2str(hs(h))));
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
% end
% set(gcf,'position',[500 0 1500 1000])

% Anm = 4*4*(cos(m*pi)-1)*(cos(n*pi)-1)/(m^2-4)/(n^2-4)/(pi^2);
% u(i+(j-1)*N,k) = u(i+(j-1)*N,k)+Anm*cos(n*pi*x0(i))*cos(m*pi*y0(j))*exp(-((n*pi)^2+(m*pi)^2)*t0(k));

%% 3(d)
%%excat solution
L = 1;%length
H = 1;%hight
T = 0.02;%time
deltax = 0.05;%space step x
deltay = 0.05;%space step y
deltat = 0.001;%time step
N = L/deltax+1;
M = H/deltay+1;
K = T/deltat+1;
for k = 1:1:K
    t0(k) = (k-1)*deltat;
    for j = 1:1:M
        y0(j) = deltay*(j-1);
        for i = 1:1:N
            x0(i) = deltax*(i-1);
            u(i+(j-1)*N,k) = sin(2*pi*x0(i))*sin(2*pi*y0(j))*exp(-((2*pi)^2+(2*pi)^2)*t0(k));
        end
    end
end

[X,Y] = meshgrid(x0,y0);

figure(8)
for i = 1:1:9
    subplot(3,3,i)
    Z = reshape(u(:,2*i-1),N,M)';
    contourf(X,Y,Z,20,'LineWidth',1,'EdgeColor','none')
    colorbar
    hold on
    grid on
    xlabel({'x'},'FontSize',20,'Interpreter','latex');
    ylabel({'y'},'FontSize',20,'Interpreter','latex');
    axis square
    title(strcat('Exact,t=',num2str(t0(2*i-1))));
    set(gca, 'FontName','Times New Roman','FontSize', 20);
end
set(gcf,'position',[500 0 1000 1000])

%%numerical solution
L = 1;%length
H = 1;%hight
T = 0.1;%time
figure(9)
ts = [0.1,0.05,0.02,0.01];
hs = [0.1,0.05,0.02,0.01];
for h = 1:1:size(hs,2)
    deltax = hs(h);%space step x
    deltay = hs(h);%space step y
    deltat = ts(h);%time step
    N = L/deltax+1;
    M = H/deltay+1;
    K = T/deltat+1;
    
    A = zeros(N*M,N*M);
    B = zeros(N*M,N*M);
    I = eye(N*M);
    for j = 1:1:M
        y(j) = deltay*(j-1);
        for i = 1:1:N
            x(i) = deltax*(i-1);
            A(i+(j-1)*N,i+(j-1)*N) = -2;
            B(i+(j-1)*N,i+(j-1)*N) = -2;
            if i == 1
                A(i+(j-1)*N,i+(j-1)*N+1) = 1;
            else if i == N
                    A(i+(j-1)*N,i+(j-1)*N-1) = 1;
                else
                    A(i+(j-1)*N,i+(j-1)*N+1) = 1;
                    A(i+(j-1)*N,i+(j-1)*N-1) = 1;
                end
            end
            if j == 1
                B(i+(j-1)*N,i+j*N) = 1;
            else if j == N
                    B(i+(j-1)*N,i+(j-2)*N) = 1;
                else
                    B(i+(j-1)*N,i+j*N) = 1;
                    B(i+(j-1)*N,i+(j-2)*N) = 1;
                end
            end
            %%initial condition
            U0(i+(j-1)*N,1) = sin(2*pi*x(i))*sin(2*pi*y(j));
        end
    end
    for i = 1:1:N
        A(i,:) = 0;
        A(i+(M-1)*N,:) = 0;
        B(i,:) = 0;
        B(i+(M-1)*N,:) = 0;
    end
    for j = 1:1:M
        A(1+(j-1)*N,:) = 0;
        A(N+(j-1)*N,:) = 0;
        B(1+(j-1)*N,:) = 0;
        B(N+(j-1)*N,:) = 0;
    end
    A = A*deltat/deltax^2/2;
    B = B*deltat/deltay^2/2;
    t0 = (K-1)*deltat;
    for j = 1:1:M
        for i = 1:1:N
            u0(i+(j-1)*N,1) = 0;
            u0(i+(j-1)*N,1) = sin(2*pi*x(i))*sin(2*pi*y(j))*exp(-((2*pi)^2+(2*pi)^2)*t0);
        end
    end
    t(1) = 0;
    for j = 2:1:K
        t(j) = (j-1)*deltat;
        Us = (I-A)\(I+B)*U0;
        U = (I-B)\(Us+A*Us);
        U0=U;
    end
    
    [X,Y] = meshgrid(x,y);
    subplot(2,size(hs,2),h)
    Z = reshape(U,N,M)';
    contourf(X,Y,Z,20,'LineWidth',0)
    colorbar
    hold on
    grid on
    xlabel({'x'},'FontSize',20,'Interpreter','latex');
    ylabel({'y'},'FontSize',20,'Interpreter','latex');
    axis square
    title(strcat('h=',num2str(hs(h))));
    set(gca, 'FontName','Times New Roman','FontSize', 20);
    
    E = abs(U-u0);
    err(h) = norm(E,'inf');
    
    subplot(2,size(hs,2),h+size(hs,2))
    Z = reshape(E,N,M)';
    contourf(X,Y,Z,20,'LineWidth',0)
    colorbar
    hold on
    grid on
    xlabel({'x'},'FontSize',20,'Interpreter','latex');
    ylabel({'y'},'FontSize',20,'Interpreter','latex');
    axis square
    title(strcat('Error'));
    set(gca, 'FontName','Times New Roman','FontSize', 20);
end
set(gcf,'position',[500 0 1500 1000])

figure(6)
plot(ts.^2,err,'LineWidth',1)
hold on;
xlabel({'k^{2}'},'FontSize',20,'Interpreter','latex');
set(gca, 'XScale', 'log')
ylabel({'error'},'FontSize',20,'Interpreter','latex');
set(gca, 'YScale', 'log')
axis square
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[1000 100 1500 600])

X = zeros(size(hs,2),2);
X(:,1) = 1;
X(:,2) = log(ts.^2)';

b = log(err)';

k = X\b;

K = k(1);
p = k(2);