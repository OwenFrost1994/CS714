clear
clc
close all

% %%2(b,c)
% h = 0.02;
% n = 1/h+1;
% h0 = 0.001;
% n0 = 1/h0+1;
% for j=1:1:n0
%     for i=1:1:n0
%     x0(i,j) = i*h0-h0;
%     y0(i,j) = j*h0-h0;
%     u0(i,j) = sin(pi*(x0(i,j)-1))*sin(2*pi*(y0(i,j)-1));
%     end
% end
% 
% A = zeros(n^2,n^2);
% f = zeros(n^2,1);
% for j=1:1:n
%     for i=1:1:n
%         x(i,j) = i*h-h;
%         y(i,j) = j*h-h;
%         if i<n && i>1 && j<n && j>1
%             A(i+(j-1)*n,i+(j-1)*n) = -20;
%             A(i+(j-1)*n,i+1+(j-1)*n) = 4;
%             A(i+(j-1)*n,i-1+(j-1)*n) = 4;
%             A(i+(j-1)*n,i+(j)*n) = 4;
%             A(i+(j-1)*n,i+(j-2)*n) = 4; 
%             A(i+(j-1)*n,i-1+(j)*n) = 1;
%             A(i+(j-1)*n,i-1+(j-2)*n) = 1;
%             A(i+(j-1)*n,i+1+(j)*n) = 1;
%             A(i+(j-1)*n,i+1+(j-2)*n) = 1;
%         else
%             A(i+(j-1)*n,i+(j-1)*n) = 1;
%         end
%         if i<n && i>1 && j<n && j>1
%             f(i+(j-1)*n,1) = F(i*h-h,j*h-h)+(F(i*h-2*h,j*h-h)+F(i*h,j*h-h)-4*F(i*h-h,j*h-h)+F(i*h-h,j*h-2*h)+F(i*h-h,j*h))/12;
%         else
%             f(i+(j-1)*n,1) = 0;
%         end
%     end
% end
% A = -A/(6*h^2);
% 
% uk = zeros(n^2,1);
% rk = f-A*uk;
% pk = rk;
% k = 1;
% us(:,1) = uk;
% while norm(rk) > 1E-6
%     err(k) = norm(rk);
%     alphak = rk'*rk/(pk'*A*pk);
%     uk1 = uk + alphak*pk;
%     rk1 = rk - alphak*A*pk;
%     betak = rk1'*rk1/(rk'*rk);
%     pk1 = rk1 + betak*pk;
%     k = k + 1;
%     us(:,k) = uk1;
%     uk = uk1;
%     rk = rk1;
%     pk = pk1;
% end
% err(k) = norm(rk);
% 
% figure(1)
% subplot(1,2,1)
% contourf(x0,y0,u0)
% hold on
% grid on
% xlabel({'x';'(a)'},'FontSize',20,'Interpreter','latex');
% xlim([0,1])
% ylabel({'y'},'FontSize',20,'Interpreter','latex');
% ylim([0,1])
% axis square
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% subplot(1,2,2)
% contourf(x,y,reshape(us(:,k),[n,n]))
% hold on
% grid on
% xlabel({'x';'(b)'},'FontSize',20,'Interpreter','latex');
% xlim([0,1])
% ylabel({'y'},'FontSize',20,'Interpreter','latex');
% ylim([0,1])
% axis square
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[100 100 1600 800])
% 
% figure(2)
% plot(1:1:k,err,'o-','LineWidth',3)
% hold on
% grid on
% xlabel({'k'},'FontSize',20,'Interpreter','latex');
% set(gca, 'YScale', 'log')
% ylabel({'error'},'FontSize',20,'Interpreter','latex');
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[1000 200 1000 800])
% 
% hs = [0.5,0.2,0.1,0.05,0.02,0.01];
% for m = 1:1:size(hs,2)
%     
%     h = hs(m);
%     n = 1/h+1;
%     
%     A = zeros(n^2,n^2);
%     f = zeros(n^2,1);
%     for j=1:1:n
%         for i=1:1:n
%             x(i,j) = i*h-h;
%             y(i,j) = j*h-h;
%             if i<n && i>1 && j<n && j>1
%                 A(i+(j-1)*n,i+(j-1)*n) = -20;
%                 A(i+(j-1)*n,i+1+(j-1)*n) = 4;
%                 A(i+(j-1)*n,i-1+(j-1)*n) = 4;
%                 A(i+(j-1)*n,i+(j)*n) = 4;
%                 A(i+(j-1)*n,i+(j-2)*n) = 4; 
%                 A(i+(j-1)*n,i-1+(j)*n) = 1;
%                 A(i+(j-1)*n,i-1+(j-2)*n) = 1;
%                 A(i+(j-1)*n,i+1+(j)*n) = 1;
%                 A(i+(j-1)*n,i+1+(j-2)*n) = 1;
%             else
%                 A(i+(j-1)*n,i+(j-1)*n) = 1;
%             end
%             if i<n && i>1 && j<n && j>1
%                 f(i+(j-1)*n,1) = F(i*h-h,j*h-h)+(F(i*h-2*h,j*h-h)+F(i*h,j*h-h)-4*F(i*h-h,j*h-h)+F(i*h-h,j*h-2*h)+F(i*h-h,j*h))/12;
%             else
%                 f(i+(j-1)*n,1) = 0;
%             end
%         end
%     end
%     A = -A/(6*h^2);
%     
%     tic;
%     uk = zeros(n^2,1);
%     rk = f-A*uk;
%     pk = rk;
%     k = 1;
% %     us(:,1) = uk;
%     while norm(rk) > 1E-6
%         err(k) = norm(rk);
%         alphak = rk'*rk/(pk'*A*pk);
%         uk1 = uk + alphak*pk;
%         rk1 = rk - alphak*A*pk;
%         betak = rk1'*rk1/(rk'*rk);
%         pk1 = rk1 + betak*pk;
%         k = k + 1;
% %         us(:,k) = uk1;
%         uk = uk1;
%         rk = rk1;
%         pk = pk1;
%     end
%     err(k) = norm(rk);
%     ks(m) = k;
%     toc
%     ts(m) = toc;
% end
% 
% figure(3)
% plot(hs,ks,'o-','LineWidth',3)
% hold on
% grid on
% xlabel({'h'},'FontSize',20,'Interpreter','latex');
% ylabel({'k'},'FontSize',20,'Interpreter','latex');
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[1000 200 1000 800])
% figure(4)
% plot(hs,ts,'o-','LineWidth',3)
% hold on
% grid on
% xlabel({'h'},'FontSize',20,'Interpreter','latex');
% set(gca, 'XScale', 'log')
% ylabel({'$\Delta t$'},'FontSize',20,'Interpreter','latex');
% set(gca, 'YScale', 'log')
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[1000 200 1000 800])

h0 = 0.0001;
n0 = (10-0.1)/h0+1;
for i=1:1:n0
    r0(i) = i*h0-h0+0.1;
    u0(i,1) = -0.2*r0(i)^-1+25.02;
end
% hs = [0.1,0.05,0.01,0.001];
% hs = [0.001,0.0005];
hs = [0.0004];
figure(1)
for m = 1:1:size(hs,2)
    h = hs(m);
    n = floor((10-0.1)/h)+1;
    
    A = zeros(n,n);
    A1 = zeros(n,n);
    A2 = zeros(n,n);
    f = zeros(n,1);
    e = ones(n,1);
    A1 = spdiags([e -2*e e],-1:1,n,n);
    A1(1,:) = 0;
    A1(n,:) = 0;
    for i=1:1:n
        r(i) = i*h-h+0.1;
        if i>1 && i<n
            A2(i,i-1) = -1/r(i);
            A2(i,i+1) = 1/r(i);
        end
        U(i,1) = -0.2*r(i)^-1+25.02;
    end
    A2(1,1) = -3/2;
    A2(1,2) = 2;
    A2(1,3) = -1/2;
    A2(n,n) = h;
    A2 = A2;
    A = A1/h^2+A2/h;
    f(1,1) = 20;
    f(n,1) = 25;
    u = A\f;
    subplot(2,2,m)
    plot(r,u,'bo','LineWidth',1)
    hold on
    plot(r0,u0,'r-','LineWidth',2)
    hold on
    grid on
    xlabel({'x'},'FontSize',20,'Interpreter','latex');
    xlim([0,10])
    ylabel({'u'},'FontSize',20,'Interpreter','latex');
    legend(strcat('u(t) at h=',num2str(hs(m))),'Exact solution u(x)','Location','northeast')
    set(gca, 'FontName','Times New Roman','FontSize', 20);
    E = u-U;
    err(m) = norm(E,'inf');
end
set(gcf,'position',[1000 200 1000 800])

X = zeros(4,2);
X(:,1) = 1;
X(:,2) = log(hs)';

b = log(err)';

k = X\b;

K = k(1);
p = k(2);

figure(2)
plot(hs,err,'o-','LineWidth',3)
hold on
grid on
xlabel({'h'},'FontSize',20,'Interpreter','latex');
set(gca, 'XScale', 'log')
ylabel({'error'},'FontSize',20,'Interpreter','latex');
set(gca, 'YScale', 'log')
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[1000 200 1000 800])