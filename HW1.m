clear
clc
close all

% %%(b)
% Hs = logspace(-2.5,0,500);
% 
% udd = exp(1);
% err = zeros(1,500);
% 
% for i=1:1:500
%     H = Hs(i);
%     for j = 1:1:10000
%         h1 = H*rand(1);
%         h2 = H*rand(1);
%         h3 = H*rand(1);
%         a = 2*(2*h2+h3)/(h1*(h1+h2)*(h1+h2+h3));
%         b = 2*(h1-2*h2-h3)/(h1*h2*(h2+h3));
%         c = 2*(-h1+h2+h3)/(h2*(h1+h2)*h3);
%         d = 2*(h1-h2)/(h3*(h2+h3)*(h1+h2+h3));
%         Udd(i) = a*exp(1-h1)+b*exp(1)+c*exp(1+h2)+d*exp(1+h2+h3);
%         err(i) = err(i)+abs(Udd(i)-udd);
%     end
% end
% err = err/10000;
% 
% figure(1)
% % subplot(1,2,1)
% % plot(Hs,Udd,'r-','Markersize',6,'LineWidth',3)
% % hold on
% % plot(Hs,udd*ones(1,500),'b-','Markersize',6,'LineWidth',3)
% % hold on
% % grid on
% % xlabel({'H'},'FontSize',20,'Interpreter','latex');
% % set(gca, 'XScale', 'log')
% % ylabel({'value'},'FontSize',20,'Interpreter','latex');
% % set(gca, 'YScale', 'log')
% % legend('u''at x=1','Approximated u''at x=1','Location','northeast')
% % subplot(1,2,2)
% plot(Hs,err,'b-','Markersize',6,'LineWidth',3)
% hold on
% grid on
% xlabel({'H'},'FontSize',20,'Interpreter','latex');
% set(gca, 'XScale', 'log')
% ylabel({'err'},'FontSize',20,'Interpreter','latex');
% set(gca, 'YScale', 'log')
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[1000 200 1000 800])
% 
% %%(c)
% A = zeros(500,2);
% A(:,1) = 1;
% A(:,2) = log(Hs)';
% 
% b = log(err)';
% 
% k = A'*A\(A'*b);
% 
% K = k(1);
% p = k(2);
% 
% figure(2)
% plot(Hs,err,'b-','Markersize',6,'LineWidth',3)
% hold on
% plot(Hs,exp(p*log(Hs))*exp(K),'r-','Markersize',6,'LineWidth',3)
% hold on
% grid on
% xlabel({'H'},'FontSize',20,'Interpreter','latex');
% set(gca, 'XScale', 'log')
% ylabel({'err'},'FontSize',20,'Interpreter','latex');
% set(gca, 'YScale', 'log')
% legend('Error','Fitted','Location','northeast')
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[1000 200 1000 800])

% %%(c)
% h = [0.05,0.01,0.002,0.001];
% h0 = 0.001;
% n0 = 10/h0+1;
%     for i=1:1:n0
%         x0(i) = i*h0-h0;
%         u0(i,1) = exp(10)*(sin(x0(i))+cos(x0(i)))/(2*cos(10))-exp(x0(i))/2;
%     end
%     figure(1)
% for j=1:1:4
%     n = 10/h(j)+1;
%     e = ones(n,1);
%     A = spdiags([e -2*e+h(j)^2*e e],-1:1,n,n)/h(j)^2;
% %     A(1,1) = -3/h(j)/2-1;
% %     A(1,2) = 2/h(j);
% %     A(1,3) = -1/h(j)/2;
% %     A(n,n) = 3/h(j)/2+1;
% %     A(n,n-1) = -2/h(j);
% %     A(n,n-2) = 1/h(j)/2;    
%     A(1,1)=-(3/2+h(j))*h(j); A(1,2)=2*h(j); A(1,3) = -h(j)/2; 
%     A(n,end)=(3/2+h(j))*h(j); A(n,end-1)=-2*h(j); A(n,end-2) = h(j)/2;
%     for i=1:1:n
%         x(i) = i*h(j)-h(j);
%         u(i,1) = exp(10)*(sin(x(i))+cos(x(i)))/(2*cos(10))-exp(x(i))/2;
%         F(i,1) = -exp(x(i));
%     end
%     F(1,1) = 0;
%     F(n,1) = 0;
%     U = A\F;
%     E = abs(U-u);
%     err(j) = norm(E,'inf');
%     subplot(2,2,j)
%     plot(x,U,'bo','LineWidth',1)
%     hold on
%     plot(x0,u0,'r-','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     ylabel({'u'},'FontSize',20,'Interpreter','latex');
%     legend(strcat('U(x) at h=',num2str(h(j))),'Exact solution u(x)','Location','northeast')
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
% end
% set(gcf,'position',[100 100 1000 800])
% 
% X = zeros(4,2);
% X(:,1) = 1;
% X(:,2) = log(h)';
% 
% b = log(err)';
% 
% k = X\b;
% 
% K = k(1);
% p = k(2);
% 
% 
% figure(2)
% plot(h,err,'-','LineWidth',3)
% hold on
% grid on
% xlabel({'h'},'FontSize',20,'Interpreter','latex');
% ylabel({'error'},'FontSize',20,'Interpreter','latex');
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[1000 200 1000 800])

%%(c)
h = [0.05,0.01,0.002,0.001];
h0 = 0.001;
n0 = 1/h0+1;
for i=1:1:n0
    x0(i) = i*h0-h0;
    u0(i,1) = sin(4*pi*x0(i))/(16*pi^2+1);
end
figure(1)
for j=1:1:4
%     n = 1/h(j)+1;
%     e = ones(n,1);
%     A = spdiags([e/12 -4*e/3 5*e/2 -4*e/3 e/12],-2:2,n,n);
%     A = A/h(j)^2+eye(n);
%     A0 = A/h(j)^2;
%     A(1,n-2) = 1/h(j)^2/12;
%     A(1,n-1) = -4/h(j)^2/3;
%     A(2,n-1) = 1/h(j)^2/12;
%     A(n-1,2) = 1/h(j)^2/12;
%     A(n,:) = 0;
%     A(n,1) = -1;
%     A(n,n) = 1;
%     A0(1,n-2) = 1/h(j)^2/12;
%     A0(1,n-1) = -4/h(j)^2/3;
%     A0(2,n-1) = 1/h(j)^2/12;
%     A0(n-1,2) = 1/h(j)^2/12;
%     A0(n,:) = 0;
%     A0(n,1) = -1;
%     A0(n,n) = 1;
    n = 1/h(j);
    e = ones(n,1);
    A = spdiags([e/12 -4*e/3 5*e/2 -4*e/3 e/12],-2:2,n,n);
    A0 = A/h(j)^2;
    A = A/h(j)^2+eye(n);
    A(1,n-1) = 1/h(j)^2/12;
    A(1,n) = -4/h(j)^2/3;
    A(2,n) = 1/h(j)^2/12;
    A(n,1) = -4/h(j)^2/3;
    A(n,2) = 1/h(j)^2/12;
    A(n-1,1) = 1/h(j)^2/12;
    A0(1,n-1) = 1/h(j)^2/12;
    A0(1,n) = -4/h(j)^2/3;
    A0(2,n) = 1/h(j)^2/12;
    A0(n,1) = -4/h(j)^2/3;
    A0(n,2) = 1/h(j)^2/12;
    A0(n-1,1) = 1/h(j)^2/12;

    [V,D] = eigs(A0);
    eig_min(j)=min(sort(diag(D)));
    [V,D] = eig(A);
    eig_min1(j)=min(sort(diag(D)));
    for i=1:1:n
        x(i) = i*h(j)-h(j);
        u(i,1) = sin(4*pi*x(i))/(16*pi^2+1);
        F(i,1) = sin(4*pi*x(i));
    end
%     F(1,1) = 0;
%     F(n,1) = 0;
    U = A\F;
    E = abs(U-u);
    err(j) = norm(E,'inf');
    subplot(2,2,j)
    plot(x,U,'bo','LineWidth',1)
    hold on
    plot(x0,u0,'r-','LineWidth',2)
    hold on
    grid on
    xlabel({'x'},'FontSize',20,'Interpreter','latex');
    xlim([0,1])
    ylabel({'u'},'FontSize',20,'Interpreter','latex');
    legend(strcat('U(x) at h=',num2str(h(j))),'Exact solution u(x)','Location','northeast')
    set(gca, 'FontName','Times New Roman','FontSize', 20);
    t = abs(A*(U-u));
    tn(j) = norm(t,'inf');
end
set(gcf,'position',[100 100 1000 800])

X = zeros(4,2);
X(:,1) = 1;
X(:,2) = log(h)';

b = log(err)';

k = X\b;

K = k(1);
p = k(2);

X0 = zeros(4,2);
X0(:,1) = 1;
X0(:,2) = log(h)';

b = log(tn)';

k = X0\b;

K = k(1);
p2 = k(2);
figure(2)
plot(h,err,'bo-','LineWidth',3)
hold on
grid on
xlabel({'h'},'FontSize',20,'Interpreter','latex');
ylabel({'error'},'FontSize',20,'Interpreter','latex');
set(gca, 'YScale', 'log')
set(gca, 'FontName','Times New Roman','FontSize', 20);
set(gcf,'position',[1000 200 1000 800])
