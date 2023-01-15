clear
clc
close all
% %% 1(b)
% r = [0.05 0.1 0.2 0.4 0.6 0.8 1];
% xh = 0:2*pi/100:2*pi;
% for j = 1:1:size(r,2)
%     for i = 1:1:size(xh,2)
%         gr(j,i)=1-r(j)*(3-4*cos(xh(i))+cos(2*xh(i)))/2+r(j)^2*(1-2*cos(xh(i))+cos(2*xh(i)))/2;
%         gi(j,i)=-r(j)*(4*sin(xh(i))-sin(2*xh(i)))/2+r(j)^2*(2*sin(xh(i))-sin(2*xh(i)))/2;
%     end
% end
% theta = 0:2*pi/100:2*pi;
% figure(1)
% legends = {};
% plot(sin(theta),cos(theta),'-r','LineWidth',3)
% hold on
% legends = [legends 'Unit circle'];
% for j = 1:1:size(r,2)
%     plot(gr(j,:),gi(j,:),'-o','LineWidth',1)
%     hold on
%     legends = [legends strcat('a*k/h=',num2str(r(j)))];
% end
% grid on
% xlabel({'Re(g)'},'FontSize',20,'Interpreter','latex');
% xlim([-1,1])
% ylabel({'Im(g)'},'FontSize',20,'Interpreter','latex');
% ylim([-1,1])
% axis square
% legend(legends)
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[100 100 800 800])
% %% 2(c)
% r = [0.05 0.1 0.2 0.4 0.6 0.8 1];
% xh = 0:2*pi/100:2*pi;
% for j = 1:1:size(r,2)
%     for i = 1:1:size(xh,2)
%         gr(j,i)=1-r(j)/2-r(j)^2+r(j)^3/2+(-r(j)/3+r(j)^2/2-r(j)^3/6)*cos(xh(i))+(r(j)+r(j)^2/2-r(j)^3/2)*cos(xh(i))...
%             +(-r(j)/6+r(j)^3/6)*cos(2*xh(i));
%         gi(j,i)=(-r(j)/3+r(j)^2/2-r(j)^3/6)*sin(xh(i))-(r(j)+r(j)^2/2-r(j)^3/2)*sin(xh(i))-(-r(j)/6+r(j)^3/6)*sin(2*xh(i));
%     end
% end
% theta = 0:2*pi/100:2*pi;
% figure(2)
% legends = {};
% plot(sin(theta),cos(theta),'-r','LineWidth',3)
% hold on
% legends = [legends 'Unit circle'];
% for j = 1:1:size(r,2)
%     plot(gr(j,:),gi(j,:),'-o','LineWidth',1)
%     hold on
%     legends = [legends strcat('a*k/h=',num2str(r(j)))];
% end
% grid on
% xlabel({'Re(g)'},'FontSize',20,'Interpreter','latex');
% xlim([-1,1])
% ylabel({'Im(g)'},'FontSize',20,'Interpreter','latex');
% ylim([-1,1])
% axis square
% legend(legends)
% set(gca, 'FontName','Times New Roman','FontSize', 20);
% set(gcf,'position',[100 100 800 800])

% %% 3(a)
% L = 1;%length
% T = 1;%time
% a = 1;
% deltax = 0.05;%space step x
% deltat = 0.01;%time step
% N = L/deltax+1;
% K = T/deltat+1;
% 
% r = a*deltat/deltax;
% 
% A = zeros(N,N);
% B = zeros(N,N);
% I = eye(N);
% for j = 1:1:N
%     x(j) = deltax*(j-1);
%     %%initial condition
%     U1(j,1) = 2*exp(-(100*x(j)-50)^2/10);
%     U2(j,1) = exp(-100*(x(j)-0.5)^2)*sin(80*pi*x(j));
%     U3(j,1) = 2*exp(-(100*x(j)-50)^2/10);
%     U4(j,1) = exp(-100*(x(j)-0.5)^2)*sin(80*pi*x(j));
% end
% 
% t(1) = 0;
% for n = 2:1:K
%     t(n) = (n-1)*deltat;
%     for j = 1:1:N
%         if j == 1
%             U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(end,n-1)+U1(end-1,n-1))/2+...
%                 +r^2*(U1(j,n-1)-2*U1(end,n-1)+U1(end-1,n-1))/2;
%             U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(end,n-1)+U2(end-1,n-1))/2+...
%                 +r^2*(U2(j,n-1)-2*U2(end,n-1)+U2(end-1,n-1))/2;
%             U3(j,n) = (-r/6+r^3/6)*U3(end-1,n-1)+(r+r^2/2-r^3/2)*U3(end,n-1)+...
%                 (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(j+1,n-1);
%             U4(j,n) = (-r/6+r^3/6)*U4(end-1,n-1)+(r+r^2/2-r^3/2)*U4(end,n-1)+...
%                 (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(j+1,n-1);
%         else if j == 2
%                 U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(j-1,n-1)+U1(end,n-1))/2+...
%                     +r^2*(U1(j,n-1)-2*U1(j-1,n-1)+U1(end,n-1))/2;
%                 U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(j-1,n-1)+U2(end,n-1))/2+...
%                     +r^2*(U2(j,n-1)-2*U2(j-1,n-1)+U2(end,n-1))/2;
%                 U3(j,n) = (-r/6+r^3/6)*U3(end,n-1)+(r+r^2/2-r^3/2)*U3(j-1,n-1)+...
%                     (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(j+1,n-1);
%                 U4(j,n) = (-r/6+r^3/6)*U4(end,n-1)+(r+r^2/2-r^3/2)*U4(j-1,n-1)+...
%                     (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(j+1,n-1);
%             else
%                 if j == N
%                     U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(j-1,n-1)+U1(j-2,n-1))/2+...
%                         +r^2*(U1(j,n-1)-2*U1(j-1,n-1)+U1(j-2,n-1))/2;
%                     U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(j-1,n-1)+U2(j-2,n-1))/2+...
%                         +r^2*(U2(j,n-1)-2*U2(j-1,n-1)+U2(j-2,n-1))/2;
%                     U3(j,n) = (-r/6+r^3/6)*U3(j-2,n-1)+(r+r^2/2-r^3/2)*U3(j-1,n-1)+...
%                         (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(1,n-1);
%                     U4(j,n) = (-r/6+r^3/6)*U4(j-2,n-1)+(r+r^2/2-r^3/2)*U4(j-1,n-1)+...
%                         (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(1,n-1);
%                 else
%                     U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(j-1,n-1)+U1(j-2,n-1))/2+...
%                         +r^2*(U1(j,n-1)-2*U1(j-1,n-1)+U1(j-2,n-1))/2;
%                     U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(j-1,n-1)+U2(j-2,n-1))/2+...
%                         +r^2*(U2(j,n-1)-2*U2(j-1,n-1)+U2(j-2,n-1))/2;
%                     U3(j,n) = (-r/6+r^3/6)*U3(j-2,n-1)+(r+r^2/2-r^3/2)*U3(j-1,n-1)+...
%                         (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(j+1,n-1);
%                     U4(j,n) = (-r/6+r^3/6)*U4(j-2,n-1)+(r+r^2/2-r^3/2)*U4(j-1,n-1)+...
%                         (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(j+1,n-1);
%                 end
%             end
%         end
%     end
% end
% 
% figure(1)
% Ts = [1,11,21,31,41,51,61,71,81];
% for j = 1:1:9
%     subplot(3,3,j)
%     plot(x,U1(:,Ts(j)),'r-o','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     xlim([0,1])
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
%     ylim([-0.5,2])
%     axis square
%     title(strcat('t=',num2str(t(Ts(j))),'s'));
%     legend('IC1')
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
% end
% set(gcf,'position',[500 0 1200 1200])
% figure(2)
% Ts = [1,11,21,31,41,51,61,71,81];
% for j = 1:1:9
%     subplot(3,3,j)
%     plot(x,U2(:,Ts(j)),'b-+','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     xlim([0,1])
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
%     ylim([-5E-15,10E-15])
%     axis square
%     title(strcat('t=',num2str(t(Ts(j))),'s'));
%     legend('IC2')
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
% end
% set(gcf,'position',[500 0 1200 1200])
% figure(3)
% Ts = [1,11,21,31,41,51,61,71,81];
% for j = 1:1:9
%     subplot(3,3,j)
%     plot(x,U3(:,Ts(j)),'r-o','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     xlim([0,1])
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
%     ylim([-0.5,2])
%     axis square
%     title(strcat('t=',num2str(t(Ts(j))),'s'));
%     legend('IC1')
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
% end
% set(gcf,'position',[500 0 1200 1200])
% figure(4)
% Ts = [1,11,21,31,41,51,61,71,81];
% for j = 1:1:9
%     subplot(3,3,j)
%     plot(x,U4(:,Ts(j)),'b-+','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     xlim([0,1])
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
%     ylim([-5E-15,10E-15])
%     axis square
%     title(strcat('t=',num2str(t(Ts(j))),'s'));
%     legend('IC2')
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
% end
% set(gcf,'position',[500 0 1200 1200])

% %% 3(b)
% L = 1;%length
% T = 0.2;%time
% a = 1;
% hs = [0.05,0.02,0.01,0.005,0.002];
% ts = [0.02,0.01,0.005,0.002,0.001];
% figure(6)
% for h = 1:1:size(hs,2)
%     deltax = hs(h);%space step x
%     deltat = ts(h);%time step0.01
%     N = L/deltax+1;
%     K = T/deltat+1;
%     
%     r = a*deltat/deltax;
%     for j = 1:1:N
%         x(j) = deltax*(j-1);
%         %%initial condition
%         U1(j,1) = 2*exp(-(100*x(j)-50)^2/10);
%         U2(j,1) = exp(-100*(x(j)-0.5)^2)*sin(80*pi*x(j));
%         U3(j,1) = 2*exp(-(100*x(j)-50)^2/10);
%         U4(j,1) = exp(-100*(x(j)-0.5)^2)*sin(80*pi*x(j));
%     end
%     
%     t(1) = 0;
%     for n = 2:1:K
%         t(n) = (n-1)*deltat;
%         for j = 1:1:N
%             if j == 1
%                 U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(end,n-1)+U1(end-1,n-1))/2+...
%                     +r^2*(U1(j,n-1)-2*U1(end,n-1)+U1(end-1,n-1))/2;
%                 U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(end,n-1)+U2(end-1,n-1))/2+...
%                     +r^2*(U2(j,n-1)-2*U2(end,n-1)+U2(end-1,n-1))/2;
%                 U3(j,n) = (-r/6+r^3/6)*U3(end-1,n-1)+(r+r^2/2-r^3/2)*U3(end,n-1)+...
%                     (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(j+1,n-1);
%                 U4(j,n) = (-r/6+r^3/6)*U4(end-1,n-1)+(r+r^2/2-r^3/2)*U4(end,n-1)+...
%                     (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(j+1,n-1);
%             else if j == 2
%                     U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(j-1,n-1)+U1(end,n-1))/2+...
%                         +r^2*(U1(j,n-1)-2*U1(j-1,n-1)+U1(end,n-1))/2;
%                     U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(j-1,n-1)+U2(end,n-1))/2+...
%                         +r^2*(U2(j,n-1)-2*U2(j-1,n-1)+U2(end,n-1))/2;
%                     U3(j,n) = (-r/6+r^3/6)*U3(end,n-1)+(r+r^2/2-r^3/2)*U3(j-1,n-1)+...
%                         (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(j+1,n-1);
%                     U4(j,n) = (-r/6+r^3/6)*U4(end,n-1)+(r+r^2/2-r^3/2)*U4(j-1,n-1)+...
%                         (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(j+1,n-1);
%                 else
%                     if j == N
%                         U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(j-1,n-1)+U1(j-2,n-1))/2+...
%                             +r^2*(U1(j,n-1)-2*U1(j-1,n-1)+U1(j-2,n-1))/2;
%                         U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(j-1,n-1)+U2(j-2,n-1))/2+...
%                             +r^2*(U2(j,n-1)-2*U2(j-1,n-1)+U2(j-2,n-1))/2;
%                         U3(j,n) = (-r/6+r^3/6)*U3(j-2,n-1)+(r+r^2/2-r^3/2)*U3(j-1,n-1)+...
%                             (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(1,n-1);
%                         U4(j,n) = (-r/6+r^3/6)*U4(j-2,n-1)+(r+r^2/2-r^3/2)*U4(j-1,n-1)+...
%                             (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(1,n-1);
%                     else
%                         U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(j-1,n-1)+U1(j-2,n-1))/2+...
%                             +r^2*(U1(j,n-1)-2*U1(j-1,n-1)+U1(j-2,n-1))/2;
%                         U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(j-1,n-1)+U2(j-2,n-1))/2+...
%                             +r^2*(U2(j,n-1)-2*U2(j-1,n-1)+U2(j-2,n-1))/2;
%                         U3(j,n) = (-r/6+r^3/6)*U3(j-2,n-1)+(r+r^2/2-r^3/2)*U3(j-1,n-1)+...
%                             (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(j+1,n-1);
%                         U4(j,n) = (-r/6+r^3/6)*U4(j-2,n-1)+(r+r^2/2-r^3/2)*U4(j-1,n-1)+...
%                             (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(j+1,n-1);
%                     end
%                 end
%             end
%         end
%     end
%     
%     subplot(4,size(hs,2),h)
%     plot(x,U1(:,end),'r-o','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     xlim([0,1])
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
%     ylim([-0.5,2])
%     axis square
%     if h == 3
%         title(strcat('Beam-Warming IC1,','t=',num2str(T),'s'));
%     end
%     legend(strcat('h=',num2str(hs(h))))
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
%     
%     subplot(4,size(hs,2),h+size(hs,2))
%     plot(x,U2(:,end),'r-+','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     xlim([0,1])
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
% %     ylim([-5E-15,10E-15])
%     axis square
%     if h == 3
%         title(strcat('Beam-Warming IC1,','t=',num2str(T),'s'));
%     end
%     legend(strcat('h=',num2str(hs(h))))
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
%     
%     subplot(4,size(hs,2),h+2*size(hs,2))
%     plot(x,U3(:,end),'b-o','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     xlim([0,1])
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
%     ylim([-0.5,2])
%     axis square
%     if h == 3
%         title(strcat('Third-Order Lax-Wendroff IC1,','t=',num2str(T),'s'));
%     end
%     legend(strcat('h=',num2str(hs(h))))
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
%     
%     subplot(4,size(hs,2),h+3*size(hs,2))
%     plot(x,U4(:,end),'b-+','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     xlim([0,1])
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
% %     ylim([-5E-15,10E-15])
%     axis square
%     if h == 3
%         title(strcat('Third-Order Lax-Wendroff IC1,','t=',num2str(T),'s'));
%     end
%     legend(strcat('h=',num2str(hs(h))))
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
% end
% set(gcf,'position',[500 0 2000 1500])
% % figure(3)
% % Ts = [1,2,3,4,5,10,15,18,21];
% % for j = 1:1:9
% %     subplot(3,3,j)
% %     plot(x,U1(:,Ts(j)),'r-o','LineWidth',2)
% %     hold on
% %     grid on
% %     xlabel({'x'},'FontSize',20,'Interpreter','latex');
% %     xlim([0,1])
% %     ylabel({'U'},'FontSize',20,'Interpreter','latex');
% %     ylim([-0.5,2])
% %     axis square
% %     title(strcat('t=',num2str(t(Ts(j))),'s'));
% %     legend('IC1')
% %     set(gca, 'FontName','Times New Roman','FontSize', 20);
% % end
% % set(gcf,'position',[500 0 1200 1200])
% figure(5)
% for h = 1:1:size(hs,2)
%     deltax = hs(h);%space step x
%     N = L/deltax+1;
%     
%     for j = 1:1:N
%         x0(j) = deltax*(j-1);
%         %%initial condition
%         U10(j,1) = 2*exp(-(100*x0(j)-50)^2/10);
%         U20(j,1) = exp(-100*(x0(j)-0.5)^2)*sin(80*pi*x0(j));
%     end
%     subplot(2,size(hs,2),h)
%     plot(x0,U10(:,1),'r-o','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     xlim([0,1])
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
%     ylim([-0.5,2])
%     axis square
%     if h == 3
%         title(strcat('IC1,','t=',num2str(0),'s'));
%     end
%     legend(strcat('h=',num2str(hs(h))))
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
%     
%     subplot(2,size(hs,2),h+size(hs,2))
%     plot(x0,U20(:,1),'b-+','LineWidth',2)
%     hold on
%     grid on
%     xlabel({'x'},'FontSize',20,'Interpreter','latex');
%     xlim([0,1])
%     ylabel({'U'},'FontSize',20,'Interpreter','latex');
% %     ylim([-0.5,2])
%     axis square
%     if h == 3
%         title(strcat('IC2,','t=',num2str(0),'s'));
%     end
%     legend(strcat('h=',num2str(hs(h))))
%     set(gca, 'FontName','Times New Roman','FontSize', 20);
% end
% set(gcf,'position',[500 0 2000 1000])

%% 3(c)
L = 1;%length
T = 1;%time
a = 1;
hs = [0.05,0.02,0.01,0.002,0.001];
ts = [0.02,0.01,0.005,0.001,0.0005];
figure(7)
for h = 1:1:size(hs,2)
    deltax = hs(h);%space step x
    deltat = ts(h);%time step
    N = L/deltax+1;
    K = T/deltat+1;
    
    r = a*deltat/deltax;
    for j = 1:1:N
        x(j) = deltax*(j-1);
        %%initial condition
        U1(j,1) = 2*exp(-(100*x(j)-50)^2/10);
        U2(j,1) = exp(-100*(x(j)-0.5)^2)*sin(80*pi*x(j));
        U3(j,1) = 2*exp(-(100*x(j)-50)^2/10);
        U4(j,1) = exp(-100*(x(j)-0.5)^2)*sin(80*pi*x(j));
    end
    
    t(1) = 0;
    for n = 2:1:K
        t(n) = (n-1)*deltat;
        for j = 1:1:N-1
            if j == 1
                U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(end-1,n-1)+U1(end-2,n-1))/2+...
                    +r^2*(U1(j,n-1)-2*U1(end-1,n-1)+U1(end-2,n-1))/2;
                U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(end-1,n-1)+U2(end-2,n-1))/2+...
                    +r^2*(U2(j,n-1)-2*U2(end-1,n-1)+U2(end-2,n-1))/2;
                U3(j,n) = (-r/6+r^3/6)*U3(end-2,n-1)+(r+r^2/2-r^3/2)*U3(end-1,n-1)+...
                    (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(j+1,n-1);
                U4(j,n) = (-r/6+r^3/6)*U4(end-2,n-1)+(r+r^2/2-r^3/2)*U4(end-1,n-1)+...
                (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(j+1,n-1);
            else if j == 2
                    U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(j-1,n-1)+U1(end-1,n-1))/2+...
                        +r^2*(U1(j,n-1)-2*U1(j-1,n-1)+U1(end-1,n-1))/2;
                    U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(j-1,n-1)+U2(end-1,n-1))/2+...
                        +r^2*(U2(j,n-1)-2*U2(j-1,n-1)+U2(end-1,n-1))/2;
                    U3(j,n) = (-r/6+r^3/6)*U3(end-1,n-1)+(r+r^2/2-r^3/2)*U3(j-1,n-1)+...
                        (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(j+1,n-1);
                    U4(j,n) = (-r/6+r^3/6)*U4(end-1,n-1)+(r+r^2/2-r^3/2)*U4(j-1,n-1)+...
                        (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(j+1,n-1);
                else
%                     if j == N-1
%                         U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(j-1,n-1)+U1(j-2,n-1))/2+...
%                             +r^2*(U1(j,n-1)-2*U1(j-1,n-1)+U1(j-2,n-1))/2;
%                         U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(j-1,n-1)+U2(j-2,n-1))/2+...
%                             +r^2*(U2(j,n-1)-2*U2(j-1,n-1)+U2(j-2,n-1))/2;
%                         U3(j,n) = (-r/6+r^3/6)*U3(j-2,n-1)+(r+r^2/2-r^3/2)*U3(j-1,n-1)+...
%                             (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(1,n-1);
%                         U4(j,n) = (-r/6+r^3/6)*U4(j-2,n-1)+(r+r^2/2-r^3/2)*U4(j-1,n-1)+...
%                             (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(1,n-1);
%                     else
                        U1(j,n) = U1(j,n-1)-r*(3*U1(j,n-1)-4*U1(j-1,n-1)+U1(j-2,n-1))/2+...
                            +r^2*(U1(j,n-1)-2*U1(j-1,n-1)+U1(j-2,n-1))/2;
                        U2(j,n) = U2(j,n-1)-r*(3*U2(j,n-1)-4*U2(j-1,n-1)+U2(j-2,n-1))/2+...
                            +r^2*(U2(j,n-1)-2*U2(j-1,n-1)+U2(j-2,n-1))/2;
                        U3(j,n) = (-r/6+r^3/6)*U3(j-2,n-1)+(r+r^2/2-r^3/2)*U3(j-1,n-1)+...
                            (1-r/2-r^2+r^3/2)*U3(j,n-1)+(-r/3+r^2/2-r^3/6)*U3(j+1,n-1);
                        U4(j,n) = (-r/6+r^3/6)*U4(j-2,n-1)+(r+r^2/2-r^3/2)*U4(j-1,n-1)+...
                            (1-r/2-r^2+r^3/2)*U4(j,n-1)+(-r/3+r^2/2-r^3/6)*U4(j+1,n-1);
%                     end
                end
            end
        end
        U1(end,n) = U1(1,n);
        U2(end,n) = U2(1,n);
        U3(end,n) = U3(1,n);
        U4(end,n) = U4(1,n);
    end
    
    subplot(4,size(hs,2),h)
    plot(x,U1(:,end),'r-o','LineWidth',1)
    hold on
    plot(x,U1(:,1),'k-','LineWidth',1)
    hold on
    grid on
    xlabel({'x'},'FontSize',20,'Interpreter','latex');
    xlim([0,1])
    ylabel({'U'},'FontSize',20,'Interpreter','latex');
    ylim([-0.5,2])
    axis square
    if h == 3
        title(strcat('Beam-Warming IC1,','t=',num2str(T),'s'));
    end
    legend(strcat('h=',num2str(hs(h))),'exact')
    set(gca, 'FontName','Times New Roman','FontSize', 20);
    
    subplot(4,size(hs,2),h+size(hs,2))
    plot(x,U2(:,end),'r-+','LineWidth',1)
    hold on
    plot(x,U2(:,1),'k-','LineWidth',1)
    hold on
    grid on
    xlabel({'x'},'FontSize',20,'Interpreter','latex');
    xlim([0,1])
    ylabel({'U'},'FontSize',20,'Interpreter','latex');
%     ylim([-5E-15,10E-15])
    axis square
    if h == 3
        title(strcat('Beam-Warming IC1,','t=',num2str(T),'s'));
    end
    legend(strcat('h=',num2str(hs(h))),'exact')
    set(gca, 'FontName','Times New Roman','FontSize', 20);
    
    subplot(4,size(hs,2),h+2*size(hs,2))
    plot(x,U3(:,end),'b-o','LineWidth',1)
    hold on
    plot(x,U3(:,1),'k-','LineWidth',1)
    hold on
    grid on
    xlabel({'x'},'FontSize',20,'Interpreter','latex');
    xlim([0,1])
    ylabel({'U'},'FontSize',20,'Interpreter','latex');
    ylim([-0.5,2])
    axis square
    if h == 3
        title(strcat('Third-Order Lax-Wendroff IC1,','t=',num2str(T),'s'));
    end
    legend(strcat('h=',num2str(hs(h))),'exact')
    set(gca, 'FontName','Times New Roman','FontSize', 20);
    
    subplot(4,size(hs,2),h+3*size(hs,2))
    plot(x,U4(:,end),'b-+','LineWidth',1)
    hold on
    plot(x,U4(:,1),'k-','LineWidth',1)
    hold on
    grid on
    xlabel({'x'},'FontSize',20,'Interpreter','latex');
    xlim([0,1])
    ylabel({'U'},'FontSize',20,'Interpreter','latex');
%     ylim([-5E-15,10E-15])
    axis square
    if h == 3
        title(strcat('Third-Order Lax-Wendroff IC1,','t=',num2str(T),'s'));
    end
    legend(strcat('h=',num2str(hs(h))),'exact')
    set(gca, 'FontName','Times New Roman','FontSize', 20);
    
%     figure(3)
%     Ts = [1,31,61,91,121,151,181,191,201];
%     for j = 1:1:9
%         subplot(3,3,j)
%         plot(x,U2(:,Ts(j)),'r-o','LineWidth',2)
%         hold on
%         grid on
%         xlabel({'x'},'FontSize',20,'Interpreter','latex');
%         xlim([0,1])
%         ylabel({'U'},'FontSize',20,'Interpreter','latex');
%         ylim([-1,1])
%         axis square
%         title(strcat('t=',num2str(t(Ts(j))),'s'));
%         legend('IC2')
%         set(gca, 'FontName','Times New Roman','FontSize', 20);
%     end
%     set(gcf,'position',[500 0 1200 1200])
    err(1,h) = norm(abs(U1(:,end)-U1(:,1)));
    err(2,h) = norm(abs(U2(:,end)-U2(:,1)));
    err(3,h) = norm(abs(U3(:,end)-U3(:,1)));
    err(4,h) = norm(abs(U4(:,end)-U4(:,1)));
end
set(gcf,'position',[500 0 2000 1500])

figure(8)
for i=1:1:4
    subplot(2,2,i)
    plot(hs,err(i,:),'o-','LineWidth',3)
    hold on
    grid on
    xlabel({'h'},'FontSize',20,'Interpreter','latex');
    set(gca, 'XScale', 'log')
    ylabel({'error'},'FontSize',20,'Interpreter','latex');
    set(gca, 'YScale', 'log')
    axis square
    if i == 1
        title('Beam-Warming IC1');
    else if i == 2
            title('Beam-Warming IC2');
        else if i == 3
                title('Third-Order Lax-Wendroff IC1');
            else
                title('Third-Order Lax-Wendroff IC2');
            end
        end
    end
    legend(strcat('h=',num2str(hs(h))))
    set(gca, 'FontName','Times New Roman','FontSize', 20);
end
set(gcf,'position',[500 0 1000 1000])

for i = 1:1:4
    X = zeros(size(hs,2),2);
    X(:,1) = 1;
    X(:,2) = log(hs)';
    b = log(err(i,:))';
    k = X\b;
    p(i) = k(2);
end