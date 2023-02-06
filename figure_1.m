
figure
plot(Mrange/3,mean(sumSINR1,2)/3,'r->','LineWidth',2);
hold on
plot(Mrange/3,mean(sumSINR2,2)/3,'b-o','LineWidth',2);
plot(Mrange/3,mean(sumSINR3,2)/3,'b-','LineWidth',2);
plot(Mrange/3,mean(sumSINR4,2)/3,'k-.','LineWidth',2);
plot(Mrange/3,mean(sumSINR5,2)/3,'k-*','LineWidth',2);
%set(gca,'xLim',[160,220]);
xlabel('Number of antennas (M)','Interpreter','latex')
ylabel('Average sum SE of each subarray (bps/Hz/subarray)','Interpreter','latex')
set(gca, 'Fontname', 'Times New Roman','FontSize',12);
legend('RZF (reference SE)','CG','Jac-PCG','GS','JOR','latex');
grid on