figure
plot(Mrange/3,mean(sumSINR1,2),'r-*','LineWidth',2);
hold on
plot(Mrange/3,mean(sumSINR2,2),'b--','LineWidth',2);
plot(Mrange/3,mean(sumSINR3,2),'k-.','LineWidth',2);
plot(Mrange/3,mean(sumSINR4,2),'g-.','LineWidth',2);
plot(Mrange/3,mean(sumSINR5,2),'m--o','LineWidth',2);
%set(gca,'xLim',[160,220]);
xlabel('M')
ylabel('Sum Rate (bps/Hz)')
legend('RZF','CG','PCG','GS','JOR','Location','best')
grid on