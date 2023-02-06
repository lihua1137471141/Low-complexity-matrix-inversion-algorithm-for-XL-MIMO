figure
plot(snr,mean(sumSINR1,2)/3,'r->','LineWidth',2);
hold on
plot(snr,mean(sumSINR2,2)/3,'b-o','LineWidth',2);
plot(snr,mean(sumSINR3,2)/3,'b-','LineWidth',2);
plot(snr,mean(sumSINR4,2)/3,'k-.','LineWidth',2);
plot(snr,mean(sumSINR5,2)/3,'k-*','LineWidth',2);
set(gca,'xLim',[15,30]);
xlabel('SNR (dB)','Interpreter','latex')
ylabel('Average sum SE of each subarray (bps/Hz/subarray)','Interpreter','latex')
set(gca, 'Fontname', 'Times New Roman','FontSize',12);
legend('RZF (reference SE)','CG','Jac-PCG','GS','JOR','latex');
grid on