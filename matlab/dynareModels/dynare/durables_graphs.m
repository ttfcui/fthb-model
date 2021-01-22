% make durables graphs

 dynare rbc_nocons.mod;
 
dur_inv_nocons = i./d;

save results_dur_nocons.mat dur_inv_nocons 

 dynare rbc.mod;
 
dur_inv = i./d;

save results_dur.mat dur_inv

close all

load  results_dur_nocons.mat  
load results_dur.mat

figure(1)
subplot(1,2,1)
bar(25:200,cumsum(dur_inv_nocons(25:200)/delta - 1),'barwidth',3);
ylim([-0.05,0.4]);
xlim([25,200]);
title('One good model', 'fontsize', 12)
subplot(1,2,2)
bar(25:200,cumsum(dur_inv(25:200)/delta - 1),'barwidth',3);
ylim([-.05,.15]);
xlim([25,200]);
title('Two good model', 'fontsize', 12)
set(gca,'fontsize', 11);
set(gcf, 'PaperPosition', [0 0 7.5 6.5]); %Position the plot further to the left and down. Extend the plot to fill entire paper.
set(gcf, 'PaperSize', [7.5 6.6]); %Keep the same paper size
fig=gcf;
saveas(fig, 'rbc/Output/RBC_nocons_dur_invest.pdf');
