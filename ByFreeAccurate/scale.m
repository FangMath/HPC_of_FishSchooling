clear all; close all;

% test results, strong scaling
p = 2.^[0:6]; % # of cores

% time elapsed on crunchy1, crunchy5, crunchy6
tS1 = [14.3724, 8.9278, 4.6605, 2.8891, 1.8193, 1.1621, 1.1672];
tS5 = [16.7837, 7.8605, 4.2547, 2.2695, 1.2194, 0.9931, 1.0197];
tS6 = [16.0404, 9.9345, 4.1601, 1.9838, 1.0847, 1.0534, 1.3045];

tW1 = [0.2385, 0.3719, 0.4562, 0.5153, 0.5535, 0.6804, 1.1672];
tW5 = [0.2890, 0.2884, 0.2882, 0.3334, 0.4320, 0.6114, 1.0197];
tW6 = [0.2815, 0.3065, 0.3799, 0.4003, 0.5123, 0.9921, 1.3045];

matlab1 = [0.1481, 0.3219, 0.5981, 1.1337, 2.2544, 4.4828, 6.0459];
matlab5 = [0.1244, 0.1993, 0.5623, 0.9899, 1.9136, 4.0066, 8.3481];
matlab6 = [0.1506, 0.3083, 0.6044, 1.2910, 2.4251, 4.7892, 9.5967];

f1 = figure(1);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontname', 'Arial')

wdth = 4; hght = 3.5;
set(gcf,'position',[100,100,wdth*100,hght*100]);
set(gcf, 'PaperPosition', [0 0 wdth hght]); %Position plot at left hand corner with width 2 and height 6.5.
set(gcf, 'PaperSize', [wdth hght]); %Set the paper to have width 2 and height 6.5.
bar(log(p)/log(2), [tS1;tS5;tS6]', 'BarWidth', 1); hold on;
legend('crunchy1', 'crunchy5', 'crunchy6', 'Location', 'best');
title(['Strong scaling for ByFree, 64 wings']);
xlabel('Number of procs');
xlim([-1,7]);
set(gca,'XTick',[0:6])  
set(gca,'XTickLabel',{'1','2', '4', '8', '16', '32', '64'})
ylabel('time = 1/efficiency (s)');

f2 = figure(2);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontname', 'Arial')

wdth = 4; hght = 3.5;
set(gcf,'position',[100,100,wdth*100,hght*100]);
set(gcf, 'PaperPosition', [0 0 wdth hght]); %Position plot at left hand corner with width 2 and height 6.5.
set(gcf, 'PaperSize', [wdth hght]); %Set the paper to have width 2 and height 6.5.
bar(log(p)/log(2), [tW1;tW5;tW6]', 'BarWidth', 1); hold on;
legend('crunchy1', 'crunchy5', 'crunchy6', 'Location', 'best');
title(['Weak scaling for ByFree']);
xlabel('Number of procs');
set(gca,'XTick',[0:6])  
set(gca,'XTickLabel',{'1','2', '4', '8', '16', '32', '64'})
xlim([-1,7]);
ylabel('time = 1/efficiency (s)');

saveas(f1, ['ByFreeStrong.fig']) 
saveas(f1, ['ByFreeStrong.png']) 
saveas(f2, ['ByFreeWeak.fig']) 
saveas(f2, ['ByFreeWeak.png']) 
