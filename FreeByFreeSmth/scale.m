clear all; close all;

% test results, strong scaling
Nw = 64; % total number of wings
p = 2.^[0:6]; % # of cores

% time elapsed on crunchy1, crunchy5, crunchy6
tS1 = [5.6951, 3.0587, 1.4574, 0.7470, 0.4061, 0.2837, 0.3101];
tS5 = [7.5522, 3.5475, 1.7854, 1.1311, 0.6257, 0.3067, 0.3601];
tS6 = [6.3187, 3.6081, 2.1575, 1.2735, 0.4428, 0.5839, 0.4183];

tW1 = [0.1071, 0.1143, 0.1602, 0.2122, 0.2471, 0.3187, 0.3101];
tW5 = [0.1171, 0.1283, 0.1283, 0.1259, 0.1796, 0.2407, 0.3601];
tW6 = [0.1012, 0.1046, 0.1462, 0.1869, 0.2482, 0.3704, 0.4183];

matlab1 = [0.0642, 0.1429, 0.2353, 0.4106, 1.2119, 2.1445, 4.0383];
matlab5 = [0.0845, 0.1836, 0.2987, 0.7143, 1.2498, 2.6172, 5.3104];
matlab6 = [0.0702, 0.1705, 0.2613, 0.6452, 1.2382, 2.3314, 5.0301];

f1 = figure(1);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontname', 'Arial')

wdth = 4; hght = 3.5;
set(gcf,'position',[100,100,wdth*100,hght*100]);
set(gcf, 'PaperPosition', [0 0 wdth hght]); %Position plot at left hand corner with width 2 and height 6.5.
set(gcf, 'PaperSize', [wdth hght]); %Set the paper to have width 2 and height 6.5.
bar(log(p)/log(2), [tS1;tS5;tS6]', 'BarWidth', 1); hold on;
legend('crunchy1', 'crunchy5', 'crunchy6', 'Location', 'NorthEast');
title(['Strong scaling for FreeSmth, 64 wings']);
xlim([-1,7]);
xlabel('Number of procs');
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
legend('crunchy1', 'crunchy5', 'crunchy6', 'Location', 'NorthWest');
title(['Weak scaling for FreeSmth']);
xlabel('Number of procs');
set(gca,'XTick',[0:6])  
set(gca,'XTickLabel',{'1','2', '4', '8', '16', '32', '64'})
xlim([-1,7]);
ylabel('time = 1/efficiency (s)');

saveas(f1, ['FreeSmthStrong.fig']) 
saveas(f1, ['FreeSmthStrong.png']) 
saveas(f2, ['FreeSmthWeak.fig']) 
saveas(f2, ['FreeSmthWeak.png']) 
