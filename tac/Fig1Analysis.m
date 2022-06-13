
% Load and preprocess results
clear;
load('./results/fixed-budget-results.mat');
idxPeakExists = (minAbscissaStats.R0 > 1 & minR0Stats.R0 > 1);
disp(['Number of simulations with a peak:' num2str(sum(idxPeakExists(:)))])
 
% Histogram + stats of cumulative infections in scenarios with peaks 
cumIfxDiff = minAbscissaStats.CumIfx(idxPeakExists) - minR0Stats.CumIfx(idxPeakExists);
disp(['Minimizing R0 led to strictly fewer cumulative infections in '...
    num2str(100 * mean(cumIfxDiff(:) > 0)) '% of simulations.' num2str(sum(cumIfxDiff(:) > 0))]); 

figure;
cumIfxDiffClean = rmoutliers(cumIfxDiff(:), 'percentiles', [2.5 97.5]); 
histogram(cumIfxDiffClean, 'Normalization', 'probability', 'FaceColor', [0 0 0]);
xlabel('Excess Cumulative Infections');
ylabel('Fraction of Simulations');
figsize = [3 2.5];
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 figsize]);
set(gcf, 'PaperSize', figsize);
saveas(gcf, './figures/excess-cumulative-p1.pdf');

% Histogram + stats of peak infections in scenarios with peaks 
peakIfxDiff = minAbscissaStats.PeakIfx(idxPeakExists) - minR0Stats.PeakIfx(idxPeakExists);
disp(['Minimizing R0 led to strictly lower peaks in '...
    num2str(100 * mean(peakIfxDiff(:) > 0)) '% of simulations.' num2str(sum(peakIfxDiff(:) > 0))]);

figure;
peakIfxDiffClean = rmoutliers(peakIfxDiff(:), 'percentiles', [2.5 97.5]);
histogram(peakIfxDiffClean, 'Normalization', 'probability', 'FaceColor', [0 0 0]);
xlabel('Excess Peak Infections');
ylabel('Fraction of Simulations'); 
ylim([0, 0.16])
figsize = [3 2.5];
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 figsize]);
set(gcf, 'PaperSize', figsize);
saveas(gcf, './figures/excess-peak-p1.pdf');

% Stats on cumulative infections in models without peaks
cumIfxDiff = minAbscissaStats.CumIfx(~idxPeakExists) - minR0Stats.CumIfx(~idxPeakExists);
disp(['Minimizing R0 led to strictly fewer cumulative infections in '...
    num2str(100 * mean(cumIfxDiff(:) > 0)) '% of simulations.' num2str(sum(cumIfxDiff(:) > 0))]); 
