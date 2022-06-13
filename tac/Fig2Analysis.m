
load('./results/fixed-model-results');

i0 = 41;
figsize = [3 2.5];
 
% Plot cumulative infections
figure;  
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 figsize]);
set(gcf, 'PaperSize', figsize);

hold on;

plot(budgets(i0:end), minR0Stats.CumIfx(1, i0:end), 'k');
set(gca, 'YScale', 'log');
plot(budgets(i0:end), minAbscissaStats.CumIfx(1, i0:end), 'k--');
set(gca, 'YScale', 'log');

plot(budgets(i0:end), minR0Stats.CumIfx(2, i0:end), 'k');
set(gca, 'YScale', 'log');
plot(budgets(i0:end), minAbscissaStats.CumIfx(2, i0:end), 'k--');
set(gca, 'YScale', 'log');

plot(budgets(i0:end), minR0Stats.CumIfx(3, i0:end), 'k');
set(gca, 'YScale', 'log');
plot(budgets(i0:end), minAbscissaStats.CumIfx(3, i0:end), 'k--');
set(gca, 'YScale', 'log');

hold off;

xlabel('Budget');
ylabel('Cumulative Infections');
ylim([0.38e+05 1.7e+05]);
legend({'Min R0', 'Min Abscissa'});
saveas(gcf, './figures/fixed-model-cumulative.pdf');  

% Plot different allocations
figure;  
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 figsize]);
set(gcf, 'PaperSize', figsize); 

hold on;

plot(budgets, sum(minR0Stats.vaxAlloc(1, :, :), 3), ...
    'Color', '#0072BD', 'LineWidth', 1);
plot(budgets, sum(minAbscissaStats.vaxAlloc(1, :, :), 3), ...
    'LineStyle', '--', 'Color', '#0072BD', 'LineWidth', 1);

plot(budgets, sum(minR0Stats.vaxAlloc(2, :, :), 3), ...
    'Color', '#D95319', 'LineWidth', 1);
plot(budgets, sum(minAbscissaStats.vaxAlloc(2, :, :), 3), ...
    'Color', '#D95319', 'LineStyle', '--', 'LineWidth', 1);

plot(budgets, sum(minR0Stats.vaxAlloc(3, :, :), 3), ...
    'Color', '#EDB120', 'LineWidth', 1);
plot(budgets, sum(minAbscissaStats.vaxAlloc(3, :, :), 3), ...
    'Color', '#EDB120', 'LineStyle', '--', 'LineWidth', 1);

leg = legend({
    'Model 1, Min R_0', 'Model 1, Min \alpha', ...
    'Model 2, Min R_0', 'Model 2, Min \alpha', ...
    'Model 3, Min R_0', 'Model 3, Min \alpha'}, ...
    'Location', 'southoutside', ...
    'Orientation', 'vertical', 'NumColumns', 3);
xlabel('Budget');
ylabel('Total Vaccine Allocation'); 

hold off;

saveas(gcf, './figures/legend.pdf');  
