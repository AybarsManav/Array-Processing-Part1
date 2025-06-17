function plot_singular_values(singular_values)
figure;
plot(singular_values, 'o-', 'LineWidth', 1.5);
title('Singular Values', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Index', 'FontSize', 12);
ylabel('Singular Value', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
end

