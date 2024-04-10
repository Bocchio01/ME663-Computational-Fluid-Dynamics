function figure_results_analysis = plot_states(rho, u, p, E, a, t, x_vec)

figure_results_analysis = figure('Name', 'Results analysis', ...
    'NumberTitle', 'off', ...
    'Position', [100 100 1000 600]);
tiledlayout(2, 2);

% rho
nexttile;
hold on
grid on

plot(x_vec, rho, 'k-', 'LineWidth', 2);

title(['rho @ t = ', num2str(t, '%.3f'), ' s']);
xlabel('x [m]');
ylabel('rho [kg/m^3]');

% u
nexttile;
hold on
grid on

plot(x_vec, u, 'k-', 'LineWidth', 2);

title(['u @ t = ', num2str(t, '%.3f'), ' s']);
xlabel('x [m]');
ylabel('u [m/s]');

% p
nexttile;
hold on
grid on

plot(x_vec, p, 'k-', 'LineWidth', 2);

title(['p @ t = ', num2str(t, '%.3f'), ' s']);
xlabel('x [m]');
ylabel('p [MPa]');

% E
nexttile;
hold on
grid on

plot(x_vec, E, 'k-', 'LineWidth', 2);

title(['E @ t = ', num2str(t, '%.3f'), ' s']);
xlabel('x [m]');
ylabel('E [J/kg]');

end
