% ME663 - Computational Fluid Dynamics
% Assignment 2
% Tommaso Bocchietti
% A.Y. 2023/24 - W24

clc
clear variables
close all

%% Parameters

fluid = struct( ...
    "gamma", 1.4, ...
    "BCs", struct( ...
    "rho", [1 0.125], ...
    "u", [0 0], ...
    "p", [100000 10000]));

sim = struct( ...
    "CFL", 1/2, ...
    "time_end", 0.01, ...
    "x", struct( ...
    "domain", [-10 10], ...
    "N_cells", 50), ...
    "flux_vector_splitting_type", "SW", ...
    "flux_vector_splitting_order", "UDS");

sim.dx = (sim.x.domain(2) - sim.x.domain(1)) / sim.x.N_cells;
sim.x_vec = linspace(sim.x.domain(1), sim.x.domain(2), sim.x.N_cells);

output = struct( ...
    'plots', true, ...
    'movie', false);

% latex_img_path = 'latex/img';


%% Solver initialization

solver = struct( ...
    "rho", zeros(1, sim.x.N_cells), ...
    "u", zeros(1, sim.x.N_cells), ...
    "p", zeros(1, sim.x.N_cells), ...
    "E", zeros(1, sim.x.N_cells), ...
    "a", zeros(1, sim.x.N_cells), ...
    "U", zeros(3, sim.x.N_cells));


%% Anonymous functions

funcs = struct();

funcs.dt = @(dx, u, a, CFL) CFL * dx / max(abs(u) + a);

funcs.p = @(rho, u, E, gamma) (gamma - 1) * rho .* (E - 1/2 * u.^2);
funcs.E = @(rho, p, gamma, u) p ./ ((gamma - 1) * rho) + 1/2 * u.^2;
funcs.a = @(rho, p, gamma)    sqrt(gamma * p ./ rho);

funcs.U = @(rho, u, E) [
    rho;
    rho .* u;
    rho .* E
    ];


%% Initial condition

half_idx = floor(length(sim.x_vec) / 2);

[U_initial, U_final] = compute_BCs_U(funcs, fluid.BCs, fluid.gamma);
solver.U(:, 1:half_idx) = repmat(U_initial, 1, half_idx);
solver.U(:, half_idx + 1:end) = repmat(U_final, 1, sim.x.N_cells - half_idx);

clear U_initial U_final half_idx


%% Main loop

iteration = 0;
dt = [0];
t = 0;

N = sim.x.N_cells;
Fe = zeros(3, N);
Fw = zeros(3, N);

[solver.rho, solver.u, solver.p, solver.E, solver.a] = compute_solver_vectors(solver.U, funcs, fluid.gamma);

while t < sim.time_end
    
    iteration = iteration + 1;
    dt(iteration) = [funcs.dt(sim.dx, solver.u, solver.a, sim.CFL)];
    t = t + dt(iteration);
    
    % Compute fluxes at the center of each cell
    switch sim.flux_vector_splitting_type
        
        case "SW"
            [F_plus, F_minus] = compute_flux_SW(solver.rho, solver.u, solver.p, solver.E, solver.a, fluid.gamma);
            
        case "VL"
            [F_plus, F_minus] = compute_flux_VL(solver.rho, solver.u, solver.p, solver.E, solver.a, fluid.gamma);
            
        otherwise
            error('Unknown flux vector splitting type');
            
    end
    
    % Compute fluxes at the cell interfaces using an interpolation of the fluxes at the cell centers
    for i = 1:N
        
        switch sim.flux_vector_splitting_order
            
            case "UDS"
                weight = 0;
                
            case "LUDS"
                weight = 1/2;
                
                if(iteration == 1 && i == 1)
                    disp('LUDS splitting order may not be stable');
                end
                
            otherwise
                error('Unknown flux vector splitting order');
                
        end
        
        WW = max(i - 2, 1);
        W = max(i - 1, 1);
        P = i;
        E = min(i + 1, N);
        EE = min(i + 2, N);
        
        Fe(:, i) = ...
            (F_plus(:, P) + weight * (F_plus(:, P) - F_plus(:, W))) + ...
            (F_minus(:, E) + weight * (F_minus(:, E) - F_minus(:, EE)));
        
        Fw(:, i) = ...
            (F_plus(:, W) + weight * (F_plus(:, W) - F_plus(:, WW))) + ...
            (F_minus(:, P) + weight * (F_minus(:, P) - F_minus(:, E)));
        
    end
    
    % Update the solution and apply boundary conditions
    solver.U(:, :) = solver.U(:, :) - dt(end) / sim.dx * (Fe - Fw);
    [solver.U(:, 1), solver.U(:, end)] = compute_BCs_U(funcs, fluid.BCs, fluid.gamma);
    
    % Compute the new state vectors
    [solver.rho, solver.u, solver.p, solver.E, solver.a] = compute_solver_vectors(solver.U, funcs, fluid.gamma);
    
    if any(solver.p < 0)
        error(['Negative pressure found @iteration=' num2str(iteration-1)]);
    end
    
    if (output.movie)
        set(0,'DefaultFigureVisible','off');
        plot_figures(iteration) = getframe(plot_states(solver.rho, solver.u, solver.p, solver.E, solver.a, t, sim.x_vec));
    end
    
end


%% Load solutions

solution = struct( ...
    x_vec = [], ...
    data = struct());

[solution.x_vec, solution.data.rho] = read_file_solution('solutions/density.exa');
[~, solution.data.S] = read_file_solution('solutions/entropy.exa');
[~, solution.data.M] = read_file_solution('solutions/mach.exa');
[~, solution.data.rho_u] = read_file_solution('solutions/massflux.exa');
[~, solution.data.p] = read_file_solution('solutions/pressure.exa');
[~, solution.data.u] = read_file_solution('solutions/velocity.exa');


%% Movie

if (output.movie)
    
    video_object = VideoWriter([char(sim.flux_vector_splitting_type) '_' char(sim.flux_vector_splitting_order) '_' num2str(sim.x.N_cells) '.avi'], 'Motion JPEG AVI');
    video_object.Quality = 95;
    video_object.FrameRate = 10;
    
    open(video_object);
    writeVideo(video_object, plot_figures);
    close(video_object);
    
end


%% Plots

set(0,'DefaultFigureVisible','on');
if output.plots
    
    % State vectors (computed vs. exact)
    figure_results_analysis = plot_states(solver.rho, solver.u, solver.p, solver.E, solver.a, t, sim.x_vec);
    
    % rho
    nexttile(1);
    plot(solution.x_vec, solution.data.rho, 'r-');
    
    % u
    nexttile(2);
    plot(solution.x_vec, solution.data.u, 'r-');
    
    % p
    nexttile(3);
    plot(solution.x_vec, solution.data.p, 'r-');
    
    % E
    nexttile(4);
    plot(solution.x_vec, funcs.E(solution.data.rho, solution.data.p, fluid.gamma, solution.data.u), 'r-');
    
    
    % Colormap with the u velocity (computed vs. exact)
    figure_velocity_gradient = figure('Name', 'Velocity analysis', 'NumberTitle', 'off');
    tiledlayout(2, 1);
    
    nexttile(1);
    hold on
    grid on
    
    colormap('parula');
    imagesc(sim.x_vec, 0, solver.u);
    colorbar;
    
    axis tight;
    xlabel('x [m]');
    
    nexttile(2);
    hold on
    grid on
    
    colormap('parula');
    imagesc(sim.x_vec, 0, solution.data.u');
    colorbar;
    
    axis tight;
    xlabel('x [m]');
    
    if(exist('latex_img_path', 'var') == 1)
        saveas(figure_results_analysis, [latex_img_path '/states/' char(sim.flux_vector_splitting_type) num2str(sim.x.N_cells) '.png']);
        saveas(figure_velocity_gradient, [latex_img_path '/velocity/' char(sim.flux_vector_splitting_type) num2str(sim.x.N_cells) '.png']);
    end
    
end

