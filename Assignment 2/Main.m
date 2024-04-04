% ME663 - Computational Fluid Dynamics
% Assignment 2
% Tommaso Bocchietti
% A.Y. 2023/24 - W24

clc
clear variables
close all

PLUS = 1;

%% Parameters

flow = struct( ...
    "gamma", 1.4, ...
    "rho", [1 0.125], ...
    "u", [0 0], ...
    "p", [100000 10000]);

sim = struct( ...
    "CFL", 1/2, ...
    "time_end", 0.01, ...
    "N_cells", 400, ...
    "x", [-10 10]);

output = struct('plots', true);


sim.dx = (sim.x(2) - sim.x(1)) / sim.N_cells;
sim.x_vec = linspace(sim.x(1), sim.x(2), sim.N_cells);

solver = struct( ...
    "rho", zeros(1, sim.N_cells), ...
    "u", zeros(1, sim.N_cells), ...
    "p", zeros(1, sim.N_cells), ...
    "E", zeros(1, sim.N_cells), ...
    "a", zeros(1, sim.N_cells), ...
    "U", zeros(3, sim.N_cells));

anonymous = struct();

anonymous.a = @(p, rho, gamma) sqrt(gamma * p ./ rho);
anonymous.p = @(rho, u, E, gamma) (gamma - 1) * rho .* (E - 1/2 * u.^2);
anonymous.E = @(p, rho, gamma, u) p ./ ((gamma - 1) * rho) + 1/2 * u.^2;

anonymous.U = @(rho, u, E) [rho'; rho' .* u'; rho' .* E'];

% anonymous.F = @(gamma, rho, a, u, eigenvalues) 1 / gamma * anonymous.Qa(gamma, rho, a, u) * [
%     (gamma - 1) * rho * eigenvalues(1);
%     +a * eigenvalues(2);
%     -a * eigenvalues(3)
%     ];

anonymous.F = @(gamma, u, a, lamda, rho) rho / (2 * gamma) * [
    2 * (gamma - 1) * lamda(1) + lamda(2) + lamda(3);
    2 * (gamma - 1) * lamda(1) * u + lamda(2) * (u + a) + lamda(3) * (u - a);
    (gamma - 1) * lamda(1) * u.^2 + 1/2 * lamda(2) * (u + a).^2 + 1/2 * lamda(3) * (u - a).^2 + (3 - gamma) / (2 * (gamma - 1)) * (lamda(2) + lamda(3)) * a.^2
    ];

anonymous.SW_plus = @(a, u) 1/2 * ([u; u + a; u - a] + abs([u; u + a; u - a]));
anonymous.SW_minus = @(a, u) 1/2 * ([u; u + a; u - a] - abs([u; u + a; u - a]));

anonymous.VL_plus = @(a, M, gamma) 1/4 * a * (M + 1)^2 * [
    1 + (M - 1)^2 / (gamma + 1);
    3 - M + (gamma - 1) / (gamma + 1) * (M - 1)^2;
    2 * (M - 1)/(gamma + 1) * (1 + (gamma - 1)/2 * M)
    ];

anonymous.VL_minus = @(a, M, gamma) - VL_plus(a, -M, gamma);

%% Initial condition

half_idx = floor(length(sim.x_vec) / 2);
solver.rho(1:half_idx) = flow.rho(1);
solver.u(1:half_idx) = flow.u(1);
solver.p(1:half_idx) = flow.p(1);

solver.rho(half_idx + 1:end) = flow.rho(2);
solver.u(half_idx + 1:end) = flow.u(2);
solver.p(half_idx + 1:end) = flow.p(2);

solver.E(:) = anonymous.E(solver.p(:), solver.rho(:), flow.gamma, solver.u(:));
solver.a(:) = anonymous.a(solver.p(:), solver.rho(:), flow.gamma);
solver.U(:, :) = anonymous.U(solver.rho(:), solver.u(:), solver.E(:));


% Solver loop

t = 0;
it = 0;
dt = sim.CFL * sim.dx / max(abs(solver.u(:)) + solver.a(:));

flux_vector_splitting_type = "SW";
flux_vector_splitting_order = "UDS";

while t < sim.time_end
    
    
    lamda_vec = compute_lambda(solver, flux_vector_splitting_type);
    
    for cell_idx = 2:sim.N_cells-1
        
        cell_idx_WW = max(cell_idx - 2, 1);
        cell_idx_W = max(cell_idx - 1, 1);
        cell_idx_P = cell_idx;
        cell_idx_E = min(cell_idx + 1, sim.N_cells);
        cell_idx_EE = min(cell_idx + 2, sim.N_cells);
        
        switch flux_vector_splitting_order
            case "UDS"
                U = [
                    solver.U(:, cell_idx_P);
                    solver.U(:, cell_idx_P);
                    solver.U(:, cell_idx_E);
                    ];
                
                % case "LUDS"
                
                
            otherwise
                
        end
        
        Fe = anonymous.F(flow.gamma, u_P, a_P, lamda_plus(:, cell_idx_P), rho_P) + ...
            anonymous.F(flow.gamma, u_E, a_E, lamda_minus(:, cell_idx_E), rho_E);
        
        Fw = anonymous.F(flow.gamma, u_W, a_W, lamda_plus(:, cell_idx_W), rho_W) + ...
            anonymous.F(flow.gamma, u_P, a_P, lamda_minus(:, cell_idx_P), rho_P);
        
        solver.U(:, cell_idx_P) = solver.U(:, cell_idx_P) - dt / sim.dx * (Fe - Fw);
        
    end
    
    
    solver.U(:, 1) = anonymous.U(flow.rho(1), flow.u(1), anonymous.E(flow.p(1), flow.rho(1), flow.gamma, flow.u(1)));
    solver.U(:, end) = anonymous.U(flow.rho(2), flow.u(2), anonymous.E(flow.p(2), flow.rho(2), flow.gamma, flow.u(2)));
    
    solver.rho(:) = solver.U(1, :);
    solver.u(:) = solver.U(2, :) ./ solver.rho(:)';
    solver.E(:) = solver.U(3, :) ./ solver.rho(:)';
    solver.p(:) = anonymous.p(solver.rho(:), solver.u(:), solver.E(:), flow.gamma);
    solver.a(:) = anonymous.a(solver.p(:), solver.rho(:), flow.gamma);
    
    if min(solver.p(:)) < 0
        disp(it)
        disp('negative pressure found!');
    end
    
    dt = sim.CFL * sim.dx / max(abs(solver.u(:)) + solver.a(:));
    
    t = t + dt;
    it = it + 1;
end


%% Plots

if output.plots
    
    figure_results_analysis = figure('Name', 'Results analysis', 'NumberTitle', 'off');
    tiledlayout(4, 1);
    
    % rho
    nexttile;
    hold on
    grid on
    
    plot(sim.x_vec(:), solver.rho(:), 'k-');
    
    ylabel('x [m]');
    ylabel('rho');
    
    % u
    nexttile;
    hold on
    grid on
    
    plot(sim.x_vec(:), solver.u(:), 'r-');
    
    ylabel('x [m]');
    ylabel('u');
    
    % p
    nexttile;
    hold on
    grid on
    
    plot(sim.x_vec(:), solver.p(:), 'b-');
    
    ylabel('x [m]');
    ylabel('p');
    
    % E
    nexttile;
    hold on
    grid on
    
    plot(sim.x_vec(:), solver.E(:), 'g-');
    
    ylabel('x [m]');
    ylabel('E');
    
end



%% Functions

function lambda_vec = compute_lambda(solver, flux_vector_splitting_type)

anonymous.VL_plus = @(a, M, gamma) 1/4 * a * (M + 1)^2 * [
    1 + (M - 1)^2 / (gamma + 1);
    3 - M + (gamma - 1) / (gamma + 1) * (M - 1)^2;
    2 * (M - 1)/(gamma + 1) * (1 + (gamma - 1)/2 * M)
    ];

anonymous.VL_minus = @(a, M, gamma) - VL_plus(a, -M, gamma);


lambda_vec = zeros(2, 3, length(solver.rho));

switch flux_vector_splitting_type
    
    
    case "SW"
        eigenvalues = [
            solver.u;
            solver.u + solver.a;
            solver.u - solver.a
            ];
        lambda_vec(1, :, :) = max(eigenvalues, 0);
        lambda_vec(2, :, :) = min(eigenvalues, 0);
        
    case "VL"
        M(:) = solver.u(:) ./ solver.a(:);
        
        lambda_vec(1, :, :) = max(eigenvalues, 0);
        
        lamda_plus = anonymous.VL_plus(solver.a, M, solver.u);
        lamda_minus = anonymous.VL_minus(solver.a, M, solver.u);
end

end
