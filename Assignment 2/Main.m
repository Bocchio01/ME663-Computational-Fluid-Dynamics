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
    "N_cells", 400), ...
    "flux_vector_splitting_type", ["VL"], ...
    "flux_vector_splitting_order", ["UDS"]);

sim.dx = (sim.x.domain(2) - sim.x.domain(1)) / sim.x.N_cells;
sim.x_vec = linspace(sim.x.domain(1), sim.x.domain(2), sim.x.N_cells);

output = struct('plots', true);


%% Solver initialization

solver = struct( ...
    "rho", zeros(1, sim.x.N_cells), ...
    "u", zeros(1, sim.x.N_cells), ...
    "p", zeros(1, sim.x.N_cells), ...
    "E", zeros(1, sim.x.N_cells), ...
    "a", zeros(1, sim.x.N_cells), ...
    "U", zeros(3, sim.x.N_cells));


%% Anonymous functions

fun = struct();

fun.dt = @(dx, u, a, CFL) CFL * dx / max(abs(u) + a);

fun.p = @(rho, u, E, gamma) (gamma - 1) * rho .* (E - 1/2 * u.^2);
fun.E = @(rho, p, gamma, u) p ./ ((gamma - 1) * rho) + 1/2 * u.^2;
fun.a = @(rho, p, gamma)    sqrt(gamma * p ./ rho);

fun.U = @(rho, u, E)               [(rho)'; (rho .* u)'; (rho .* E)'];
fun.F = @(rho, u, a, gamma, lamda) rho / (2 * gamma) * [
    2 * (gamma - 1) * lamda(1) + lamda(2) + lamda(3);
    2 * (gamma - 1) * lamda(1) * u + lamda(2) * (u + a) + lamda(3) * (u - a);
    (gamma - 1) * lamda(1) * u.^2 + 1/2 * lamda(2) * (u + a).^2 + 1/2 * lamda(3) * (u - a).^2 + (3 - gamma) / (2 * (gamma - 1)) * (lamda(2) + lamda(3)) * a.^2
    ];


%% Initial condition

half_idx = floor(length(sim.x_vec) / 2);

[U_initial, U_final] = compute_BCs_U(fun, fluid.BCs, fluid.gamma);
solver.U(:, 1:half_idx) = repmat(U_initial, 1, half_idx);
solver.U(:, half_idx + 1:end) = repmat(U_final, 1, sim.x.N_cells - half_idx);

[solver.rho, solver.u, solver.p, solver.E, solver.a] = compute_solver_vectors(solver.U, fun, fluid.gamma);

clear U_initial U_final half_idx


%% Main loop

iteration = 0;
dt = [0];

while max(cumsum(dt)) < sim.time_end
    
    iteration = iteration + 1;
    dt(iteration) = [fun.dt(sim.dx, solver.u, solver.a, sim.CFL)];
    
    [rho, u, ~, ~, a] = compute_solver_vectors(solver.U, fun, fluid.gamma);
    lambdas = compute_lambda(fluid, solver.U, fun, sim.flux_vector_splitting_type);
    
    for cell_idx = 1:sim.x.N_cells
        
        idx_WW = max(cell_idx - 2, 1);
        idx_W = max(cell_idx - 1, 1);
        idx_P = cell_idx;
        idx_E = min(cell_idx + 1, sim.x.N_cells);
        idx_EE = min(cell_idx + 2, sim.x.N_cells);
        
        
        switch sim.flux_vector_splitting_order
            case "UDS"
                
                weight = 0;
                
            case "LUDS"
                
                weight = 1/2;
                
            otherwise
                
                error('Unknown flux vector splitting order');
                
        end
        
        Fe = fun.F(...
            rho(idx_P) + weight * (rho(idx_P) - rho(idx_W)), ...
            u(idx_P) + weight * (u(idx_P) - u(idx_W)), ...
            a(idx_P) + weight * (a(idx_P) - a(idx_W)), ...
            fluid.gamma, ...
            lambdas(1, :, idx_P) + weight * (lambdas(1, :, idx_P) - lambdas(1, :, idx_W))) + ...
            fun.F(...
            rho(idx_E) + weight * (rho(idx_E) - rho(idx_EE)), ...
            u(idx_E) + weight * (u(idx_E) - u(idx_EE)), ...
            a(idx_E) + weight * (a(idx_E) - a(idx_EE)), ...
            fluid.gamma, ...
            lambdas(2, :, idx_E) + weight * (lambdas(2, :, idx_E) - lambdas(2, :, idx_EE)));
        
        Fw = fun.F(...
            rho(idx_W) + weight * (rho(idx_W) - rho(idx_WW)), ...
            u(idx_W) + weight * (u(idx_W) - u(idx_WW)), ...
            a(idx_W) + weight * (a(idx_W) - a(idx_WW)), ...
            fluid.gamma, ...
            lambdas(1, :, idx_W) + weight * (lambdas(1, :, idx_W) - lambdas(1, :, idx_WW))) + ...
            fun.F(...
            rho(idx_P) + weight * (rho(idx_P) - rho(idx_E)), ...
            u(idx_P) + weight * (u(idx_P) - u(idx_E)), ...
            a(idx_P) + weight * (a(idx_P) - a(idx_E)), ...
            fluid.gamma, ...
            lambdas(2, :, idx_P) + weight * (lambdas(2, :, idx_E) - lambdas(2, :, idx_E)));
        
        % Fe = fun.F(rho(idx_P), u(idx_P), a(idx_P), fluid.gamma, lambdas(1, :, idx_P)) + ...
        %     fun.F(rho(idx_E), u(idx_E), a(idx_E), fluid.gamma, lambdas(2, :, idx_E));
        
        % Fw = fun.F(rho(idx_W), u(idx_W), a(idx_W), fluid.gamma, lambdas(1, :, idx_W)) + ...
        %     fun.F(rho(idx_P), u(idx_P), a(idx_P), fluid.gamma, lambdas(2, :, idx_P));
        
        solver.U(:, idx_P) = solver.U(:, idx_P) - dt(end) / sim.dx * (Fe - Fw);
        
    end
    
    [solver.U(:, 1), solver.U(:, end)] = compute_BCs_U(fun, fluid.BCs, fluid.gamma);
    [solver.rho, solver.u, solver.p, solver.E, solver.a] = compute_solver_vectors(solver.U, fun, fluid.gamma);
    
    if min(solver.p) < 0
        disp(iteration)
        disp('negative pressure found!');
    end
    
end


%% Plots

if output.plots
    
    figure_results_analysis = figure('Name', 'Results analysis', 'NumberTitle', 'off');
    tiledlayout(4, 1);
    
    % rho
    nexttile;
    hold on
    grid on
    
    plot(sim.x_vec, solver.rho, 'k-');
    
    ylabel('x [m]');
    ylabel('rho');
    
    % u
    nexttile;
    hold on
    grid on
    
    plot(sim.x_vec, solver.u, 'r-');
    
    ylabel('x [m]');
    ylabel('u');
    
    % p
    nexttile;
    hold on
    grid on
    
    plot(sim.x_vec, solver.p, 'b-');
    
    ylabel('x [m]');
    ylabel('p');
    
    % E
    nexttile;
    hold on
    grid on
    
    plot(sim.x_vec, solver.E, 'g-');
    
    ylabel('x [m]');
    ylabel('E');
    
end



%% Functions

function [rho, u, a, lambdas] = interpolate_flux_splitting(flux_vector_splitting_order, fluid, solver, cell_idx)

decide = @(east, west, lambdas) east * (lambdas(1) > 0) + west * (lambdas(2) < 0);

idx_WW = max(cell_idx - 2, 1);
idx_W = max(cell_idx - 1, 1);
idx_P = cell_idx;
idx_E = min(cell_idx + 1, sim.x.N_cells);
idx_EE = min(cell_idx + 2, sim.x.N_cells);

U = solver.U(:, [idx_WW, idx_W, idx_P, idx_E, idx_EE]);

[rho, u, ~, ~, a] = compute_solver_vectors(U, fun, fluid.gamma);
lambdas = compute_lambda(fluid, U, flux_vector_splitting_type);

switch flux_vector_splitting_order
    case "UDS"
        
        U_interpolated_e.rho = decide(rho(3), rho(4), lambdas(1:2, 3));
        U_interpolated_e.u = decide(u(3), u(4), lambdas(1:2, 3));
        U_interpolated_e.a = decide(a(3), a(4), lambdas(1:2, 3));
        
        U_interpolated_w.rho = decide(rho(2), rho(3), lambdas(1:2, 2));
        U_interpolated_w.u = decide(u(2), u(3), lambdas(1:2, 2));
        U_interpolated_w.a = decide(a(2), a(3), lambdas(1:2, 2));
        
    case "LUDS"
        
        U_interpolated_e.rho = decide(rho(3) + 1/2 * (rho(3) - rho(2)), rho(4) - 1/2 * (rho(4) - rho(5)), lambdas(1:2, 3));
        U_interpolated_e.u = decide(u(3) + 1/2 * (u(3) - u(2)), u(4) - 1/2 * (u(4) - u(5)), lambdas(1:2, 3));
        U_interpolated_e.a = decide(a(3) + 1/2 * (a(3) - a(2)), a(4) - 1/2 * (a(4) - a(5)), lambdas(1:2, 3));
        
        U_interpolated_w.rho = decide(rho(2) + 1/2 * (rho(2) - rho(1)), rho(3) - 1/2 * (rho(3) - rho(4)), lambdas(1:2, 2));
        U_interpolated_w.u = decide(u(2) + 1/2 * (u(2) - u(1)), u(3) - 1/2 * (u(3) - u(4)), lambdas(1:2, 2));
        U_interpolated_w.a = decide(a(2) + 1/2 * (a(2) - a(1)), a(3) - 1/2 * (a(3) - a(4)), lambdas(1:2, 2));
        
    otherwise
        error('Unknown flux vector splitting order');
end

rho = [U_interpolated_w.rho, U_interpolated_e.rho];
u = [U_interpolated_w.u, U_interpolated_e.u];
a = [U_interpolated_w.a, U_interpolated_e.a];
lambdas = [lambdas(1, 2), lambdas(2, 3)];

end


function lambda_vec = compute_lambda(fluid, U, fun, flux_vector_splitting_type)
% Compute the eigenvalues of the flux vector splitting

lambda_vec = zeros(2, 3, size(U, 2));
[~, u, ~, ~, a] = compute_solver_vectors(U, fun, fluid.gamma);

switch flux_vector_splitting_type
    
    
    case "SW"
        eigenvalues = @(u, a) [
            u;
            u + a;
            u - a
            ];
        
        lambda_vec(1, :, :) = max(eigenvalues(u, a), 0);
        lambda_vec(2, :, :) = min(eigenvalues(u, a), 0);
        
    case "VL"
        eigenvalues = @(a, M, gamma) [
            1/4 * a .* (M + 1).^2 .* (1 - (M - 1).^2 / (gamma + 1));
            1/4 * a .* (M + 1).^2 .* (3 - M + (gamma - 1) / (gamma + 1) * (M - 1).^2);
            1/2 * a .* (M + 1).^2 .* (M - 1)/(gamma + 1) .* (1 + (gamma - 1)/2 * M)
            ];
        
        M = u ./ a;
        lambda_vec(1, :, :) = + eigenvalues(a, M, fluid.gamma);
        lambda_vec(2, :, :) = - eigenvalues(a, - M, fluid.gamma);
        
    otherwise
        error('Unknown flux vector splitting type');
        
end

end


function [U_initial, U_final] = compute_BCs_U(fun, BCs, gamma)

U_initial = fun.U(...
    BCs.rho(1), ...
    BCs.u(1), ...
    fun.E(BCs.rho(1), BCs.p(1), gamma, BCs.u(1)));

U_final = fun.U(...
    BCs.rho(2), ...
    BCs.u(2), ...
    fun.E(BCs.rho(2), BCs.p(2), gamma, BCs.u(2)));

end


function [rho, u, p, E, a] = compute_solver_vectors(U, fun, gamma)

rho = U(1, :);
u = U(2, :) ./ rho;
E = U(3, :) ./ rho;
p = fun.p(rho, u, E, gamma);
a = fun.a(rho, p, gamma);

end