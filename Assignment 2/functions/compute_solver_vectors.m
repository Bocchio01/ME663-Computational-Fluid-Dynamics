function [rho, u, p, E, a] = compute_solver_vectors(U, funcs, gamma)

rho = U(1, :);
u = U(2, :) ./ rho;
E = U(3, :) ./ rho;
p = funcs.p(rho, u, E, gamma);
a = funcs.a(rho, p, gamma);

end