function [F_plus, F_minus] = compute_flux_SW(rho, u, ~, ~, a, gamma)

SW = @(u, a) [
    u;
    u + a;
    u - a
    ];

F = @(rho, u, a, gamma, lamda) rho / (2 * gamma) * [
    2 * (gamma - 1) * lamda(1) + lamda(2) + lamda(3);
    2 * (gamma - 1) * lamda(1) * u + lamda(2) * (u + a) + lamda(3) * (u - a);
    (gamma - 1) * lamda(1) * u.^2 + 1/2 * lamda(2) * (u + a).^2 + 1/2 * lamda(3) * (u - a).^2 + (3 - gamma) / (2 * (gamma - 1)) * (lamda(2) + lamda(3)) * a.^2
    ];

N = length(rho);
F_plus = zeros(3, N);
F_minus = zeros(3, N);

lambda = SW(u, a);
for i = 1:N
    F_plus(:, i) = F(rho(i), u(i), a(i), gamma, max(lambda(:, i), 0));
    F_minus(:, i) = F(rho(i), u(i), a(i), gamma, min(lambda(:, i), 0));
end

end
