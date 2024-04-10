function [F_plus, F_minus] = compute_flux_VL(rho, u, p, E, a, gamma)

F = @(rho, u, p, E, ~, ~) [
    rho * u;
    rho * u^2 + p;
    (rho * E + p) * u
    ];

Fmod_plus = @(rho, u, p, E, a, gamma) rho*a/4 * (u/a+1)^2 * [
    1;
    ((gamma-1) * u + 2*a) / gamma;
    ((gamma-1) * u + 2*a)^2 / (2*(gamma+1)*(gamma-1))
    ];

Fmod_minus = @(rho, u, p, E, a, gamma) -rho*a/4 * (u/a-1)^2 * [
    1;
    ((gamma-1) * u - 2*a) / gamma;
    ((gamma-1) * u - 2*a)^2 / (2*(gamma+1)*(gamma-1))
    ];

% F = @(rho, u, a, gamma, lamda) rho / (2 * gamma) * [
%     2 * (gamma - 1) * lamda(1) + lamda(2) + lamda(3);
%     2 * (gamma - 1) * lamda(1) * u + lamda(2) * (u + a) + lamda(3) * (u - a);
%     (gamma - 1) * lamda(1) * u.^2 + 1/2 * lamda(2) * (u + a).^2 + 1/2 * lamda(3) * (u - a).^2 + (3 - gamma) / (2 * (gamma - 1)) * (lamda(2) + lamda(3)) * a.^2
%     ];

% VL = @(a, M, gamma) [
%     1/4 * a .* (M + 1).^2 .* (1 - (M - 1).^2 / (gamma + 1));
%     1/4 * a .* (M + 1).^2 .* (3 - M + (gamma - 1) / (gamma + 1) * (M - 1).^2);
%     1/2 * a .* (M + 1).^2 .* (M - 1)/(gamma + 1) .* (1 + (gamma - 1)/2 * M)
%     ];

N = length(rho);
F_plus = zeros(3, N);
F_minus = zeros(3, N);

M = u ./ a;
for i = 1:N
    F_plus(:, i) = F(rho(i), u(i), a(i), gamma) * (abs(M(i)) >= 1) + Fmod_plus(rho(i), u(i), p(i), E(i), a(i), gamma) * (abs(M(i)) < 1);
    F_minus(:, i) = F(rho(i), u(i), a(i), gamma) * (abs(M(i)) >= 1) + Fmod_minus(rho(i), u(i), p(i), E(i), a(i), gamma) * (abs(M(i)) < 1);
end

end
