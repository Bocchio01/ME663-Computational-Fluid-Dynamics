%% Rocket physics modelling

funcs.rho       = @(s)    1 / ((8.314 / 28.96 * 1000) * (273.15 + 15)) * 101325 * exp(-(9.81 * 28.96 / 1000 * (s)) / (8.31 * (273.15 + 15)));
funcs.mass      = @(t)    ((rocket.m(end) - rocket.m(1)) / (thrust_curve(end).Time - thrust_curve(1).Time) * t + rocket.m(1)) .* (t <= thrust_curve(end).Time) + rocket.m(end) * (t > thrust_curve(end).Time);
funcs.F_gravity = @(t)    9.81 * funcs.mass(t);
funcs.F_drag    = @(s, v) 1/2 * rocket.Cd * funcs.rho(s) .* rocket.A .* v.^2;
funcs.F_thrust  = @(t)    interp1([thrust_curve.Time 100], [thrust_curve.Thrust 0] * thrust_reduction_coefficient, t, "pchip");


%% Solution of the ODE for the Equation of Motion of the system

[t, z] = ode45(@(t, z) EquationOfMotion(t, z, funcs), ...
    [0 timestamp(end)], ...
    [funcs.velocity(0); funcs.altitude(0)], ...
    odeset('RelTol', 1e-6, 'AbsTol', 1e-6));


%% Functions

function [z_dot] = EquationOfMotion(t, z, funcs)

velocity = z(1,1);
altitude = z(2,1);

acceleration = (funcs.F_thrust(t) - funcs.F_gravity(t) - funcs.F_drag(altitude, velocity)) / funcs.mass(t);

z_dot(1,1) = acceleration;
z_dot(2,1) = velocity;

end

