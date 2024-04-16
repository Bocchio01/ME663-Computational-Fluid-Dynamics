clc
clear variables
close all

%% Section 0: Settings & Parameters

rocket = struct( ...
    'A', (34e-3/2)^2 * pi, ...
    'm', [0.316 0.224], ...
    'Cd', 0.86034355);

output = struct('plot', true);

funcs = struct();

thrust_reduction_coefficient = 1.0;


%% Section 1: Data import

flight_data = table2struct(readtable('Flight data.xlsx'));
thrust_curve = table2struct(readtable('Thrust curve.xlsx'));
openrocket = table2struct(readtable('OpenRocket.csv'));


%% Section 2: Flight data interpolation

liftoff_idx = 3;
[~, apogee_idx] = max([flight_data.Altitude]);

timestamp = [flight_data(liftoff_idx:apogee_idx).Time_s_] - flight_data(liftoff_idx).Time_s_;
altitude = [flight_data(liftoff_idx:apogee_idx).Altitude];

altitude_poly_coefficient = polyfit(timestamp, altitude, 7);
velocity_poly_coefficient = polyder(altitude_poly_coefficient);
acceleration_poly_coefficient = polyder(velocity_poly_coefficient);

funcs.altitude = @(t) polyval(altitude_poly_coefficient, t);
funcs.velocity = @(t) polyval(velocity_poly_coefficient, t);
funcs.acceleration = @(t) polyval(acceleration_poly_coefficient, t);


%% Section 2: Rocket physics modelling & Thrust data interpolation

funcs.rho       = @(s)    1 / ((8.314 / 28.96 * 1000) * (273.15 + 15)) * 101325 * exp(-(9.81 * 28.96 / 1000 * (s)) / (8.31 * (273.15 + 15)));
funcs.mass      = @(t)    ((rocket.m(end) - rocket.m(1)) / (thrust_curve(end).Time - thrust_curve(1).Time) * t + rocket.m(1)) .* (t <= thrust_curve(end).Time) + ...
    rocket.m(end) * (t > thrust_curve(end).Time);
funcs.F_gravity = @(t)    9.81 * funcs.mass(t);
funcs.F_drag    = @(s, v) 1/2 * rocket.Cd * funcs.rho(s) .* rocket.A .* v.^2;
funcs.F_thrust  = @(t)    interp1([thrust_curve.Time 100], [thrust_curve.Thrust 0] * thrust_reduction_coefficient, t, "pchip");


%% Section 3: Solution of the ODE for the Equation of Motion of the system
% The EoM is based on the anonymus functions defined above

[t, z] = ode45(@(t, z) EquationOfMotion(t, z, funcs), ...
    [0 timestamp(end)], ...
    [0; funcs.altitude(0)], ...
    odeset('RelTol', 1e-6, 'AbsTol', 1e-6));


%% Plots

[~, apogee_idx] = max([openrocket.Altitude_m_]);
OR_time = [openrocket(1:apogee_idx).x_Time_s_];
OR_altitude = [openrocket(1:apogee_idx).Altitude_m_] + altitude(1);
OR_velocity = [openrocket(1:apogee_idx).VerticalVelocity_m_s_];

if (output.plot)

    figure_data_validation = figure('Name', 'Data validation', 'NumberTitle', 'off', 'Position', [100 200 1200 600]);
    tile = tiledlayout(1, 2);
    title(tile, ['C_d = ' num2str(rocket.Cd) ' | Thrust @' num2str(thrust_reduction_coefficient * 100) '% nominal']);

    % Altitude vs. time
    nexttile(1);
    hold on
    grid on

    plot(t, z(:, 2), ...
        OR_time, OR_altitude, ...
        t, funcs.altitude(t), ...
        'LineWidth', 2);

    xlabel('Time [s]');
    ylabel('Height above mean sea level [m]')
    title('Altitude analysis');
    legend('From Ansys CFD', 'From OpenRocket', 'Interpolated flight data', 'Location', 'best');


    % Velocity vs. time
    nexttile(2);
    hold on
    grid on

    plot(t, z(:, 1), ...
        OR_time, OR_velocity, ...
        t, funcs.velocity(t), ...
        'LineWidth', 2);

    xlabel('Time [s]');
    ylabel('Velocity [m/s]')
    title('Velocity analysis');
    legend('From Ansys CFD', 'From OpenRocket', 'Interpolated flight data', 'Location', 'best');

end



%% Functions

function [z_dot] = EquationOfMotion(t, z, funcs)

velocity = z(1,1);
altitude = z(2,1);

acceleration = (funcs.F_thrust(t) - funcs.F_gravity(t) - funcs.F_drag(altitude, velocity)) / funcs.mass(t);

z_dot(1,1) = acceleration;
z_dot(2,1) = velocity;

end

