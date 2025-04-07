function nonisothermalpfr
    clc; clear; close all;

    % Create UI figure
    fig = figure('Name', 'HOTSPOT01 Reactor Simulation', 'NumberTitle', 'off', 'Position', [100, 100, 900, 500]);

    % UI Input Fields
    uicontrol('Style', 'text', 'Position', [20, 440, 200, 20], 'String', 'Feed Temperature (K):');
    Tf_input = uicontrol('Style', 'edit', 'Position', [220, 440, 100, 20], 'String', '600');

    uicontrol('Style', 'text', 'Position', [20, 410, 200, 20], 'String', 'Initial Molar Flow Rate (mol/s):');
    F_Af_input = uicontrol('Style', 'edit', 'Position', [220, 410, 100, 20], 'String', '1.5');

    % Button to run the simulation
    uicontrol('Style', 'pushbutton', 'String', 'Run Simulation', 'Position', [100, 370, 150, 30], 'Callback', @run_simulation);

    % Axes for plotting results
    ax1 = axes('Position', [0.4, 0.55, 0.55, 0.35]);
    ax2 = axes('Position', [0.4, 0.10, 0.55, 0.35]);

    % Text panel for results
    results_panel = uicontrol('Style', 'text', 'Position', [20, 300, 300, 50], 'String', '', 'FontSize', 10, 'HorizontalAlignment', 'left');

    function run_simulation(~, ~)
        % Get user inputs
        Tf = str2double(get(Tf_input, 'String'));
        F_Af = str2double(get(F_Af_input, 'String'));

        % Validate inputs
        if isnan(Tf) || isnan(F_Af)
            errordlg('Please enter valid numeric values for all parameters.', 'Input Error');
            return;
        end

        % Simulation parameters
        lspan = [0, 1.5]; % Reactor length in meters
        xo = [F_Af, Tf]; % Initial conditions

        % Solve ODE system
        [l, x] = ode15s(@(l, x) HOTSPOT01(l, x), lspan, xo);

        % Extract outlet values
        outlet_temperature = x(end, 2);
        outlet_molar_flow = x(end, 1);

        % Display results in UI
        results_text = sprintf('Outlet Reactor Temperature: %.0f K\nOutlet Molar Flow Rate: %.2f mol/s', outlet_temperature, outlet_molar_flow);
        set(results_panel, 'String', results_text);

        % Plot Molar Flow Rate vs Length
        axes(ax1);
        plot(l, x(:, 1), 'k.-', 'LineWidth', 1.5);
        xlabel('Length (m)');
        ylabel('F (mol/s)');
        title('Molar Flow Rate of A vs Length');
        grid on;

        % Plot Temperature vs Length
        axes(ax2);
        plot(l, x(:, 2), 'b.-', 'LineWidth', 1.5);
        xlabel('Length (m)');
        ylabel('T (K)');
        title('Temperature vs Length');
        grid on;
    end

    % Define the ODE system for HOTSPOT01
    function dxdl = HOTSPOT01(~, x)
        % Parameters
        km = 1923;  % 1/s
        E = 13636;  % K (E/R)
        Tm = 625;  % K
        DELTAH = -1361;  % J/mol

        L = 1.5;  % m
        RT = 0.0125;  % m
        A = pi * RT^2;  % Cross-sectional area (m²)

        h = 373;  % W/(m².K)
        Tc = 625;  % K

        Mdot = 2.6371e-3;  % kg/s
        CP = 992;  % J/(kg.K)

        P = 1e5;  % Pa
        R = 8.314;  % J/(mol.K)
        y_Af = 0.019;

        F_Af = 1.5;  % mol/s
        F_T = F_Af / y_Af;  % mol/s

        % Reaction rate constant
        k = km * exp(-E * (1 / x(2) - 1 / Tm));

        % Reaction rate
        RATE = k * (x(1) * P / (F_T * R * x(2)));

        % Model equations
        dNdl = -A * RATE;
        dTdl = -(DELTAH * A * RATE / (Mdot * CP)) - (h * 2 * pi * RT * (x(2) - Tc) / (Mdot * CP));

        dxdl = [dNdl; dTdl];
    end
end
