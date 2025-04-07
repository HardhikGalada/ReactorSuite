function adiabaticcstr
    clc; clear; close all;

    % Create the main figure window
    fig = figure('Name', 'Adiabatic CSTR Simulation', 'NumberTitle', 'off', 'Position', [100, 100, 900, 500]);

    % UI Input Fields
    uicontrol('Style', 'text', 'Position', [20, 440, 200, 20], 'String', 'Initial Concentration (mol/m³):');
    C0_input = uicontrol('Style', 'edit', 'Position', [220, 440, 100, 20], 'String', '2000');

    uicontrol('Style', 'text', 'Position', [20, 410, 200, 20], 'String', 'Density (kg/m³):');
    rho_input = uicontrol('Style', 'edit', 'Position', [220, 410, 100, 20], 'String', '1000');

    uicontrol('Style', 'text', 'Position', [20, 380, 200, 20], 'String', 'Heat Capacity (J/(kg.K)):');
    cp_input = uicontrol('Style', 'edit', 'Position', [220, 380, 100, 20], 'String', '4000');

    uicontrol('Style', 'text', 'Position', [20, 350, 200, 20], 'String', 'Heat of Reaction (J/mol):');
    DELTAH_input = uicontrol('Style', 'edit', 'Position', [220, 350, 100, 20], 'String', '-3e5');

    uicontrol('Style', 'text', 'Position', [20, 320, 200, 20], 'String', 'Residence Time (sec):');
    tau_input = uicontrol('Style', 'edit', 'Position', [220, 320, 100, 20], 'String', '1000');

    uicontrol('Style', 'text', 'Position', [20, 290, 200, 20], 'String', 'Reference Temperature (K):');
    Tm_input = uicontrol('Style', 'edit', 'Position', [220, 290, 100, 20], 'String', '298');

    uicontrol('Style', 'text', 'Position', [20, 260, 200, 20], 'String', 'Rate Constant at Tm (1/sec):');
    km_input = uicontrol('Style', 'edit', 'Position', [220, 260, 100, 20], 'String', '0.001/60');

    uicontrol('Style', 'text', 'Position', [20, 230, 200, 20], 'String', 'Activation Energy (K):');
    E_input = uicontrol('Style', 'edit', 'Position', [220, 230, 100, 20], 'String', '8000');

    uicontrol('Style', 'text', 'Position', [20, 200, 200, 20], 'String', 'Feed Temperature (K):');
    Tf_input = uicontrol('Style', 'edit', 'Position', [220, 200, 100, 20], 'String', '298');

    % Button to run the simulation
    uicontrol('Style', 'pushbutton', 'String', 'Run Simulation', 'Position', [100, 160, 150, 30], 'Callback', @run_simulation);

    % Axes for plotting results
    ax1 = axes('Position', [0.4, 0.55, 0.55, 0.35]);
    ax2 = axes('Position', [0.4, 0.10, 0.55, 0.35]);

    % Text panel for results
    results_panel = uicontrol('Style', 'text', 'Position', [20, 100, 300, 40], 'String', '', 'FontSize', 10, 'HorizontalAlignment', 'left');

    function run_simulation(~, ~)
        % Get user inputs
        C0 = str2double(get(C0_input, 'String'));
        rho = str2double(get(rho_input, 'String'));
        cp = str2double(get(cp_input, 'String'));
        DELTAH = str2double(get(DELTAH_input, 'String'));
        tau = str2double(get(tau_input, 'String'));
        Tm = str2double(get(Tm_input, 'String'));
        km = str2double(get(km_input, 'String'));
        E = str2double(get(E_input, 'String'));
        Tf = str2double(get(Tf_input, 'String'));

        % Validate inputs
        if isnan(C0) || isnan(rho) || isnan(cp) || isnan(DELTAH) || isnan(tau) || isnan(Tm) || isnan(km) || isnan(E) || isnan(Tf)
            errordlg('Please enter valid numeric values for all parameters.', 'Input Error');
            return;
        end

        % Simulation parameters
        tspan = [0 20000];
        y0 = [1000, 340]; % Initial conditions: [C (mol/m^3), T (K)]

        % Solve ODE
        [t, y] = ode15s(@(t, y) CSTR_ODE_MODEL(t, y, C0, tau, km, E, Tm, Tf, rho, cp, DELTAH), tspan, y0);

        % Extract steady-state values
        C_ss = y(end, 1);
        T_ss = y(end, 2);

        % Display results in UI
        results_text = sprintf('Steady-State Concentration: %.1f mol/m³\nSteady-State Temperature: %.0f K', C_ss, T_ss);
        set(results_panel, 'String', results_text);

        % Plot Concentration
        axes(ax1);
        plot(t / 3600, y(:, 1), 'r.-', 'LineWidth', 1.5);
        xlabel('Time (hrs)');
        ylabel('C (mol/m³)');
        title('Concentration vs Time');
        grid on;

        % Plot Temperature
        axes(ax2);
        plot(t / 3600, y(:, 2), 'b.-', 'LineWidth', 1.5);
        xlabel('Time (hrs)');
        ylabel('T (K)');
        title('Temperature vs Time');
        grid on;
    end

    % Define the ODE system for the adiabatic CSTR
    function dydt = CSTR_ODE_MODEL(~, y, C0, tau, km, E, Tm, Tf, rho, cp, DELTAH)
        % Extract variables
        C = y(1);
        T = y(2);

        % Rate constant
        k = km * exp(-E * (1 / T - 1 / Tm)); % 1/s

        % Mass and energy balances
        dCdt = (1 / tau) * (C0 - C) - k * C;
        dTdt = (1 / tau) * (Tf - T) + (k / (rho * cp)) * C * (-DELTAH);

        dydt = [dCdt; dTdt];
    end
end
