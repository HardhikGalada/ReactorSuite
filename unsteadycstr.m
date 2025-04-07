function unsteadycstr
    clc; clear; close all;

    % Create UI figure
    fig = figure('Name', 'Isothermal CSTR Simulation', 'NumberTitle', 'off', 'Position', [100, 100, 900, 500]);

    % UI Input Fields
    uicontrol('Style', 'text', 'Position', [20, 440, 200, 20], 'String', 'Reactor Volume (L):');
    V_input = uicontrol('Style', 'edit', 'Position', [220, 440, 100, 20], 'String', '100');

    uicontrol('Style', 'text', 'Position', [20, 410, 200, 20], 'String', 'Inlet Flow Rate (L/s):');
    Q_input = uicontrol('Style', 'edit', 'Position', [220, 410, 100, 20], 'String', '1');

    uicontrol('Style', 'text', 'Position', [20, 380, 200, 20], 'String', 'Inlet Concentration (mol/L):');
    C_in_input = uicontrol('Style', 'edit', 'Position', [220, 380, 100, 20], 'String', '1');

    uicontrol('Style', 'text', 'Position', [20, 350, 200, 20], 'String', 'Reaction Rate Constant (1/s):');
    k_input = uicontrol('Style', 'edit', 'Position', [220, 350, 100, 20], 'String', '0.1');

    uicontrol('Style', 'text', 'Position', [20, 320, 200, 20], 'String', 'Initial Concentration (mol/L):');
    C0_input = uicontrol('Style', 'edit', 'Position', [220, 320, 100, 20], 'String', '1');

    % Button to run the simulation
    uicontrol('Style', 'pushbutton', 'String', 'Run Simulation', 'Position', [100, 280, 150, 30], 'Callback', @run_simulation);

    % Axes for plotting results
    ax = axes('Position', [0.4, 0.15, 0.55, 0.75]);

    % Text panel for results
    results_panel = uicontrol('Style', 'text', 'Position', [20, 100, 300, 100], 'String', '', 'FontSize', 10, 'HorizontalAlignment', 'left');

    function run_simulation(~, ~)
        % Get user inputs
        V = str2double(get(V_input, 'String'));
        Q = str2double(get(Q_input, 'String'));
        C_in = str2double(get(C_in_input, 'String'));
        k = str2double(get(k_input, 'String'));
        C0 = str2double(get(C0_input, 'String'));

        % Check if Q = 0 (Batch Reactor Case)
        if Q == 0
            errordlg('Flow rate cannot be zero for a CSTR. Try a small nonzero value.', 'Input Error');
            return;
        end

        % Residence Time Calculation
        tau = V / Q; % Space time/Residence time (s)

        % Solve ODE
        t_span = [0 200]; % Simulation time span
        [t, C] = ode45(@(t, C) cstr_mass_balance(t, C, V, Q, C_in, k), t_span, C0);
        
        % Compute steady-state outlet concentration
        C_out = C(end);

        % Display results
        results_text = sprintf(['Steady-State Outlet Concentration: %.2f mol/L\n', ...
                                'Residence Time: %.1f s'], C_out, tau);
        set(results_panel, 'String', results_text);

        % Plot the concentration profile over time
        axes(ax);
        plot(t, C, 'b.-', 'LineWidth', 2);
        xlabel('Time (s)');
        ylabel('Concentration (mol/L)');
        title('Dynamic Simulation of an Isothermal CSTR');
        grid on;
        legend('C_A (mol/L)');
    end

    % Function defining the ODE for Isothermal CSTR
    function dCdt = cstr_mass_balance(~, C, V, Q, C_in, k)
        dCdt = (Q / V) * (C_in - C) - k * C;
    end
end
