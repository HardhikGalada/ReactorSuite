function steadypfr
    clc; clear; close all;

    % Create UI figure
    fig = figure('Name', 'Isothermal PFR Simulation', 'NumberTitle', 'off', 'Position', [100, 100, 900, 500]);

    % UI Input Fields
    uicontrol('Style', 'text', 'Position', [20, 440, 200, 20], 'String', 'Flow Rate (m³/s):');
    Q_input = uicontrol('Style', 'edit', 'Position', [220, 440, 100, 20], 'String', '0.001');

    uicontrol('Style', 'text', 'Position', [20, 410, 200, 20], 'String', 'Reaction Rate Constant (1/s):');
    k_input = uicontrol('Style', 'edit', 'Position', [220, 410, 100, 20], 'String', '0.1');

    uicontrol('Style', 'text', 'Position', [20, 380, 200, 20], 'String', 'Inlet Concentration (mol/L):');
    C_in_input = uicontrol('Style', 'edit', 'Position', [220, 380, 100, 20], 'String', '1');

    uicontrol('Style', 'text', 'Position', [20, 350, 200, 20], 'String', 'Reactor Length (m):');
    L_input = uicontrol('Style', 'edit', 'Position', [220, 350, 100, 20], 'String', '15');

    uicontrol('Style', 'text', 'Position', [20, 320, 200, 20], 'String', 'Reactor Diameter (m):');
    d_input = uicontrol('Style', 'edit', 'Position', [220, 320, 100, 20], 'String', '0.05');

    % Button to run the simulation
    uicontrol('Style', 'pushbutton', 'String', 'Run Simulation', 'Position', [100, 280, 150, 30], 'Callback', @run_simulation);

    % Axes for plotting results
    ax = axes('Position', [0.4, 0.15, 0.55, 0.75]);

    % Text panel for results
    results_panel = uicontrol('Style', 'text', 'Position', [20, 100, 300, 100], 'String', '', 'FontSize', 10, 'HorizontalAlignment', 'left');

    function run_simulation(~, ~)
        % Get user inputs
        Q = str2double(get(Q_input, 'String'));
        k = str2double(get(k_input, 'String'));
        C_in = str2double(get(C_in_input, 'String'));
        L = str2double(get(L_input, 'String'));
        d = str2double(get(d_input, 'String'));

        % Compute reactor properties
        A = (pi * d^2) / 4;  % Cross-sectional area of reactor (m²)
        V = A * L * 1000;  % Reactor volume (L)
        tau = (V / 1000) / Q;  % Residence time (s)

        % Solve ODE
        z_span = [0 L];  % Reactor length range
        C0 = C_in;  % Initial concentration
        [z, C] = ode45(@(z, C) pfr_mole_balance(z, C, k, Q, A), z_span, C0);
        
        % Compute outlet concentration
        C_out = C(end);

        % Display results
        results_text = sprintf(['Reactor Volume: %.0f L\n', ...
                                'Residence Time: %.0f s\n', ...
                                'Outlet Concentration: %.2f mol/L'], V, tau, C_out);
        set(results_panel, 'String', results_text);

        % Plot the concentration profile over reactor length
        axes(ax);
        plot(z, C, 'b.-', 'LineWidth', 2);
        xlabel('Reactor Length (m)');
        ylabel('Concentration of A (mol/L)');
        title('Steady-State Concentration Profile in a PFR');
        grid on;
        legend('C_A (mol/L)');
    end

    % Function defining the ODE for Isothermal PFR
    function dCdz = pfr_mole_balance(~, C, k, Q, A)
        dCdz = -k * C / (Q / A);
    end
end
