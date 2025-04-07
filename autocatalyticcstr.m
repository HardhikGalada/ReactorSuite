function autocatalyticcstr
    clc; clear; close all;

    % Create UI figure
    fig = figure('Name', 'Autocatalytic CSTR Simulation', 'NumberTitle', 'off', 'Position', [100, 100, 900, 500]);

    % UI Input Fields
    uicontrol('Style', 'text', 'Position', [20, 440, 200, 20], 'String', 'Flow Rate (L/min):');
    Q_input = uicontrol('Style', 'edit', 'Position', [220, 440, 100, 20], 'String', '1.0');

    uicontrol('Style', 'text', 'Position', [20, 410, 200, 20], 'String', 'Reactor Volume (L):');
    V_input = uicontrol('Style', 'edit', 'Position', [220, 410, 100, 20], 'String', '3.0');

    uicontrol('Style', 'text', 'Position', [20, 380, 200, 20], 'String', 'Rate Constant (L/mol/min):');
    k_input = uicontrol('Style', 'edit', 'Position', [220, 380, 100, 20], 'String', '1.0');

    uicontrol('Style', 'text', 'Position', [20, 350, 200, 20], 'String', 'Inlet Concentration of A (mol/L):');
    CA_in_input = uicontrol('Style', 'edit', 'Position', [220, 350, 100, 20], 'String', '1.0');

    uicontrol('Style', 'text', 'Position', [20, 320, 200, 20], 'String', 'Inlet Concentration of R (mol/L):');
    CR_in_input = uicontrol('Style', 'edit', 'Position', [220, 320, 100, 20], 'String', '0.0');

    uicontrol('Style', 'text', 'Position', [20, 290, 200, 20], 'String', 'Initial Concentration of A (mol/L):');
    CA0_input = uicontrol('Style', 'edit', 'Position', [220, 290, 100, 20], 'String', '0.5');

    uicontrol('Style', 'text', 'Position', [20, 260, 200, 20], 'String', 'Initial Concentration of R (mol/L):');
    CR0_input = uicontrol('Style', 'edit', 'Position', [220, 260, 100, 20], 'String', '0.1');

    % Button to run the simulation
    uicontrol('Style', 'pushbutton', 'String', 'Run Simulation', 'Position', [100, 220, 150, 30], 'Callback', @run_simulation);

    % Axes for plotting results
    ax1 = axes('Position', [0.4, 0.55, 0.55, 0.35]);
    ax2 = axes('Position', [0.4, 0.10, 0.55, 0.35]);

    % Text panel for results
    results_panel = uicontrol('Style', 'text', 'Position', [20, 150, 300, 50], 'String', '', 'FontSize', 10, 'HorizontalAlignment', 'left');

    function run_simulation(~, ~)
        % Get user inputs
        Q = str2double(get(Q_input, 'String'));
        V = str2double(get(V_input, 'String'));
        k = str2double(get(k_input, 'String'));
        CA_in = str2double(get(CA_in_input, 'String'));
        CR_in = str2double(get(CR_in_input, 'String'));
        CA0 = str2double(get(CA0_input, 'String'));
        CR0 = str2double(get(CR0_input, 'String'));

        % Validate inputs
        if any(isnan([Q, V, k, CA_in, CR_in, CA0, CR0]))
            errordlg('Please enter valid numeric values for all parameters.', 'Input Error');
            return;
        end

        % Simulation time span
        tspan = [0, 50];
        y0 = [CA0, CR0];  % Initial conditions

        % Solve ODE system
        [t, y] = ode15s(@(t, y) cstr_autocatalytic(t, y, Q, V, k, CA_in, CR_in), tspan, y0);

        % Extract outlet values
        CA_steady_state = y(end, 1);
        CR_steady_state = y(end, 2);

        % Display results in UI
        results_text = sprintf('Steady-state Concentration of A: %.2f mol/L\nSteady-state Concentration of R: %.2f mol/L', CA_steady_state, CR_steady_state);
        set(results_panel, 'String', results_text);

        % Plot Concentration of A vs Time
        axes(ax1);
        plot(t, y(:, 1), 'r.-', 'LineWidth', 1.5);
        xlabel('Time (min)');
        ylabel('C_A (mol/L)');
        title('Concentration of A vs Time');
        grid on;

        % Plot Concentration of R vs Time
        axes(ax2);
        plot(t, y(:, 2), 'b.-', 'LineWidth', 1.5);
        xlabel('Time (min)');
        ylabel('C_R (mol/L)');
        title('Concentration of R vs Time');
        grid on;
    end

    % Define the ODE system for CSTR with Autocatalytic Reaction
    function dydt = cstr_autocatalytic(~, y, F, V, k, CA_in, CR_in)
        CA = y(1);
        CR = y(2);
        dCA_dt = (F / V) * (CA_in - CA) - k * CA * CR;
        dCR_dt = (F / V) * (CR_in - CR) + k * CA * CR;
        dydt = [dCA_dt; dCR_dt];
    end
end
