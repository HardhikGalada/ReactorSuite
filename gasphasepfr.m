function gasphasepfr
    clc; clear; close all;

    % Create UI figure
    fig = figure('Name', 'Gas Phase PFR Simulation', 'NumberTitle', 'off', 'Position', [100, 100, 900, 500]);

    % UI Input Fields
    uicontrol('Style', 'text', 'Position', [20, 440, 200, 20], 'String', 'Reactor Length (m):');
    L_input = uicontrol('Style', 'edit', 'Position', [220, 440, 100, 20], 'String', '15');

    uicontrol('Style', 'text', 'Position', [20, 410, 200, 20], 'String', 'Reactor Diameter (m):');
    d_input = uicontrol('Style', 'edit', 'Position', [220, 410, 100, 20], 'String', '0.025');

    uicontrol('Style', 'text', 'Position', [20, 380, 200, 20], 'String', 'Reactor Pressure (Pa):');
    P_input = uicontrol('Style', 'edit', 'Position', [220, 380, 100, 20], 'String', '200000');

    uicontrol('Style', 'text', 'Position', [20, 350, 200, 20], 'String', 'Temperature (K):');
    T_input = uicontrol('Style', 'edit', 'Position', [220, 350, 100, 20], 'String', '750');

    uicontrol('Style', 'text', 'Position', [20, 320, 200, 20], 'String', 'Volumetric Flow Rate (m³/s):');
    Q_input = uicontrol('Style', 'edit', 'Position', [220, 320, 100, 20], 'String', '0.0005');

    % Button to run the simulation
    uicontrol('Style', 'pushbutton', 'String', 'Run Simulation', 'Position', [100, 280, 150, 30], 'Callback', @run_simulation);

    % Axes for plotting results
    ax = axes('Position', [0.4, 0.15, 0.55, 0.75]);

    % Text panel for results
    results_panel = uicontrol('Style', 'text', 'Position', [20, 100, 300, 100], 'String', '', 'FontSize', 10, 'HorizontalAlignment', 'left');

    function run_simulation(~, ~)
        % Get user inputs
        L = str2double(get(L_input, 'String'));
        d = str2double(get(d_input, 'String'));
        P = str2double(get(P_input, 'String'));
        T = str2double(get(T_input, 'String'));
        Q = str2double(get(Q_input, 'String'));

        % Constants
        R = 8.314;  % Ideal gas constant in J/(mol·K)
        A_c = pi * (d^2) / 4;  % Cross-sectional area (m²)
        A = 1.0e15;  % Pre-exponential factor in Arrhenius equation (s⁻¹)
        E_a = 3.0e5; % Activation energy (J/mol)
        k = A * T * exp(-E_a / (R * T)); % Rate constant (1/s): Modified Arrhenius expression

        % Initial molar flow rate of A using the ideal gas law
        F_A_in = P * Q / (R * T);  % Initial molar flow rate of A (mol/s)

        % Solve ODE
        L_span = [0, L];
        initial_conditions = F_A_in;

        [L_vals, solution] = ode15s(@(z, F_A) gas_phase_PFR(z, F_A, k, A_c, P, R, T, F_A_in), L_span, initial_conditions);
        F_A_values = solution;

        % Compute output parameters
        F_A_out = F_A_values(end);
        F_T_out = 2 * F_A_in - F_A_out; % Total molar flow rate at reactor outlet
        X_A = (F_A_in - F_A_out) / F_A_in; % Conversion of A
        Q_in = Q * 1000 * 60; % Inlet flow rate in L/min
        Q_out = (F_T_out * (R * T) / P) * 1000 * 60; % Outlet flow rate in L/min

        % Display results
        results_text = sprintf(['Inlet Molar Flow Rate of A: %.4f mol/s\n', ...
                                'Outlet Molar Flow Rate of A: %.4f mol/s\n', ...
                                'Conversion of A: %.2f\n', ...
                                'Inlet Volumetric Flow Rate: %.1f L/min\n', ...
                                'Outlet Volumetric Flow Rate: %.1f L/min'], ...
                                F_A_in, F_A_out, X_A, Q_in, Q_out);
        set(results_panel, 'String', results_text);

        % Plot the molar flow rate along the reactor length
        axes(ax);
        plot(L_vals, F_A_values, 'b-', 'LineWidth', 2, 'Marker', '.');
        xlabel('Reactor Length (m)');
        ylabel('Molar Flow Rate of A (mol/s)');
        title('Isothermal Gas Phase PFR with Volume Change');
        grid on;
        legend('F_A');
    end

    % Function defining the ODE for Gas Phase PFR
    function dF_A_dz = gas_phase_PFR(z, F_A, k, A_c, P, R, T, F_A_in)
        F_T = 2 * F_A_in - F_A; % Total molar flow rate
        dF_A_dz = -(A_c * k * P / (R * T)) * (F_A / F_T);
    end
end
