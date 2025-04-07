function recyclepfr
    clc; clear; close all;

    % Create UI figure
    fig = figure('Name', 'PFR with Recycle', 'NumberTitle', 'off', 'Position', [100, 100, 900, 500]);

    % UI Input Fields
    uicontrol('Style', 'text', 'Position', [20, 440, 200, 20], 'String', 'Rate Constant (L/mol/min):');
    k_input = uicontrol('Style', 'edit', 'Position', [220, 440, 100, 20], 'String', '1.0');

    uicontrol('Style', 'text', 'Position', [20, 410, 200, 20], 'String', 'Reactor Volume (L):');
    V_input = uicontrol('Style', 'edit', 'Position', [220, 410, 100, 20], 'String', '3.0');

    uicontrol('Style', 'text', 'Position', [20, 380, 200, 20], 'String', 'Fresh Feed Flow Rate (L/min):');
    Q0_input = uicontrol('Style', 'edit', 'Position', [220, 380, 100, 20], 'String', '1.0');

    uicontrol('Style', 'text', 'Position', [20, 350, 200, 20], 'String', 'Inlet Concentration of A (mol/L):');
    CA0_input = uicontrol('Style', 'edit', 'Position', [220, 350, 100, 20], 'String', '1.0');

    uicontrol('Style', 'text', 'Position', [20, 320, 200, 20], 'String', 'Inlet Concentration of R (mol/L):');
    CR0_input = uicontrol('Style', 'edit', 'Position', [220, 320, 100, 20], 'String', '0.0');

    % Button to run the simulation
    uicontrol('Style', 'pushbutton', 'String', 'Run Simulation', 'Position', [100, 280, 150, 30], 'Callback', @run_simulation);

    % Axes for plotting results
    ax = axes('Position', [0.4, 0.15, 0.55, 0.75]);

    % Text panel for results
    results_panel = uicontrol('Style', 'text', 'Position', [20, 150, 300, 100], 'String', '', 'FontSize', 10, 'HorizontalAlignment', 'left');

    function run_simulation(~, ~)
        % Get user inputs
        k = str2double(get(k_input, 'String'));
        V = str2double(get(V_input, 'String'));
        Q0 = str2double(get(Q0_input, 'String'));
        CA0 = str2double(get(CA0_input, 'String'));
        CR0 = str2double(get(CR0_input, 'String'));

        % Validate inputs
        if any(isnan([k, V, Q0, CA0, CR0]))
            errordlg('Please enter valid numeric values for all parameters.', 'Input Error');
            return;
        end

        % Define range of recycle ratios
        recycle_ratios = linspace(0, 20, 50);
        CA_out_values = zeros(size(recycle_ratios));
        CR_out_values = zeros(size(recycle_ratios));

        % Compute outlet concentrations for each recycle ratio
        for i = 1:length(recycle_ratios)
            Rr = recycle_ratios(i);
            [CA_out_values(i), CR_out_values(i)] = compute_outlet_concentrations(Rr, k, V, Q0, CA0, CR0);
        end

        % Display results for selected recycle ratios
        recycle_ratios_to_print = [5, 10, 20];
        results_text = "Steady-State Outlet Concentrations:\n";
        for Rr = recycle_ratios_to_print
            [CA_out, CR_out] = compute_outlet_concentrations(Rr, k, V, Q0, CA0, CR0);
            results_text = results_text + sprintf('Rr = %d: A = %.2f mol/L, R = %.2f mol/L\n', Rr, CA_out, CR_out);
        end
        set(results_panel, 'String', results_text);

        % Plot outlet concentrations vs. recycle ratio
        axes(ax);
        plot(recycle_ratios, CA_out_values, 'bo-', 'LineWidth', 2, 'MarkerSize', 5);
        hold on;
        plot(recycle_ratios, CR_out_values, 'ro-', 'LineWidth', 2, 'MarkerSize', 5);
        hold off;
        xlabel('Recycle Ratio (Rr)');
        ylabel('Outlet Concentration (mol/L)');
        title('Effect of Recycle Ratio on PFR Outlet Concentrations');
        legend('A (Reactant)', 'R (Autocatalyst/Product)');
        grid on;
    end

    % Function to compute outlet concentrations for a given recycle ratio
    function [CA_out, CR_out] = compute_outlet_concentrations(Rr, k, V, Q0, CA0, CR0)
        tol = 1e-6;    % Convergence tolerance
        max_iter = 500; % Max iterations for steady-state convergence
        A_guess = 1.0;  % Initial guess for A concentration
        R_guess = 1.0;  % Initial guess for R concentration

        Qr = Rr * Q0;  % Recycle flow rate
        Qin = Q0 + Qr; % Total flow rate entering the reactor

        % Initial guess for outlet concentrations
        CA_out = A_guess;
        CR_out = R_guess;

        for iter = 1:max_iter
            % Compute new inlet concentrations using a mole balance at the merging node
            CA_in = (Q0 * CA0 + Qr * CA_out) / Qin;
            CR_in = (Q0 * CR0 + Qr * CR_out) / Qin;

            % Solve the PFR equations using ode15s
            [~, sol] = ode15s(@(V, C) pfr_reactor(V, C, k, Qin), [0, V], [CA_in, CR_in]);

            % New outlet concentrations
            CA_new_out = sol(end, 1);
            CR_new_out = sol(end, 2);

            % Check for convergence
            if abs(CA_new_out - CA_out) < tol && abs(CR_new_out - CR_out) < tol
                CA_out = CA_new_out;
                CR_out = CR_new_out;
                break;
            end

            % Update concentrations
            CA_out = CA_new_out;
            CR_out = CR_new_out;
        end
    end

    % Function defining the PFR differential equations
    function dC_dV = pfr_reactor(~, C, k, Q)
        CA = C(1);
        CR = C(2);
        r = k * CA * CR;  % Autocatalytic reaction rate
        dCA_dV = -r / Q;
        dCR_dV = r / Q;
        dC_dV = [dCA_dV; dCR_dV];
    end
end
