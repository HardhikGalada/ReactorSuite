function seriescstrs
    clc; clear; close all;

    % Create UI figure
    fig = figure('Name', 'CSTRs in Series vs. Ideal PFR', 'NumberTitle', 'off', 'Position', [100, 100, 900, 500]);

    % UI Input Fields
    uicontrol('Style', 'text', 'Position', [20, 440, 200, 20], 'String', 'Number of CSTRs:');
    n_cstrs_input = uicontrol('Style', 'edit', 'Position', [220, 440, 100, 20], 'String', '5');

    uicontrol('Style', 'text', 'Position', [20, 410, 200, 20], 'String', 'Reactor Volume (L):');
    V_input = uicontrol('Style', 'edit', 'Position', [220, 410, 100, 20], 'String', '10.0');

    uicontrol('Style', 'text', 'Position', [20, 380, 200, 20], 'String', 'Volumetric Flow Rate (L/s):');
    Q_input = uicontrol('Style', 'edit', 'Position', [220, 380, 100, 20], 'String', '1.0');

    uicontrol('Style', 'text', 'Position', [20, 350, 200, 20], 'String', 'Reaction Rate Constant (L/(mol.s)):');
    k_input = uicontrol('Style', 'edit', 'Position', [220, 350, 100, 20], 'String', '0.5');

    uicontrol('Style', 'text', 'Position', [20, 320, 200, 20], 'String', 'Initial Concentration of A (mol/L):');
    CA0_input = uicontrol('Style', 'edit', 'Position', [220, 320, 100, 20], 'String', '1.0');

    uicontrol('Style', 'text', 'Position', [20, 290, 200, 20], 'String', 'Initial Concentration of B (mol/L):');
    CB0_input = uicontrol('Style', 'edit', 'Position', [220, 290, 100, 20], 'String', '1.0');

    uicontrol('Style', 'text', 'Position', [20, 260, 200, 20], 'String', 'Initial Concentration of C (mol/L):');
    CC0_input = uicontrol('Style', 'edit', 'Position', [220, 260, 100, 20], 'String', '0.0');

    % Button to run the simulation
    uicontrol('Style', 'pushbutton', 'String', 'Run Simulation', 'Position', [100, 220, 150, 30], 'Callback', @run_simulation);

    % Axes for plotting results
    ax = axes('Position', [0.4, 0.15, 0.55, 0.75]);

    % Text panel for results
    results_panel = uicontrol('Style', 'text', 'Position', [20, 100, 300, 100], 'String', '', 'FontSize', 10, 'HorizontalAlignment', 'left');

    function run_simulation(~, ~)
        % Get user inputs
        n_cstrs = str2double(get(n_cstrs_input, 'String'));
        V = str2double(get(V_input, 'String'));
        Q = str2double(get(Q_input, 'String'));
        k = str2double(get(k_input, 'String'));
        CA0 = str2double(get(CA0_input, 'String'));
        CB0 = str2double(get(CB0_input, 'String'));
        CC0 = str2double(get(CC0_input, 'String'));

        % Validate inputs
        if any(isnan([n_cstrs, V, Q, k, CA0, CB0, CC0])) || n_cstrs < 1
            errordlg('Please enter valid numeric values for all parameters.', 'Input Error');
            return;
        end

        % Simulate the n CSTRs in series
        tspan = [0 200];
        initial_conditions = [ones(n_cstrs, 1) * CA0; ones(n_cstrs, 1) * CB0; ones(n_cstrs, 1) * CC0];

        [~, solution] = ode15s(@(t, C) cstr_odes(t, C, V, Q, k, n_cstrs, CA0, CB0, CC0), tspan, initial_conditions);

        CA_steady_state = solution(end, 1:n_cstrs);
        CC_steady_state = solution(end, 2 * n_cstrs + 1:end);

        % Simulate the ideal PFR
        V_span = [0, V * n_cstrs];
        initial_conditions_pfr = [CA0, CB0, CC0];
        [~, pfr_solution] = ode15s(@(V, C) pfr_odes(V, C, k, Q), V_span, initial_conditions_pfr);
        pfr_outlet_concentrations = pfr_solution(end, :);

        % Display results for CSTRs and PFR
        results_text = "Steady-State Outlet Concentrations:\n";
        for i = 1:n_cstrs
            results_text = results_text + sprintf('CSTR %d: A = %.2f mol/L, C = %.2f mol/L\n', i, CA_steady_state(i), CC_steady_state(i));
        end
        results_text = results_text + sprintf('\nIdeal PFR: A = %.2f mol/L, C = %.2f mol/L\n', pfr_outlet_concentrations(1), pfr_outlet_concentrations(3));
        set(results_panel, 'String', results_text);

        % Plot outlet concentrations
        axes(ax);
        plot(1:n_cstrs, CA_steady_state, 'bo-', 'LineWidth', 2, 'MarkerSize', 5);
        hold on;
        plot(1:n_cstrs, CC_steady_state, 'ro-', 'LineWidth', 2, 'MarkerSize', 5);
        yline(pfr_outlet_concentrations(1), '--r', 'PFR A');
        yline(pfr_outlet_concentrations(3), '--b', 'PFR C');
        hold off;
        xlabel('CSTR Number');
        ylabel('Outlet Concentration (mol/L)');
        title('CSTRs in Series vs. PFR');
        legend('CSTR [A]', 'CSTR [C]', 'PFR [A]', 'PFR [C]');
        grid on;
    end

    % Function defining the ODEs for CSTRs
    function dCdt = cstr_odes(~, concentrations, V, Q, k, n_cstrs, CA0, CB0, CC0)
        CA = concentrations(1:n_cstrs);
        CB = concentrations(n_cstrs+1:2*n_cstrs);
        CC = concentrations(2*n_cstrs+1:end);

        dCAdt = zeros(n_cstrs, 1);
        dCBdt = zeros(n_cstrs, 1);
        dCCdt = zeros(n_cstrs, 1);

        for i = 1:n_cstrs
            if i == 1
                inflow_A = CA0;
                inflow_B = CB0;
                inflow_C = CC0;
            else
                inflow_A = CA(i-1);
                inflow_B = CB(i-1);
                inflow_C = CC(i-1);
            end

            r = k * CA(i) * CB(i);
            dCAdt(i) = Q/V * (inflow_A - CA(i)) - r;
            dCBdt(i) = Q/V * (inflow_B - CB(i)) - r;
            dCCdt(i) = Q/V * (inflow_C - CC(i)) + r;
        end

        dCdt = [dCAdt; dCBdt; dCCdt];
    end

    % Function defining the ODE for the PFR
    function dCdt = pfr_odes(~, C, k, Q)
        r = k * C(1) * C(2);
        dCdt = [-r / Q; -r / Q; r / Q];
    end
end
