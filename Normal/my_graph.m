classdef my_graph
    methods(Static)
        %% Plot wealth distribution by skill type (Partial Equilibrium).
        function [] = plot_dist(par, sol, sim)
            figure(1)
            set(gcf, 'Position', [100, 100, 1200, 400])

            asim_middle = sim.asim(2, :); % Middle-aged savings
            skill_assign = par.skill_assign;

            subplot(1, 2, 1)
            histogram(asim_middle(skill_assign == 1), 20, 'Normalization', 'probability')
            xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 12)
            ylabel('Probability', 'Interpreter', 'latex', 'FontSize', 12)
            title('Wealth Distribution (Low-Skill, Middle-Aged)', 'Interpreter', 'latex', 'FontSize', 14)

            subplot(1, 2, 2)
            histogram(asim_middle(skill_assign == 2), 20, 'Normalization', 'probability')
            xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 12)
            ylabel('Probability', 'Interpreter', 'latex', 'FontSize', 12)
            title('Wealth Distribution (High-Skill, Middle-Aged)', 'Interpreter', 'latex', 'FontSize', 14)

            saveas(gcf, 'wealth_distribution_pe.png')
        end

        %% Plot wealth distribution by skill type (General Equilibrium).
        function [] = plot_dist_ge(sim_ge, par_ge, sol_ge)
            figure(2)
            set(gcf, 'Position', [100, 100, 1200, 400])

            asim_middle = sim_ge.asim(2, :); % Middle-aged savings
            skill_assign = par_ge.skill_assign;

            subplot(1, 2, 1)
            histogram(asim_middle(skill_assign == 1), 20, 'Normalization', 'probability')
            xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 12)
            ylabel('Probability', 'Interpreter', 'latex', 'FontSize', 12)
            title('Wealth Distribution (Low-Skill, Middle-Aged)', 'Interpreter', 'latex', 'FontSize', 14)

            subplot(1, 2, 2)
            histogram(asim_middle(skill_assign == 2), 20, 'Normalization', 'probability')
            xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 12)
            ylabel('Probability', 'Interpreter', 'latex', 'FontSize', 12)
            title('Wealth Distribution (High-Skill, Middle-Aged)', 'Interpreter', 'latex', 'FontSize', 14)

            saveas(gcf, 'wealth_distribution_pe.png')
        end

        %% Plot consumption policy functions by age and skill (Partial Equilibrium).
        function [] = cfun(par, sol, sim) 
            figure(3)
            set(gcf, 'Position', [100, 100, 1200, 800])

            agrid = par.agrid;
            zgrid = par.zgrid;
            z_median_ind = ceil(length(zgrid) / 2); % Median productivity state
            cpol_young = squeeze(sol.c2(:, z_median_ind, 1, :)); % Young consumption (alen, slen)
            cpol_middle = squeeze(sol.c(:, z_median_ind, 1, :)); % Middle-aged consumption
            cpol_old = squeeze(sol.c3(:, z_median_ind, 1, :)); % Old consumption

            for s = 1:2 % 1 = low-skill, 2 = high-skill
                % Young
                subplot(3, 2, (1 - 1) * 2 + s)
                plot(agrid, cpol_young(:, s), 'LineWidth', 2)
                xlabel('$a_t$', 'Interpreter', 'latex', 'FontSize', 12)
                ylabel('$c_t$', 'Interpreter', 'latex', 'FontSize', 12)
                title(sprintf('Young Consumption (Skill: %s)', ...
                    my_graph.ifelse(s == 1, 'Low', 'High')), 'Interpreter', 'latex', 'FontSize', 14)
                grid on

                % Middle-aged
                subplot(3, 2, (2 - 1) * 2 + s)
                plot(agrid, cpol_middle(:, s), 'LineWidth', 2)
                xlabel('$a_t$', 'Interpreter', 'latex', 'FontSize', 12)
                ylabel('$c_t$', 'Interpreter', 'latex', 'FontSize', 12)
                title(sprintf('Middle-Aged Consumption (Skill: %s)', ...
                    my_graph.ifelse(s == 1, 'Low', 'High')), 'Interpreter', 'latex', 'FontSize', 14)
                grid on

                % Old
                subplot(3, 2, (3 - 1) * 2 + s)
                plot(agrid, cpol_old(:, s), 'LineWidth', 2)
                xlabel('$a_t$', 'Interpreter', 'latex', 'FontSize', 12)
                ylabel('$c_t$', 'Interpreter', 'latex', 'FontSize', 12)
                title(sprintf('Old Consumption (Skill: %s)', ...
                    my_graph.ifelse(s == 1, 'Low', 'High')), 'Interpreter', 'latex', 'FontSize', 14)
                grid on
            end

            saveas(gcf, 'consumption_policy_functions_pe.png')
        end

        %% Plot value functions by age and skill (Partial Equilibrium).
        function [] = vfun(par, sol, sim)
            figure(4)
            set(gcf, 'Position', [100, 100, 1200, 800])

            agrid = par.agrid;
            zgrid = par.zgrid;
            z_median_ind = ceil(length(zgrid) / 2); % Median productivity state
            vpol = squeeze(sol.v(:, z_median_ind, :, :)); % Value function (alen, nage, slen)

            for s = 1:2 % 1 = low-skill, 2 = high-skill
                for age = 1:par.nage
                    subplot(par.nage, 2, (age - 1) * 2 + s)
                    plot(agrid, vpol(:, age, s), 'LineWidth', 2)
                    xlabel('$a_t$', 'Interpreter', 'latex', 'FontSize', 12)
                    ylabel('$V_t$', 'Interpreter', 'latex', 'FontSize', 12)
                    title(sprintf('Value Function (Age: %d, Skill: %s)', ...
                        age, my_graph.ifelse(s == 1, 'Low', 'High')), 'Interpreter', 'latex', 'FontSize', 14)
                    grid on
                end
            end

            saveas(gcf, 'value_functions_pe.png')
        end

        %% Plot value functions for different scenarios (No Tax/UBI, Tax Only, Tax+UBI).
        function [] = plot_value_functions_scenarios(par, sol_no_tax_no_ubi, sol_tax_only, sol_tax_ubi)
            figure(5)
            set(gcf, 'Position', [100, 100, 1600, 1200])

            agrid = par.agrid;
            zgrid = par.zgrid;
            z_median_ind = ceil(length(zgrid) / 2); % Median productivity state

            % Extract value functions for each scenario
            v_no_tax_no_ubi = squeeze(sol_no_tax_no_ubi.v(:, z_median_ind, :, :)); % (alen, nage, slen)
            v_tax_only = squeeze(sol_tax_only.v(:, z_median_ind, :, :));
            v_tax_ubi = squeeze(sol_tax_ubi.v(:, z_median_ind, :, :));

            % Plot for each age and skill type
            for s = 1:2 % 1 = low-skill, 2 = high-skill
                for age = 1:par.nage
                    subplot(par.nage, 2, (age - 1) * 2 + s)
                    hold on
                    plot(agrid, v_no_tax_no_ubi(:, age, s), 'b-', 'LineWidth', 2, 'DisplayName', 'No Tax, No UBI')
                    plot(agrid, v_tax_only(:, age, s), 'r--', 'LineWidth', 2, 'DisplayName', 'Tax Only')
                    plot(agrid, v_tax_ubi(:, age, s), 'g-', 'LineWidth', 2, 'DisplayName', 'Tax + UBI')
                    hold off
                    xlabel('$a_t$', 'Interpreter', 'latex', 'FontSize', 12)
                    ylabel('$V_t$', 'Interpreter', 'latex', 'FontSize', 12)
                    title(sprintf('Value Function (Age: %d, Skill: %s)', ...
                        age, my_graph.ifelse(s == 1, 'Low', 'High')), 'Interpreter', 'latex', 'FontSize', 14)
                    legend('show', 'Interpreter', 'latex', 'FontSize', 10)
                    grid on
                end
            end

            saveas(gcf, 'value_functions_scenarios_pe.png')
        end

        %% Plot labor supply distribution by skill type (Partial Equilibrium).
        function [] = plot_labor_dist(par, sol, sim)
            figure(6)
            set(gcf, 'Position', [100, 100, 1200, 400])

            lsim_middle = sim.lsim(2, :); % Middle-aged labor supply
            skill_assign = par.skill_assign;

            subplot(1, 2, 1)
            histogram(lsim_middle(skill_assign == 1), 20, 'Normalization', 'probability')
            xlabel('$l$', 'Interpreter', 'latex', 'FontSize', 12)
            ylabel('Probability', 'Interpreter', 'latex', 'FontSize', 12)
            title('Labor Supply Distribution (Low-Skill, Middle-Aged)', 'Interpreter', 'latex', 'FontSize', 14)

            subplot(1, 2, 2)
            histogram(lsim_middle(skill_assign == 2), 20, 'Normalization', 'probability')
            xlabel('$l$', 'Interpreter', 'latex', 'FontSize', 12)
            ylabel('Probability', 'Interpreter', 'latex', 'FontSize', 12)
            title('Labor Supply Distribution (High-Skill, Middle-Aged)', 'Interpreter', 'latex', 'FontSize', 14)

            saveas(gcf, 'labor_distribution_pe.png')
        end

        %% Helper function for conditional string
        function str = ifelse(condition, true_str, false_str)
            if condition
                str = true_str;
            else
                str = false_str;
            end
        end
    end
end