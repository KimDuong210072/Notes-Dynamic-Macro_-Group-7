%% File Info.
%{
    solve.m
    -------
    This code solves the model with UBI, linear taxes, and computes welfare using backward induction.
%}

%% Solve class.
classdef solvewithtax
    methods(Static)
        function sol = hh_problem(par, sol)
            beta = par.beta;
            agrid = par.agrid;
            alen = par.alen;
            zgrid = par.zgrid;
            zlen = par.zlen;
            pmat = par.pmat;
            r = par.r;
            w = par.w;
            ubi = par.ubi;
            phi = 0; % No borrowing

            % Initialize arrays
            v = nan(alen, zlen, par.nage); % Value function for each age
            c = nan(alen, zlen, par.nage); % Consumption for each age
            a_next = nan(alen, zlen, par.nage); % Savings (a_{t+1}) for each age

            fprintf('------------Beginning Backward Induction.------------\n\n')

            % Compute tax revenue and adjust tau for budget balance
            income = w * par.e(2) * zgrid; % Middle-aged labor income
            avg_income = mean(income); % Average labor income
            total_ubi = par.ubi * par.nage; % Total UBI payments (3 age groups)
            par.tau = total_ubi / avg_income; % Tax rate to balance budget
            tau = par.tau;
            fprintf('Adjusted tax rate for budget balance: tau = %.4f\n', tau)

            % Backward induction: Loop over ages in reverse order (3, 2, 1)
            for age = par.nage:-1:1
                fprintf('Solving for age %d\n', age)
                for i = 1:alen
                    if agrid(i) >= phi
                        for j = 1:zlen
                            % Current resources
                            income = (1 - tau) * w * par.e(age) * zgrid(j); % Labor income (0 for young/old)
                            resources = (1 + r) * agrid(i) + income + ubi;

                            if age == par.nage % Old age (age 3)
                                % Consume all resources (no savings)
                                c(i, j, age) = resources;
                                a_next(i, j, age) = 0; % No savings
                                v(i, j, age) = modelwithtax.utility(c(i, j, age), par);
                            else % Young (age 1) or Middle-aged (age 2)
                                % Asset choices for next period
                                asset_choices = agrid;
                                c_candidates = resources - asset_choices;
                                c_candidates(c_candidates <= 0) = -inf; % Ensure positive consumption

                                % Expected value for next age
                                if age == 1
                                    % Young: Transition to middle-aged, average over z'
                                    ev = mean(v(:, :, age + 1), 2); % Average over productivity states
                                else
                                    % Middle-aged: Transition to old, use pmat
                                    ev = v(:, :, age + 1) * pmat(j, :)';
                                end

                                % Utility + discounted expected value
                                vall = modelwithtax.utility(c_candidates, par) + beta * ev;
                                vall(c_candidates <= 0) = -inf;

                                % Maximize
                                [vmax, ind] = max(vall);
                                v(i, j, age) = vmax;
                                c(i, j, age) = c_candidates(ind);
                                a_next(i, j, age) = asset_choices(ind);
                            end
                        end
                    end
                end
            end

            fprintf('------------End of Backward Induction.------------\n')

            % Compute welfare (average utility of consumption)
            util = nan(alen, zlen, par.nage);
            for age = 1:par.nage
                util(:, :, age) = modelwithtax.utility(c(:, :, age), par);
            end
            avg_utility = mean(util(:)); % Average over all ages, assets, productivity
            fprintf('Average utility of consumption: %.4f\n', avg_utility)

            % Store results
            sol.a = a_next(:, :, 2); % Savings from middle-aged to old
            sol.c = c(:, :, 2); % Middle-aged consumption
            sol.c2 = c(:, :, 1); % Young consumption
            sol.c3 = c(:, :, 3); % Old consumption
            sol.v = v;
            sol.welfare = avg_utility;
            sol.tau = tau;
        end

        function [par, sol] = firm_problem(par)
            delta = par.delta;
            alpha = par.alpha;
            r = par.r;
            k = ((r + delta) / alpha) ^ (1 / (alpha - 1));

            sol = struct();
            sol.k = k;
            par.w = (1 - alpha) * k ^ alpha;
        end
    end
end