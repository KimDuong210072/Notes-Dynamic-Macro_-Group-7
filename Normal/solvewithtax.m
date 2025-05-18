classdef solvewithtax
    methods(Static)
        function [par, sol] = firm_problem(par)
            alpha = par.alpha;
            delta = par.delta;
            r = par.r;
            e = par.e;

            k = ((r + delta) / alpha)^(1 / (alpha - 1)); % Capital stock (from FOC)
            w = (1 - alpha) * (k^alpha); % Wage (from FOC)

            par.w = w;
            sol.k = k;
            sol.wage = w * e(2); % Wage for middle-aged (e(2) = 1)
        end

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
            wage_mult = par.wage_mult;
            slen = par.slen;
            tau0 = par.tau0; % Base tax rate
            tau1 = par.tau1; % Marginal tax rate
            phi = 0; % No borrowing
            psi = par.psi; % Disutility of labor
            eta = par.eta; % Inverse Frisch elasticity

            % Initialize arrays
            v = nan(alen, zlen, par.nage, slen); % Value function
            c = nan(alen, zlen, par.nage, slen); % Consumption
            a_next = nan(alen, zlen, par.nage, slen); % Savings
            l = nan(alen, zlen, par.nage, slen); % Labor supply

            fprintf('------------Beginning Backward Induction.------------\n\n')

            % Compute tax revenue and adjust tau0 for budget balance
            total_tax_revenue = 0;
            for s = 1:slen
                for j = 1:zlen
                    y = w * wage_mult(s) * zgrid(j); % Middle-aged income (assume l=1 for initial estimate)
                    tax_rate = tau0 + tau1 * y;
                    total_tax_revenue = total_tax_revenue + par.skill_prob(s) * tax_rate * y / zlen;
                end
            end
            total_ubi = par.ubi * par.nage; % Total UBI payments
            if total_tax_revenue < total_ubi
                scale = total_ubi / total_tax_revenue;
                par.tau0 = par.tau0 * scale; % Adjust only tau0 to balance budget
                % par.tau1 = par.tau1 * scale; % Do not scale tau1
            end
            tau0 = par.tau0;
            tau1 = par.tau1;
            fprintf('Adjusted progressive tax: tau0 = %.4f, tau1 = %.4f\n', tau0, tau1)

            % Backward induction
            for age = par.nage:-1:1
                fprintf('Solving for age %d\n', age)
                for s = 1:slen
                    for i = 1:alen
                        if agrid(i) >= phi
                            for j = 1:zlen
                                % Current resources (without labor income)
                                resources = (1 + r) * agrid(i) + ubi;

                                if age == par.nage % Old age
                                    c(i, j, age, s) = resources;
                                    a_next(i, j, age, s) = 0;
                                    l(i, j, age, s) = 0;
                                    v(i, j, age, s) = modelwithtax.utility(c(i, j, age, s), par);
                                else % Young or Middle-aged
                                    % Asset choices
                                    asset_choices = agrid(agrid >= phi);
                                    if age == 2 % Middle-aged: Optimize labor supply
                                        l_grid = linspace(0, 1, 50); % Labor supply grid
                                        vall = -inf * ones(length(asset_choices), length(l_grid));
                                        for il = 1:length(l_grid)
                                            labor = l_grid(il);
                                            y = w * wage_mult(s) * zgrid(j) * labor;
                                            tax_rate = tau0 + tau1 * y;
                                            income = (1 - tax_rate) * y;
                                            c_candidates = resources + income - asset_choices;
                                            c_candidates(c_candidates <= 0) = -inf;
                                            u = modelwithtax.utility(c_candidates, par) - psi * (l_grid(il).^(1+eta))/(1+eta); % Use element-wise power .^
                                            ev = v(:, :, age + 1, s) * pmat(j, :)';
                                            vall(:, il) = u + beta * ev;
                                            vall(c_candidates <= 0, il) = -inf;
                                        end
                                        [vmax, ind] = max(vall(:));
                                        [ia, il] = ind2sub([length(asset_choices), length(l_grid)], ind);
                                        v(i, j, age, s) = vmax;
                                        y = w * wage_mult(s) * zgrid(j) * l_grid(il);
                                        c(i, j, age, s) = resources + (1 - (tau0 + tau1 * y)) * y - asset_choices(ia);
                                        a_next(i, j, age, s) = asset_choices(ia);
                                        l(i, j, age, s) = l_grid(il);
                                    else % Young: No labor
                                        c_candidates = resources - asset_choices;
                                        c_candidates(c_candidates <= 0) = -inf;
                                        ev = mean(v(:, :, age + 1, s), 2);
                                        vall = modelwithtax.utility(c_candidates, par) + beta * ev;
                                        vall(c_candidates <= 0) = -inf;
                                        [vmax, ind] = max(vall);
                                        v(i, j, age, s) = vmax;
                                        c(i, j, age, s) = c_candidates(ind);
                                        a_next(i, j, age, s) = asset_choices(ind);
                                        l(i, j, age, s) = 0;
                                    end
                                end
                            end
                        end
                    end
                end
            end

            % Store solutions
            sol.c = c(:, :, 2, :); % Middle-aged consumption
            sol.c2 = c(:, :, 1, :); % Young consumption
            sol.c3 = c(:, :, 3, :); % Old consumption
            sol.v = v;
            sol.a = a_next(:, :, 2, :); % Middle-aged savings
            sol.l = l(:, :, 2, :); % Middle-aged labor supply
            sol.tau0 = tau0;
            sol.tau1 = tau1;

            % Compute welfare
            util = nan(alen, zlen, par.nage, slen);
            for age = 1:par.nage
                for s = 1:slen
                    util(:, :, age, s) = modelwithtax.utility(c(:, :, age, s), par) - psi * (l(:, :, age, s).^(1+eta))/(1+eta); % Use element-wise power .^
                end
            end
            avg_utility_low = mean(util(:, :, :, 1), 'all', 'omitnan');
            avg_utility_high = mean(util(:, :, :, 2), 'all', 'omitnan');
            avg_utility = mean(util, 'all', 'omitnan');
            sol.welfare_low = avg_utility_low;
            sol.welfare_high = avg_utility_high;
            sol.welfare = avg_utility;

            fprintf('------------Backward Induction Complete.------------\n\n')
            fprintf('Welfare (overall): %.4f\n', sol.welfare)
            fprintf('Welfare (low-skill): %.4f\n', sol.welfare_low)
            fprintf('Welfare (high-skill): %.4f\n', sol.welfare_high)
        end
    end
end