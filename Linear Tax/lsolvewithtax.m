classdef lsolvewithtax
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

        function sol = hh_problem(par, sol, use_tax, use_ubi)
            if nargin < 3
                use_tax = true; % Default: include tax
            end
            if nargin < 4
                use_ubi = true; % Default: include UBI
            end

            beta = par.beta;
            agrid = par.agrid;
            alen = par.alen;
            zgrid = par.zgrid;
            zlen = par.zlen;
            pmat = par.pmat;
            r = par.r;
            w = par.w;
            ubi = par.ubi * use_ubi; % Apply UBI if use_ubi is true
            wage_mult = par.wage_mult;
            slen = par.slen;
            theta = par.theta; % Initialize theta
            phi = -0.1; % Allow small borrowing
            psi = par.psi; % Disutility of labor
            eta = par.eta; % Inverse Frisch elasticity
            l_min = 0.2; % Minimum labor supply to ensure positive income

            % Initialize arrays
            v = nan(alen, zlen, par.nage, slen); % Value function
            c = nan(alen, zlen, par.nage, slen); % Consumption
            a_next = nan(alen, zlen, par.nage, slen); % Savings
            l = nan(alen, zlen, par.nage, slen); % Labor supply

            fprintf('------------Beginning Backward Induction (Tax: %d, UBI: %d).------------\n\n', use_tax, use_ubi)

            % Compute tax rate theta to balance budget (if using tax and UBI)
            if use_tax && use_ubi
                total_tax_revenue = 0;
                for s = 1:slen
                    for j = 1:zlen
                        y = w * wage_mult(s) * zgrid(j); % Middle-aged income (assume l=1 for initial estimate)
                        tax_paid = theta * y; % Flat tax
                        total_tax_revenue = total_tax_revenue + par.skill_prob(s) * tax_paid / zlen;
                    end
                end
                total_ubi = par.ubi * par.nage * par.N; % Total UBI payments across all agents and periods
                if total_tax_revenue > 0
                    theta = total_ubi / total_tax_revenue; % Compute theta to balance budget
                    par.theta = theta; % Update par.theta
                    fprintf('Computed flat tax rate theta: %.4f\n', theta)
                else
                    fprintf('Warning: Zero tax revenue. Setting theta to 0.\n')
                    theta = 0;
                end
            elseif use_tax && ~use_ubi
                theta = 0.1; % Default flat tax rate (10%) when no UBI
                par.theta = theta;
                fprintf('Using default flat tax rate theta: %.4f\n', theta)
            else
                theta = 0; % No tax
            end

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
                                    if resources <= 0
                                        c(i, j, age, s) = 0.01; % Ensure positive consumption
                                    else
                                        c(i, j, age, s) = resources;
                                    end
                                    a_next(i, j, age, s) = 0;
                                    l(i, j, age, s) = 0;
                                    v(i, j, age, s) = lmodelwithtax.utility(c(i, j, age, s), par);
                                else % Young or Middle-aged
                                    % Asset choices
                                    asset_choices = agrid(agrid >= phi);
                                    if age == 2 % Middle-aged: Optimize labor supply
                                        l_grid = linspace(l_min, 1, 50); % Labor supply grid with minimum
                                        vall = -inf * ones(length(asset_choices), length(l_grid));
                                        valid_choices = false(length(asset_choices), length(l_grid));
                                        for ia = 1:length(asset_choices)
                                            a_next_val = asset_choices(ia);
                                            [~, idx] = min(abs(agrid - a_next_val)); % Find closest asset index for interpolation
                                            ev = v(idx, :, age + 1, s) * pmat(j, :)'; % Scalar expected value
                                            for il = 1:length(l_grid)
                                                labor = l_grid(il);
                                                y = w * wage_mult(s) * zgrid(j) * labor;
                                                if use_tax
                                                    income = (1 - theta) * y; % Flat tax
                                                else
                                                    income = y; % No tax
                                                end
                                                c_candidate = resources + income - a_next_val;
                                                if c_candidate > 0
                                                    valid_choices(ia, il) = true;
                                                    u = lmodelwithtax.utility(c_candidate, par) - psi * (labor^(1+eta))/(1+eta);
                                                    vall(ia, il) = u + beta * ev; % Both u and ev are scalars
                                                end
                                            end
                                        end
                                        if any(valid_choices(:))
                                            [vmax, idx] = max(vall(:));
                                            [ia, il] = ind2sub(size(vall), idx);
                                            v(i, j, age, s) = vmax;
                                            a_next_val = asset_choices(ia);
                                            labor = l_grid(il);
                                            y = w * wage_mult(s) * zgrid(j) * labor;
                                            if use_tax
                                                income = (1 - theta) * y;
                                            else
                                                income = y;
                                            end
                                            c(i, j, age, s) = resources + income - a_next_val;
                                            a_next(i, j, age, s) = a_next_val;
                                            l(i, j, age, s) = labor;
                                        else
                                            fprintf('Warning: No valid choice found at age %d, skill %d, state (a=%d, z=%d)\n', age, s, i, j);
                                            c(i, j, age, s) = max(0.01, resources); % Ensure positive consumption
                                            a_next(i, j, age, s) = 0;
                                            l(i, j, age, s) = l_min;
                                            v(i, j, age, s) = lmodelwithtax.utility(c(i, j, age, s), par) - psi * (l(i, j, age, s)^(1+eta))/(1+eta);
                                        end
                                    else % Young: No labor
                                        c_candidates = resources - asset_choices;
                                        valid_choices = c_candidates > 0;
                                        if any(valid_choices)
                                            valid_indices = find(valid_choices);
                                            ev = zeros(length(valid_indices), 1);
                                            for k = 1:length(valid_indices)
                                                a_next_val = asset_choices(valid_indices(k));
                                                [~, idx] = min(abs(agrid - a_next_val));
                                                ev(k) = mean(v(idx, :, age + 1, s)); % Scalar expected value
                                            end
                                            vall = lmodelwithtax.utility(c_candidates(valid_choices), par) + beta * ev;
                                            [vmax, idx] = max(vall);
                                            v(i, j, age, s) = vmax;
                                            c(i, j, age, s) = c_candidates(valid_indices(idx));
                                            a_next(i, j, age, s) = asset_choices(valid_indices(idx));
                                            l(i, j, age, s) = 0;
                                        else
                                            fprintf('Warning: No valid choice found at age %d, skill %d, state (a=%d, z=%d)\n', age, s, i, j);
                                            c(i, j, age, s) = max(0.01, resources); % Ensure positive consumption
                                            a_next(i, j, age, s) = 0;
                                            l(i, j, age, s) = 0;
                                            v(i, j, age, s) = lmodelwithtax.utility(c(i, j, age, s), par);
                                        end
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
            sol.theta = theta; % Store computed tax rate
            sol.ubi = ubi;

            % Compute welfare
            util = nan(alen, zlen, par.nage, slen);
            for age = 1:par.nage
                for s = 1:slen
                    util(:, :, age, s) = lmodelwithtax.utility(c(:, :, age, s), par) - psi * (l(:, :, age, s).^(1+eta))/(1+eta);
                end
            end
            avg_utility_low = mean(util(:, :, :, 1), 'all', 'omitnan');
            avg_utility_high = mean(util(:, :, :, 2), 'all', 'omitnan');
            avg_utility = mean(util, 'all', 'omitnan');
            sol.welfare_low = avg_utility_low;
            sol.welfare_high = avg_utility_high;
            sol.welfare = avg_utility;

            fprintf('------------Backward Induction Complete (Tax: %d, UBI: %d).------------\n\n', use_tax, use_ubi)
            fprintf('Welfare (overall): %.4f\n', sol.welfare)
            fprintf('Welfare (low-skill): %.4f\n', sol.welfare_low)
            fprintf('Welfare (high-skill): %.4f\n', sol.welfare_high)
        end
    end
end