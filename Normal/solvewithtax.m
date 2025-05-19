classdef solvewithtax
    methods(Static)
        function [par, sol] = firm_problem(par)
            alpha = par.alpha;
            delta = par.delta;
            r = par.r;

            % Capital from FOC
            k = ((r + delta) / alpha)^(1 / (alpha - 1)); 

            % Base wage (per efficiency unit of labor)
            w_base = (1 - alpha) * k^alpha;

            % Compute skill-specific wages
            par.w_base = w_base;
            par.w_skill = w_base * par.wage_mult;  % [low-skill, high-skill]

            sol.k = k;
            sol.wage_low = par.w_skill(1);  % for reporting
            sol.wage_high = par.w_skill(2);
        end

        function sol = hh_problem(par, sol, use_tax, use_ubi)
            if nargin < 3, use_tax = true; end
            if nargin < 4, use_ubi = true; end

            % Parameters
            beta = par.beta;
            agrid = par.agrid; alen = par.alen;
            zgrid = par.zgrid; zlen = par.zlen;
            pmat = par.pmat;
            r = par.r;
            ubi = par.ubi * use_ubi;
            wage_mult = par.wage_mult;
            slen = par.slen;
            lambda = par.lambda * use_tax;
            tau = par.tau * use_tax;
            phi = -0.1;
            psi = par.psi;
            eta = par.eta;
            l_min = 0.2;

            % Value and policy function containers
            v = nan(alen, zlen, par.nage, slen);
            c = nan(alen, zlen, par.nage, slen);
            a_next = nan(alen, zlen, par.nage, slen);
            l = nan(alen, zlen, par.nage, slen);

            fprintf('------------Beginning Backward Induction (Tax: %d, UBI: %d).------------\n\n', use_tax, use_ubi)

            % --- Optional tax revenue calculation ---
            total_tax_revenue = 0;
            if use_tax && use_ubi
                for s = 1:slen
                    for j = 1:zlen
                        wage = par.w_skill(s);
                        y = wage * zgrid(j); % assume l = 1
                        net_income = lambda * y^(1 - tau);
                        tax_paid = y - net_income;
                        total_tax_revenue = total_tax_revenue + par.skill_prob(s) * tax_paid / zlen;
                    end
                end
                total_ubi = par.ubi * par.nage;
                if total_tax_revenue < total_ubi
                    fprintf('Warning: Tax revenue (%.4f) less than UBI cost (%.4f)\n', total_tax_revenue, total_ubi)
                end
            end

            % --- Backward induction ---
            for age = par.nage:-1:1
                fprintf('Solving for age %d\n', age)
                for s = 1:slen
                    wage = par.w_skill(s);
                    for i = 1:alen
                        if agrid(i) >= phi
                            for j = 1:zlen
                                resources = (1 + r) * agrid(i) + ubi;

                                if age == par.nage  % Old
                                    c(i, j, age, s) = max(0.01, resources);
                                    a_next(i, j, age, s) = 0;
                                    l(i, j, age, s) = 0;
                                    v(i, j, age, s) = modelwithtax.utility(c(i, j, age, s), par);

                                else
                                    asset_choices = agrid(agrid >= phi);

                                    if age == 2  % Middle-aged (labor decision)
                                        l_grid = linspace(l_min, 1, 50);
                                        vall = -inf * ones(length(asset_choices), length(l_grid));
                                        valid = false(size(vall));

                                        for ia = 1:length(asset_choices)
                                            a1 = asset_choices(ia);
                                            [~, idx] = min(abs(agrid - a1));
                                            ev = v(idx, :, age + 1, s) * pmat(j, :)';

                                            for il = 1:length(l_grid)
                                                labor = l_grid(il);
                                                y = wage * zgrid(j) * labor;
                                                income = use_tax * lambda * y^(1 - tau) + ~use_tax * y;
                                                c_try = resources + income - a1;

                                                if c_try > 0
                                                    valid(ia, il) = true;
                                                    u = modelwithtax.utility(c_try, par) - psi * (labor^(1 + eta)) / (1 + eta);
                                                    vall(ia, il) = u + beta * ev;
                                                end
                                            end
                                        end

                                        if any(valid(:))
                                            [vmax, idx] = max(vall(:));
                                            [ia, il] = ind2sub(size(vall), idx);
                                            a1 = asset_choices(ia);
                                            labor = l_grid(il);
                                            y = wage * zgrid(j) * labor;
                                            income = use_tax * lambda * y^(1 - tau) + ~use_tax * y;
                                            c(i, j, age, s) = resources + income - a1;
                                            a_next(i, j, age, s) = a1;
                                            l(i, j, age, s) = labor;
                                            v(i, j, age, s) = vmax;
                                        else
                                            c(i, j, age, s) = max(0.01, resources);
                                            a_next(i, j, age, s) = 0;
                                            l(i, j, age, s) = l_min;
                                            v(i, j, age, s) = modelwithtax.utility(c(i, j, age, s), par) - psi * (l_min^(1 + eta)) / (1 + eta);
                                        end
                                    else  % Young (no labor)
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
                                            vall = modelwithtax.utility(c_candidates(valid_choices), par) + beta * ev;
                                            [vmax, idx] = max(vall);
                                            v(i, j, age, s) = vmax;
                                            c(i, j, age, s) = c_candidates(valid_indices(idx));  % <-- FIXED
                                            a_next(i, j, age, s) = asset_choices(valid_indices(idx));  % <-- FIXED
                                            l(i, j, age, s) = 0;
                                        else
                                            fprintf('Warning: No valid choice found at age %d, skill %d, state (a=%d, z=%d)\n', age, s, i, j);
                                            c(i, j, age, s) = max(0.01, resources); % Ensure positive consumption
                                            a_next(i, j, age, s) = 0;
                                            l(i, j, age, s) = 0;
                                            v(i, j, age, s) = modelwithtax.utility(c(i, j, age, s), par);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            % Output results
            sol.c = c(:, :, 2, :); sol.c2 = c(:, :, 1, :); sol.c3 = c(:, :, 3, :);
            sol.a = a_next(:, :, 2, :); sol.l = l(:, :, 2, :);
            sol.v = v;
            sol.lambda = lambda; sol.tau = tau; sol.ubi = ubi;

            % Compute welfare
            util = nan(size(c));
            for age = 1:par.nage
                for s = 1:slen
                    util(:, :, age, s) = modelwithtax.utility(c(:, :, age, s), par) - psi * (l(:, :, age, s).^(1+eta)) / (1+eta);
                end
            end
            sol.welfare_low = mean(util(:, :, :, 1), 'all', 'omitnan');
            sol.welfare_high = mean(util(:, :, :, 2), 'all', 'omitnan');
            sol.welfare = mean(util, 'all', 'omitnan');

            fprintf('------------Backward Induction Complete (Tax: %d, UBI: %d).------------\n\n', use_tax, use_ubi)
            fprintf('Welfare (overall): %.4f\n', sol.welfare)
            fprintf('Welfare (low-skill): %.4f\n', sol.welfare_low)
            fprintf('Welfare (high-skill): %.4f\n', sol.welfare_high)
        end
    end
end
