classdef lsolvewithtax
    methods(Static)
        function [par, sol] = firm_problem(par)
            alpha = par.alpha;
            delta = par.delta;
            r = par.r;
            e = par.e;

            k = ((r + delta) / alpha)^(1 / (alpha - 1));
            w = (1 - alpha) * (k^alpha);

            par.w = w;
            sol.k = k;
            sol.wage = w * e(2); % wage for middle-aged
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
            tau = par.tau; % flat tax rate
            phi = 0;
            psi = par.psi;
            eta = par.eta;

            v = nan(alen, zlen, par.nage, slen);
            c = nan(alen, zlen, par.nage, slen);
            a_next = nan(alen, zlen, par.nage, slen);
            l = nan(alen, zlen, par.nage, slen);

            fprintf('------------Beginning Backward Induction.------------\n\n')

            % Compute tax revenue and adjust tau for budget balance
            total_tax_revenue = 0;
            for s = 1:slen
                for j = 1:zlen
                    y = w * wage_mult(s) * zgrid(j); % assume l = 1
                    total_tax_revenue = total_tax_revenue + par.skill_prob(s) * tau * y / zlen;
                end
            end
            total_ubi = par.ubi * par.nage;
            if total_tax_revenue < total_ubi
                scale = total_ubi / total_tax_revenue;
                par.tau = par.tau * scale;
            end
            tau = par.tau;
            fprintf('Adjusted flat tax: tau = %.4f\n', tau)

            for age = par.nage:-1:1
                fprintf('Solving for age %d\n', age)
                for s = 1:slen
                    for i = 1:alen
                        if agrid(i) >= phi
                            for j = 1:zlen
                                resources = (1 + r) * agrid(i) + ubi;

                                if age == par.nage % old
                                    c(i, j, age, s) = resources;
                                    a_next(i, j, age, s) = 0;
                                    l(i, j, age, s) = 0;
                                    v(i, j, age, s) = lmodelwithtax.utility(c(i, j, age, s), par);
                                else
                                    asset_choices = agrid(agrid >= phi);

                                    if age == 2 % middle-aged
                                        l_grid = linspace(0, 1, 50);
                                        vall = -inf * ones(length(asset_choices), length(l_grid));
                                        for il = 1:length(l_grid)
                                            labor = l_grid(il);
                                            y = w * wage_mult(s) * zgrid(j) * labor;
                                            income = (1 - tau) * y;
                                            c_candidates = resources + income - asset_choices;
                                            u = lmodelwithtax.utility(c_candidates, par) - psi * (labor.^(1+eta)) / (1+eta);
                                            ev = v(:, :, age + 1, s) * pmat(j, :)';
                                            vall(:, il) = u + beta * ev;
                                            vall(c_candidates <= 0, il) = -inf;
                                        end
                                        [vmax, ind] = max(vall(:));
                                        [ia, il] = ind2sub([length(asset_choices), length(l_grid)], ind);
                                        v(i, j, age, s) = vmax;
                                        labor = l_grid(il);
                                        y = w * wage_mult(s) * zgrid(j) * labor;
                                        c(i, j, age, s) = resources + (1 - tau) * y - asset_choices(ia);
                                        a_next(i, j, age, s) = asset_choices(ia);
                                        l(i, j, age, s) = labor;
                                    else % young
                                        c_candidates = resources - asset_choices;
                                        ev = mean(v(:, :, age + 1, s), 2);
                                        vall = lmodelwithtax.utility(c_candidates, par) + beta * ev;
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

            sol.c = c(:, :, 2, :);
            sol.c2 = c(:, :, 1, :);
            sol.c3 = c(:, :, 3, :);
            sol.v = v;
            sol.a = a_next(:, :, 2, :);
            sol.l = l(:, :, 2, :);
            sol.tau = tau;

            util = nan(alen, zlen, par.nage, slen);
            for age = 1:par.nage
                for s = 1:slen
                    util(:, :, age, s) = lmodelwithtax.utility(c(:, :, age, s), par) - psi * (l(:, :, age, s).^(1+eta)) / (1+eta);
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
