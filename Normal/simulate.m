classdef simulate
    methods(Static)
        function sim = economy(par, sol)
            %% Set up.
            agrid = par.agrid;
            zgrid = par.zgrid;
            pmat = par.pmat;
            nage = par.nage;
            N = par.N;
            skill_assign = par.skill_assign;

            cpol = sol.c;         % Middle-aged consumption
            cpol_young = sol.c2;  % Young consumption
            cpol_old = sol.c3;    % Old consumption
            apol = sol.a;         % Middle-aged savings
            lpol = sol.l;         % Labor supply

            T = nage; % 3 periods

            zsim = nan(T, N);
            asim = nan(T, N);
            csim = nan(T, N);
            usim = nan(T, N);
            lsim = nan(T, N);

            %% Initial conditions
            rng(par.seed);
            pmat0 = pmat^100;
            cmat = cumsum(pmat, 2);
            z0_ind = randsample(par.zlen, N, true, pmat0(1, :))';
            a0_ind = ones(N, 1);

            %% Simulation loop
            for i = 1:N
                s = skill_assign(i);

                % Young
                at_ind = a0_ind(i);
                zt_ind = z0_ind(i);
                csim(1, i) = cpol_young(at_ind, zt_ind, 1, s);
                asim(1, i) = agrid(at_ind);
                usim(1, i) = modelwithtax.utility(csim(1, i), par);
                lsim(1, i) = 0;
                a_next = 0;
                at_ind = find(agrid == a_next);

                % Middle-aged
                zt_ind = z0_ind(i);
                csim(2, i) = cpol(at_ind, zt_ind, 1, s);
                asim(2, i) = apol(at_ind, zt_ind, 1, s);
                lsim(2, i) = lpol(at_ind, zt_ind, 1, s);
                usim(2, i) = modelwithtax.utility(csim(2, i), par);
                zsim(2, i) = zgrid(zt_ind);

                at_ind = find(agrid <= asim(2, i), 1, 'last');
                z1_ind = find(rand <= cmat(zt_ind, :), 1);
                zt_ind = z1_ind;

                % Old
                csim(3, i) = cpol_old(at_ind, zt_ind, 1, s);
                asim(3, i) = 0;
                usim(3, i) = modelwithtax.utility(csim(3, i), par);
                lsim(3, i) = 0;
            end

            %% Aggregates and welfare
            sim = struct();
            sim.zsim = zsim;
            sim.asim = asim;
            sim.csim = csim;
            sim.usim = usim;
            sim.lsim = lsim;

            sim.asup = mean(asim(2, :));
            sim.csup = mean(csim, 'all');
            sim.lsup = mean(lsim(2, :));

            sim.welfare = mean(usim, 'all', 'omitnan');
            sim.welfare_low = mean(usim(:, skill_assign == 1), 'all', 'omitnan');
            sim.welfare_high = mean(usim(:, skill_assign == 2), 'all', 'omitnan');

            fprintf('Simulation completed.\n')
            fprintf('Avg savings (middle-aged): %.4f\n', sim.asup)
            fprintf('Avg consumption: %.4f\n', sim.csup)
            fprintf('Avg labor (middle-aged): %.4f\n', sim.lsup)
            fprintf('Welfare (overall): %.4f\n', sim.welfare)
            fprintf('Welfare (low-skill): %.4f\n', sim.welfare_low)
            fprintf('Welfare (high-skill): %.4f\n', sim.welfare_high)
        end
    end
end
