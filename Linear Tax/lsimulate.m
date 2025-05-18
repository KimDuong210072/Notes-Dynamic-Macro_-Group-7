classdef lsimulate
    methods(Static)
        function sim = economy(par, sol)
            %% Set up.
            agrid = par.agrid;
            zgrid = par.zgrid;
            pmat = par.pmat;
            nage = par.nage;
            N = par.N;
            skill_assign = par.skill_assign;
            psi = par.psi;
            eta = par.eta;

            cpol = sol.c; % Middle-aged consumption
            cpol_young = sol.c2; % Young consumption
            cpol_old = sol.c3; % Old consumption
            apol = sol.a; % Middle-aged savings
            lpol = sol.l; % Labor supply

            T = nage; % 3 periods (young, middle-aged, old)

            zsim = nan(T, N);
            asim = nan(T, N);
            csim = nan(T, N);
            usim = nan(T, N);
            lsim = nan(T, N);

            %% Begin simulation.
            rng(par.seed);

            % Initial conditions
            pmat0 = pmat^100; % Stationary distribution
            cmat = cumsum(pmat, 2); % CDF matrix
            z0_ind = randsample(par.zlen, N, true, pmat0(1, :))'; % Initial productivity indices
            a0_ind = ones(N, 1); % Initial assets = 0 (agrid(1) = 0)

            % Debug skill assignment
            fprintf('Number of low-skill agents: %d\n', sum(skill_assign == 1));
            fprintf('Number of high-skill agents: %d\n', sum(skill_assign == 2));

            %% Simulate life-cycle
            for i = 1:N
                s = skill_assign(i); % Skill type (1 = low, 2 = high)

                % Age 1: Young
                at_ind = a0_ind(i);
                zt_ind = z0_ind(i);
                csim(1, i) = cpol_young(at_ind, zt_ind, 1, s);
                asim(1, i) = agrid(at_ind); % Assets = 0
                usim(1, i) = lmodelwithtax.utility(csim(1, i), par);
                lsim(1, i) = 0;
                a_next = 0; % Young don't save
                at_ind = find(agrid == a_next);

                % Age 2: Middle-aged
                zt_ind = z0_ind(i);
                csim(2, i) = cpol(at_ind, zt_ind, 1, s);
                asim(2, i) = apol(at_ind, zt_ind, 1, s);
                lsim(2, i) = lpol(at_ind, zt_ind, 1, s);
                fprintf('Agent %d, Skill %d, Savings: %.4f, Labor: %.4f\n', i, s, asim(2, i), lsim(2, i)); % Debug
                usim(2, i) = lmodelwithtax.utility(csim(2, i), par) - psi * (lsim(2, i)^(1+eta))/(1+eta);
                zsim(2, i) = zgrid(zt_ind);
                at_ind = find(agrid <= asim(2, i), 1, 'last');
                z1_ind = find(rand <= cmat(zt_ind, :));
                zt_ind = z1_ind(1);

                % Age 3: Old
                csim(3, i) = cpol_old(at_ind, zt_ind, 1, s);
                asim(3, i) = 0;
                usim(3, i) = lmodelwithtax.utility(csim(3, i), par);
                lsim(3, i) = 0;
            end

            %% Compute aggregates and welfare
            sim = struct();
            sim.zsim = zsim;
            sim.asim = asim;
            sim.csim = csim;
            sim.usim = usim;
            sim.lsim = lsim;
            sim.asup = mean(asim(2, :));
            sim.csup = mean(csim, 'all');
            sim.lsup = mean(lsim(2, :));

            welfare_low = mean(usim(:, skill_assign == 1), 'all', 'omitnan');
            welfare_high = mean(usim(:, skill_assign == 2), 'all', 'omitnan');
            welfare = mean(usim, 'all', 'omitnan');

            sim.welfare = welfare;
            sim.welfare_low = welfare_low;
            sim.welfare_high = welfare_high;

            fprintf('Simulation completed.\n')
            fprintf('Average savings (middle-aged): %.4f\n', sim.asup)
            fprintf('Average consumption: %.4f\n', sim.csup)
            fprintf('Average labor supply (middle-aged): %.4f\n', sim.lsup)
            fprintf('Welfare (overall): %.4f\n', sim.welfare)
            fprintf('Welfare (low-skill): %.4f\n', sim.welfare_low)
            fprintf('Welfare (high-skill): %.4f\n', sim.welfare_high)
        end
    end
end