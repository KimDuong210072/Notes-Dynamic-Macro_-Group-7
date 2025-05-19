classdef lmodelwithtax
    methods(Static)
        function par = setup()            
            par = struct();

            %% Preferences.
            par.beta = 0.94; % Discount factor
            par.sigma = 1.5; % Risk aversion (CRRA utility)
            par.psi = 0.3; % Disutility of labor
            par.eta = 0.3; % Inverse Frisch elasticity
            assert(par.beta > 0 && par.beta < 1.00, 'Beta must be in (0,1)')
            assert(par.sigma > 0, 'Sigma must be positive')
            assert(par.psi > 0, 'Psi must be positive')
            assert(par.eta >= 0, 'Eta must be non-negative')

            %% Technology.
            par.alpha = 0.33; % Capital share in Cobb-Douglas production
            par.delta = 0.03; % Depreciation rate
            assert(par.alpha > 0 && par.alpha < 1.00, 'Alpha must be in (0,1)')
            assert(par.delta >= 0 && par.delta <= 1.00, 'Delta must be in [0,1]')

            par.sigma_eps = 0.07; % Standard deviation of productivity shocks
            par.rho = 0.85; % Persistence of productivity shocks
            par.mu = 0.0; % Mean of productivity shocks
            assert(par.sigma_eps > 0, 'Sigma_eps must be positive')
            assert(abs(par.rho) < 1, 'Rho must be in (-1,1)')

            par.age = [1, 2, 3]; % 1 = young, 2 = middle-aged, 3 = old
            par.nage = 3; % Number of age groups
            par.e = [0.0, 1.0, 0.0]; % Productivity profile (only middle-aged work)

            %% Skills.
            par.skill_types = [1, 2]; % 1 = low-skill, 2 = high-skill
            par.slen = 2; % Number of skill types
            par.skill_prob = [0.5, 0.5]; % Probability of low/high skill
            par.wage_mult = [1.0, 1.5]; % Wage multipliers (low = 1.0, high = 1.5)
            par.N = 100; % Number of agents
            par.seed = 2025; % Random seed
            rng(par.seed); % Set random seed for reproducibility
            par.skill_assign = rand(1, par.N) < par.skill_prob(2); % Random skill assignment
            par.skill_assign = par.skill_assign + 1; % Convert to 1 (low), 2 (high)
            assert(sum(par.skill_prob) == 1, 'Skill probabilities must sum to 1')
            assert(all(par.wage_mult > 0), 'Wage multipliers must be positive')

            %% Prices and Policy.
            par.r = 0.06; % Interest rate (fixed for PE)
            par.w = 0.2; % Base wage (updated by firm_problem)
            par.theta = 0.0; % Fixed tax rate (to be computed)
            par.ubi = 0.05;
            assert(par.r >= -1, 'Interest rate must allow non-negative returns')
            assert(par.w > 0, 'Wage must be positive')
            assert(par.theta >= 0, 'Tax rate must be non-negative')
            assert(par.ubi >= 0, 'UBI must be non-negative')

            %% Simulation parameters.
            par.T = 10; % Simulation periods (for extensions)
            assert(par.N > 0, 'Number of agents must be positive')
        end

        function u = utility(c, par)
            sigma = par.sigma;
            if c <= 0
                u = -inf;
            else
                if sigma == 1
                    u = log(c);
                else
                    u = (c.^(1 - sigma)) / (1 - sigma); % Use element-wise power .^
                end
            end
        end

        function par = gen_grids(par)
            %% Asset grid
            par.amin = 0;
            par.amax = 10;
            par.alen = 500;
            par.agrid = linspace(par.amin, par.amax, par.alen)';

            %% Productivity grid
            par.zlen = 5;
            [zgrid, pmat] = lmodelwithtax.tauchen(par.mu, par.rho, par.sigma_eps, par.zlen, 3);
            par.zgrid = exp(zgrid);
            par.pmat = pmat;
        end

        function [zgrid, pmat] = tauchen(mu, rho, sigma, nz, nstd)
            % Discretize an AR(1) process using Tauchen's method
            zgrid = linspace(mu - nstd * sigma / sqrt(1 - rho^2), ...
                            mu + nstd * sigma / sqrt(1 - rho^2), nz)';
            dz = zgrid(2) - zgrid(1);
            pmat = zeros(nz, nz);
            
            for i = 1:nz
                for j = 1:nz
                    if j == 1
                        pmat(i, j) = normcdf((zgrid(j) + dz/2 - mu - rho * (zgrid(i) - mu)) / sigma);
                    elseif j == nz
                        pmat(i, j) = 1 - normcdf((zgrid(j) - dz/2 - mu - rho * (zgrid(i) - mu)) / sigma);
                    else
                        pmat(i, j) = normcdf((zgrid(j) + dz/2 - mu - rho * (zgrid(i) - mu)) / sigma) - ...
                                     normcdf((zgrid(j) - dz/2 - mu - rho * (zgrid(i) - mu)) / sigma);
                    end
                end
            end
            % Ensure rows sum to 1
            pmat = pmat ./ sum(pmat, 2);
        end
    end
end