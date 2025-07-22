function erdosrenyiexperiment(dimgt, dimopt, n, ps, num_repeats, generategt, initscheme, options, noiselevel)
% Perform SNL experiments on Erdős–Rényi graphs.
%
% Inputs:
%   dimgt        - Dimension of the ground truth configuration (e.g., 2 or 3).
%   dimopt       - Embedding dimension used during optimization.
%   n            - Number of nodes (points) in the network.
%   ps           - Array of Erdős–Rényi graph edge probabilities (densities).
%   num_repeats  - Number of random trials per density.
%   generategt   - Function handle to generate ground truth configurations (size n x dimgt).
%   initscheme   - Function handle to generate initializations for optimization (size n x dimopt).
%   options      - Struct of options for the optimization algorithm (e.g., tolerance, max iterations).
%   noiselevel   - Variance of Gaussian noise added to the squared distance measurements.
%
% This function runs multiple SNL optimization problems using random Erdős–Rényi graphs
% with varying edge densities, and evaluates performance under different initialization schemes.
%
    tic

    fid = fopen("/dev/random");
    rndval = fread(fid, 1, "uint32");
    fclose(fid);
    rng(rndval);

    % If you don't have the /dev/random file, comment the above and use the following line:
    % rng(milliseconds(datetime("now") - datetime("2025-07-15"))/100);

    num_ps = numel(ps);
    chunk_size = 100;
    num_outer_repeats = num_repeats / chunk_size;

    fprintf("k = 0/%d", num_repeats);
    fprintf(" (" + string(toc) + "s)\n");
    for k = 1:num_outer_repeats
        idstring = string(datetime("now", "Format", "MMM-d-y-HH:mm:SS"));

        ground_truths = nan(chunk_size, num_ps, n, dimgt);
        masks = nan(n, n, num_ps);
        inits = nan(chunk_size, num_ps, n, dimopt);

        connectedness = nan(chunk_size, num_ps);

        costs_ground_truths = nan(chunk_size, num_ps);
        costs_ground_truths_reopt = nan(chunk_size, num_ps);

        costs = nan(chunk_size, num_ps);
        gradnorms = nan(chunk_size, num_ps);
        mineigvals = nan(chunk_size, num_ps);
        niters = nan(chunk_size, num_ps);
        normprocrustdiff = nan(chunk_size, num_ps);
        solutions = nan(chunk_size, num_ps, n, dimopt);

        costs_proj = nan(chunk_size, num_ps);
        gradnorms_proj = nan(chunk_size, num_ps);
        mineigvals_proj = nan(chunk_size, num_ps);
        normprocrustdiff_proj = nan(chunk_size, num_ps);
        solutions_proj = nan(chunk_size, num_ps, n, dimgt);

        costs_reopt = nan(chunk_size, num_ps);
        gradnorms_reopt = nan(chunk_size, num_ps);
        mineigvals_reopt = nan(chunk_size, num_ps);
        niters_reopt = nan(chunk_size, num_ps);
        normprocrustdiff_reopt = nan(chunk_size, num_ps);
        normprocrustdiffysreopt_reopt = nan(chunk_size, num_ps);
        solutions_reopt = nan(chunk_size, num_ps, n, dimgt);

        for i = 1:chunk_size
            for j = 1:num_ps
                inits(i, j, :, :) = initscheme(n, dimopt);
                ground_truths(i, j, :, :) = generategt(n, dimgt);
                masks(:, :, j) = erdosrenyimask(n, ps(j));
            end
            parfor j = 1:num_ps
                Ys = squeeze(ground_truths(i, j, :, :));
                D = delta(Ys * Ys');
                mask = masks(:, :, j);

                rdmat = randn(n, n);
                noisemat = noiselevel * (rdmat + rdmat')/2;

                Dnoise = D + noisemat;

                connectedness(i, j) = isconnected(mask);

                % Reoptimize ground truth for noisy settings
                [Ys_reopt, Ys_reopt_cost, ~, ~, ~] = findminimum(n, dimgt, Dnoise, mask, Ys, options);
                costs_ground_truths_reopt(i, j) = Ys_reopt_cost;

                % Base problem
                [Y, cost, info, ~, problem] = findminimum(n, dimopt, Dnoise, mask, squeeze(inits(i, j, :, :)), options);

                costs_ground_truths(i, j) = getCost(problem, Ys);
                costs(i, j) = cost;
                gradnorms(i, j) = norm(gradf(Y, Dnoise, mask), "fro");
                mineigvals(i, j) = min(eig(hessianmatrix(problem, Y)));
                niters(i, j) = max([info.iter]);
                [~, ~, r] = procrust(Y, [Ys zeros(n, dimopt - dimgt)]);
                normprocrustdiff(i, j) = r;
                solutions(i, j, :, :) = Y;

                % Project to dimgt
                [U, S, ~] = svd(Y, "econ");
                Y_proj = U(:, 1:dimgt) * S(1:dimgt, 1:dimgt);

                % Reoptimize in dimgt
                if dimopt > dimgt
                    [Y_reopt, cost_reopt, info_reopt, ~, problem_reopt] = findminimum(n, dimgt, Dnoise, mask, Y_proj, options);

                    costs_reopt(i, j) = cost_reopt;
                    gradnorms_reopt(i, j) = norm(gradf(Y_reopt, Dnoise, mask), "fro");
                    mineigvals_reopt(i, j) = min(eig(hessianmatrix(problem_reopt, Y_reopt)));
                    niters_reopt(i, j) = max([info_reopt.iter]);
                    [~, ~, r] = procrust(Y_reopt, Ys);
                    normprocrustdiff_reopt(i, j) = r;
                    solutions_reopt(i, j, :, :) = Y_reopt;
                else % No reoptimization
                    Y_reopt = Y;
                    costs_reopt(i, j) = cost;
                    gradnorms_reopt(i, j) = gradnorms(i, j);
                    mineigvals_reopt(i, j) = mineigvals(i, j);
                    niters_reopt(i, j) = niters(i, j);
                    normprocrustdiff_reopt(i, j) = normprocrustdiff(i, j);
                    solutions_reopt(i, j, :, :) = solutions(i, j, :, :);
                    problem_reopt = problem;
                end

                % Procrustes Y_reopt vs Ys_reopt
                [~, ~, r] = procrust(Y_reopt, Ys_reopt);
                normprocrustdiffysreopt_reopt(i, j) = r;

                % Data for Y_proj
                costs_proj(i, j) = f(Y_proj, Dnoise, mask);
                gradnorms_proj(i, j) = norm(gradf(Y_proj, Dnoise, mask), "fro");
                mineigvals_proj(i, j) = min(eig(hessianmatrix(problem_reopt, Y_proj)));
                [~, ~, r] = procrust(Y_proj, Ys);
                normprocrustdiff_proj(i, j) = r;
                solutions_proj(i, j, :, :) = Y_proj;
            end
            if mod(i, 10) == 0
                fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
                fprintf("k = %d/%d", i + (k - 1) * chunk_size, num_repeats);
                fprintf(" (" + string(toc) + "s)\n");

                save(idstring + "_l=" + dimgt + "_k=" + dimopt + "_n=" + n, ...
                    "dimgt", "dimopt", "n", "ps", "num_ps", "chunk_size", "options", "noiselevel", ...
                    "ground_truths", "inits", "connectedness", "costs_ground_truths", "costs_ground_truths_reopt", ...
                    "costs", "gradnorms", "mineigvals", "niters", "normprocrustdiff", "solutions", ...
                    "costs_proj", "gradnorms_proj", "mineigvals_proj", "normprocrustdiff_proj", "solutions_proj", ...
                    "costs_reopt", "gradnorms_reopt", "mineigvals_reopt", "niters_reopt", "normprocrustdiff_reopt", "normprocrustdiffysreopt_reopt", "solutions_reopt")
            end
        end
    end
end
