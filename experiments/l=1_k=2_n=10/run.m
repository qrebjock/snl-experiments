clear; clc;

addpath("../..")

generategt = @(n, dimgt) center(randn(n, dimgt));

initscheme = @(n, dimopt) center(randn(n, dimopt));

options.verbosity = 0;
options.maxiter = 2000;
options.tolgradnorm = 1e-8;

dimgt = 1;
dimopt = 2;
n = 10;
ps = linspace(0, 1, 41);
num_repeats = 10000;
noiselevel = 0;

erdosrenyiexperiment(dimgt, dimopt, n, ps, num_repeats, generategt, initscheme, options, noiselevel);
