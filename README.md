# Sensor network localization experiments

This repository contains the code to reproduce the experiments from the paper:

> *Sensor network localization has a benign landscape after low-dimensional relaxation*
>
> Christopher Criscitiello, Andrew D. McRae, Quentin Rebjock and Nicolas Boumal
>
> <https://arxiv.org/abs/2507.15662>

It requires [Manopt](https://www.manopt.org/) version 8.0.0 or later.

The subfolders in `experiments/` contain the experimental setups.
Run the `run.m` files in those folders to generate `.mat` files containing the results.

Use the `plot_results.ipynb` notebook to visualize the results.
Figures will be saved in the `figures/` folder.

**Note:** Accurately reproducing the figures requires repeating the experiments many times, which may take a significant amount of time.

The file `/dev/random` is used to set the RNG seed.
If your system does not provide this file, see the comments in `erdosrenyiexperiment.m` for instructions on how to adapt the code accordingly.

