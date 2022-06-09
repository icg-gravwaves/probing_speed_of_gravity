# Probing the speed of gravity with LVK, LISA and joint observations - Data Release

This is the data release for the paper "Probing the speed of gravity with LVK, LISA, and joint observations", which can be found
\[ADD ARXIV LINK WHEN ON ARXIV\]. We release the data behind each of the figures and the code and configuration files that were used to produce each of them.

As always with the python landscape, software dependencies will change in the coming years, and this code may be incompatible with future releases of any dependency. This code was checked against the v2.0.3 PyCBC Docker image, which can be found here https://zenodo.org/record/6583784. The docker image provides a full suite of software and it should be possible to "just use this", with the changes and additions documented below.

## Supplementary data and images

We release the following supplementary plots. The data used to create each of the plots shown in the manuscript is discussed in the sections below.

* An animated version of Figure 2B [Can be seen here](https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2B/figure2b.mp4)
* Figure 2A, produced using superluminal data [Can be seen here](https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2A/figure2a_sl.png)
* Figure 3, produced using superluminal data [Can be seen here](https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_3/figure3_sl.png)

We highlight that the superluminal versions of figure 2A and figure 3 are shown to demonstrate that these are identical to the subluminal versions of these plots
shown in the manuscript.

## Figure 1 (and Figure 4 in appendix)

Figure 1 (and Figure 4 in the appendix) is generated by the python script which can be downloaded by clicking here:

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_1/make_plots.py>

It also requires the matplotlib style files, which can be downloaded from here:

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/matplotlib_style_files/paper.mplstyle>

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/matplotlib_style_files/presentation.mplstyle>


## Figure 2A

Figure 2A is generated in a few steps. First we generate the data by running the following two scripts:

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2A/make_figure_data_LIGO.py>

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2A/make_figure_data_ET.py>

This produces `.npz` data files, which are also made available here, as these scripts can take hours to run.

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2A/figure2_dataLIGO.npz>

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2A/figure2_dataET.npz>

The plot is then produced by running the script in

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2A/make_plot.py>

This plot uses the subluminal results, but you can see that the superluminal results are identical by running the script here

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2A/make_plot_sl.py>

to generate the figure with the superluminal data points.

## Figure 2B

Figure 2B is created using an edited version of PyCBC Inference.

Please see here:

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2B>

for details of this figure.

## Figure 3

Figure 3 uses a similar script to Figure 2A but must examine a lot more points. The scripts for generating the LIGO and LISA data respectively are here:

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_3/make_data_ligo.py>

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_3/make_data_lisa.py>

These files are split into chunks so you must run

```
python make_data_ligo.py 0
python make_data_ligo.py 1
...
python make_data_ligo.py 24
```

and similar for the lisa script. In my case I ran these in parallel on a single machine.

The outputs can then be combined using `cat` to produce `LIGO_output.txt` and `LISA_output.txt`. These two files are available here:

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_3/LIGO_output.txt>

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_3/LISA_output.txt>

Figure 3 is then made from these input files by running:

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_3/make_plots.py>

To generate the same plot using superluminal data, you can use:

<https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_3/make_plots_sl.py>

