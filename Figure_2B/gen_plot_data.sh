./pycbc_inference_plot_posterior --input-file test_superluminal.hdf --plot-scatter --output-file scatter_2param_mk2.pdf --z-arg snr --parameters non_gr_fac_n "log10(non_gr_fac_a_at_30)"
mv figure2_data.npz figure2_data_subluminal.npz

./pycbc_inference_plot_posterior --input-file test.hdf.bkup --plot-scatter --output-file scatter_2param_mk2.pdf --z-arg snr --parameters non_gr_fac_n "log10(non_gr_fac_a_at_30)"

mv figure2_data.npz figure2_data_subluminal.npz

