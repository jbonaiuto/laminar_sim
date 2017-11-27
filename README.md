MEG laminar simulations
Simulation and analysis code for Bonaiuto, et al
"Non-invasive laminar inference with MEG: Comparison of methods and source
inversion algorithms"

Requirements:
* McNemarextest: http://www.mathworks.com/matlabcentral/fileexchange/6297
* mengz: http://uk.mathworks.com/matlabcentral/fileexchange/37867
* myBinomTest: https://uk.mathworks.com/matlabcentral/fileexchange/24813/
* MEGSurfer: https://github.com/jbonaiuto/MEGsurfer


%%

% To run simulations:

%%

% Create subject structure
  subjects=create_subject_structure();

% Run cross validation whole brain simulations
  simlayer_cross_val(subjects(1), 1, [10 30], -20);

% Run free energy whole brain simulations
  simlayer_free_energy(subjects(1), 1, [10 30], -20);

% Run ROI simulations
  simlayer_roi(subjects(1), 1, [10 30], -20);

% Run free energy - patch size simulations
  simlayer_free_energy_patch_size(subjects(1), 1, [10 30], -20);

% Run ROI - patch size simulations
  simlayer_roi_patch_size(subjects(1), 1, [10 30], -20);

%% 

% To analyze results

%%

% Plot whole brain and ROI stats for each simulation
  plot_wholebrain_roi_sim_results(subjects(1),1,[10 30], -20);

% Plot metric correlations - EBB
  plot_metric_correlations(subjects(1), 1, [10 30], -20, 1);

% Plot metric correlations - MSP
  plot_metric_correlations(subjects(1), 1, [10 30], -20, 4);

% Plot SNR - % pial for whole brain ad ROI, EBB and MSP algorithms
  plot_snr_perc_pial(subjects(1), 1, [10 30], [-100 -50 -20 -5 0 5])

% Plot SNR - % correct for whole brain and ROI, EBB and MSP algorithms
  plot_snr_perc_correct(subjects(1), 1, [10 30], [-100 -50 -20 -5 0 5]);

% Compare classification accuracy of whole brain/ROI and EBB/MSP algorithms
  compare_classifier_accuracy(subjects(1), 1, [10 30], -100);
  compare_classifier_accuracy(subjects(1), 1, [10 30], -50);
  compare_classifier_accuracy(subjects(1), 1, [10 30], -20);
  compare_classifier_accuracy(subjects(1), 1, [10 30], -5);
  compare_classifier_accuracy(subjects(1), 1, [10 30], 0);
  compare_classifier_accuracy(subjects(1), 1, [10 30], 5);

% Plot whole brain analysis - patch size results
  plot_free_energy_patch_size_results(subjects(1), 1, [10 30], [-100 -50 -20 -5 0 5]);

% Plot patch size results as in Troebinger (2014)
  plot_patch_size_luzia(subjects(1), 1, [10 30], -20)

% Plot ROI analysis - patch size results
  plot_roi_patch_size_results(subjects(1), 1, [10 30], [-100 -50 -20 -5 0 5]);

% Plot surface statistics - free energy relationship
  plot_surface_statistic_free_energy('depth', subjects(1), 1, [10 30])
  plot_surface_statistic_free_energy('curvature', subjects(1), 1, [10 30])
  plot_surface_statistic_free_energy('thickness', subjects(1), 1, [10 30])
  plot_surface_statistic_free_energy('lead_field_norm', subjects(1), 1, [10 30])

% Compare surface statistic - free energy correlations
  compare_surface_statistic_free_energy(subjects(1), 1, [10 30]);
