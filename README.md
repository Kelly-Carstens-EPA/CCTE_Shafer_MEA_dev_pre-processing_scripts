# Microelectrode Array Network Formation Assay Pre-Processing Scripts

These scripts can be used to process the experimental data from microelectrode array recordings.

For a step-by-step guide, see `documentation/guide_to_prepare_mea_nfa_level0_for_tcpl.html`. For an example script that runs the functions, see `template/run_me_Template.Rmd`. This guide was designed for pre-processing data from the Network Formation Assay done in the EPA Shafer Lab, but it can be adapted for others.

## Acknowledgements

These scripts are derived from the packages `sjemea` and `meadq`. These packages are availabe on GitHub and contain most of the functions used to calculate the feature values in these scripts. (Links to pages on GitHub: [`sjemea`](https://github.com/sje30/sjemea), [`meadq`](https://github.com/dianaransomhall/meadq) )

## Narrative of the process

### Raw data: spike list files
The inital raw data used is the spike list csv files created by Axion's Spike Detector from microelectrode array recordings. One spike list file is created for each plate and for each day of recording, days in vitro (DIV) 5, 7, 9, and 12. The spike list files are csv files with the columns:

| Time (s) | Electrode | Amplitude(mV) |
| ----------- | ----------- | ----------- |
| ... | ... | ... |

where the first column records the time of each spike, the second column records the ID of the electrode that spiked, and the third column records the amplitude of the spike. Each recording spans ~900 seconds (15 minutes).

### h5 files
The scripts `h5_conversion.R` and `spike_list_functions.R` convert the raw data into the Hierarchical Data Format .h5. This file type is designed to handle large amounts of data. One h5 file is created for each spike list file.

### 16 features values

The scripts `create_burst_ont_Data.R`, `local.corr.all.ont.ae.filter.R`, and `create_ont_csv.R` calculate 16 features from the spike list files. The table below summarizes the 16 features.

Note that these are all well-level values. 

A "burst" on an electrode is a set of spikes that occur in rapid succession. Using the max-interval method with the default set of parameters, a set of spikes must satisfy these conditions in order to be considered a burst (adapted from [Brown *et al*., 2016](https://academic.oup.com/toxsci/article/154/1/126/2422066)):

* The time in between the first 2 spikes in the burst <= 0.1s<br>
* The maximum time in between any 2 spikes in the burst <= 0.25s (This is called the "inter-spike interval," or ISI. There can be an ISI up to 0.8s in a burst if those spikes are followed by another set of spikes with ISI less than 0.1)<br>
* The number of spikes in the burst >= 5 spikes<br>
* The duration of the burst >= 0.05s<br>
* The amount of time in between bursts >= 0.8s

A "network spike" is a group of spikes that occur on several electrodes at the same time. The definition for a network spike was inspired by [Eytan & Marom, 2006](https://www.jneurosci.org/content/26/33/8465). In our scripts, a network spike is calculated as follows:

The entire recording is divided into 0.05 s time bins. A network spike occurs when at least 5 electrodes fire in 1 time bin. The peak of a network spike is the maximum number of electrodes involved in a spike. The duration of a network spike is the length of time between when half of the peak number of electrodes spiked before and after the time at the peak. See `documentation/feature_calcuation_notes.html` for more details.

| Name | Description | Abbreviation | TCPL acnm |
| ----------- | ----------- | ----------- | ----------- |
| Number of Active Electrodes (AE) | # of electrodes where mean firing rate >= 5 spikes/min | nAE | CCTE_Shafer_MEA_dev_active_electrodes_number |
| Mean Firing Rate | # spikes per second, averaged over all AE in each well | meanfiringrate | CCTE_Shafer_MEA_dev_firing_rate_mean |
| Burst Rate | # bursts per minute, averaged over all AE in each well | burst.per.min | CCTE_Shafer_MEA_dev_burst_rate |
| Number of Actively Bursting Electrodes (ABE) | # of electrodes where burst rate >= 0.5 bursts/min | nABE | CCTE_Shafer_MEA_dev_bursting_electrodes_number |
| Mean Burst Duration | mean duration of bursts (s), averaged over all ABE in each well | mean.dur | CCTE_Shafer_MEA_dev_burst_duration_mean |
| Mean Interburst Interval | mean interval between bursts (s), averaged over all ABE in each well | mean.IBIs | CCTE_Shafer_MEA_dev_interburst_interval_mean |
| Interspike Interval in a Burst | mean interspike interval (s) within a burst, averaged over all ABE in each well | mean.isis | CCTE_Shafer_MEA_dev_per_burst_interspike_interval |
| Percent of Spikes in Burst | # of spikes within bursts divided by total spike count, averaged over all ABE in each well | per.spikes.in.burst | CCTE_Shafer_MEA_dev_per_burst_spike_percent |
| Number of Network Spikes | # of network spikes in each well during the 15 minute recording | ns.n | CCTE_Shafer_MEA_dev_network_spike_number |
| Network Spike Peak | max # of electrodes particpating in a network spike, averaged over all network spikes during recording | ns.peak.m | CCTE_Shafer_MEA_dev_network_spike_peak |
| Mean Network Spike Duration | mean duration (s) of all network spikes in each well | ns.durn.m | CCTE_Shafer_MEA_dev_spike_duration_mean |
| Standard Deviation of Network Spike Duration | standard deviation of duration of all network spikes in well | ns.durn.sd | CCTE_Shafer_MEA_dev_network_spike_duration_std |
| Mean correlation | mean Pearson correlations between pairs of AE, averaged on AE | r | CCTE_Shafer_MEA_dev_correlation_coefficient_mean |
| Percent of Spikes in Network Spike | total # of spikes within 0.05 s of the peak of a network spike / total # of spikes in well during recording | ns.percent.of.spikes.in.ns | CCTE_Shafer_MEA_dev_per_network_spike_spike_percent |
| Mean Number of Spikes in Network Spikes |  total # of spikes within 0.05 s of the peak of a network spike / # of network spikes | ns.mean.spikes.in.ns | CCTE_Shafer_MEA_dev_per_network_spike_spike_number_mean |
| Inter-Network Spike Interval | mean time between peaks of consecutive network spikes (s) | ns.mean.insis | CCTE_Shafer_MEA_dev_inter_network_spike_interval_mean |

The scripts will create one csv file containing these feature values will for each plate. Values from each DIV recording will be in separate rows.


### Mutual Information

The mutual information is a robust feature that describes both the global synchrony and level of activity in a network. See *A multivariate extension of mutual information for growing neural networks* by Ball et al., 2017 ([doi](https://doi.org/10.1016/j.neunet.2017.07.009), [link](https://www.sciencedirect.com/science/article/abs/pii/S0893608017301612?via%3Dihub) for more information. The scripts `spikeLoadRoutines.R`, `nmi_wrapper.R`, and `nmi2_final.R` contain the functions used to calculate the mutual information.

| Name | Description | Abbreviation | TCPL acsn |
| ----------- | ----------- | ----------- | ----------- |
| Normalized Mutual Information | concurrently measures synchrony and activity of the neural network | mi | CCTE_Shafer_MEA_dev_mutual_information_norm |

The calculation of the mutual information is computationally intensive. Therefore, this feature is calculated separately from the rest of the features. The script `MI_script_all.R` is designed to calculate the mutual information for all plates. The scripts will create one csv file containing the mutual information for each plate. 


### Area Under the Curve

We want to quantify the alterations to development over time in treated wells compared to control wells. In order to "sum up" the overall changes in a feature value, we calculate the trapezoidal Area Under the Curve (AUC). See this [example](https://github.com/amycarpenter/CCTE_Shafer_MEA_dev_pre-processing_scripts/blob/master/images/meanfiringrate_development_example.jpeg) image of the development of the mean firing rate in a given well over time. The AUC value will be used to compare the overall increase or decrease of a feature in treated wells versus control wells.

The script `parameter_values_to_AUC.R` uses the `trapz` function from the `pracma` package to calculate the trapezoidal AUC for each feature.

Notes:

* When there are no bursts or network spikes in a well, many features that measure some aspect of bursts or network spikes are NA. In order to calculate the area under the curve, these NA values are set to 0. Below is a list of the features that are sometimes NA. Setting NA values to 0 might be reconsidered in the future for some of these endpoints. 

  * Mean Burst Duration
  * Network Spike Peak<br>
  * Network Spike Duration<br>
  * Mean Number of Spikes in Network Spikes<br>
  * Interspike Interval<br>
  * Mean Interburst interval<br>
  * Interspike Interval in Network Spikes<br>
  * Standard Deviation of Network Spike Duration<br>

* Optionally, data points at DIV 2 (of 0 value) are added before calculating the AUC for consistency with data historically produced in the Shafer lab. 


### Cytotoxicity data

After the cells are grown on the plates for 12 days, two assays are used to assess to the cell viability - CellTiter-Blue and Lactate Dehydrogenase.

The CellTiter-Blue assay (also called Alamar Blue) measures the amount of reagent metabolized by living cells in each well. First, reagent is added to each well. Then, some media from each well is transferred to an opaque 96-well plate. The fluorescence of resazurin, a metabolite of the reagent, is measured. Three blank wells that contain only reagent are used as a baseline fluorescence values. The average of the fluorescence in the three blank-corrected wells is subtracted from the raw fluorescence values in the remaining 48 wells. (See for example [Brown *et al.* (2016)](https://academic.oup.com/toxsci/article/154/1/126/2422066) for more information). In the ToxCast Pipeline, the blank-corrected values in the treated wells are normalized to the median value in control wells. 

The total lactate dehydrogenase (LDH) assay is also used to quantify the number of living cells in treated wells versus control wells. First, all of the media is removed from each well in order to remove any LDH already released from dying cells. Then, a lysis solution is added to lyse all living cells. Next, the lysis solution is transferred to another plate with a solution containing tetrazolium salt. The LDH in the solution is allowed to transform the tetrazolium salt in red formazan for 30 minutes until a stop solution is added. The optical density of red formazan is measured. The amount of red formazan produced reflects the amount of LDH released from the cells that were living at the start of the assay. Three blank wells containing only lysis solution are used as a baseline optical density value. This information was synthesized from [Frank *et al.* (2017)](https://academic.oup.com/toxsci/article/160/1/121/4083261), the product description for [CytoTox 96 Non-Radioactive Cytotoxicity Assay](https://www.promega.com/products/cell-health-assays/cell-viability-and-cytotoxicity-assays/cytotox-96-non_radioactive-cytotoxicity-assay/?catNum=G1780), and additional summary data on the LDH assay as done in the Shafer lab. In the ToxCast Pipeline, the blank-corrected values in the treated wells are normalized to the median value in control wells.

For both assays, the script `run_cytotox_functions.R` handles the steps to extract the blank-corrected values from excel sheets. Any negative blank-corrected values are set to zero.

### Data transformations and checking

Finally, the `run_me_Template.Rmd` script contains many suggested steps to transform and clean the output .csv files input a "level 0" table that is ready for hit-calling and curve-fitting with the ToxCast Pipeline R package.
