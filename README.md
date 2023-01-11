---
output:
  pdf_document: default
  word_document: default
  html_document: default
---
# Microelectrode Array Network Formation Assay Pre-Processing Scripts

## Purpose

These scripts are designed to process the raw data for the microelectrode array network formation concentration-response assay.

Briefly, cortical cells are grown on 48-well microelectrode-containing plates. On each plate, 6 compounds are tested at 7 concentrations, plus 6 control wells. Each group of compounds and concentrations is replicated on 3 plates. The test compounds are added to the wells on days 0, 5, and 9. The electrical activity of the neurons are recorded as the network develops, on days 5, 7, 9, and 12. The recordings provide a lot of information concerning the general activity, organization, and connectivity of the neurons in each well. These scripts calculate several features from the recordings in order to describe these 3 aspects of the developing neuronal networks.

We want to plot the dose response data for each compound and each endpoint, calculate an EC50 value, and determine if the compound is a positive hit. These scripts transform the raw the recording data into a long file format. Then, the data can be processed with the functions in the ToxCast Pipeline to fit either a hill, gain-loss, or constant model to the dose-response plot and determine the hit calls for each compound.

## Acknowledgements

These scripts are derived from the packages `sjemea` and `meadq`. These packages are availabe on GitHub and contain most of the functions used to calculate the feature values in these scripts. (Links to pages on GitHub: [`sjemea`](https://github.com/sje30/sjemea), [`meadq`](https://github.com/dianaransomhall/meadq) )

## How to use these scripts

For a step-by-step guide, see the document *Step-by-Step_Guide.docx*. Checkout this [diagram](https://ncct-bitbucket.epa.gov/projects/NSLTM/repos/nfa-spike-list-to-mc0-r-scripts/browse/images/SpikeList_to_mc0_overview.jpg) to visualize the steps. Raw data files are shown in blue, intermediate output files are in purple, and scripts are in orange.

## Narrative of the process

### Raw data: spike list files
The inital raw data used is the spike list csv files created by Axion's Spike Detector from a micro-electrode array recording. One spike list file is created for each plate and for each day of recording, days in vitro (DIV) 5, 7, 9, and 12. The spike list files are csv files with the columns:

| Time (s) | Electrode | Amplitude(mV) |
| ----------- | ----------- | ----------- |
where the first column records the time of each spike, the second column records the ID of the electrode that spiked, and the third column records the amplitude of the spike. Each recording spans ~900 seconds (15 minutes). (If a recording goes over 900 seconds, the script `spike_list_functions.R` will truncate the spike list file at 900 seconds).

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

The entire recording is divided into 0.05 s time bins. A network spike occurs when at least 5 electrodes fire in 1 time bin. The peak of a network spike is the maximum number of electrodes involved in a spike. The duration of a network spike is the length of time between when half of the peak number of electrodes spiked before and after the time at the peak. See *feature_calcuation_notes.md* for more details.

| Name | Description | Abbreviation | TCPL acsn |
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

One csv file containing these feature values will be created for each plate. Values from each DIV recording will be in separate rows.

Note:

- The scripts currently calculate 2 additional features: cv.time and cv.network. However, there are concerns about how these features are calculated. These features should be ignored. These values are removed from the data in the next script, `burst_parameter_to_AUC.R`.


### Mutual Information

The mutual information is a robust feature that desribes both the global synchrony and level of activity in a network. See [*A multivariate extension of mutual information for growing neural networks*](https://www.sciencedirect.com/science/article/abs/pii/S0893608017301612?via%3Dihub) by K. Ball et al. for more information. The scripts `spikeLoadRoutines.R`, `nmi_wrapper.R`, and `nmi2_final.R` contain the functions used to calculate the mutual information.

| Name | Description | Abbreviation | TCPL acsn |
| ----------- | ----------- | ----------- | ----------- |
| Normalized Mutual Information | concurrently measures synchrony and activity of the neural network | mi | CCTE_Shafer_MEA_dev_mutual_information_norm |

The calculation of the mutual information is computationally intensive. Therefore, this feature is calculated separatley from the rest of the features. The script `MI_script_all.R` is designed to calculate the mutual information for all plates, so that the task could be done overnight or remotely for the entire dataset. One csv file containing the mutual information values will be created for all plates from same culture date.


### Area Under the Curve

We want to quatify the alterations to development from DIV 0 - 12 in treated wells compared to control wells. In order to "sum up" the overall changes in a feature value, we calculate the trapezoidal area under the curve. See this [example](https://ncct-bitbucket.epa.gov/projects/NSLTM/repos/nfa-spike-list-to-mc0-r-scripts/browse/images/meanfiringrate_development_example.jpeg) image of the development of the mean firing rate in a given well over time. This value will be used to compare the overall increase or decrease of a feature in treated wells versus control wells.

The script `burst_parameter_to_AUC.R` uses the `trapz` function from the `pracma` package to calculate the trapezoidal area under the curve (AUC) for each feature. One csv file will contain the AUC values for all plates and features.

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

* Historically, activity was recorded on DIV 2. There was usually very little activity in these recordings. At some point, it was decided that all DIV 2 feature values should be set to 0 before calculating the AUC in `burst_parameter_to_AUC.R`. Now, even though activity is no longer recorded on DIV 2, we still calculate the area under the curve with the first point at DIV = 2, feature value = 0.0 (instead of starting at DIV 5).


### Cytotoxicity data

After the cells are grown on the plates for 12 days, two assays are used to assess to the cell viability - CellTiter-Blue and Lactate Dehydrogenase.

The CellTiter-Blue assay (also called Alamar Blue) measures the amount of reagent metabolized by living cells in each well. First, reagent is added to each well. Then, some media from each well is transferred to an opaque 96-well plate. The fluorescense of resazurin, a metabolite of the reagent, is measured. Three blank wells that contain only reagent are used as a baseline fluorescense values. The average of the fluorescence in the three blank-corrected wells is substracted from the raw fluorescence values in the remaining 48 wells. (See [Brown *et al.* (2016)](https://academic.oup.com/toxsci/article/154/1/126/2422066) for more information). The blank-corrected values in the treated wells will be normalized to the median value in control wells in the ToxCast Pipeline. 

The total lactate dehydrogenase (LDH) assay is also used to quantify the number of living cells in treated wells versus control wells. First, all of the media is removed from each well in order to remove any LDH already released from dying cells. Then, a lysis solution is added to lyse all living cells. Next, the lysis solution is transfered to another plate with a solution containing tetrazolium salt. The LDH in the solution is allowed to transfrom the tetrazolium salt in red formazan for 30 minutes until a stop solution is added. The optical density of red formazan is measured. The amount of red formazan produced reflects the amount of LDH released from the cells that were living at the start of the assay. Three blank wells containing only lysis solution are used as a baseline optical density value. This information was synthesized from [Frank *et al.* (2017)](https://academic.oup.com/toxsci/article/160/1/121/4083261), the product description for [CytoTox 96 Non-Radioactive Cytotoxicity Assay](https://www.promega.com/products/cell-health-assays/cell-viability-and-cytotoxicity-assays/cytotox-96-non_radioactive-cytotoxicity-assay/?catNum=G1780), and from the LDH Assay Summary file *CCTE_Shafer_MWP_LDHo_dn_final.docx*. The blank-corrected values in the treated wells will be normalized to the median value in control wells in the ToxCast Pipeline.

For both assays, the script `cytotox_prep06.R` extracts these blank-corrected values from excel sheets created by the lab technicians. Any negative blank-corrected values are set to zero.

### Format data into long file

The script `tcpl_MEA_dev_AUC.R` formats all of the AUC and cytotoxicity data into one long file with the columns needed for "level 0" data in the ToxCast Pipeline. 
This script will also:<br>
- Set wllt ("well type") to 't' for all treated wells and 'n' for all control wells (where conc = 0)
- Set the treatment column to the corresponding vehicle control for control wells
- Check if any plate ID was used in multiple culture dates. An alphabetical suffix will be added to re-used plate ID's to differentiate the data.

### Make adjustments to the data set if needed

Once all of the cytotoxicity and ontogeny data is combined, you can make adjustments to the entire data set as needed. For example,<br>
- If there was an issue with any wells (e.g. contamination, misdose, etc.), set wllq = 0 for those wells. These data rows will be removed in level 2 in the ToxCast Pipeline.
- Remove any values that should not be included (e.g. unregistered compounds, or rows of plates that were already pipelined with another data set).
These tasks are currently not integrated into any of the scripts in this repository, so custom scripts must be made to make these adjustments.

### Replace treatment column with sample ID

The script `spid_mapping.R` maps the treatment names to the corresponding sample IDs. The sample ID is the unique ID corresponding to registered samples of each compound in the ToxCast database. In the future, the chemical names in the Master Chemical Lists may be replaced with the spid's in order to remove this step.

The output from `spid_mapping.R` should be ready for the ToxCast Pipeline.

