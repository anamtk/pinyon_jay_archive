# Title: Data and code for models of relationships between pinyon jay abundance and cone availability using eBird

# Abstract

In this project, we examined the relationship between pinyon jay abundance and two-needle pinyon pine cone availability across the SW United States. We used data from eBird and published data on cone availability and stochastic antecedent modeling to explore the temporal aspects of this species interaction. This repository contains all the data and code to reproduce our data cleaning, analysis, and visualizations.

# Creators

| First Name | Middle Initial | Last Name | Organization | email |
|----|----|----|----|----|
| Ana | T | Miller-ter Kuile | Northern Arizona University | ana.miller-ter-kuile\@nau.edu |
| Andreas | P | Wion | Forest Stewards Guild | awion\@forestguild.org |
| Kyle | C | Rodman | Ecological Restoration Institute | kyle.rodman\@nau.edu |
| Kiona |  | Ogle | Northern Arizona University | kiona.ogle\@nau.edu |
| Jamie |  | Sanderlin | Rocky Mountain Research Station | jamie.l.sanderlin\@usda.gov |

# License

CC-By. Please cite this data package if using code.

# Keywords

*Gymnorhinus cyanocephalus*, occupancy model, *Pinus edulis*, resource availability, seed production, species interactions, temporal legacies, time lag

# Funding of this work

| PI | Title of Grant | Funding Agency |
|----|----|----|
| Jamie Sanderlin, Kiona Ogle | Application and development of data integration methods for focal wildlife species’ distribution models for management applications in ecology | USDA Forest Service |

# Timeframe

-   begin date: 2010
-   end date: 2023
-   Data collection ongoing/completed: Completed

# Geographic location

-   Verbal description: Southwest Region of the United States of America
-   North bounding coordinate: 41.975007
-   South bounding coordinate: 31.249435
-   East bounding coordinate: -101.896744
-   West bounding coordinate: -114.047623

# Taxonomic species or groups

pinyon jay (*Gymnorhinus cyanocephalus*); two-needle pinyon pine (*Pinus edulis*)

# Methods

In this study, we quantify the relationship between pinyon jay and pinyon pine seed (cone) production by leveraging two broad-scale, long-term datasets from the pinyon pine-pinyon jay system. We combine a large (\~800,000 observations) eBird dataset and estimates of pinyon pine cone production; these datasets cover more than half the range of pinyon jays (Figure 1), the full range of pinyon pine, and capture over a decade (2010-2022) of bird abundance and cone production. This observation period is also characterized by high variation in climate and habitat conditions . We explicitly model the temporal dynamics of this mutualism using a stochastic antecedent modeling (SAM) framework to examine how populations of pinyon jay respond to resource availability (cones) and climate across space and time. We address two key questions: 1) Do pinyon jays exhibit anticipatory, immediate, or delayed relationships with pinyon pine cone availability, and 2) How do interactions between pinyon jays and pinyon pine vary across a range of habitat quality and climate conditions? Addressing these questions is not only necessary for conservation of this system but also bolsters our understanding of how species interactions shape population trends.

To quantify coupling between pinyon jay abundance and pinyon cone (seed) production across various timescales and seasonal periods, we fit an N-mixture (occupancy) model to the jay abundance data. The model explicitly accounts for the effects of spatially and temporally varying pinyon seed availability, climate conditions (e.g., seasonal temperature and precipitation), and pinyon basal area, and site-specific monsoonality, which describes the contribution of summer rainfall to total annual precipitation and shapes alternative food availability for jays (e.g., invertebrates). Importantly, we included interactions between cone production and the other predictor variables (precipitation, temperature, basal area, and monsoonality) to understand how these climate and habit factors modulate the coupling between pinyon jays and cone production. Importantly, we simultaneously estimated the most influential seasonal and lagged periods for which jay abundance is most strongly coupled to climate and cone production, thus quantifying the relative importance of anticipatory, immediate, and delayed responses.

More methods will be found in the forthcoming manuscript associated with this dataset.

# Data Provenance

| Dataset title | Dataset DOI or URL | Creator | Contact |
|----|----|----|----|
| eBird Basic Dataset - May 2024 | <https://science.ebird.org/en/use-ebird-data/download-ebird-data-products> | Cornell Lab of Ornithology |  |
| eBird Basic Dataset - June 2024 | <https://science.ebird.org/en/use-ebird-data/download-ebird-data-products> | Cornell Lab of Ornithology |  |
| Two-needle pinyon pine yearly cone production data | <https://datadryad.org/dataset/doi:10.5061/dryad.7h44j103t> | Andreas Wion | awion\@forestguild.org |
|  |  |  |  |

# Data Tables

This repository includes many (\~75-100) data files. Many are produced as downstream products of data processing steps. Rather than include metadata on all of these data files, I include general descriptions of all of the data folders and their contents.

## Data Folder

## Folder: 01_ebird_data

This folder contains both raw and processed data from eBird. The subfolders that begin with **ebd** are the datasets directly downloaded from eBird; the **pinjay_range_2023** folder contains processed range maps from EBird Status and Trends (<https://science.ebird.org/en/status-and-trends>).

## Subfolder: cleaned_data

### Subfolder: 01_auk_filtering

This subfolder contains outputs from filtering eBird data with the auk package in R

### Subfolder: 02_all_auk_filtered

This folder contains a dataset that combines the results of each state-level filtering process in the above folder

### Subfolder: 03_subsampled

This folder contains data (both .csv and spatial) of eBird data that was processed using eBird standard practices (<https://ebird.github.io/ebird-best-practices/>) to reduce spatial and temporal sampling bias. **These are the data that we used in our regression model exploring relationships between jays and cones**.

### Subfolder: 04_JAGS_indexIDs

This folder contains indexed data to match eBird checklist IDs to the IDs used to index in the JAGS Bayesian regression model.

### Subfolder: 05_oos

This folder contains the dataset that we used to test our model's predictive accuracy using out-of-sample data. These datasets are the products of the exact same processing that the subsampled data in **03_subsampled** went through, and explicitly exclude those checklists that were used to train the model.

## Folder: 02_spatial_data

This folder contains spatial data for all covariates and also cleaned datasets that link spatial information from those datasets to the eBird dataset

### Subfolder: cleaned_data

This subfolder contains all the cleaned data and cleaning steps to get final covariate data.

#### Sub-subfolder: 00_metadata

This folder contains a dataset that has all the unique grid cell IDs for rasters to be able to link raster cell/grid IDs across eBird and covariate datasets

#### Sub-subfolder: 01_get_gridIDs

These datasets are the product of processing steps and link original spatial data for each covariate to the cell IDs present for eBird checklist locations.

#### Sub-subfolder: 02_weighted_blob_covariates

These datasets are processed from the **01_get_gridIDs** datasets and include weighted average covariate values at the level of the multi-polygon surrounding each set of eBird checklists in a location and their sampling location uncertainty buffer.

#### Sub-subfolder: 03_oos

These datasets are the repeat of the **00_metadata** and **02_weighted_blob_covariates** datasets for the out of sample data. These datasets were produced in the same way as those used for the training data.

### Subfolder: masting_data

This folder contains both the original (quantile_pied_predictions_scaled.tif) and processed (ConePredictions_final.tif) datasets used for cone abundance in the model. The processed dataset is a product of updating dates to match the yearly values in other datasets (e.g., eBird, basal area).

### Subfolder: monsoon

This folder contains the spatial data on monsoon moisture percent for the SW region.

### Subfolder: pinyonBA

This folder contains the pinyon basal area data that was combined from two sources and then annualized for use in analyses (final dataset used: PinyonBA_4km_sqmPerHa.tif).

### Subfolder: prism_monthly_ppt

This folder contains many subfolders that include gridded monthly precipitation data downloaded from PRISM.

### Subfolder: prism_monthly_tmax

This folder contains many subfolders that include gridded monthly maximum temperature data downloaded from PRISM.

### Subfolder: prism_monthly_tmean

This folder contains many subfolders that include gridded monthly mean temperature data downloaded from PRISM.

## Folder: 03_jags_input_data

These are data list objects necessary to run a JAGS model, including the input data list and the list of initial conditions for the model

### Subfolder: oos

These are the same data as the data list for the JAGS model, but prepped for the out-of-sample dataset

## Folder: 04_jags_output_data

This folder contains covariate effect means for covariates corrected in a script to account for the weighting of covariates in the SAM model

## Folder: 05_cross_validation

This folder contains the out-of-sample RMSE values we used to compare to the test data RMSE values for model validation.

## Folder: 06_variable_importance

This folder contains data derived from computing the variable importance for each individual variable using permuted variable importance.

## Monsoon Folder

## Folder: inputs

The subfolder in this folder **00_data_lists** contains replicates of the JAGS input data list and the initial values list described above

## Folder: outputs

This folder contains outputs from running the JAGS model on a computing cluster ("Monsoon") at Northern Arizona University. The objects include:

-   ebird_abund_model2_Rhat.RDS - convergence statistics for the model

-   ebird_abund_model2_summary.RDS - summary of covariate effects for the model

-   ebird_abund_model_yyrepr2.RDS - samples of R2 between y and replicated data for the full model

-   ebird_abund_model_yrep.RDS - replicated data summary

-   ebird_abund_model_RMSE_samples.RDS - RMSE samples to evaluate model accuracy

-   ebird_abund_model_residuals.RDS - residuals summary to test for spatial autocorrelation

-   ebird_abund_model_covariate_effect_samples.RDS - samples of covariate effects used to standardize covariate effects based on weights for covariate effect visualization

# Spatial Data Objects

Several of the datasets listed above are spatial datasets, but we list the spatial datasets here explicitly

#### 01_ebird_data -\> cleaned_data -\> 03_subsampled

**all_ebird_data_buffercellIDs.shp**

The eBird checklists including the cellIDs for all raster grid cells in the spatial uncertainty in sampling location, based on the travel distance for each checklist

**all_ebird_data_conefiltered.shp**

Legacy dataset that doesn't buffer around each checklist point

#### 01_ebird_data -\> cleaned_data -\> 05_oos

**all_oos_ebird_data_buffercellIDs.shp**

The same dataset as that for the test data, but processed for the out-of-sample dataset

#### 02_spatial_data -\> masting_data

These are spatial yearly tif files for cone production (processed data is in **ConePredictions_final.tif**

#### 02_spatial_data -\> monsoon

These are spatial data (.tif) on monsoon percent contribution across the SW region

#### 02_spatial_data -\> pinyonBA

These are spatial data on pinyon basal area in the SW, including a final processed dataset used in our analyses (**PinyonBA_4km_sqmPerHa.tif**)

#### 02_spatial_data -\> prism_monthly_ppt

These are data downloaded using the prism package in R on monthly spatial precipitation values

#### 02_spatial_data -\> prism_monthly_tmax

These are data downloaded using the prism package in R on monthly spatial maximum temperature values

#### 02_spatial_data -\> prism_monthly_tmean

These are data downloaded using the prism package in R on monthly spatial mean temperature values

# Scripts/code

| Name | Location | What it does | Language |
|----|----|----|----|
| ebird.filtering.R | code/01_data_prep/01_ebird | uses auk to filter eBird checklists based on eBird best practices | R |
| 00_FormattingSpeciesBA_Data.R | code/01_data_prep/02_spatial | annualizes pinyon basal area data and updates year IDs in cone dataset | R |
| 00_importprism.R | code/01_data_prep/02_spatial | imports climate data from PRISM | R |
| 01_subsample_ebird_and_buffer.R | code/01_data_prep/02_spatial | spatial and temporally subsample eBird checklists based on best practices and make spatial buffers based on sampling distances | R |
| 02_buffered_ebird_gridIDs.R | code/01_data_prep/02_spatial | get metadata on grid cell IDs in each eBird sampling multi-polygon | R |
| 03_get_gridIDs_covariates.R | code/01_data_prep/02_spatial | Link grid and multi-polygon IDs from eBird to the covariate datasets | R |
| 04_blob_weighted_covariates.R | code/01_data_prep/02_spatial | get weighted average covariate values for the raster grid cells covered by each eBird multi-polygon | R |
| ebird_jags_data_prep_checkN.R | code/01_data_prep/03_model | preps the data list with data structures required to run analyses in JAGS | R |
| 01_subsample_oos_ebird_and_buffer.R | code/01_data_prep/04_out_of_sample | Perform subsampling for out of sample dataset | R |
| 02_buffered_oos_ebird_gridIDs.R | code/01_data_prep/04_out_of_sample | Buffers and grid IDs for out of sample dataset | R |
| 03_oos_blob_weighted_covariates.R | code/01_data_prep/04_out_of_sample | get weighted covariate values linked to grid IDs for the out of sample dataset | R |
| 04_oos_data_list_prep_checkN.R | code/01_data_prep/04_out_of_sample | Prep data list for the out of sample function to validate model | R |
| 00_check_cov_distributions.R | code/02_analyses | check that covariate distributions of background and checklists are \~equal to be able to use eBird data in standardized N-mixture model | R |
| 01_check_convergence.R | code/02_analyses | check convergence of JAGS model | R |
| 02_test_spatial_autocorrelation.R | code/02_analyses | test that residuals aren't spatially autocorrelated | R |
| 03_oos_predict_function.R | code/02_analyses | Function for model validation which uses model outputs and out of sample dataset | R |
| 04_permuted_variable_importance_function.R | code/02_analyses | Function to determine variable importance via permutation | R |
| figx_cone_weights.R | code/03_visualization | figure of yearly cone weights | R |
| figx_covariates_effects.R | code/03_visualization | figure of all covariate effects | R |
| figx_interaction_plots.R | code/03_visualization | figure of interactions between cones and other variables | R |
| figx_maps.R | code/03_visualization | map figure with 3 panels | R |
| FigX_R2_contributions.R | code/03_visualization | contributions of different covariate groups to R2 | R |
| sifigx_oos_test_RMSE_graphs.R | code/03_visualization | graph of test and oos RMSE | R |
| sifigx_weights.R | code/03_visualization | figure of seasonal weights for temperature and precipitation | R |
| SIfigx_y_yrep_linear_R2.R | code/03_visualization | figure of goodness of fit based on relationships between observed and replicated data in the model | R |
| ebird_abund_JAGS_checkN.R | code/jags_models | JAGS Bayesian model to test relationship between jay abundance and cone abundance | JAGS |
| ebird_abund_JAGS_checkN.R | monsoon/inputs/01_run_model | replicate JAGS Bayesian model | JAGS |
| script_checkN.R | monsoon/inputs/01_run_model | JAGS wrapper script to run the JAGS model | R |
| script2_checkN.R | monsoon/inputs/01_run_model | JAGS wrapper script to include initials from first script and re-run JAGS model | R |
| GOF.R | monsoon/inputs/02_coverged_model | get replicated data, RMSE, covariate samples from converged model | R |
| summary.R | monsoon/inputs/02_coverged_model | get summary of covariate effects from converged model | R |
| y_yrep_r2.R | monsoon/inputs/02_coverged_model | get R2 value between observed and replicated data for the converged model | R |
| ebird_abund_JAGS_checkN_climateonly.R | monsoon/inputs/03_variable_importance | JAGS model including only climate main effects | JAGS |
| ebird_abund_JAGS_checkN_habitatonly.R | monsoon/inputs/03_variable_importance | JAGS model including only habitat main effects | JAGS |
| GOF_climateonly.R | monsoon/inputs/03_variable_importance | get replicated data for the climate only model | R |
| GOF_habitatonly.R | monsoon/inputs/03_variable_importance | get replicated data for the habitat only model | R |
| script_checkN_climateonly.R | monsoon/inputs/03_variable_importance | run the climate only model | R |
| script_checkN_habitatonly.R | monsoon/inputs/03_variable_importance | run the habitat only model | R |
| script2_checkN_climateonly.R | monsoon/inputs/03_variable_importance | update the climate only model with initials from previous model run and rerun | R |
| script2_checkN_habitatonly.R | monsoon/inputs/03_variable_importance | update the habitat only model with initials from the previous model run and rerun | R |
| y_yrep_r2_climateonly.R | monsoon/inputs/03_variable_importance | R2 between replicated and observed data for the climate only model | R |
| y_yrep_r2_habitatonly.R | monsoon/inputs/03_variable_importance | R2 between replicated and observed data for the habitat only model | R |
