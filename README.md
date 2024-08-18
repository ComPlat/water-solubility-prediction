# Water Solubility Prediction Study

## Overview

This repository contains the code, data, and models used in a comprehensive study on water solubility prediction. The study focuses on training and testing various machine learning models using curated datasets from literature and benchmark datasets, comparing the performance of different descriptors, and evaluating the models against existing state-of-the-art methods.

## Datasets

The study utilizes four primary datasets for training and testing, which are sourced from the literature. These datasets are combined, preprocessed, and then split into final training and test sets. The datasets used are:

1. **water_solubility_data.csv** (900 samples)
2. **dataset-not-FA.csv** (6,154 samples)
3. **Supplementary_data.csv** (9,943 samples)
4. **data_paper.csv** (11,862 samples)
5. **dataset-E.csv** (1,282 Test set - samples)

The combined dataset results in a total of 28,859 datapoints.

### Preprocessing

Preprocessing of the datasets is conducted in the `data_preprocess_3_24.ipynb` notebook. The steps include:
- **Removing duplicates**: Ensuring that no duplicate entries exist within the combined dataset.
- **Canonical SMILES**: Generating canonical SMILES strings to identify and remove matching data points between the training and test datasets.
- **Final Dataset**: After preprocessing, the final datasets are:
  - Training data: **unique_train4_new24.csv** (16,626 samples)
  - Test data: **unique_test_new24.csv** (1,282 samples)

## Model Training

Once the training and test datasets were finalized, various models were trained and evaluated. The training process and model evaluation are detailed in the `4_24.ipynb` notebook.

### Descriptors and Feature Engineering

The study explored a wide range of descriptors to improve model performance:
- **Basic to Advanced Descriptors**: The `utilities.py` script generates combinations of descriptors ranging from 4 basic descriptors to 123 advanced descriptors, including fingerprints of varying lengths (from 128 bits to 1024 bits).
- **Feature Engineering**: The `feature_fg7_fe_38.ipynb` notebook introduces 38 feature-engineered descriptors and 7 functional group descriptors.
## Complete Descriptors

The complete list of descriptors used in this project is provided in the [descriptors_detail.md](descriptors_detail.md) file.

### Model Evaluation

The best-performing model and its parameters were selected based on minimizing the Mean Absolute Error (MAE). The comparative results are stored in the `model_results.csv` file. Additionally, the performance of the Message Passing Neural Network (MPNN) model is evaluated and stored in the `mpnn.ipynb` notebook.

## Comparison with Sorkun's Work

The study also includes a detailed comparison with Sorkun et al.'s work, which used the same benchmark dataset (`dataset-E.csv`). A crucial finding was the presence of overlapping data between the training and test sets in Sorkun's preprocessing. After removing these overlaps, the MAE increased from 0.35 to 0.54, indicating the importance of proper data preprocessing.

This analysis is documented in the `Sorkundata_improve_Preprocess.ipynb` notebook, with the overlapping compounds listed in the `overlap_data_new.csv` file.

## External Evaluation

The model was further evaluated using other online prediction tools, such as VCC Lab, and compared against Sorkun's model using self-experimented compound solubility data. These compounds were experimentally tested in the lab to obtain their solubility values. The experimental details and results are provided in the `Sol_exp` folder which has excel file for more detail abaout the Solubility experimet.The results of this comparison are saved in the `results/compare_results.csv` file.

- **Experimental Procedure**: Step-by-step process available in the `Steps_for_solubility_experiment.md` file.

This comparison highlights the accuracy and reliability of our model.

## Repository Structure

- **data/**: Contains all datasets used in the study.
  - `water_solubility_data.csv`
  - `dataset-not-FA.csv`
  - `Supplementary_data.csv`
  - `data_paper.csv`
  - `dataset-E.csv`
  - `unique_train4_new24.csv`
  - `unique_test_new24.csv`
  - `overlap_data_new.csv`

- **notebooks/**: Jupyter notebooks containing the analysis, preprocessing, and model training code.
  - `data_preprocess_3_24.ipynb`
  - `4_24.ipynb`
  - `feature_fg7_fe_38.ipynb`
  - `Sorkundata_improve_Preprocess.ipynb`
  - `mpnn.ipynb`

- **scripts/**: Python scripts used in the study.
  - `utilities.py`

- **results/**: Contains the comparative model results.
  - `model_results.csv`
  - `compare_results.csv`

- **Sol_exp/**: Contains the Experimental solubility values for five specific compounds..
  - `EXP40_823.xlsx`
  - `EXP42_827.xlsx`
  - `EXP56_1562.xlsx`
  - `EXP76_1593.xlsx`
  - 'EXP8_260.xlsx'
  - `Steps_for_solubility_experiment.md`
## Conclusion

This study presents a thorough investigation into the prediction of water solubility using machine learning models. By carefully curating datasets, engineering features, and comparing with state-of-the-art methods, the study provides insights into the challenges and opportunities in this field. The findings highlight the importance of data preprocessing and feature selection in building robust predictive models.

For any questions or further information, please feel free to open an issue or contact me directly over mail id mushtaq.ali@kit.edu.
