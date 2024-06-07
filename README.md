# st-rep

This repository contains the code of the paper "Representation learning for multi-modal spatially resolved transcriptomics data".

**Authors**: Kalin Nonchev, Sonali Andani, Joanna Ficek-Pascual, Marta Nowak, Bettina Sobottka, Tumor Profiler Consortium, Viktor Hendrik Koelzer, and Gunnar Rätsch

The preprint is available [here](https://www.medrxiv.org/content/10.1101/2024.06.04.24308256v1).

You can find the AESTETIK code and tutorial on how to use it [here](https://github.com/ratschlab/aestetik).

## Snakemake overview

We provide the Snakemake pipeline we used to generate our paper's results. In the [workflows folder](workflows), we share the code base used for producing the cluster assignments on each dataset. 
Briefly:
1) The dataset-specific preprocessing can be found [here](workflows/preprocess).
2) The model scripts can be found [here](workflows/models).
3) The evaluation scripts can be found [here](workflows/analysis).


## How to start

### 1. Install conda environments

We start by installing the [conda environments](environments/) required for the different rules.

```
conda env create -f=path/to/yaml
```

### 2. Config files

In the [Snakemake_info.yaml](Snakemake_info.yaml) we specify the general rule requirements and resources. Please adjust based on your setup.

In each dataset folder, there should be an `info.yaml` file (e.g., [`maynard_human_brain_analysis/info.yaml`](maynard_human_brain_analysis/info.yaml)), where we specify the reversed leave-one-out cross-validation folds along with dataset-specific information and the models to use. This file is used as input for the Snakemake pipeline. 

### 3. Data structuring and preprocessing

The datasets can be downloaded from:
1) LIBD Human DLPFC dataset is available at https://github.com/LieberInstitute/HumanPilot and http://research.libd.org/spatialLIBD;
2) Human Breast Cancer - Zenodo https://doi.org/10.5281/zenodo.4739739,
3) Human Liver Normal and Cancer - https://nanostring.com/products/cosmx-spatial-molecular-imager/human-liver-rna-ffpe-dataset/.
4) Metastatic melanoma dataset - tba.

The downloaded raw spatial transcriptomics data has a different structure so it has to be unified. For the discussed datasets, we provide examples [here](workflows/preprocess). For the Human Liver Normal and Cancer datasets, we start by [grouping and preprocessing the individual FOVs](prepareData_cosMx_human_liver_fov.ipynb), and then create the RGB images (e.g., [cosMx_human_liver_normal/create_rgb_image.ipynb](cosMx_human_liver_normal/create_rgb_image.ipynb).

### Execute preprocessing pipeline

Navigate to one of the dataset folders:
  - 10x_TuPro_v2
  - cosMx_human_liver_cancer
  - cosMx_human_liver_normal
  - maynard_human_brain_analysis
  - human_breast_cancers

#### Variant 1: You can unify the raw spatial transcriptomics datasets yourself

```
conda activate st_rep_python_snakemake
snakemake -s ../Snakefile.preprocess -k --use-conda --rerun-incomplete --rerun-triggers mtime --cluster "sbatch --mem={resources.mem_mb} --cpus-per-task={threads} -t {resources.time} --gres={resources.gpu} -p {resources.p} -o {resources.log} -J {resources.jobname} --tmp {resources.tmp}" -j 10
```

#### Variant 2: You can download the already unified transcriptomics datasets from [here](https://zenodo.org/records/10658804)
It has the following structure: {dataset}/{out_folder}/{data}/:
  - h5ad (transcriptomic h5ad)
  - rds (transcriptomic rds)
  - image (morphology)
  - meta (spot size, etc)

```
wget -c https://zenodo.org/records/10658804/files/st_data.tar.gz?download=1 # 21.2G 
tar -xvzf st_data.tar.gz
```

### 4. Execute evaluation pipeline

Within the folder execute the following command in the terminal:

```python
conda activate st_rep_python_snakemake
snakemake -s ../Snakefile.evaluate -k --use-conda --rerun-incomplete --rerun-triggers mtime --cluster "sbatch --mem={resources.mem_mb} --cpus-per-task={threads} -t {resources.time} --gres={resources.gpu} -p {resources.p} -o {resources.log} -J {resources.jobname} --tmp {resources.tmp}" -j 10
```

#### NB: Computational data analysis was performed at Leonhard Med (https://sis.id.ethz.ch/services/sensitiveresearchdata/) secure trusted research environment at ETH Zurich. Our pipeline aligns with the specific cluster requirements and resources.

## Ablation study

For the ablation study, we specify the fixed hyperparameters for each model [here](workflows/configs) and then we create `info.yaml` file with the model and the fixed hyperparameter value (e.g., [`maynard_human_brain_analysis/info_ablation.yaml`](maynard_human_brain_analysis/info_ablation.yaml)).

## Simulated data

We provide the code for simulating spatial transcriptomics data as a [separate module linked to this repo](https://github.com/ratschlab/simulate_spatial_transcriptomics_tool). In summary:

> We adapted the simulation approach suggested in [5] by introducing spatial structure in the experiment. Briefly, relying on simulated ground truth labels, we simulate transcriptomics and morphology modalities, allowing partial observation of true clusters within each modality individually. However, combining both modalities enables the identification of all clusters. Spatial coordinates are incorporated by sorting the ground truth in spatial space.

First, we have to generate the simulated data by running the [notebook](simulate_data.ipynb). It creates 3 datasets with 2500 cells each with 5, 10, and 15 clusters together with the corresponding `info.yaml` files.
Please note that for the simulated data we use directly [Snakefile.evaluate](Snakefile.evaluate) file, because we start from the already structured files.

Navigate to one of the dataset folders:
  - simulated_data_5_clusters
  - simulated_data_10_clusters
  - simulated_data_15_clusters

Within the folder execute the following command in the terminal:

```
conda activate st_rep_python_snakemake
snakemake -s ../Snakefile.evaluate -k --use-conda --rerun-incomplete --rerun-triggers mtime --cluster "sbatch --mem={resources.mem_mb} --cpus-per-task={threads} -t {resources.time} --gres={resources.gpu} -p {resources.p} -o {resources.log} -J {resources.jobname} --tmp {resources.tmp}" -j 10
```

## Test run

We also provide an example dataset to quickly test your pipeline setup before running the full pipeline on the real data

```
cd test_data
conda activate st_rep_python_snakemake
snakemake -s ../Snakefile.evaluate -k --use-conda --rerun-incomplete --rerun-triggers mtime --cluster "sbatch --mem={resources.mem_mb} --cpus-per-task={threads} -t {resources.time} --gres={resources.gpu} -p {resources.p} -o {resources.log} -J {resources.jobname} --tmp {resources.tmp}" -j 10
```

## Citation

In case you found our work useful, please consider citing us:

```
@article{nonchev2024representation,
  title={Representation learning for multi-modal spatially resolved transcriptomics data},
  author={Nonchev, Kalin and Andani, Sonali and Ficek-Pascual, Joanna and Nowak, Marta and Sobottka, Bettina and Tumor Profiler Consortium and Koelzer, Viktor Hendrik and Raetsch, Gunnar},
  journal={medRxiv},
  pages={2024--06},
  year={2024},
  publisher={Cold Spring Harbor Laboratory Press}
}
```

The code for reproducing the paper results can be found [here](https://github.com/ratschlab/st-rep).

## Contact

In case, you have questions, please get in touch with [Kalin Nonchev](https://bmi.inf.ethz.ch/people/person/kalin-nonchev).
