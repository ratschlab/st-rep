from aestetik.utils.utils_transcriptomics import preprocess_adata
from aestetik.utils.utils_clustering import clustering
from sklearn.metrics import adjusted_rand_score
import matplotlib.pyplot as plt
from plotnine_prism import *
import plotnine as p9
import squidpy as sq
import scanpy as sc
import pandas as pd
import numpy as np
import glob
import yaml


def get_st_modality_plot(dataset_path, sample, model, dot_size=10, img_alpha=0.8,
                         ground_truth=False, out_path=False, ncols=4,
                         n_components=15, pad=0, legend_loc="right margin",
                         out_folder="out_benchmark", grount_truth_ordered=None):
    adata = sc.read(f"../{dataset_path}/{out_folder}/data/h5ad/{sample}.h5ad")
    adata = preprocess_adata(adata)
    adata.obsm["Transcriptomics"] = adata.obsm["X_pca"][:, 0:n_components]

    img_features_path = f"../{dataset_path}/{out_folder}/data/image_features/{sample}_inception.npy"
    img_features = np.load(img_features_path)
    adata.obsm["Morphology"] = img_features[:, 0:n_components]

    clustering(adata, adata.obs.ground_truth.unique().size, "Transcriptomics", "bgm")
    adata.obs["Transcriptomics"] = adata.obs["Transcriptomics_cluster"]
    clustering(adata, adata.obs.ground_truth.unique().size, "Morphology", "bgm")
    adata.obs["Morphology"] = adata.obs["Morphology_cluster"]

    adata.var_names_make_unique()
    best_split_sample = pd.read_csv(f"../{dataset_path}/{out_folder}/summary/summary_best_split_sample.csv")
    best_split_sample["sample"] = best_split_sample["sample"].astype(str)

    best_split_sample = best_split_sample[best_split_sample.model.isin([model])]
    best_split_sample = best_split_sample.query(f"sample == '{sample}'")
    best_split_sample = best_split_sample.sort_values("ari")
    cluster_paths = best_split_sample.path.values
    adata.obs["Barcode"] = adata.obs.index

    bounds = (adata.obsm["spatial"][:, 0].min() - pad,
              adata.obsm["spatial"][:, 1].min() - pad,
              adata.obsm["spatial"][:, 0].max() + pad,
              adata.obsm["spatial"][:, 1].max() + pad)

    for cluster_path in cluster_paths:
        df = pd.read_csv(f"../{dataset_path}/{cluster_path}")
        cluster_label_dict = pd.Series(df[df.columns[1]].values, index=df[df.columns[0]].values).to_dict()
        adata.obs[cluster_path.split("/")[2].replace("_evaluate",
                                                     "")] = adata.obs.Barcode.apply(lambda x: cluster_label_dict[x]).astype(str)
    models = ["Transcriptomics", "Morphology", model]

    if grount_truth_ordered is None:
        grount_truth_ordered = list(adata.obs.ground_truth.unique())
    elif isinstance(grount_truth_ordered, dict):
        adata.obs["ground_truth"] = adata.obs["ground_truth"].apply(lambda x: grount_truth_ordered[x])
        grount_truth_ordered = list(grount_truth_ordered.values())

    adata.obs["ground_truth"] = pd.Categorical(adata.obs["ground_truth"], grount_truth_ordered)
    sq.pl._color_utils._maybe_set_colors(adata, adata, "ground_truth")
    for model in models:

        result_df = adata.obs[["ground_truth", model]].groupby(
            "ground_truth")[model].apply(lambda x: x.mode().iloc[0]).reset_index()
        result_df.set_index("ground_truth", inplace=True)
        result_df = result_df.loc[grount_truth_ordered]
        result_df = result_df[~result_df[model].duplicated()]

        model_ordered = list(result_df[model].values)

        for label in adata.obs[model].unique():
            if label not in model_ordered:
                model_ordered.append(label)

        model_ordered_dict = dict(zip(model_ordered, range(1, len(model_ordered) + 1)))

        adata.obs[model] = adata.obs[model].apply(lambda x: model_ordered_dict[x])
        adata.obs[model] = pd.Categorical(adata.obs[model], model_ordered_dict.values())
        adata.uns[f"{model}_colors"] = adata.uns[f"ground_truth_colors"][:len(model_ordered)]

    if ground_truth:
        color = ["ground_truth", *models]
        title = ["Manual annotation", *
                 [f"{model} ARI: {adjusted_rand_score(adata.obs['ground_truth'].values, adata.obs[model].values):.2f}" for model in models]]
    else:
        color = models
        title = [
            *[f"{model} ARI: {adjusted_rand_score(adata.obs['ground_truth'].values, adata.obs[model].values):.2f}" for model in models]]
    if out_path:
        out_path = f"{dataset_path}_{sample}_{'_'.join(models)}.png"
    return sq.pl.spatial_scatter(adata, img_alpha=img_alpha, crop_coord=bounds, wspace=0, color=color, size=dot_size,
                                 ncols=ncols, title=title, save=out_path, dpi=300,
                                 frameon=False, legend_loc=legend_loc, figsize=(4, 8))


def get_st_plot(dataset_path, sample, models, dot_size=10, ground_truth=False,
                out_path=False, wspace=0, ncols=2, pad=0, legend_loc="right margin",
                figsize=(4, 5), img_alpha=0.8, alpha=1, ground_truth_name="Manual annotation", grount_truth_ordered=None):
    adata = sc.read(f"../{dataset_path}/out_benchmark/data/h5ad/{sample}.h5ad")
    adata.var_names_make_unique()
    best_split_sample = pd.read_csv(f"../{dataset_path}/out_benchmark/summary/summary_best_split_sample.csv")
    best_split_sample["sample"] = best_split_sample["sample"].astype(str)
    best_split_sample = best_split_sample[best_split_sample.model.isin(models)]
    best_split_sample = best_split_sample.query(f"sample == '{sample}'")
    best_split_sample = best_split_sample.sort_values("ari")
    cluster_paths = best_split_sample.path.values
    adata.obs["Barcode"] = adata.obs.index

    bounds = (adata.obsm["spatial"][:, 0].min() - pad,
              adata.obsm["spatial"][:, 1].min() - pad,
              adata.obsm["spatial"][:, 0].max() + pad,
              adata.obsm["spatial"][:, 1].max() + pad)

    for cluster_path in cluster_paths:
        df = pd.read_csv(f"../{dataset_path}/{cluster_path}")
        if min(df[df.columns[1]]) == 0:
            df[df.columns[1]] = df[df.columns[1]] + 1
        cluster_label_dict = pd.Series(df[df.columns[1]].values, index=df[df.columns[0]].values).to_dict()
        adata.obs[cluster_path.split("/")[2].replace("_evaluate",
                                                     "")] = adata.obs.Barcode.apply(lambda x: cluster_label_dict[x]).astype(str)

    if grount_truth_ordered is None:
        grount_truth_ordered = list(adata.obs.ground_truth.unique())
    elif isinstance(grount_truth_ordered, dict):
        adata.obs["ground_truth"] = adata.obs["ground_truth"].apply(lambda x: grount_truth_ordered[x])
        grount_truth_ordered = list(grount_truth_ordered.values())

    adata.obs["ground_truth"] = pd.Categorical(adata.obs["ground_truth"], grount_truth_ordered)
    sq.pl._color_utils._maybe_set_colors(adata, adata, "ground_truth")
    for model in best_split_sample.model:

        result_df = adata.obs[["ground_truth", model]].groupby(
            "ground_truth")[model].apply(lambda x: x.mode().iloc[0]).reset_index()
        result_df.set_index("ground_truth", inplace=True)
        result_df = result_df.loc[grount_truth_ordered]
        result_df = result_df[~result_df[model].duplicated()]

        model_ordered = list(result_df[model].values)

        for label in adata.obs[model].unique():
            if label not in model_ordered:
                model_ordered.append(label)

        model_ordered_dict = dict(zip(model_ordered, range(1, len(model_ordered) + 1)))

        adata.obs[model] = adata.obs[model].apply(lambda x: model_ordered_dict[x])
        adata.obs[model] = pd.Categorical(adata.obs[model], model_ordered_dict.values())
        adata.uns[f"{model}_colors"] = adata.uns[f"ground_truth_colors"][:len(model_ordered)]

    if ground_truth:
        color = ["ground_truth", *best_split_sample.model]
        title = [ground_truth_name, *
                 [f"{model} ARI: {adjusted_rand_score(adata.obs['ground_truth'].values, adata.obs[model].values):.2f}" for model in best_split_sample.model]]
    else:
        color = best_split_sample.model
        title = [
            *[f"{model} ARI: {adjusted_rand_score(adata.obs['ground_truth'].values, adata.obs[model].values):.2f}" for model in best_split_sample.model]]
    if out_path and not isinstance(out_path, str):
        out_path = f"{dataset_path}_{sample}_{'_'.join(models)}.png"

    return sq.pl.spatial_scatter(adata, img_alpha=img_alpha, crop_coord=bounds, wspace=wspace, color=color,
                                 size=dot_size, ncols=ncols, title=title, save=out_path, alpha=alpha,
                                 dpi=300, frameon=False, legend_loc=legend_loc, figsize=figsize)
