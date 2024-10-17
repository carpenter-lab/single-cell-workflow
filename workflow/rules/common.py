import textwrap
from collections.abc import MutableMapping
from itertools import chain, product
from typing import Any, Callable, Generator, Optional, Self

from snakemake import exceptions as smk


def valid_dict_key(d: dict, key: str) -> bool:
    try:
        d[key]
    except AttributeError:
        return False
    except KeyError:
        return False

    if d[key] == "None":
        return False

    return True


class DefaultLabels(dict):
    title: Optional[str]
    grouping: Optional[str]
    split: Optional[str]
    assay: Optional[str]
    reduction: Optional[str]
    type: Optional[str]

    @classmethod
    def create(
        cls,
        title="No Title",
        grouping="No Grouping",
        split="No Split",
        assay="Undefined",
        reduction="Undefined",
        type="Undefined"
    ) -> Self:
        return DefaultLabels(
            title=title,
            grouping=grouping,
            split=split,
            assay=assay,
            reduction=reduction,
            type=type
        )


def report_plot_labels(wildcards: dict, params: dict) -> DefaultLabels:
    labels = DefaultLabels()

    if valid_dict_key(params, "plot"):
        labels["type"] = PlotTitle.human_readable_plot_type(params.get("plot"))
        if valid_dict_key(wildcards, "type"):
            labels["type"] = wildcards["type"] + "-" + labels["type"]

    if valid_dict_key(wildcards, "grouping"):
        labels["grouping"] = wildcards.get("grouping")

    if valid_dict_key(wildcards, "group_by"):
        labels["grouping"] = wildcards.get("group_by")
        
        if valid_dict_key(wildcards, "split_by"):
            labels["split"] = wildcards.get("split_by")
            
            if labels["split"] == labels["grouping"]:
                labels["split"] = "grouping"
        
    
    if valid_dict_key(wildcards, "assay"):
        labels["assay"] = wildcards.get("assay")
    
    if valid_dict_key(params, "reduction"):
        labels["reduction"] = params.get("reduction")
        
        if valid_dict_key(labels, "assay"):
            if labels.get("assay") in labels.get("reduction"):
                labels["reduction"].replace(f"{labels.get('assay')}_", "")
        
        if valid_dict_key(labels, "grouping"):
            if "clusters" in labels.get("grouping"):
                if "label" in labels.get("grouping"):
                    labels["grouping"] = "cluster labels"
                else:
                    labels["grouping"] = "clusters"

    for k, v in labels.items():
        labels[k] = v.replace("_", " ")

    return labels


def get_plot_type(wildcards: dict) -> str:
    try:
        wildcards["split_by"]
    except KeyError:
        return "cluster"

    if (
        "predicted.ann" in wildcards["split_by"]
        and "predicted.ann" in wildcards["group_by"]
    ):
        return "cell_type"
    elif wildcards["split_by"] is None:
        return "cluster"
    else:
        return "cluster"


class PlotTitle:
    def __init__(self, plot_type):
        self.plot_type = plot_type
        
    @staticmethod
    def human_readable_plot_type(plot_type):
        match plot_type:
            case "dot":
                return "Dot Plot"
            case "cluster" | "split" | "cell_type":
                return "DimPlot"
            case "heatmap":
                return "Heatmap"
            case "de_plot":
                return "Group wise DE Plot"
            case "cluster_cor":
                return "Cluster by Cluster Correlation Matrix"
            case "bar":
                return "Bar Plot"
            case "qc":
                return "Quality Control Metrics"
            case _:
                return "No Title Can Be Created"
                

    def _plot_title_type(self, wildcards) -> str:
        if isinstance(self.plot_type, Callable):
            self.plot_type = self.plot_type(wildcards)

        return self.human_readable_plot_type(self.plot_type)

    def make_title(self, wildcards):
        assay = wildcards["assay"]
        if wildcards["assay"] == "SCT":
            assay = "RNA with SCTransform"
        elif wildcards["assay"] == "RNA":
            assay = "Log Normalized RNA"

        return f"{wildcards.subset} | {self._plot_title_type(wildcards)} on {assay}"


def make_plot_subtitle(wildcards: dict) -> str:
    subtitle = f"Groups: {wildcards['group_by']}"
    if valid_dict_key(wildcards, "split_by"):
        subtitle = subtitle + f" | Split: {wildcards['split_by']}"

    return subtitle


def validate_de_method(method: str) -> str:
    valid_methods = [
        "wilcox",
        "wilcox_limma",
        "bimod",
        "roc",
        "t",
        "negbinom",
        "poisson",
        "LR",
        "MAST",
        "DESeq2",
    ]

    if method not in valid_methods:
        raise smk.WorkflowError(
            f"DE_METHOD must be one of: {', '.join(valid_methods)} not {method}. "
            "Please note, these are case sensitive."
        )

    return method


def get_dot_plot_features(config) -> Callable[[dict], dict]:
    def _get_dot_plot_features(wildcards: dict) -> dict:
        features = config["plotting"]["dot_plot"]["features"].get(wildcards["assay"])
        return features

    return _get_dot_plot_features


def get_category_name(wildcards: dict) -> str:
    subset: str = wildcards["subset"]
    return subset.replace("_", " ")


def get_proper_clustering_output(config, rules) -> Callable:
    def _get_proper_clustering_output(wildcards):
        if wildcards["subset"] == config["cluster"].get("all_data_key"):
            return "results/gliph/labelled_tcr.rds"
        else:
            return "results/clustering/{subset}.rds"

    return _get_proper_clustering_output


def _flatten_dict_gen(d: MutableMapping) -> Generator:
    for key, value in d.items():
        if isinstance(value, MutableMapping):
            yield from flatten_dict(value).items()
        else:
            yield key, value


def flatten_dict(d: MutableMapping) -> dict:
    return dict(_flatten_dict_gen(d))


def get_subcluster_params(config: dict) -> Optional[list[dict]]:
    return [flatten_dict(d) for d in config["cluster"]["subclusters"]]


def unzip_dict(d: dict) -> Generator[dict, None, None]:
    keys = [k for k, v in d.items() if isinstance(v, list)]
    values = [d[k] for k in keys]

    for values in product(*values):
        yield {**d, **dict(zip(keys, values))}


def get_umap_plot_output_files(config: dict) -> list[str]:
    umap = config["plotting"]["umap_plot"]
    subsets = umap["subsets"]
    assays = umap["assays"]
    reductions = umap["reductions"]
    groups = umap["group_by"] + umap["group_by_yte_only"]
    params = list(chain.from_iterable([unzip_dict(flatten_dict(d)) for d in groups]))

    path_list = list()

    for param in params:
        group, split = param.values()
        for subset, assay, reduction in product(subsets, assays, reductions):
            path_list.append(
                f"results/clustering/{subset}/plots/{assay}/{reduction}/{group}_split_{split}.png".format(
                    subset=subset,
                    assay=assay,
                    reduction=reduction,
                    group=group,
                    split=split,
                )
            )

    return path_list


class WorkflowResults:
    def __init__(self, config, file_path, ext=None):
        self.config = config
        self.file_path = file_path

        self.do_not_run = config is None

        if not self.do_not_run:
            self.group_by = config.get("group_by", [])
            self.group_by_yte_only = config.get("group_by_yte_only", [])
            self.assays = config.get("assays", [])
            self.reductions = config.get("reductions", [])
            self.subsets = config.get("subsets", [])
            self._extensions = ext

    def __str__(self):
        string_rep = f"""
        A WorkflowResults object on {len(self.subsets)} subset{"" if len(self.subsets) == 1 else "s"} with:
        {len(self.group_by) + len(self.group_by_yte_only)} group by params
        {len(self.assays)} assay{"" if len(self.assays) == 1 else "s"}
        {len(self.reductions)} reduction{"" if len(self.reductions) == 1 else "s"}
        And is {"not " if self.do_not_run else ""}set to run.
        """
        return textwrap.dedent(string_rep)

    @property
    def extensions(self):
        return self._extensions

    @extensions.setter
    def extensions(self, value):
        self._extensions = value

    def combine_group_by_keys(self):
        return self.group_by + self.group_by_yte_only

    @staticmethod
    def unlist_dict_keys(input_list: list) -> list:
        return list(
            chain.from_iterable([unzip_dict(flatten_dict(d)) for d in input_list])
        )

    @staticmethod
    def _curly_braces_to_parens_string(string: str) -> str:
        return string.replace("{", "%(").replace("}", ")s")

    @classmethod
    def _create_path_list(
        cls,
        path_list,
        subsets,
        assays,
        file_path,
        group,
        reductions: Optional[list],
        split=None,
        extensions=None,
    ):
        if reductions:
            for subset, assay, reduction in product(subsets, assays, reductions):
                init_path = file_path.format(
                    **{
                        "subset": subset,
                        "assay": assay,
                        "reduction": reduction,
                        "group": group,
                        "split": split,
                    }
                )
                avail_fields = string.Formatter().parse(init_path)
                avail_fields = [ x for _, x, _, _ in avail_fields]
                if any([f == "reduction" for f in avail_fields]):
                    updated_path = init_path.format(
                        **{
                            "reduction": reduction,
                        }
                    )
                else:
                    updated_path = init_path
                path_list.append(updated_path)
        else:
            for subset, assay in product(subsets, assays):
                if extensions:
                    for extension in extensions:
                        init_path = file_path.format(
                            **{
                                "subset": subset,
                                "assay": assay,
                                "group": group,
                                "split": split,
                                "ext": extension,
                            }
                        )
                        path_list.append(init_path)
                else:
                    init_path = file_path.format(
                        **{
                            "subset": subset,
                            "assay": assay,
                            "group": group,
                            "split": split,
                        }
                    )
                    path_list.append(init_path)

    def create_path_list(self):
        if self.do_not_run:
            return []

        params = self.combine_group_by_keys()
        path_list = list()

        def _safe_list_keys(obj: Any) -> list[str]:
            try:
                return obj.keys()
            except AttributeError:
                return [""]

        if any(["split" in _safe_list_keys(d) for d in params]):
            params = self.unlist_dict_keys(self.combine_group_by_keys())
            for param in params:
                group, split = param.values()
                self._create_path_list(
                    path_list,
                    subsets=self.subsets,
                    file_path=self.file_path,
                    assays=self.assays,
                    reductions=self.reductions,
                    group=group,
                    split=split,
                )
        else:
            for group in params:
                self._create_path_list(
                    path_list,
                    file_path=self.file_path,
                    subsets=self.subsets,
                    assays=self.assays,
                    reductions=self.reductions,
                    group=group,
                    extensions=self._extensions,
                )

        return path_list
