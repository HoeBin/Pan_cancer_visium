"""STopover ligand-receptor co-localization runner.

Consolidates STo_global_All.py and STo_global_3region.py (STopver_run folder) into
one parameterized script. Runs STopover's LR topological-similarity (Jaccard)
analysis over every pre-processed Visium slide under `output_root`, either on the
whole slide ("global" scope) or separately within each Malignant/Boundary/Normal
compartment ("3region" scope).

All STopover call parameters match exactly what was actually run (transcribed from
the non-commented code in both scripts, which were identical).

Usage:
    python run_stopover.py --scope global
    python run_stopover.py --scope 3region
    python run_stopover.py --scope global --start-index 30   # resume a partial run
"""

import argparse
import gc
import os
import time

from STopover import STopover_visium

# "Paths per scope" -----------------------------------------------------------
# The paths below are example directory names. Edit them for your environment
# before running.
SCOPE_CONFIG = {
    "global": dict(
        output_root="/path/to/input/stopover_slides/global_all",
        new_output_root="/path/to/output/stopover_lr/global_all",
        compartments=None,
    ),
    "3region": dict(
        output_root="/path/to/input/stopover_slides/global_3region",
        new_output_root="/path/to/output/stopover_lr/global_3region",
        compartments=["Mal", "Bdy", "Normal"],
    ),
}

# "STopover parameters" -------------------------------------------------------
TOPOLOGICAL_SIMILARITY_PARAMS = dict(
    use_lr_db=True,
    lr_db_species="human",
    db_name="CellTalk",
    jaccard_type="default",
    J_result_name="result",
    num_workers=10,
)
# recorded as metadata columns on the Jaccard result table (matches original code;
# these are not passed into topological_similarity(), which uses its own defaults)
JACCARD_RESULT_METADATA = dict(min_size=20, fwhm=2.5, percent=30)


def run_slide(sp_load_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    spatial_model = STopover_visium(sp_load_path=sp_load_path, lognorm=False, save_path=out_dir)
    spatial_model.topological_similarity(**TOPOLOGICAL_SIMILARITY_PARAMS)

    df_mod = spatial_model.uns["J_result_1"]
    for key, value in JACCARD_RESULT_METADATA.items():
        df_mod[key] = value
    df_mod.to_csv(os.path.join(out_dir, "jaccard_composite_lr.csv"), sep=",", header=True, index=False)

    spatial_model.save_connected_loc_data(save_format="h5ad", filename="lr_pair_cc")

    del spatial_model, df_mod
    gc.collect()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", choices=sorted(SCOPE_CONFIG), required=True)
    parser.add_argument(
        "--start-index", type=int, default=0,
        help="Slide index to resume from (folder_list order is os.listdir order)",
    )
    args = parser.parse_args()

    config = SCOPE_CONFIG[args.scope]
    output_root = config["output_root"]
    new_output_root = config["new_output_root"]
    compartments = config["compartments"]
    os.makedirs(new_output_root, exist_ok=True)

    folder_list = os.listdir(output_root)
    print(f"## Scope: {args.scope}, {len(folder_list)} slides")

    start_total = time.time()
    for i in range(args.start_index, len(folder_list)):
        loop_start = time.time()
        tmp_folder = folder_list[i]
        print(f"[{i + 1}/{len(folder_list)}] Processing folder: {tmp_folder}")

        out_dir = os.path.join(new_output_root, tmp_folder)
        os.makedirs(out_dir, exist_ok=True)

        if compartments is None:
            sp_load_path = os.path.join(output_root, tmp_folder, "sp_visium_celltype_adata.h5ad")
            run_slide(sp_load_path, out_dir)
        else:
            for compartment in compartments:
                sp_load_path = os.path.join(
                    output_root, tmp_folder, compartment, "sp_visium_celltype_adata.h5ad"
                )
                out_dir2 = os.path.join(out_dir, compartment)
                run_slide(sp_load_path, out_dir2)

        print(f"--> Done [{i + 1}/{len(folder_list)}], time taken: {time.time() - loop_start:.2f} seconds\n")

    print(f"All done! Total time: {time.time() - start_total:.2f} seconds")


if __name__ == "__main__":
    main()
