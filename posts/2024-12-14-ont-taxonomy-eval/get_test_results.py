#! /usr/bin/env python3

import subprocess
import os

mgs_results_dir = "/Users/simongrimm/code/simons-notebook/posts/2024-12-14-ont-taxonomy-eval/mgs-results"


run_types = ["test-ont-ww"]


def get_test_results(mgs_results_dir):
    for run_type in run_types:
        s3_results = f"s3://nao-mgs-simon/{run_type}/output/"
        results_dir = os.path.join(mgs_results_dir, run_type, "output")
        os.makedirs(results_dir, exist_ok=True)
        subprocess.run(
            [
                "aws",
                "s3",
                "sync",
                s3_results,
                results_dir,
            ]
        )


# def get_test_work_dirs(mgs_results_dir):
#     for run_type in run_types:
#         s3_work = f"s3://nao-mgs-simon/{run_type}/work/"
#         results_dir = os.path.join(mgs_results_dir, run_type, "work")
#         os.makedirs(results_dir, exist_ok=True)
#         subprocess.run(["aws", "s3", "sync", s3_work, results_dir])


if __name__ == "__main__":
    get_test_results(mgs_results_dir)
    # get_test_work_dirs(mgs_results_dir)
