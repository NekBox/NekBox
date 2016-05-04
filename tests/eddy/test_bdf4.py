import json
from nekpy.dask.subgraph import series
from nekpy.dask import run_all
from nekpy.dask.tasks import configure
from nekpy.dask.utils import work_name, outer_product
from nekpy.tools.log import grep_log
from os.path import join, dirname
from os import getcwd
import numpy as np
import pytest

with open(join(dirname(__file__), "eddy_uv_f90.tusr"), "r") as f:
    tusr = f.read()

base_dir = join(getcwd(), "scratch")

@pytest.mark.pr
def test_bdf4():
    with open(join(dirname(__file__), "eddy_uv.json")) as f:
        base = json.load(f)
    base["prefix"] = "test_bdf4"
    base["torder"] = 4
    sweep = {"courant" : [0.5, 0.25, 0.125]}
    overrides = list(outer_product(sweep))
    for ov in overrides:
        ov["name"] = work_name(base["prefix"], ov)

    workdir = join(base_dir, base["prefix"])
    configs = [configure(base, ov, join(base_dir, ov["name"])) for ov in overrides]
    res = [series(config, tusr) for config in configs]
    run_all(res, base)

    errs = {}
    for config in configs:
        with open(join(config["workdir"], "{}-0.stdout".format(config["name"]))) as f:
            test = f.readlines()    
        errs[config["courant"]] =  np.max(np.abs(grep_log(test, "X err", pos=2)[10:]))

    assert errs[.5] / errs[.25] > 12
    assert errs[.5] / errs[.25] < 24
    assert errs[.25] / errs[.125] > 12
    assert errs[.25] / errs[.125] < 24

    return
