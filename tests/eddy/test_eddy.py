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

def test_single():
    with open(join(dirname(__file__), "eddy_uv.json")) as f:
        base = json.load(f)
    base["prefix"] = "test_single"

    workdir = join(base_dir, base["prefix"])
    config = configure(base, {'name': base["prefix"]}, workdir)
    res = series(config, tusr)
    run_all([res,], base)[0]

    with open(join(workdir, "{}-0.stdout".format(config["name"]))) as f:
        test = f.readlines()    

    errs = grep_log(test, "X err", pos=2)[10:]
    assert np.max(np.abs(errs)) < 2.0e-08

    return
