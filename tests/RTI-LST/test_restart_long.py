import json
from nekpy.dask.subgraph import series
from nekpy.dask import run_all
from nekpy.dask.tasks import configure
from os.path import join, dirname
from os import getcwd
from nekpy.tools.log import grep_log
import numpy as np
import pytest

base_dir = join(getcwd(), "scratch")

with open(join(dirname(__file__), "LST_f90.tusr"), "r") as f:
    tusr = f.read()

def test_dummy():
    assert 1 = 1

@pytest.mark.pr
def test_restart_long():
    with open(join(dirname(__file__), "LST.json")) as f:
        base = json.load(f)
    base["prefix"] = "test_restart_long"

    workdir = join(base_dir, base["prefix"])
    config = configure(base, {'name': base["prefix"]}, workdir)
    res = series(config, tusr, job_time = 0.25)
    run_all([res,], base)[0]

    with open(join(workdir, "{}-3.stdout".format(config["name"]))) as f:
        test = f.readlines()    
    with open(join(dirname(__file__), "LST_ref.out")) as f:
        ref = f.readlines()

    tests = grep_log(test, "Maximum scalar")
    refs  = grep_log(ref,  "Maximum scalar")
    max_scalar_diff = np.max(np.abs(tests[-1]-refs[-1]))
    assert max_scalar_diff < 1.e-9

    tests = grep_log(test, "Maximum velocity")
    refs  = grep_log(ref, "Maximum velocity")
    max_velocity_diff = np.max(np.abs(tests[-1]-refs[-1]))
    assert max_velocity_diff < 1.e-9

    return
