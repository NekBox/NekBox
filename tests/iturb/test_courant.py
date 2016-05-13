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

with open(join(dirname(__file__), "iturb_f90.tusr"), "r") as f:
    tusr = f.read()

def ndiff(test, ref):
    tests = grep_log(test, "Maximum scalar") 
    refs  = grep_log(ref, "Maximum scalar") 
    assert tests.shape[0] == refs.shape[0]
    max_scalar_diff = np.max(np.abs(tests-refs))
    assert max_scalar_diff < 1.e-9

    tests = grep_log(test, "Maximum velocity")
    refs  = grep_log(ref, "Maximum velocity")
    assert tests.shape[0] == refs.shape[0]
    max_velocity_diff = np.max(np.abs(tests-refs))
    assert max_velocity_diff < 1.e-9
    return

def old_test(name):
    with open(join(dirname(__file__), "{}.json".format(name))) as f:
        base = json.load(f)
    base["prefix"] = "test_{}".format(name)
    workdir = join(base_dir, base["prefix"])
    config = configure(base, {'name': base["prefix"]}, workdir)
    res = series(config, tusr)
    run_all([res,], base)[0]

    with open(join(workdir, "{}-0.stdout".format(base["prefix"]))) as f:
        test = f.readlines()    

    with open(join(dirname(__file__), "{}_ref.out".format(name))) as f:
        ref = f.readlines()

    ndiff(test, ref)

    return

def test_courant():
    old_test("iturb_courant")
    return
