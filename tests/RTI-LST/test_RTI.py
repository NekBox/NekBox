import json
from nekpy.dask.subgraph import series
from nekpy.dask import run_all
from nekpy.dask.tasks import configure
from os.path import join, dirname
from os import getcwd
import numpy as np

with open(join(dirname(__file__), "LST_f90.tusr"), "r") as f:
    tusr = f.read()

def scalars(lines):
    res = []
    for line in lines:
        if "Maximum scalar" in line:
            res.append(float(line.split()[2]))
    return np.array(res)

def velocities(lines):
    res = []
    for line in lines:
        if "Maximum velocity" in line:
            toks = line.split()
            res.append([float(toks[2]), float(toks[3]), float(toks[4])])
    return np.array(res)

def test_periodic():
    with open(join(dirname(__file__), "LST.json")) as f:
        base = json.load(f)
    workdir = join(getcwd(), "test_RTI")
    base["prefix"] = "test_RTI"
    config = configure(base, {'name': 'test_RTI'}, workdir)
    res = series(config, tusr)
    run_all([res,], base)[0]

    with open(join(workdir, "test_RTI-0.stdout")) as f:
        test = f.readlines()    

    with open(join(dirname(__file__), "LST_ref.out")) as f:
        ref = f.readlines()

    tests = scalars(test)
    refs  = scalars(ref)
    assert tests.shape[0] == refs.shape[0]
    max_scalar_diff = np.max(np.abs(tests-refs))
    assert max_scalar_diff < 1.e-9

    tests = velocities(test)
    refs  = velocities(ref)
    assert tests.shape[0] == refs.shape[0]
    max_velocity_diff = np.max(np.abs(tests-refs))
    assert max_velocity_diff < 1.e-9

    return
