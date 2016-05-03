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

@pytest.mark.pr
def test_bdf2():
    with open(join(dirname(__file__), "eddy_uv.json")) as f:
        base = json.load(f)
    base["prefix"] = "test_bdf2"
    base["torder"] = 2
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

    assert errs[.5] / errs[.25] > 3
    assert errs[.5] / errs[.25] < 6
    assert errs[.25] / errs[.125] > 3
    assert errs[.25] / errs[.125] < 6

    return

@pytest.mark.pr
def test_bdf3():
    with open(join(dirname(__file__), "eddy_uv.json")) as f:
        base = json.load(f)
    base["prefix"] = "test_bdf3"
    base["torder"] = 3
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

    assert errs[.5] / errs[.25] > 6
    assert errs[.5] / errs[.25] < 12
    assert errs[.25] / errs[.125] > 6
    assert errs[.25] / errs[.125] < 12

    return

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
