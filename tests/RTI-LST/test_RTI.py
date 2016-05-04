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

@pytest.mark.pr
def test_periodic():
    old_test("LST")
    return

@pytest.mark.pr
def test_fft():
    old_test("LST_fft")
    return

@pytest.mark.pr
def test_asym():
    old_test("LST_asym")
    return

def test_nonuniform():
    old_test("LST_nonuniform")
    return

@pytest.mark.pr
def test_all():
    old_test("LST_all")
    return

@pytest.mark.pr
def test_sym():
    old_test("LST_sym")
    return

def test_sym_half():
    old_test("LST_sym_half")
    return

def test_restart():
    with open(join(dirname(__file__), "LST.json")) as f:
        base = json.load(f)
    base["prefix"] = "test_restart"

    workdir = join(base_dir, base["prefix"])
    config = configure(base, {'name': base["prefix"]}, workdir)
    res = series(config, tusr, job_time = 0.5)
    run_all([res,], base)[0]

    with open(join(workdir, "{}-1.stdout".format(config["name"]))) as f:
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
