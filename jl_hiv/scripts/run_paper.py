#!/usr/bin/env python
"""Reproducible script to download data and run analysis from Li et al paper (PLoS One in review).
"""
import os
import subprocess

def main():
    make_directory_structure()
    config_file = get_data()
    cmds = install_code()
    run_analysis(cmds, config_file)

def run_analysis(cmds, config_file):
    subprocess.check_call([cmds["python"], cmds["variant_script"], config_file])

def get_data():
    """Retrieve data and configuration files for running the analysis.
    """
    pname = "A5262"
    for sample in ["S1", "S2", "S3", "S4", "S5", "Control"]:
        for read in ["1", "2"]:
            fname = os.path.join("input", "%s-%s-%s.fq.gz" % (pname, sample, read))
            if not os.path.exists(fname):
                raise NotImplementedError("Get from SRA when available: %s" % fname)
    config_file = os.path.join(os.getcwd(), "config", "20120111.yaml")
    if not os.path.exists(config_file):
        subprocess.check_call(["wget", "-O", config_file,
                               "https://raw.github.com/hbc/projects/master/jl_hiv/config/20120111.yaml"])
    if not os.path.exists("refinfo"):
        subprocess.check_call(["wget", "https://s3.amazonaws.com/hbc_data/li_refinfo.tar.gz"])
        subprocess.check_call(["tar", "-xzvpf", "li_refinfo.tar.gz"])
        os.remove("li_refinfo.tar.gz")
    return config_file

def install_code():
    """Install the required code for running the analysis.

    Requires a basic python installation to bootstrap installation.
    """
    python_cmd = os.path.join(os.getcwd(), "local", "anaconda", "bin", "python")
    if not os.path.exists(python_cmd):
        subprocess.check_call(["wget", "https://raw.github.com/chapmanb/bcbio-nextgen/master/"
                               "scripts/bcbio_nextgen_install.py"])
        subprocess.check_call(["python", "bcbio_nextgen_install.py", "local", "--nodata", "--nosudo", "--isolate"])
        os.remove("bcbio_nextgen_install.py")
    orig_dir = os.getcwd()
    os.chdir("build")
    dirname = "projects-master/jl_hiv"
    if not os.path.exists(dirname):
        subprocess.check_call(["wget", "https://github.com/hbc/projects/archive/master.zip"])
        subprocess.check_call(["unzip", "master.zip"])
        os.chdir(dirname)
        subprocess.check_call([python_cmd, "setup.py", "install"])
    os.chdir(orig_dir)
    return {"python": python_cmd,
            "variant_script": os.path.join(os.getcwd(), "build", dirname, "scripts", "variant_identify.py")}

def make_directory_structure():
    for dname in ["input", "config", "local", "build"]:
        if not os.path.exists(dname):
            os.makedirs(dname)

if __name__ == "__main__":
    main()
