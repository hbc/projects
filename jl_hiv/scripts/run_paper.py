#!/usr/bin/env python
"""Reproducible script to download data and run analysis from Li et al paper (PLoS One in review).

Bootstraps local third party tools, processing code and data to fully reproduce
work in the paper. Requires a minimal unix system with:

- Python 2.7
- Java
- Ruby
- Git
"""
import os
import subprocess
import sys

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
        fname = os.path.join("input", "%s-%s.fq.gz" % (pname, sample))
        if not os.path.exists(fname):
            raise NotImplementedError("Get from SRA when available: %s" % fname)
    config_file = os.path.join(os.getcwd(), "config", "20120111-paper.yaml")
    if not os.path.exists(config_file):
        subprocess.check_call(["wget", "-O", config_file,
                               "https://raw.github.com/hbc/projects/master/jl_hiv/config/20120111-paper.yaml"])
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
    install_dir = os.path.join(orig_dir, "local")
    install_thirdparty(install_dir,
                       os.path.join(orig_dir, "build", dirname, "cbl_install"))
    install_supporting(python_cmd, install_dir)
    os.chdir(orig_dir)
    return {"python": python_cmd,
            "variant_script": os.path.join(os.getcwd(), "build", dirname, "scripts", "variant_identify.py")}

def install_thirdparty(install_dir, flavor):
    """Install third party software needed by the pipeline.
    """
    cbl_dir = os.path.join(os.getcwd(), "cloudbiolinux")
    if not os.path.exists(cbl_dir):
        subprocess.check_call(["git", "clone", "https://github.com/chapmanb/cloudbiolinux.git"])
    s = {"fabricrc_overrides": {"system_install": install_dir,
                                "local_install": install_dir,
                                "use_sudo": False, "edition":
                                "minimal"},
         "flavor": flavor,
         "vm_provider": "novm",
         "hostname": "localhost",
         "actions": ["install_biolinux"]}
    if not os.path.exists(os.path.join(install_dir, "bin", "novoalign")):
        sys.path.insert(0, cbl_dir)
        cbl_deploy = __import__("cloudbio.deploy", fromlist=["deploy"])
        cbl_deploy.deploy(s)

def install_supporting(python_cmd, install_dir):
    """Install supporting scripts: bloom uniquification from @brentp
    """
    support_script = os.path.join(install_dir, "bin", "fastq-unique-bloom.py")
    if not os.path.exists(support_script):
        pip_cmd = os.path.join(os.path.dirname(python_cmd), "pip")
        subprocess.check_call([pip_cmd, "install", "git+https://github.com/brentp/pybloomfaster.git"])
        subprocess.check_call(["wget",
                               "https://raw.github.com/brentp/pybloomfaster/master/examples/fastq-unique.py"])
        with open(support_script, "w") as out_handle:
            with open("fastq-unique.py") as in_handle:
                out_handle.write("#!%s\n" % python_cmd)
                for line in in_handle:
                    out_handle.write(line)
        subprocess.check_call(["chmod", "a+x", support_script])
        os.remove("fastq-unique.py")

def make_directory_structure():
    for dname in ["input", "config", "local", "local/bin", "build"]:
        if not os.path.exists(dname):
            os.makedirs(dname)

if __name__ == "__main__":
    main()
