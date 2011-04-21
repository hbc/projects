"""Automated upload and installation of ISA-tab and Galaxy data
"""
from fabric.api import *
from fabric.contrib.files import *

import yaml
import boto

from bcbio.deploy import environ, java
from bcbio.deploy.shared import (_fetch_and_unpack)

# ## Install targets

def scde_data():
    config = _read_config()
    _setup_environment(config)
    config = _install_prereqs(config)
    _add_isatab_to_bii(config)

# ## Add ISA-tab prepared data to a BII instance

def _add_isatab_to_bii(config):
    dirs = _setup_bii_directories(config)
    dirs = _get_isatab_config(dirs, config)
    for upload_dir, perm_type in _get_isatab_uploads(dirs, config):
        _upload_isatab_to_bii(upload_dir, perm_type, dirs, config)
    _run_bii_cl("org.isatools.isatab.commandline.ReindexShellCommand",
                [], dirs, config)

def _upload_isatab_to_bii(upload_dir, perm_type, dirs, config):
    """Use SimpleManager target from Eamonn to upload to database.
    """
    if _check_if_uploaded(upload_dir, dirs):
        return
    args = ["load", upload_dir, dirs["config"]]
    out = _run_bii_cl("org.isatools.isatab.manager.SimpleManager",
                      args, dirs, config)
    if out.find("Loading failed") >= 0 or out.find("FATAL:") >= 0:
        raise ValueError("Upload to BII failed")
    for study in _get_studies(upload_dir):
        if perm_type == "public":
            args = ["-p", study]
        elif perm_type == "private":
            user = _add_bii_user(dirs, config)
            args = ["-v", study, "-o", user]
        _run_bii_cl("org.isatools.isatab.commandline.PermModShellCommand",
                    args, dirs, config)

def _add_bii_user(dirs, config):
    args = ["-e", config["galaxy_adminusers"].split(",")[0],
            "-p", config["db_pass"],
            "-f", "Administrator",
            "-r", "curator",
            "-n", "Private",
            "-s", "Administrator",
            config["db_user"]]
    _run_bii_cl("org.isatools.isatab.commandline.UserAddShellCommand",
                args, dirs, config)
    return config["db_user"]

def _run_bii_cl(main, args, dirs, config):
    common = "-Xms256m -Xmx1024m -XX:PermSize=64m -XX:MaxPermSize=128m"
    isatools_jar = os.path.join(dirs["mgr"], "isatools_deps.jar")
    classpath = "%s:%s" % (isatools_jar, config["jdbc_driver"])
    with settings(warn_only=True):
        with cd(os.path.dirname(isatools_jar)):
            out = run("java %s -cp %s %s %s" % (common, classpath, main,
                                                " ".join(args)))
    return out

def _check_if_uploaded(upload_dir, dirs):
    """Look at existing loaded metadata to determine if experiment is loaded.
    """
    study_name = _get_studies(upload_dir)[0]
    check_glob = os.path.join(dirs["bii_data"], "meta_data", "study_%s_*" % study_name)
    with settings(hide("everything"), warn_only=True):
        result = run("ls -1d %s" % check_glob)
    if result.find("No such file") >= 0:
        return False
    else:
        return True

def _get_studies(isatab_dir):
    """Retrieve a list of studies from an ISA-tab directory.
    """
    study_text = "Study Identifier"
    file_to_check = "i_Investigation.txt"
    with settings(hide("everything")):
        result = run("grep '%s' %s/%s" % (study_text, isatab_dir, file_to_check))
    out = []
    for line in result.split("\n"):
        out.append(line.replace(study_text, "").replace('"', "").strip())
    return out

def _get_isatab_uploads(dirs, config):
    """Retrieve directories of ISA-tab files to upload from S3.
    """
    isatab_cfile = os.path.join(config["fab_config_dir"], "isatab_studies.yaml")
    conn = boto.connect_s3()
    bucket = conn.get_bucket(config["isatab_s3_bucket"])
    with open(isatab_cfile) as in_handle:
        isatab_config = yaml.load(in_handle)
    for permission_type in ["public", "private"]:
        for isatab_name in isatab_config.get(permission_type, []):
            isatab_dir = _download_isatab_s3(isatab_name, permission_type, bucket, dirs)
            yield isatab_dir, permission_type

def _download_isatab_s3(name, perm_type, bucket, dirs):
    """Download an ISA-tab data file from the S3 bucket.
    """
    s3_key = bucket.get_key("%s.tar.gz" % name)
    assert s3_key is not None, name
    url = s3_key.generate_url(expires_in=60)
    out_dir = os.path.join(dirs["isatab"], perm_type)
    if not exists(out_dir):
        run("mkdir -p %s" % out_dir)
    with cd(out_dir):
        cur_dir = _fetch_and_unpack(url)
    return os.path.join(out_dir, cur_dir)

def _setup_bii_directories(config):
    base_dir = os.path.join(config["base_install"], config["bii_dirname"])
    mgr_glob = os.path.join(base_dir, "ISA*", "*manager*", "target", "BII-data-mgr-*")
    with settings(hide("everything"), warn_only=True):
        result = run("ls -d1 %s" % mgr_glob)
    mgr_dir = result.split()[0].strip()
    index_dir = os.path.join(base_dir, config["bii_lucene_index"])
    sudo("chown -R %s %s" % (env.user, index_dir))
    isatab_data_dir = os.path.join(base_dir, config["isatab_data_dirname"])
    if not exists(isatab_data_dir):
        run("mkdir -p %s" % isatab_data_dir)
    bii_data_dir = os.path.join(base_dir, config["bii_data_dirname"])
    return dict(base=base_dir, mgr=mgr_dir, isatab=isatab_data_dir,
                bii_data=bii_data_dir)

def _get_isatab_config(dirs, config):
    """Retrieve ISA-tab configuration files.

    Need to add SCDE configurations to GitHub repo.
    """
    config_git = "git clone git://github.com/ISA-tools/Configuration-Files.git"
    with cd(dirs["base"]):
        config_dir = _fetch_and_unpack(config_git)
    our_config_dir = os.path.join(dirs["base"], config_dir, config["isatab_config"])
    assert exists(our_config_dir), our_config_dir
    dirs["config"] = our_config_dir
    return dirs

# ## Setup and system pre-requesites

def _install_prereqs(config):
    config["jdbc_driver"] = java.install_postgresql_jdbc(config["base_install"])
    return config

def _setup_environment(config):
    if config["distribution"] == "centos":
        environ.setup_centos()
    elif config["distribution"] == "ubuntu":
        environ.setup_ubuntu()
    else:
        raise NotImplementedError
    if env.hosts == ["vagrant"]:
        environ.setup_vagrant_environment()
    result = run("echo $HOSTNAME")
    config["host"] = result.strip()
    run("chmod a+rx $HOME")
    shell_config = None
    for test in [".bash_profile", ".bashrc"]:
        if exists(test):
            shell_config = test
            break
    assert shell_config is not None
    config["shell_config"] = shell_config

def _read_config():
    config_dir = os.path.join(os.path.dirname(__file__), "config")
    config_file = os.path.join(config_dir, "scde.yaml")
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    config["fab_base_dir"] = os.path.dirname(config_dir)
    config["fab_config_dir"] = config_dir
    return config
