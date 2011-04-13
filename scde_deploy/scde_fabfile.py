"""Deployment file for the Stem Cell Discovery Engine.
"""
from fabric.api import *
from fabric.contrib.files import *
import yaml

from bcbio.deploy import environ, java
from bcbio.deploy.shared import (_fetch_and_unpack)

# ## Install targets

def install_scde():
    config = _read_config()
    _setup_environment(config)
    config = _install_prereqs(config)
    _install_bii(config)

# ## BII installation

def _install_bii(config):
    """Main target for installing the bioinvestigator index.
    """
    dirs = _install_bii_tools(config)
    print dirs
    _configure_bii(dirs, config)
    _configure_manager(dirs, config)

def _install_bii_tools(config):
    mgr_version = "1.2"
    base_dir = os.path.join(config["base_install"], config["bii_dirname"])
    if not exists(base_dir):
        sudo("mkdir -p %s" % base_dir)
        sudo("chown %s %s" % (env.user, base_dir))
    with cd(base_dir):
        git_cmd = "git clone git://github.com/ISA-tools/BioInvIndex.git"
        bii_dir = _fetch_and_unpack(git_cmd)
        bii_dir = os.path.join(base_dir, bii_dir)
    url = "https://github.com/downloads/ISA-tools/" \
          "ISAvalidator-ISAconverter-BIImanager/BII-data-mgr-%s.zip" % mgr_version
    with cd(base_dir):
        mgr_dir = _fetch_and_unpack(url)
        mgr_dir = os.path.join(base_dir, mgr_dir)
    index_dir = os.path.join(base_dir, config["bii_lucene_index"])
    if not exists(index_dir):
        run("mkdir -p %s" % index_dir)
    bii_data_dir = os.path.join(base_dir, config["bii_data_dirname"])
    if not exists(bii_data_dir):
        run("mkdir -p %s" % bii_data_dir)
    return dict(base=base_dir, bii=bii_dir, mgr=mgr_dir, index=index_dir,
                bii_data=bii_data_dir)

def _configure_bii(dirs, config):
    """Site specific configuration for BII.
    """
    pass

def _configure_manager(dirs, config):
    """Site specific configuration for the upload manager.
    """
    _configure_manager_jdbc(dirs, config)
    _configure_manager_hibernate(dirs, config)
    _configure_manager_datalocation(dirs, config)

def _configure_manager_jdbc(dirs, config):
    run_file = os.path.join(dirs["mgr"], "run.sh")
    old_str = "JDBCPATH=/path/to/jdbc_driver.jar"
    jdbc_str = "JDBCPATH=%s" % config["jdbc_driver"]
    if contains(run_file, old_str):
        sed(run_file, old_str, jdbc_str)

def _configure_manager_hibernate(dirs, config):
    hib_file = os.path.join(dirs["mgr"], "config", "hibernate.properties")
    props = ["hibernate.connection.username", "hibernate.connection.password",
             "hibernate.search.default.indexBase"]
    defaults = ["CHANGEME", "CHANGEME", "/tmp/bii/testdb/lucene"]
    new_vals = [config["db_user"], config["db_pass"], dirs["index"]]
    test_val = "%s=%s" % (props[0], new_vals[0])
    if not contains(hib_file, test_val):
        run("cp %s.postgresql %s" % (hib_file, hib_file))
        for prop, new_val, default in zip(props, new_vals, defaults):
            old = "%s=%s" % (prop, default)
            new = "%s=%s" % (prop, new_val)
            if not contains(hib_file, new):
                sed(hib_file, old, new)

def _configure_manager_datalocation(dirs, config):
    def dl_str(prop, info):
        return '%s="%s' % (prop, info)
    dl_file = os.path.join(dirs["mgr"], "config", "data_locations.xml")
    props = ["filesystem_path", "web_url"]
    orig_vals = ["/tmp/bii/test_repo", "file:///tmp/bii/test_repo/"]
    new_vals = [dirs["bii_data"], "http://%s/%s/" % (config["host"], config["bii_data_url"])]
    if contains(dl_file, orig_vals[1]):
        for prop, orig_val, new_val in zip(props, orig_vals, new_vals):
            sed(dl_file, dl_str(prop, orig_val), dl_str(prop, new_val))

# ## Setup and system pre-requesites

def _install_prereqs(config):
    java.install_maven(config["base_install"])
    jboss_dir = java.install_jboss(config["base_install"])
    config["jboss"] = jboss_dir
    jdbc_driver = java.install_postgresql_jdbc(config["base_install"])
    config["jdbc_driver"] = jdbc_driver
    return config

def _setup_environment(config):
    if config["distribution"] == "centos":
        environ.setup_centos()
    else:
        raise NotImplementedError
    if env.hosts == ["vagrant"]:
        environ.setup_vagrant_environment()
    result = run("echo $HOSTNAME")
    config["host"] = result.strip()

def _read_config():
    config_dir = os.path.join(os.path.dirname(__file__), "config")
    config_file = os.path.join(config_dir, "scde.yaml")
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    return config
