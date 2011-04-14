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
    _configure_system(config)
    #_install_bii(config)
    _start_servers(config)

# ## Deployed servers

def _start_servers(config):
    jboss_bin = os.path.join(config["jboss"], "bin")
    jboss_run = os.path.join(jboss_bin, "run.sh")
    jboss_stop = os.path.join(jboss_bin, "shutdown.sh")
    jboss_args = "-Djboss.service.binding.set=ports-01"
    with settings(hide="everything", warn_only=True):
        sudo("%s -S" % jboss_stop)
    sudo("screen -d -m %s %s >/dev/null 2>&1" % (jboss_run, jboss_args))

# ## BII installation

def _install_bii(config):
    """Main target for installing the bioinvestigator index.
    """
    dirs = _install_bii_tools(config)
    _configure_bii(dirs, config)
    _configure_manager(dirs, config)
    _deploy_bii(dirs, config)

def _deploy_bii(dirs, config):
    with cd(dirs["bii"]):
        run("mvn clean package install -Pdeploy,postgresql,index_local " \
            "-Dmaven.test.skip=true")
        ear_reg = os.path.join("ear", "target", "bii*ear")
        ear_file = run("ls %s" % ear_reg).strip().split()[0]
        jboss_dir = os.path.join(config["jboss"], "server", "default", "deploy")
        run("cp %s %s" % (ear_file, jboss_dir))

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
    _configure_bii_webapp(dirs, config)
    _configure_bii_profile(dirs, config)

def _configure_bii_profile(dirs, config):
    pro_file = os.path.join(dirs["bii"], "profiles.xml")
    props = ["jdbc.username", "jdbc.password", "hibernate.search.default.indexBase"]
    orig_vals = ["CHANGEME", "CHANGEME", "/tmp/bii/luceneindex"]
    new_vals = [config["db_user"], config["db_pass"], dirs["index"]]
    def make_str(prop, val):
        return "<%s>%s" % (prop, val)
    if contains(pro_file, orig_vals[0]):
        for prop, orig, new in zip(props, orig_vals, new_vals):
            sed(pro_file, make_str(prop, orig), make_str(prop, new))

def _configure_bii_webapp(dirs, config):
    web_dir = os.path.join(dirs["bii"], "web", "src", "main", "webapp")
    with settings(hide('everything')):
        to_mod = [f.strip() for f in run("ls -1 %s/*.xhtml" % web_dir).split()]
    ignore = ("debug",)
    def tpl_str(name):
        return "layout/%s.xhtml" % name
    for fname in (f for f in to_mod if not os.path.basename(f).startswith(ignore)):
        if not contains(fname, tpl_str(config["bii_layout_template"])):
            sed(fname, tpl_str("local-template"),
                tpl_str(config["bii_layout_template"]))

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

# ## Configuration for standard system utilities -- nginx, postgresql

def _configure_system(config):
    _configure_postgres(config)
    _configure_nginx(config)

def _configure_nginx(config):
    config_file = "/etc/nginx/nginx.conf"
    local_config = os.path.join(config["fab_base_dir"], "install_files",
                                "nginx.conf")
    if not contains(config_file, "SCDE"):
        sudo("mv %s %s.orig" % (config_file, config_file))
        put(local_config, config_file, use_sudo=True)
        sudo("/etc/init.d/nginx restart")

def _configure_postgres(config):
    """Setup required databases and access permissions in postgresql.
    """
    pg_base = "/var/lib/pgsql"
    _configure_postgres_access(pg_base)
    _create_postgres_db("bioinvindex", pg_base, config)

def _create_postgres_db(db_name, pg_base, config):
    with cd(pg_base):
        with settings(warn_only=True):
            result = sudo("psql --list | grep %s" % db_name, user="postgres")
        if not result.strip():
            user_sql = "CREATE USER %s WITH PASSWORD '%s'" % (config["db_user"],
                                                              config["db_pass"])
            sudo('psql -c "%s"' % user_sql, user="postgres")
            sudo('psql -c "CREATE DATABASE %s WITH OWNER %s"' %
                 (db_name, config["db_user"]), user="postgres")

def _configure_postgres_access(pg_base):
    def restart():
        with settings(warn_only=True):
            sudo("/etc/init.d/postgresql restart")
    conf_file = os.path.join(pg_base, "data", "postgresql.conf")
    hba_file = os.path.join(pg_base, "data", "pg_hba.conf")
    orig_auth = "host    all         all         127.0.0.1/32          ident sameuser"
    new_auth = "host    all         all         127.0.0.1/32          md5"
    if (not exists(hba_file) and
          not (contains(hba_file, new_auth, use_sudo=True).strip() == new_auth)):
        restart()
        uncomment(conf_file, "^#listen_addresses = ", use_sudo=True)
        uncomment(conf_file, "^#port = 5432", use_sudo=True)
        comment(hba_file, orig_auth, use_sudo=True)
        append(hba_file, new_auth, use_sudo=True)
        restart()

# ## Setup and system pre-requesites

def _install_prereqs(config):
    java.install_maven(config["base_install"])
    config["jboss"] = java.install_jboss(config["base_install"])
    config["jdbc_driver"] = java.install_postgresql_jdbc(config["base_install"])
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
    run("chmod a+rx $HOME")
    run("chmod a+r .bash_profile")

def _read_config():
    config_dir = os.path.join(os.path.dirname(__file__), "config")
    config_file = os.path.join(config_dir, "scde.yaml")
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    config["fab_base_dir"] = os.path.dirname(config_dir)
    return config
