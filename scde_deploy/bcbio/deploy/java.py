"""Install targets for java packages.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from bcbio.deploy.shared import (_if_not_installed, _fetch_and_unpack)

def _java_dir(base_dir):
    java_dir = os.path.join(base_dir, "java")
    if not exists(java_dir):
        sudo("mkdir -p %s" % java_dir)
        sudo("chown %s %s" % (env.user, java_dir))
    return java_dir

@_if_not_installed("mvn")
def install_maven(base_dir):
    version = "2.2.1"
    mirror = "mirror.cc.columbia.edu/pub/software/apache"
    java_dir = _java_dir(base_dir)
    url = "http://%s/maven/binaries/apache-maven-%s-bin.tar.gz" % ( mirror, version)
    with cd(java_dir):
        mvn_dir = _fetch_and_unpack(url)
    install_mvn = os.path.join(base_dir, "bin", "mvn")
    dl_mvn = os.path.join(java_dir, mvn_dir, "bin", "mvn")
    sudo("ln -s %s %s" % (dl_mvn, install_mvn))

def install_jboss(base_dir):
    version = "5.1.0"
    url = "http://downloads.sourceforge.net/project/jboss/JBoss/JBoss-%s.GA/jboss-%s.GA.zip" % \
          (version, version)
    java_dir = _java_dir(base_dir)
    with cd(java_dir):
        jboss_dir = _fetch_and_unpack(url)
    return os.path.join(java_dir, jboss_dir)

def install_postgresql_jdbc(base_dir):
    version = "9.0-801"
    url = "http://jdbc.postgresql.org/download/postgresql-%s.jdbc3.jar" % version
    java_dir = _java_dir(base_dir)
    jdbc_dir = os.path.join(java_dir, "postgresql_jdbc")
    out_jar = os.path.basename(url)
    if not exists(jdbc_dir):
        run("mkdir -p %s" % jdbc_dir)
    with cd(jdbc_dir):
        if not exists(out_jar):
            run("wget %s" % url)
    return os.path.join(jdbc_dir, out_jar)
