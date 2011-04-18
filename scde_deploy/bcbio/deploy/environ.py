"""Functions to help with setting up Fabric build environments.
"""
import subprocess

from fabric.api import env

def setup_vagrant_environment():
    """Use vagrant commands to get connection information.
    https://gist.github.com/1d4f7c3e98efdf860b7e
    """
    raw_ssh_config = subprocess.Popen(["vagrant", "ssh-config"],
                                      stdout=subprocess.PIPE).communicate()[0]
    ssh_config = dict([l.strip().split() for l in raw_ssh_config.split("\n") if l])
    env.user = ssh_config["User"]
    env.hosts = [ssh_config["HostName"]]
    env.port = ssh_config["Port"]
    env.host_string = "%s@%s:%s" % (env.user, env.hosts[0], env.port)
    env.key_filename = ssh_config["IdentityFile"]

def setup_centos():
    env.python_version_ext = "2.6"
    if not env.has_key("java_home"):
        env.java_home = "/etc/alternatives/java_sdk"

def setup_ubuntu():
    default_dist = "maverick"
    default_version = "10.10"
    env.sources_file = "/etc/apt/sources.list"
    version = (env.get("dist_name", default_dist),
               env.get("dist_version", default_version))
    sources = [
      "deb http://us.archive.ubuntu.com/ubuntu/ %s universe",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s-updates universe",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s-updates universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s multiverse",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s multiverse",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s-updates multiverse",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s-updates multiverse",
      "ppa:sun-java-community-team/sun-java6", # sun-java
      "deb http://downloads.mongodb.org/distros/ubuntu % s 10gen", # mongodb
      "deb http://cran.stat.ucla.edu/bin/linux/ubuntu %s/", # lastest R versions
      "deb http://nebc.nox.ac.uk/bio-linux/ unstable bio-linux", # Bio-Linux
      "deb http://archive.cloudera.com/debian %s-cdh3 contrib", # Hadoop
      "deb http://ppa.launchpad.net/freenx-team/ppa/ubuntu lucid main", # FreeNX PPA
      "deb http://download.virtualbox.org/virtualbox/debian %s contrib", # virtualbox
    ]
    env.std_sources = _add_source_versions(version, sources)
    env.python_version_ext = ""
    if not env.has_key("java_home"):
        # XXX look for a way to find JAVA_HOME automatically
        env.java_home = "/usr/lib/jvm/java-6-openjdk"

def _add_source_versions(version, sources):
    name, num = version
    final = []
    for s in sources:
        if s.find("%s") > 0:
            s = s % name
        elif s.find("% s") > 0:
            s = s % num
        final.append(s)
    return final


