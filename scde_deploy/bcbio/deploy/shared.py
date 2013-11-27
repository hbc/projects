"""Reusable decorators and functions for custom installations.
"""
import os
import time

from fabric.api import *
from fabric.contrib.files import *

# --- Running servers and daemons

def _is_running(cmd):
    """Check if a given command is currently running.
    """
    with settings(hide('everything')):
        result = run("ps ax | grep '%s'" % cmd)
    is_running = False
    for line in result.split("\n"):
        if line.find(cmd) >= 0 and line.find("grep") == -1:
            is_running = True
            break
    return is_running

def _run_in_screen(name, run_dir, cmd, check_cmd=None, use_sudo=False):
    """Run the given command in a named screen session in the background.

    check_cmd is optional and used to check if the command is running, in cases
    where the running script spawns a process with a different name.
    """
    if check_cmd is None: check_cmd = cmd
    do_run = sudo if use_sudo else run
    send_return = "`echo -ne '\015'`"
    stdout_redirect = ">/dev/null 2>&1"
    if not _is_running(check_cmd):
        with cd(run_dir):
            # Start a detached screen session and then send the command to it
            if use_sudo is False:
                do_run("screen -d -m -S %s %s" % (name, stdout_redirect), pty=False)
            time.sleep(5)
            do_run("screen -S %s -p0 -X stuff '%s'%s %s" % (name, cmd,
                                                            send_return, stdout_redirect))
