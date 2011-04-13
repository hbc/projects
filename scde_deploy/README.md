Automated deployment for [Galaxy][1] and the
[BioInvestigatorIndex (BII)][2] to support the [Stem Cell Discovery Engine][3].

- Install [VirtualBox][5], [Vagrant][4] and Python

- Build libraries

      python setup.py build && sudo python setup.py install

- Install CloudBioLinux box with Vagrant

      mkdir tmp/biolinux
      cd tmp/biolinux
      vagrant box add biolinux_centos_20110412 https://s3.amazonaws.com/chapmanb/biolinux_centos_20110412.box
      vagrant init biolinux_centos_20110412
      vagrant up

- Run fabric script to install BII

      cd tmp/biolinux
      fab -f ~/hsph/projects/scde_deploy/scde_fabfile.py -H vagrant install_scde

[1]: http://usegalaxy.org
[2]: http://isatab.sourceforge.net/
[3]: http://discovery.hsci.harvard.edu/
[4]: http://vagrantup.com/
[5]: http://www.virtualbox.org/
