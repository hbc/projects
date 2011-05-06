Automated deployment for [Galaxy][1] and the
[BioInvestigatorIndex (BII)][2] to support the [Stem Cell Discovery Engine][3].

## Galaxy and BII Installation

- Install [VirtualBox][5], [Vagrant][4] and Python

- Build libraries

        cd ~/hsph/projects/scde_deploy
        python setup.py build && sudo python setup.py install

- Edit `config/scde.yaml` to specify system directories and passwords

- Install CloudBioLinux box with Vagrant. Creates a `Vagrantfile` for
  configuring and running the virtual machine

        mkdir tmp/vagrant/galaxy_bii
        cd tmp/vagrant/galaxy_bii
        vagrant box add biolinux_centos_20110412 https://s3.amazonaws.com/chapmanb/biolinux_centos_20110412.box
        vagrant init biolinux_centos_20110412

- Allow web access to your machine from a port on your computer. Edit
  the `Vagrantfile` to forward port 80 from the machine to 8082 on
  your local computer:

        config.vm.forward_port "http", 80, 8082

- Start the virtual machine:

        cd tmp/vagrant/galaxy_bii
        vagrant up

- Run fabric script to install dependencies, BII and Galaxy; all
  servers are started in screen sessions:

        cd tmp/vagrant/galaxy_bii
        fab -f ~/hsph/projects/scde_deploy/scde_fabfile.py -H vagrant install_scde

- Browse BII and Galaxy:

http://localhost:8082/bioinvindex
http://localhost:8082/galaxy

[1]: http://usegalaxy.org
[2]: http://isatab.sourceforge.net/
[3]: http://discovery.hsci.harvard.edu/
[4]: http://vagrantup.com/
[5]: http://www.virtualbox.org/

## ISA-Tab record management

- Directories of ISA-Tab metadata, raw and derived data files are
  stored in Amazon S3 buckets for reliable backup

- Edit `config/isatab_studies.yaml` to specify public and private
  datasets that should be downloaded.

- Download data from S3 bucket and install into BII

       fab -f scde_data_fabfile.py -H localhost scde_data
