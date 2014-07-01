#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#

system "module load bio/BAPS";
echo $SCRIPTDIR;

echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  echo Setting up environment variables
  MCRROOT=$1
  echo ---
  MWE_ARCH="glnxa64" ;
  if [ "$MWE_ARCH" = "sol64" ] ; then
	LD_LIBRARY_PATH=.:/usr/lib/lwp:${MCRROOT}/runtime/glnxa64 ; 
  else
  	LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  fi
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
  if [ "$MWE_ARCH" = "mac" -o "$MWE_ARCH" = "maci" ]; then
	DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/System/Library/Frameworks/JavaVM.framework/JavaVM:/System/Library/Frameworks/JavaVM.framework/Libraries;
  else
	MCRJREVER=`cat ${MCRROOT}/sys/java/jre/glnxa64/jre.cfg` ; 
	echo Found MCR Java JRE version: $MCRJREVER ; 
	MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre${MCRJREVER}/lib/amd64 ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; 
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;  
  fi
  XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
  export LD_LIBRARY_PATH;
  export XAPPLRESDIR;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
  shift 1
~croucher/bin/BAPS/BAPS_package/baps6 $*
fi
exit
