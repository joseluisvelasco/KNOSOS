#!/bin/bash -l
# Initial working directory:
#SBATCH -D ./
# Memory usage [MB] of the job is required, 3800 MB per task:
#SBATCH --mem=30400
#SBATCH --time=04:30:00

ulimit -s unlimited

source $HOME/KNOSOS/SCRIPTS/clean.sh
test -e phi2d_mar*
if [ $? -eq 0 ]
then
  ln -s phi2d_mar* ph1_2d.in
fi

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
srun -n $SLURM_NTASKS $HOME/KNOSOS/SOURCES/knosos.x

test -e flux.knosos.00
if [ $? -eq 0 ]
then

    for n in $(seq 1 100)
    do
    
	if [ $n -eq 1 ]
	then
	    sdir=flux.knosos
	elif [ $n -eq 2 ]
	then
	    sdir=flux.amb
	elif [ $n -eq 3 ]
	then
	    sdir=B.map
	elif [ $n -eq 4 ]
	then
	    sdir=varphi1.map
	elif [ $n -eq 5 ]
	then
	    sdir=Er.map
	elif [ $n -eq 6 ]
	then
	    sdir=imp.knosos
        elif [ $n -eq 7 ]
        then
            sdir=knososTASK3D.flux
	elif [ $n -eq 8 ]
	then
	    sdir=knososTASK3D.ambEr
	elif [ $n -eq 9 ]
	then
	    sdir=knosos.dk
	elif [ $n -eq 9 ]
	then
	    sdir=flux.amb.comp
	elif [ $n -eq 10 ]
	then
	    sdir=varphi1.map.comp
	elif [ $n -eq 11 ]
	then
	    sdir=flux.modes
	fi
	
	dir=$sdir.00
	test -e $dir
	if [ $? -eq 0 ]
	then
           head -1 $dir > $sdir
	   for dir in $sdir.??
	   do
	       cat $dir|grep -v "\[" >> $sdir
	   done
	fi
    done
fi

test -e log.tar.gz
if [ $? -eq 1 ]
then
  mkdir LOG
  for dir in *.modes B.map ???_2d.in fort.* *.[0-9][0-9]
  do
    test -e $dir
    if [ $? -eq 0 ]
    then
	mv $dir LOG
    fi 
  done
  test -e LOG/flux.av
  if [ $? -eq 0 ]
  then
      mv LOG/flux.av .
  fi 
  touch LOG/*
  tar cfz ./log.tar.gz LOG/* && rm -r LOG
fi


