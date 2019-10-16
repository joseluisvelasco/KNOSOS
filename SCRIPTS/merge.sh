
test -e flux.knosos.00
if [ $? -eq 0 ]
then

    for n in $(seq 1 10)
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
	    sdir=knososTASK3D.flux
	elif [ $n -eq 6 ]
	then
	    sdir=knososTASK3D.ambEr
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
  mv fort.* *.[0-9][0-9] LOG
  tar cvfz log.tar.gz LOG && rm -r LOG
fi



