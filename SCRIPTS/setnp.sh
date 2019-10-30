
test -e input.surfaces 
if [ $? -eq 0 ]
then 
  ns=`cat input.surfaces|grep NS|cut -f2 -d=`
else
  test -e ../input.surfaces
  if [ $? -eq 0 ]
  then
    ns=`cat ../input.surfaces|grep NS|cut -f2 -d=`
  else
    ns=1
  fi
fi

test -e input.model
if [ $? -eq 0 ]
then
  grep TASK3D=.TRUE. input.model > /dev/null 
  if [ $? -eq 0 ]
  then
    ns=40
  fi
  grep NEOTRANSP=.TRUE. input.model > /dev/null 
  if [ $? -eq 0 ]
  then
    ns=7
  fi
else
  test -e ../input.model
  if [ $? -eq 0 ]
  then
     grep TASK3D=.TRUE. ../input.model > /dev/null
     if [ $? -eq 0 ]
     then
       ns=40
     fi
     grep NEOTRANSP=.TRUE. ../input.model > /dev/null
     if [ $? -eq 0 ]
     then
       ns=7
     fi
  fi
fi

test -e input.parameters 
if [ $? -eq 0 ] 
then 
  nerr=`cat input.parameters|grep NERR|cut -f2 -d=`
else 
  test -e ../input.parameters
  if [ $? -eq 0 ]
  then
    nerr=`cat ../input.parameters|grep NERR|cut -f2 -d=`
  fi
fi
if [ -z "$nerr" ]
then 
  nerr=1
fi



np=`awk -v f=1 -v nerr=$nerr -v ns=$ns 'BEGIN{print ( f*ns*nerr) }'`
if [ $np -gt 64 ]
then 
  np=$ns
fi



if [ $np -eq 1  ]
then
  nodes=1
  ppn=1
elif [ $np -gt 1 ] && [ $np -le 8 ]
then
  nodes=1
  ppn=8
elif [ $np -gt 8 ] && [ $np -le 16 ]
then
  nodes=2
  ppn=8
elif [ $np -gt 16 ] && [ $np -le 24 ]
then
  nodes=3
  ppn=8
elif [ $np -gt 24 ] && [ $np -le 32 ]
then
  nodes=4
  ppn=8
elif [ $np -gt 32 ] && [ $np -le 40 ]
then
  nodes=5
  ppn=8
elif [ $np -gt 40 ] && [ $np -le 48 ]
then
  nodes=6
  ppn=8
elif [ $np -gt 48 ] && [ $np -le 54 ]
then
  nodes=7
  ppn=8
elif [ $np -gt 56 ] && [ $np -le 64 ]
then
  nodes=8
  ppn=8
fi	
