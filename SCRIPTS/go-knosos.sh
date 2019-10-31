PATH_KNOSOS=/home/LHD/velasco/KNOSOS_OFFICIAL #Change this patch accordingly
PATH_KNOSOS=$HOME/KNOSOS #Change this patch accordingly

ln -s input_dkes.vmec ddkes2.data

echo "&model
 PENTA=.TRUE.
/
" > input.model

echo "&parameters
 MAL=32
 MLAMBDA=64
 EFIELD=efield
 CMUL=cmul
/
" > input.parameters

efield=`grep efield input_dkes.vmec|cut -f2 -d=|cut -f1 -d,`
cmul=`grep cmul input_dkes.vmec|cut -f2 -d=|cut -f1 -d,`
perl -p -i -e s/efield/"$efield"/ input.parameters
perl -p -i -e s/cmul/"$cmul"/ input.parameters

$PATH_KNOSOS/SOURCES/knosos.x

cp -p results.knosos results.vmec
