for dir in kn[0-9]*[ls].sh ;do cp -p kns.sh $dir ; done ; 
perl -p -i -e s/04:30/50:00/ kn*l.sh
perl -p -i -e s/PPN/1/ kn1[sl].sh
perl -p -i -e s/PPN/8/ kn*[sl].sh
perl -p -i -e s/NODES/1/ kn[18][sl].sh
perl -p -i -e s/NODES/2/ kn16[sl].sh
perl -p -i -e s/NODES/3/ kn24[sl].sh
perl -p -i -e s/NODES/4/ kn32[sl].sh
perl -p -i -e s/NODES/5/ kn40[sl].sh
perl -p -i -e s/NODES/6/ kn48[sl].sh
perl -p -i -e s/NODES/7/ kn56[sl].sh
perl -p -i -e s/NODES/8/ kn64[sl].sh

