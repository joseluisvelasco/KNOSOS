for dir in s_0*;do vals=0.`echo $dir|cut -f2-10 -d0`;echo $vals;cp -p $HOME/KNOSOS/SCRIPTS/input.surfaces $dir ; perl -p -i -es/vals/$vals/ $dir/input.surfaces ;done
