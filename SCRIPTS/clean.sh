for dir in slurm* fort.* *knosos* *.map* *.modes* *.amb* *.comp* *.av* ???_2d.in ? *.read log.tar.gz fit.log *.[eo][0-9][0-9][0-9][0-9][0-9][0-9][0-9]* STD*
do 
    test -e $dir
    if [ $? -eq 0 ]
    then
	rm $dir
    fi
done
test -e LOG
if [ $? -eq 0 ]
    then
    rm -r LOG
fi