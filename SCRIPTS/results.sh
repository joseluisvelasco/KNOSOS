s=0.`pwd|rev|cut -f1 -d_|rev|cut -f2-10 -d0`
for i in `seq 1 500`
do
 echo $s
done > s
paste s results.data >> ../results.data
paste s results.knosos >> ../results.knosos
rm s
