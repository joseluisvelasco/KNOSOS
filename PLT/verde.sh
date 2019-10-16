file=`ls -rt *.ps  | tail -1`
perl -p -i -e s/"\/LC1 \{0 1 0\} def"/"\/LC1 \{0 0.5 0\} def"/ $file
