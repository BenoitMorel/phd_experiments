
files=$1

sed -i 's/)\([0-9]*\)/)0.\1/g' ${files}
