#!/bin/bash

# tool to merge csv files from TopCount into one big csv file
# should be in the folder with all the CSV files

# user inputs experiment date
echo "Exp date DDMMYY:"
read EXP

FINALCSV="$(echo $EXP$"_all.csv")"

# get the barcode number
for i in *.* # take the first file in the folder
do
  FILE="$i"
  break 1
done

BAR="$(echo "$i" | cut -d'.' -f 1)"

# merge csv files
ulimit -Sn 10240
touch $FINALCSV # create the big csv file
paste -d , "$BAR".* > "$FINALCSV"
ulimit -Sn 256
sed -i '' '$d' "$FINALCSV"
mv "$FINALCSV" ..

echo ">> Created $FINALCSV"
