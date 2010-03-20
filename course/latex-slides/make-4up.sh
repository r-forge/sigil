#!/bin/sh

if [ $# != 1 ]
then
	echo "Usage:  make-4up.sh my_slides.pdf"
else
	outfile=`basename $1 .pdf`.4up.pdf
	echo "$1 => $outfile"
	pdfnup --nup 2x2 --paper a4paper --orient auto --trim "-8mm -8mm -8mm -8mm" --outfile "$outfile" "$1"
fi
