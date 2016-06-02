for i in */*/pp-report.tex; do
  # get dir and file names
  dn=`dirname $i`;
  fn=`basename $i`;
  # new file name
  nf=`echo $dn | sed 's/\//_/'`
  cd $dn; pdflatex $fn;
  # copy rrerport
  cp pp-report.pdf ../../rep/$nf.pdf;
  cd ../../;
done
