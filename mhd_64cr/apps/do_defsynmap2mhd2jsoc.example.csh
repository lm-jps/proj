#! /bin/csh -f
unlimit
setenv KMP_STACKSIZE 16M
setenv OMP_NUM_THREADS 1

## path(s), now this script does not need general path
set mypath=/home/keiji/cvs/JSOC/bin/linux_x86_64
# set path=(. /bin /usr/bin)
## some directory
set mytmp="./"
# set mytmp="/home/keiji/temp"

#### switch and flag : 1=do it, others=skip it (for speeding up in test mode)
## fetch data
set lgetdat=1
## run MHD
set lmhdrun=1
## send data to JSOC
set ljsoc=1
## make figures
set lmkfigs=1
## make record etc. on local disk
set lrecord=1
##

## some variables for jsoc_export_as_fits. Notice codes in this script assume Br (or Mr) be input.  
set inseries="hmi.BrBlossynop_720s"
set insegmnt="Br"
set outseries="hmi.MHDcorona"

## this script may assume two sub-directories "waste" and "done" exist.
# mkdir waste done


## do-it

# set crnum=$1
foreach crnum (2145)
 if ($lgetdat == 1) then
  /home/keiji/cvs/JSOC/bin/linux_x86_64/jsoc_export_as_fits reqid=JSOC_20090514_001 expversion=0.5 \
      rsquery=$inseries"["$crnum"]{"$insegmnt"}" path=$mytmp ffmt="{segment}" \
      method=url protocol=FITS cparms="**NONE**" JSOC_DBNAME=jsoc JSOC_DBUSER=keiji
#  /bin/mv -f $mytmp/$insegmnt.fits ./data.fits
#  echo data.fits > hmisynofits2txt_2.lst
  echo $insegmnt.fits > hmisynofits2txt_2.lst
  $mypath/hmisynofits2txt_2 < hmisynofits2txt_2.lst  > tmp1.log
  $mypath/hmitxt2wso >> tmp1.log
## calculate coefficients
#  echo 1000 180 >  gtpotco7.lst
#  echo 1000 180 >> gtpotco7.lst
#  echo 0        >> gtpotco7.lst
#  $mypath/gtpotco7 < gtpotco7.lst >> tmp1.log
##
 endif

## MHD part...
 if ($lmhdrun == 1) then
##
  setenv OMP_NUM_THREADS 8
## Iterative Laplace solver using 1-deg map
  taskset -c 0-7 $mypath/mhd_64raw > 1000.$crnum.log
## Spherical harmonics, order up to 5. (for test)  
#  taskset -c 0-7 $mypath/mhd_64pc5 > 1000.$crnum.log
## if relax-time be set longer than 40 hours, we need to create "nn.lst" file.
# echo 6000 > nn.lst 
# echo 0000 >> nn.lst
# echo 6000 >> nn.lst
# setenv OMP_NUM_THREADS 8
# taskset -c 0-7 $mypath/mhd_64raw > 1000.$crnum.log
## Then, adjust output file name for post-processing(s) that assume the 40-hour time relaxation 
# /bin/mv -f x306000.1000 x304000.1000
##
  /bin/mv -f init1000.log init.dat
  setenv OMP_NUM_THREADS 1
 else
  /bin/cp -p -f test_dat/init.dat .
  /bin/cp -p -f test_dat/x304000.$crnum ./x304000.1000
 endif
## end of MHD part ##

## post process(es)
 echo 1000 0000 2 >  xy2uv.1000_4000.lst
 echo 1000 4000 2 >> xy2uv.1000_4000.lst
 echo  -1 -1 -1   >> xy2uv.1000_4000.lst
 $mypath/xy2uv_64 <  xy2uv.1000_4000.lst > tmp3.log

## send data to JSOC
 if ($ljsoc == 1) then
  $mypath/mhd2equidist_txt4jsoc_72x64x128to72x144x80
  $mypath/mhdtxt2jsoc_64cr crdate=$crnum out=$outseries synomap=$inseries
  /bin/mv -f d3equidist144x72.txt done/d3equidist144x72.$crnum.txt
 endif

## make some files and figures
 if ($lmkfigs == 1) then
  echo $crnum > carr.txt
  echo 180   >> carr.txt
  echo 90.0 180.0 > latlon.txt

## lon-lat plot
  echo 9       >  slcalt15.mhddef.lst
  echo 3 3 9   >> slcalt15.mhddef.lst
  echo 2 1 1 1 >> slcalt15.mhddef.lst
  echo 1000 4000   0  3 1 0   2.0   >> slcalt15.mhddef.lst
  echo 1000 4000   0  1 1 0   5.0e6 >> slcalt15.mhddef.lst
  echo 1000 4000   0  6 1 0   1.0e5 >> slcalt15.mhddef.lst
  echo 1000 4000  35  3 1 0  10.0   >> slcalt15.mhddef.lst
  echo 1000 4000  35  1 1 0   2.5e4 >> slcalt15.mhddef.lst
  echo 1000 4000  35  6 1 0   3.0e3 >> slcalt15.mhddef.lst
  echo 1000 4000  48  3 1 0   5.0   >> slcalt15.mhddef.lst
  echo 1000 4000  48  1 1 0   5.0e2 >> slcalt15.mhddef.lst
  echo 1000 4000  48  6 1 0   2.5e1 >> slcalt15.mhddef.lst
  $mypath/slcalt15_64cr < slcalt15.mhddef.lst > tmp4.log
  /bin/mv -f slal1501.ps slaldef.new.ps
  /usr/bin/gs -q -dNOPAUSE -dBATCH -dGraphicsAlphaBits=4 -dTextAlphaBits=4 \
        -g840x560 -r144 \
        -sDEVICE=ppmraw -sOutputFile=slaldef.new.ppm \
        -c save pop 0 neg 330 neg translate -f slaldef.new.ps
  /home/jsoc/bin/linux_x86_64/pnmtojpeg slaldef.new.ppm > slaldef.$crnum.new.jpg
  /usr/bin/gs -q -dNOPAUSE -dBATCH -dGraphicsAlphaBits=4  -dTextAlphaBits=4 \
        -g420x280 -r72 \
        -sDEVICE=ppmraw -sOutputFile=slaldef.small.new.ppm \
        -c save pop 0 neg 330 neg translate -f slaldef.new.ps
  /home/jsoc/bin/linux_x86_64/pnmtojpeg slaldef.small.new.ppm > slaldef.$crnum.small.new.jpg
  /bin/mv -f slaldef.new.ps  slaldef.$crnum.new.ps
## make csv file
  echo 60 1000 4000 > mhdati2csv.lst
  $mypath/mhdati2csv_64 < mhdati2csv.lst > tmp5.log
  /bin/mv -f solwind_vr.csv solwind_vr.$crnum.csv
  /bin/mv -f solwind_br.csv solwind_br.$crnum.csv
## make MK4-like plot
  echo 3 >  mkcgsw9w.mhddef.lst
  echo 2 >> mkcgsw9w.mhddef.lst 
  echo 2 >> mkcgsw9w.mhddef.lst 
  echo 1 >> mkcgsw9w.mhddef.lst 
  echo 1000  4000 5.0 0 7 0 0 0   10  90.0  180.0 0 0 0 7 255 234 175  -1.0 >> mkcgsw9w.mhddef.lst
  $mypath/mkcgsw9w_64cr < mkcgsw9w.mhddef.lst > tmp6.log
  /bin/mv -f sw0001.ppm swtoday.ppm
  /home/jsoc/bin/linux_x86_64/pnmtojpeg swtoday.ppm > swtoday.$crnum.jpg
## make plain text volume data
  $mypath/mhd2equidist_64 > tmp7.log
  /bin/mv -f d3equidist.dat d3equidist0.$crnum.dat
  /usr/bin/gzip -f d3equidist0.$crnum.dat
## make fits data, this may be made indendpently of JSOC module
  $mypath/mhd2equidist4_64_fits > tmp8.log
  /bin/mv -f d3equidist.fits  d3equidist0.$crnum.fits

# endif for lmkfigs
 endif

## finalize, clean-up, rename for record
 if ($lrecord == 1) then
  /bin/mv -f stn1000.dat  stn1000.$crnum.dat
  /bin/mv -f hmi2wso.txt hmi2wso.$crnum.txt
  /bin/mv -f hminew.txt  hminew.$crnum.txt
  /bin/mv -f $insegmnt.fits $insegmnt.$crnum.fits
  /bin/mv -f carr.txt  carr.$crnum.txt
  /bin/mv -f x304000.1000 x304000.$crnum
  /bin/mv -f *.$crnum *.$crnum.* done/
  /bin/mv -f *.ppm *.ps u3*.???? *.log x3*.???? *.txt *.lst *.dat  waste/
 else
  /bin/rm -f *.ppm *.ps u3*.???? *.log x3*.???? *.txt *.lst *.dat *.fits
 endif

 echo ...def mhd for $crnum done
end

echo def-MHD loop done
## end of this file ##
