#! /usr/bin/ruby

# ruby_args.rb
#
# This script generates coverage maps for 
# the observables, sharp, and fd10 data.
# 
# Obserables are defined as: hmi.S_720s, hmi.M_720s, hmi.M_45s, hmi.V_45s
# Sharp are defined as     : hmi.sharp_720s
# fd10 are defined as      : hmi.ME_720s_fd10 and hmi.ME_720s_nrt
#
#  Usage:
#
#  > ruby maps_args.rb [/path/to/output/] [type] [YYYY] [MM]
#    where [type] equals 'fd10', 'observables', 'sharp', or 'all'
#    where [YYYY] equals a four-digit year typed as an int
#    where [MM]   equals a two-digit month typed as an int
#    where [/path/to/output/] equals the path to the directory
#          that contains all the output files (the assumption is that
#          this path is not /web/jsoc/htdocs/doc/data/hmi/coverage_maps/.
#
#  E.g. 
#
#  > ruby maps_args.rb /path/to/output all 2013 01
#  > ruby maps_args.rb /path/to/output observables 2013 01
#
# Written by Monica Bobra 31 January 2013

# determine the date range for which to generate coverage maps
require 'date'

indexes = []

puts "Year is " "#{ARGV[2]}"
puts "Month is " "#{ARGV[3]}"

type = ARGV[1]
year = ARGV[2].to_i
month = ARGV[3].to_i

next_month = month + 1
next_month = 1 if next_month > 12
next_year  = year
next_year += 1 if next_month == 1

beginning_of_time     = Time.utc(2010, 1, 1).to_i
month_start_timestamp = Time.utc(year, month, 1).to_i  - beginning_of_time
month_end_timestamp   = Time.utc(next_year, next_month, 1).to_i - beginning_of_time
  
label = "#{Date::MONTHNAMES[month]} #{year}"
month = month > 9 ? month.to_s : "0#{month}"
next_month = next_month > 9 ? next_month.to_s : "0#{next_month}"

`cat ~thailand/stats_lev0/HMI/HMI_LOOPOPEN_LIST.#{year}.#{month}* > EOI_LOOP_#{year}_#{month}.txt`

if type == 'observables' or type == 'all'
  # make all the necessary show_info and show_coverage calls
  #### OBSERVABLES ####
  # hmi.V_45s

  # figures out the missing data for dopplergrams
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.V_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss.txt`

  # figures out the unknown data for the dopplergrams
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.V_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'UNK' | ./reform_coverage.csh > #{year}_#{month}_unk.txt`

  # figures out the missing data under mask for dopplergrams 
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.V_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI mask=0xffffffff | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss_mask.txt`

  # hmi.M_45s
  
  # figures out the missing data for 45 sec mags
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss_m.txt`

  # figures out the unknown data for 45 sec mags
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'UNK' | ./reform_coverage.csh > #{year}_#{month}_unk_m.txt`

  # figures out the missing data under mask for 45 sec mags
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI mask=0xffffffff | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss_mask_m.txt`

  # hmi.S_720s 

  # figures out the missing data for 720 sec stokes
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.S_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'MISS' | ./reform_coverage_720.csh > #{year}_#{month}_miss_s720.txt`

  # figures out the unknown data for 720 sec stokes
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.S_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'UNK' | ./reform_coverage_720.csh > #{year}_#{month}_unk_s720.txt`

  # figures out the missing data under mask for 720 sec stokes
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.S_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI mask=0xffffffff | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss_mask_s.txt`

  # hmi.M_720s
  
  # figures out the missing data for 720 sec manetograms
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'MISS' | ./reform_coverage_720.csh > #{year}_#{month}_miss_s720_m.txt`

  # figures out the unknown data for 720 sec magnetograms
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'UNK' | ./reform_coverage_720.csh > #{year}_#{month}_unk_s720_m.txt`

  # figures out the missing data under mask for 720 sec magnetograms
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI mask=0xffffffff | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss_mask_s_m.txt`

  idl = <<-IDLSESSION
  
  	  .r observables.pro

	  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	  ; MAPS FOR hmi.V_45s
	  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	  print, 'Printing Maps for hmi.V_45s'
	  coveragemap, monthnumber=#{month}, junstart=#{month_start_timestamp}L, junend=#{month_end_timestamp}L, monthname='#{label}', file1='#{year}_#{month}_miss.txt',file2='EOI_LOOP_#{year}_#{month}.txt', file3='#{year}_#{month}_miss_mask.txt', file5='#{year}_#{month}_unk.txt', outfile='#{year}_#{month}_map.ps', ds='hmi.V_45s', cad=45L

	  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;
	  ; MAPS FOR hmi.M_45s
	  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;

	  print, 'Printing Maps for hmi.M_45s'
	  coveragemap, monthnumber=#{month}, junstart=#{month_start_timestamp}L, junend=#{month_end_timestamp}L, monthname='#{label}', file1='#{year}_#{month}_miss_m.txt', file2='EOI_LOOP_#{year}_#{month}.txt', file3='#{year}_#{month}_miss_mask_m.txt', file5='#{year}_#{month}_unk_m.txt',  outfile='#{year}_#{month}_map_m.ps', ds='hmi.M_45s', cad=45L

	  ; ;;;;;;;;;;;;;;;;;;;;;;;;;
	  ; MAPS FOR hmi.S_720s
	  ; ;;;;;;;;;;;;;;;;;;;;;;;;;

	  print, 'Printing Maps for hmi.S_720s'
	  coveragemap, monthnumber=#{month}, junstart=#{month_start_timestamp}L, junend=#{month_end_timestamp}L, monthname='#{label}', file1='#{year}_#{month}_miss_s720.txt', file2='EOI_LOOP_#{year}_#{month}.txt', file3='#{year}_#{month}_miss_mask_s.txt',  file5='#{year}_#{month}_unk_s720.txt',  outfile='#{year}_#{month}_map_s720s.ps', ds='hmi.S_720s', cad=720L

	  ; ;;;;;;;;;;;;;;;;;;;;;;
	  ; MAPS FOR hmi.M_720s
	  ; ;;;;;;;;;;;;;;;;;;;;;;

	  print, 'Printing Maps for hmi.M_720s'
	  coveragemap, monthnumber=#{month}, junstart=#{month_start_timestamp}L, junend=#{month_end_timestamp}L, monthname='#{label}', file1='#{year}_#{month}_miss_s720_m.txt', file2='EOI_LOOP_#{year}_#{month}.txt', file3='#{year}_#{month}_miss_mask_s_m.txt',  file5='#{year}_#{month}_unk_s720_m.txt',  outfile='#{year}_#{month}_map_m720s.ps', ds='hmi.M_720s', cad=720L

     print, 'Finished making all maps. Closing IDL session, making .pngs, and copying files over to JSOC.'

  IDLSESSION

  `/home3/ssw.new/gen/setup/ssw_idl << endofidlsession #{idl}`

  Dir.chdir ARGV[0]

  # Convert all files to png; rotate all files
  (`/bin/ls *.ps | cut -d . -f1`).split.each do |file|
    `convert -density 300 #{file}.ps -resize 100% #{file}.png`
    `convert -rotate -90 #{file}.png #{file}.png`
  end

  # Move Observables
  `ls *map*.png`
  `mv *map*.png /web/jsoc/htdocs/doc/data/hmi/coverage_maps/observables/`
end



if type == 'fd10' or type == 'all'
  # make all the necessary show_info and show_coverage calls
  #### OBSERVABLES ####
  # hmi.V_45s

  # figures out the missing data for dopplergrams
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.V_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss.txt`

  # figures out the unknown data for the dopplergrams
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.V_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'UNK' | ./reform_coverage.csh > #{year}_#{month}_unk.txt`

  # figures out the missing data under mask for dopplergrams 
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.V_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI mask=0xffffffff | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss_mask.txt`

  # hmi.M_45s
  
  # figures out the missing data for 45 sec mags
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss_m.txt`

  # figures out the unknown data for 45 sec mags
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'UNK' | ./reform_coverage.csh > #{year}_#{month}_unk_m.txt`

  # figures out the missing data under mask for 45 sec mags
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_45s[][2]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI mask=0xffffffff | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss_mask_m.txt`

  # hmi.S_720s 

  # figures out the missing data for 720 sec stokes
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.S_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'MISS' | ./reform_coverage_720.csh > #{year}_#{month}_miss_s720.txt`

  # figures out the unknown data for 720 sec stokes
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.S_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'UNK' | ./reform_coverage_720.csh > #{year}_#{month}_unk_s720.txt`

  # figures out the missing data under mask for 720 sec stokes
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.S_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI mask=0xffffffff | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss_mask_s.txt`

  # hmi.M_720s
  
  # figures out the missing data for 720 sec manetograms
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'MISS' | ./reform_coverage_720.csh > #{year}_#{month}_miss_s720_m.txt`

  # figures out the unknown data for 720 sec magnetograms
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'UNK' | ./reform_coverage_720.csh > #{year}_#{month}_unk_s720_m.txt`

  # figures out the missing data under mask for 720 sec magnetograms
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI mask=0xffffffff | grep 'MISS' | ./reform_coverage.csh > #{year}_#{month}_miss_mask_s_m.txt`

  #### Milne-Eddington Inverted Data ###
  # Query hmi.M_720s for missing values in definitive data

  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_720s[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'MISS' | ./reform_coverage_720.csh > #{year}_#{month}_miss_fd10_s720_m.txt`

  # Query hmi.M_720s_nrt for missing values in nrt data

  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_coverage -qi ds='hmi.M_720s_nrt[][1]' low=#{year}.#{month}.01_00_TAI high=#{next_month < month ? next_year : year}.#{next_month}.01_00_TAI | grep 'MISS' | ./reform_coverage_720.csh > #{year}_#{month}_miss_fd10_s720_m.txt_nrt`

  # % of pixels in hmi.S_720s  
  
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_info -q key=t_rec_index,RSUN_OBS,CDELT1 ds=hmi.S_720s'[#{year}.#{month}.01_00:00:00_TAI-#{next_month < month ? next_year : year}.#{next_month}.01_00:00:00_TAI][? QUALITY >= 0 ?]' > #{year}_#{month}_pixels.txt`

  # INVNPRCS in hmi.ME_720s_fd10 

  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_info -q key=t_rec_index,INVNPRCS ds=hmi.ME_720s_fd10'[#{year}.#{month}.01_00:00:00_TAI-#{next_month < month ? next_year : year}.#{next_month}.01_00:00:00_TAI]' > #{year}_#{month}_invnprcs.txt`

  # combine pixels and INVNPRCS for definitive data 
 
  `join -t'	'  #{year}_#{month}_invnprcs.txt #{year}_#{month}_pixels.txt > #{year}_#{month}_join_fd10.txt`

  # % of pixels in hmi.S_720s_nrt 
 
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_info -q key=t_rec_index,RSUN_OBS,CDELT1 ds=hmi.S_720s_nrt'[#{year}.#{month}.01_00:00:00_TAI-#{next_month < month ? next_year : year}.#{next_month}.01_00:00:00_TAI][? QUALITY >= 0 ?]' > #{year}_#{month}_pixels.txt_nrt`
 
  # INVNPRCS in hmi.ME_720s_fd10_nrt 
 
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_info -q key=t_rec_index,INVNPRCS ds=hmi.ME_720s_fd10_nrt'[#{year}.#{month}.01_00:00:00_TAI-#{next_month < month ? next_year : year}.#{next_month}.01_00:00:00_TAI]' > #{year}_#{month}_invnprcs.txt_nrt`

  # combine pixels and INVNPRCS for nrt data 
 
  `join -t'	'  #{year}_#{month}_invnprcs.txt_nrt #{year}_#{month}_pixels.txt_nrt > #{year}_#{month}_join_fd10.txt_nrt`
 
   idl = <<-IDLSESSION
   
   	  ; ;;;;;;;;;;;;;;;;;;;;;;
	  ; MAPS FOR hmi.ME_720s_fd10
	  ; ;;;;;;;;;;;;;;;;;;;;;;

      .r inversion.pro

      print, 'Printing Maps for hmi.ME_720s_fd10'

      coveragemap, monthnumber=#{month}, junstart=#{month_start_timestamp}L, junend=#{month_end_timestamp}L, monthname='#{label}', file1='#{year}_#{month}_miss_fd10_s720_m.txt', file4='#{year}_#{month}_join_fd10.txt', outfile='#{year}_#{month}_fd10.ps', ds='hmi.ME_720s_fd10', cad=720L

	  ; ;;;;;;;;;;;;;;;;;;;;;;
	  ; MAPS FOR hmi.ME_720s_fd10_nrt
	  ; ;;;;;;;;;;;;;;;;;;;;;;

      print, 'Printing Maps for hmi.ME_720s_fd10_nrt'

      coveragemap, monthnumber=#{month}, junstart=#{month_start_timestamp}L, junend=#{month_end_timestamp}L, monthname='#{label}', file1='#{year}_#{month}_miss_fd10_s720_m.txt_nrt', file4='#{year}_#{month}_join_fd10.txt_nrt',outfile='#{year}_#{month}_fd10_nrt.ps', ds='hmi.ME_720s_fd10_nrt', cad=720L

      print, 'Finished making all maps. Closing IDL session, making .pngs, and copying files over to JSOC.'

  IDLSESSION

  `/home3/ssw.new/gen/setup/ssw_idl << endofidlsession #{idl}`
  
  Dir.chdir ARGV[0]

  # Convert all files to png; rotate all files
  (`/bin/ls *.ps | cut -d . -f1`).split.each do |file|
    `convert -density 300 #{file}.ps -resize 100% #{file}.png`
    `convert -rotate -90 #{file}.png #{file}.png`
  end

  # Montage fd10
  files_fd10 = (`/bin/ls *fd10.png | cut -d . -f1`).split
  files_nrts = (`/bin/ls *fd10_nrt.png | cut -d . -f1`).split
  files_maps = (`/bin/ls *map_m720s.png | cut -d . -f1`).split

  files_maps.zip(files_fd10).each do |file1, file2|
      `montage #{file1}.png #{file2}.png -background none -tile 2x -geometry +0+0  #{file2}_montage.png`
  end

  files_maps.zip(files_nrts).each do |file1, file2|
      `montage #{file1}.png #{file2}.png -background none -tile 2x -geometry +0+0  #{file2}_montage.png`
  end

  # Move fd10
  `mv *montage*.png /web/jsoc/htdocs/doc/data/hmi/coverage_maps/fd10/`
end



if type == 'sharp' or type == 'all'

  #### SHARP data ###
  # Query the SHARP series:
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_info -q key=T_REC_index ds='hmi.sharp_720s[][#{year}.#{month}.00_00_TAI-#{next_month < month ? next_year : year}.#{next_month}.00_00_TAI][? (NPIX>2000) and (NPIX<400000) ?]' | sort | uniq -c > #{year}_#{month}_sharp.txt`

  # Query the MHARP series:
  `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_info -q key=T_REC_index ds='hmi.mharp_720s[][#{year}.#{month}.00_00_TAI-#{next_month < month ? next_year : year}.#{next_month}.00_00_TAI][? (NPIX>2000) and (NPIX<400000) ?]' | sort | uniq -c > #{year}_#{month}_mharp.txt`

  # Join the two series:
  # e.g. join -1 2 -2 2 tmpa.txt tmpb.txt > tmpc.txt
  `join -1 2 -2 2 #{year}_#{month}_sharp.txt #{year}_#{month}_mharp.txt > #{year}_#{month}_join_sharp.txt`

  idl = <<-IDLSESSION
  
  	  ; ;;;;;;;;;;;;;;;;;;;;;;
	  ; MAPS FOR hmi.sharp_720s
	  ; ;;;;;;;;;;;;;;;;;;;;;;

      .r sharp_coverage.pro

	  print, 'Printing Maps for hmi.sharp_720s'

      coveragemap, monthnumber=#{month}, junstart=#{month_start_timestamp}L, junend=#{month_end_timestamp}L, monthname='#{label}', file1='#{year}_#{month}_join_sharp.txt',outfile='#{year}_#{month}_sharp.ps', ds='hmi.sharp_720s', cad=720L

      print, 'Finished making all maps. Closing IDL session, making .pngs, and copying files over to JSOC.'

  IDLSESSION

  `/home3/ssw.new/gen/setup/ssw_idl << endofidlsession #{idl}`
  
  Dir.chdir ARGV[0]

  # Convert all files to png; rotate all files
  (`/bin/ls *.ps | cut -d . -f1`).split.each do |file|
    `convert -density 300 #{file}.ps -resize 100% #{file}.png`
    `convert -rotate -90 #{file}.png #{file}.png`
  end

  # Move Sharps
  `ls *sharp.png`
  `mv *sharp.png /web/jsoc/htdocs/doc/data/hmi/coverage_maps/sharp/`
  
end

# Clean up
# Cannibalize all text files
`rm *.txt`
`rm *.txt_nrt`

# Cannibalize all .ps files
`rm *.ps`

# Cannibalize all remaining .png files
`rm *.png`
