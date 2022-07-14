#! /bin/tcsh

set list=`cat`

set i=0
while ($i < $#list)
  @ i1 = $i + 1
  @ i2 = $i + 2
  @ i3 = $i + 3

  @ index0 = ($list[$i2] / 1920) * 1920
  @ index1 = $index0 + 1920
  if ($list[$i2] + $list[$i3] < $index1) then
    echo $list[$i1-$i3]
  else
    set ntot=$list[$i3]
    @ n0 = $index1 - $list[$i2]
    echo $list[$i1] $list[$i2] $n0
    @ ntot = $ntot - $n0
    @ ind = $index1
    while ($ntot > 1920)
      echo $list[$i1] $ind 1920
      @ ind = $ind + 1920
      @ ntot = $ntot - 1920
    end
    echo $list[$i1] $ind $ntot
  endif
@ i = $i + 3
end
