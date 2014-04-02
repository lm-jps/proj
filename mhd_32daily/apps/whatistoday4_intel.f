       program whatday
       implicit none
*       character*8 sdate ! for Fujitsu/SUN
       character*9 sdate ! for Intel
       integer itm(3)
       character*3 smonth
       integer idatetime(8)
       character*10 bigben(3)

       call date_and_time(bigben(1),bigben(2),bigben(3),idatetime)
       itm(3) = idatetime(3) ! DD
       itm(2) = idatetime(2) ! MM
       itm(1) = idatetime(1) ! YYYY
       itm(1) = mod(itm(1),100)

       if (itm(2) .EQ.  1) smonth='Jan'
       if (itm(2) .EQ.  2) smonth='Feb'
       if (itm(2) .EQ.  3) smonth='Mar'
       if (itm(2) .EQ.  4) smonth='Apr'
       if (itm(2) .EQ.  5) smonth='May'
       if (itm(2) .EQ.  6) smonth='Jun'
       if (itm(2) .EQ.  7) smonth='Jul'
       if (itm(2) .EQ.  8) smonth='Aug'
       if (itm(2) .EQ.  9) smonth='Sep'
       if (itm(2) .EQ. 10) smonth='Oct'
       if (itm(2) .EQ. 11) smonth='Nov'
       if (itm(2) .EQ. 12) smonth='Dec'
       write(*,'(i2.2,''-'',a3,''-'',i2.2)') itm(3),smonth,itm(1)
*       write(*,'(i2.2,''-'',a3,''-'',i2.2)') itm(1),smonth,itm(3)
       stop
       end
