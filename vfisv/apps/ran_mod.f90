module ran_mod 
! module contains three functions 
! ran1 returns a uniform random number between 0-1 
! spread returns random number between min - max 
! normal returns a normal distribution

contains 
    function ran1()  !returns random number between 0 - 1 
        use cons_param
        implicit none 
        real(dp) ran1,x 
        call random_number(x) ! built in fortran 90 random number function 
        ran1=x 
    end function ran1

    function spread(min,max)  !returns random number between min - max 
        use cons_param
        implicit none 
        real(dp) spread 
        real(dp) min,max 
        spread=(max - min) * ran1() + min 
    end function spread

    function normal(mean,sigma) !returns a normal distribution 
        use cons_param
        implicit none 
        real(dp) normal,tmp 
        real(dp) mean,sigma 
        integer flag 
        real(dp) fac,gsave,rsq,r1,r2 
        save flag,gsave 
        data flag /0/ 
        if (flag.eq.0) then 
        rsq=2.0_dp 
            do while(rsq.ge.1.0_dp.or.rsq.eq.0.0_dp) ! new from for do 
                r1=2.0_dp*ran1()-1.0_dp
                r2=2.0_dp*ran1()-1.0_dp 
                rsq=r1*r1+r2*r2 
            enddo 
            fac=sqrt(-2.0_dp*log(rsq)/rsq) 
            gsave=r1*fac 
            tmp=r2*fac 
            flag=1 
        else 
            tmp=gsave 
            flag=0 
        endif 
        normal=tmp*sigma+mean 
        return 
    end function normal

end module ran_mod 
!CVSVERSIONINFO "$Id: ran_mod.f90,v 1.4 2012/04/09 22:21:48 keiji Exp $"
