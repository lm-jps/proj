!
!  translated to f95 Nov. 2009
!===============================================================
!   NOTE THE INTERPOLATION ROUTINE STARTS AT BIN 0 NOT 1        
!								
!			Cart_to_Polar				
!								
!        	      Bradley W. Hindman				
!		      version 2: June 19, 1998		        
!                     rewritten as subroutine                   
!                     Deborah Haber (DH) Oct 22, 2007           
!                     added azimuthal average Jan. 2009 - DH      
!                     works now with ringanalysis_sngl.f  
!                             16 Sept 2009 - DH                      
!                     passes arrays to be written to files
!                         back to ringanalysis_sngl.f to be
!                         passed to c wrapper for actual output
!								
!	This routine converts the power spectra contained in	
!   powxy from Cartesian (kx,ky,nu) coordinates to polar		
!   (theta,k,nu) coordinates.  The resulting power spectra has	
!   pixels evenly spaced in k and in theta.  Cubic convolution	
!   interpolation is used to convert from the initially evenly	
!   spaced rectangular data.					
!	This routine outputs the data in array "powthtk",	
!   which holds the data in polar format,	  		
!								
!		(azimuth, wavenumber, frequency).			
!								
!  Inputs_____________________________________________________	
!
!    powxy     Input array  = the (nkx,nky,nnu) original spectrum. Real*4	
!    powthtk   Output array = the unwrapped (ntht,nk,nnu) spectrum. Real*4	
!    x0        Center of powxy in x direction in pixels. Real*4 = crpix1.
!    y0        Center of powxy in y direction in pixels. Real*4 = crpix2.
!								
!  Data Array_________________________________________________	
!
!    powbuf(nkx,nky)	Cartesian frequency slice		
!								
!  Book Keeping Variables_____________________________________	
!
!    nkx,nky		The number of data elements in the 	
!    ntht,nk	        appropriate dimension of the Carte-	
!			   sian or polar coordinate systems.	
!    nnu		Number of data elements in the frequen-	
!			   cy dimension in both Cartesian and	
!			   polar coordinates.			
!    ikx,iky,inu	Index variables for the appropriate	
!        itht,ik	   data dimensions.			
!    kx,ky		Location polar grid points in Cartesian	
!			   coordinates.				
!     
!    lnbase             Keyword specifying whether powxy is the
!                       power or log(power), passed from main routine
!                       = 0.0 data is in units of power
!                       = 2.718 data is in units of log(power)
!								
!   ROUTINES___________________________________________________ 
!
!	CCINT2(powbuf,nkx,nky,kx,ky)	Function		
!								
!===============================================================
    module cart
        contains
	SUBROUTINE cart_to_polar(powxy,powthtk,nkx,nky,nnu,ntht,nk, &
        &                        x0,y0,lnbase,verbose,ierr)
        implicit none
        integer, parameter :: single = 4
        integer :: nk,nkx,nky,nnu,ntht,ier,ierr
        integer :: ik,inu,itht,ikx,iky,ik3,inu3
        integer :: verbose
        real(kind=single) :: kx,ky,x0,y0,dtht,theta,lnbase,pi
        real(kind=single), dimension(nkx,nky) :: powbuf
        real(kind=single), dimension(nkx,nky,nnu) :: powxy
        real(kind=single), dimension(ntht,nk,nnu) :: powthtk
        real(kind=single), dimension(ntht) :: cosine, sine
        real(kind=single) :: CCINT2

!---------------------------------------------------------------
!   Initialize cosine and sine arrays for later use.	
!---------------------------------------------------------------
       pi=ACOS(-1.0)
       dtht=2.0*pi/real(ntht)
       DO itht=1,ntht
         theta=(itht-1)*dtht
         cosine(itht)=COS(theta)
         sine(itht)=SIN(theta)
       END DO 
!       if (verbose == 1) then
!         write(*,'("cosine")')
!         write(*,'(8(f8.4,2x))')(cosine(itht),itht=1,ntht)
!         write(*,'("sine")')
!         write(*,'(8(f8.4,2x))')(cosine(itht),itht=1,ntht)
!       endif


!----------------------------------------------------------
!  Initialize error catching variable
!----------------------------------------------------------

       ierr=0

!---------------------------------------------------------
!  If debugging write out some intermediary variables
!---------------------------------------------------------

       IF (verbose == 1) THEN
          WRITE(*,'(" nkx, nky, nnu")') 
          WRITE(*,*) nkx, nky, nnu    
          WRITE(*,'(" ntht, dtht, x0, y0")') 
          WRITE(*,*) ntht, dtht, x0, y0    
          WRITE(*,'(" lnbase = ",f8.4)') lnbase
       ENDIF

!---------------------------------------------------------------
!   Step through each azimuth and wave number and interpolate   
!   for each frequency value inu out of total nnu values        
!   If lnbase is not 0.0 THEN the power is already in log space
!   else take log of power before proceeding
!---------------------------------------------------------------

       IF (lnbase /= 0.0) THEN
          DO inu=1,nnu
            DO iky = 1,nky 
              DO ikx = 1,nkx
                 powbuf(ikx,iky)=powxy(ikx,iky,inu)
!                 write(*,'("ikx,iky,inu",3i8)')ikx,iky,inu
!                 write(*,'("powbuf, powxy ",2(1pe14.4,2x))')powbuf(ikx,iky), &
!                      & powxy(ikx,iky,inu)
              END DO 
            END DO 
            DO ik = 1,nk 
              DO itht = 1,ntht
                kx=real(ik)*cosine(itht)+x0
                ky=real(ik)*sine(itht)+y0
                powthtk(itht,ik,inu)=EXP(ccint2(powbuf,nkx,nky,kx,ky))
              END DO 
            END DO 
!            write(*,'("powthtk = ",4pe13.4)') ((powthtk(itht,ik,inu),itht=1,2),ik=10,11)
          END DO 
       ELSE
          DO inu=1,nnu
            DO iky = 1,nky 
              DO ikx = 1,nkx
               if (powxy(ikx,iky,inu) > 0.0) then
                 powbuf(ikx,iky)=alog(powxy(ikx,iky,inu))
               else
                 write(*,*)" Negative power in Cart_to_Polar"
                 ierr=11
                 goto 9999
               endif
              END DO 
            END DO 
            DO ik = 1,nk 
              DO itht = 1,ntht
                kx=real(ik)*cosine(itht)+x0
                ky=real(ik)*sine(itht)+y0
                powthtk(itht,ik,inu)=EXP(ccint2(powbuf,nkx,nky,kx,ky))
              END DO 
            END DO 
          END DO 
       END IF
       IF (verbose == 1)THEN
         WRITE(*,*)" cart powthtk " 
         ik3=nk/3
         inu3=nnu/3
         WRITE(*,20)((powthtk(1,ik,inu),ik=ik3,ik3+2),inu=inu3,inu3+1)
 20      format(3(2(1pe13.5,2x)))
         WRITE(*,*)" Finishing Cart_to_Polar"
       ENDIF
 9999  IF (ierr .NE. 0) THEN
         WRITE(*,'(" Problem in cart_to_polar subroutine ")')
         WRITE(*,'(" ierr = ",i4)') ierr
       ENDIF
       RETURN 
       END subroutine cart_to_polar
     end module cart 
