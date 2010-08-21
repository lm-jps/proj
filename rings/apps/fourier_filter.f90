!===============================================================
!                                                               
!                       FourierFilter                           
!                                                               
!                       Bradley W. Hindman                             
!                       Sept 23, 1998                           
!                       adapted for HMI                          
!                       Deborah A. Haber                        
!                       Oct. 31, 2007                             
!                       Sept. 15, 2009
!                       Sept. 29, 2009
!                       revised to pass filter and filtered power
!                         to ringanalysis_sngl.f to be passed
!                         back to c wrapper for writing to drms
!                       revised jan. 18, 2010 to comment out some of 
!                         debugging info                                                                 
!       Depending on what iflgs are set, this routine applies   
!   a multiplicative filter to the spectra in the powthtk array 
!   to remove low order fourier	components from the azimuthal   
!   variation (order one is not	removed) and subsamples the     
!   data to reduce	
!   the size of the data array.  This subsampling is preceded	
!   by removing high order fourier components directly from	
!   the data.     						
!	The routine creates two output arrays to return to      
!   ringanalysis.f.  The "filtarr" array contains the filter    
!   applied and "powfilt" is the filtered and subsampled     
!   data set.							
!                                                               
!   Book Keeping And Dimension Variables ________________________________ 
!
!       ntht,nk         The number of data elements in the        
!                          theta and k dimensions of the polar  
!                          coordinate system.  (input)  	
!
!       nthts       	The number of data elements in the      
!                          azimuthal dimension of the sub-	
!			   sampled data set.   (input)	
!
!       nnu             Number of data elements in the frequen- 
!                         cy dimension in both Cartesian and    
!                         polar coordinates.   (input)                  
!
!       ik,itht,inu     Index variables for the appropriate     
!                          data dimensions.    (book keeping)             
!
!       inumn            1.5 mHz in bins       (calculated)                 
!
!       inumx            5.5 mHz in bins       (calculated)                
!     
!       verbose          Flag set for debugging (input)
!                           = 1  intermediate values and comments are printed
!                           = 0  no comments printed 
!                                                                
!   Variables Used For Subsampling And Filtering -------------------------- 
! 
!       avg             Array holding the average of theta dependence in 
!                           power between 1.5 and 5.5 mHz as a function 
!                           of theta and k.
!
!       buf             Intermediate buffer arrays  Real*4
!       bufs            dimensions (ntht) and (nthts) respectively.
!
!       cbuf            Complex buffer array used for Fourier transforms.
! 
!       dnu             Delta nu in units of microHz.   (input)          
!
!       j               Number of elements to skip when	subsampling.
!                           (calculated)	
! 
!	mdata		Number of Fourier components to keep	
!			   when subsampling the data.  (input)		
!
!	mfilt		Highest fourier component to use in the 
!			   filter.   (input)			
!
!       powthtk         Unwrapped power spectrum (Real*4).  
!                        (input)   Dimensions (ntht,nk,nnu)
!
!       powfilt         Subsampled and Filtered power spectrum (Real*8).    
!                        (output)  Dimensions (nthts,nk,nnu)
! 
!       thetaavg        Power spectrum averaged in theta for each k and nu.
!                           (calculated)
!
!       filter          Array containing the filter that was applied
!                        to obtain powfilt,  (output) Dimensions(nthts,nk) 
!   
!       work            Work space for Fourier Transform routine
!                         Dimension (2*ntht)    (scratch space)
!                                                               
!   ROUTINES___________________________________________________ 
!	FFTCC(cbuf,n,sgn,cbuf,work,ier)	Subroutine		
!                                                               
!===============================================================

      module fourier
        contains
	SUBROUTINE fourier_filter(powthtk,powfilt,filter,ntht,nk,nnu, &
                  nthts,dnu,mdata,mfilt,inumn,inumx,verbose,ierr)
        IMPLICIT NONE
        integer, parameter :: single = 4
        integer, parameter :: double = 8
	INTEGER :: ntht,nk,nthts,nnu,ik,inu,j,ier,ierr,itht
        INTEGER :: mdata,mfilt,verbose,inumn,inumx,ik3,inu3
	REAL(kind=single) :: thetaavg
	REAL(kind=single), dimension(ntht,nk,nnu) :: powthtk
        REAL(kind=single), dimension(nthts) :: bufs
        REAL(kind=single), dimension(ntht) :: buf
        REAL(kind=single), dimension(2*ntht) :: work
        REAL(kind=single), dimension(nk,nnu) :: nrmlz
        REAL(kind=single), dimension(ntht,nk) :: avg
        REAL(kind=single), dimension(nthts,nk) :: filter
        REAL(kind=single), dimension(nthts,nk,nnu) :: powfilt
	COMPLEX(kind=single), dimension(ntht) :: cbuf
	REAL(kind=double) :: dnu, dnu2

!---------------------------------------------------------------

        ierr = 0

!
!  setup formats for use in verbose mode
!

!                                                               
!   Determine subsampling spacing		
!                                                               
!---------------------------------------------------------------

       j=int(ntht/nthts)
       ik3=nk/3
       inu3=nnu/3

       IF (verbose .EQ. 1) THEN 
          WRITE(*,'("Executing FourierFilter")')
          WRITE(*,'("Subsample by factor of = ",i4)')j	
          write(*,*)" nk, ntht, nthts, nnu "
          write(*,*) nk, ntht, nthts, nnu 
          write(*,'("dnu = ",1pd13.5)') dnu
          write(*,'("mdata =",i5)')mdata
       ENDIF

!---------------------------------------------------------------
!								
!   Calculate the azimuthal average for each frequency and	
!   wave number and store in avg(,).  This will be used to	
!   construct the filter.  The only frequencies used in the	
!   average are between frequency bins 150 and 550.  At the 	
!   same time remove all fourier components of order higher than 
!   mdata from the data.  Subsample the data and store it in 	
!   powfilt.							 
!								
!   Change Feb. 16, 1999 to make the program general to data of 
!   different time lengths.  Make it average data between 1.5   
!   and 5.5 mHz                                                 
!
!   Change Sept. 29, 2009
!   Get dnu in mHz from header in main routine and pass it through.
!   Use it to determine inumx and inumn (bins corresponding to 
!   5.5 mHz and 1.5 mHz respectively).          
!---------------------------------------------------------------
     
! 
!  change dnu to units of mHz 
!
        dnu2 = dnu / 1000.0d0
        inumn=nint(1.5/sngl(dnu2))
        inumx=nint(5.5/sngl(dnu2))
        if (verbose .eq. 1) then
           write(*,'(" inumn,inumx ",2I6)') inumn, inumx
           WRITE(*,10)((powthtk(1,ik,inu),ik=ik3,ik3+2),inu=inu3,inu3+1)
        endif
 
 10     format(2(3(1pe13.5,2x)))
!--------------------------------------------------------------
!   Initialize filter and average arrays
!-------------------------------------------------------------
        avg=0.0
        filter=0.0
	DO inu=1,nnu
	   DO ik=1,nk

!---------------------------------------------------------------
!   Calculate the theta average for use in normalization.		
!---------------------------------------------------------------

              DO itht=1,ntht
                buf(itht)=powthtk(itht,ik,inu) 
!                if (verbose .eq. 1)then 
!                  write(*,'("buf,powthtk,itht,ik,inu",2(1pe15.5),3i6)') &
!                    buf(itht),powthtk(itht,ik,inu),itht,ik,inu
!                endif
              ENDDO
	      thetaavg=0.0
!              write(*,*)thetaavg, ntht, buf(1200)
	      DO itht=1,ntht
		 thetaavg=thetaavg+buf(itht)
	      ENDDO
          
!              write(*,*)thetaavg
	      nrmlz(ik,inu)=thetaavg/REAL(ntht)
!              if ((verbose == 1) .and. (inu == inu3)) then
!                write(*,'("thetaavg = ", 1pe13.5)')thetaavg
!                write(*,'("nrmlz(ik3,inu3) = ", 1pe13.5)') nrmlz(ik3,inu3)
!              endif
              thetaavg=thetaavg*REAL(inumx-inumn)/REAL(ntht)
!                write(*,'("thetaavg = ", 1pe13.5)')thetaavg
!                write(*,'("inu inumn inumx ntht", 4i5)')inu, inumn, inumx, ntht
              DO itht=1,ntht
	        IF ((inu .GE. inumn) .AND. (inu .LE. inumx)) THEN
                  avg(itht,ik)=avg(itht,ik)+buf(itht)/thetaavg
                ENDIF
	        cbuf(itht)=CMPLX(buf(itht),0.0)
              ENDDO
!---------------------------------------------------------------
!   Remove the high order fourier components.			
!---------------------------------------------------------------

!              if ((verbose .eq. 1) .and. (inu .eq. 1)) then
!                write(*,'("made it to here")')
!                write(*,'("cbuf before fft ",2(4(1pe14.4,2x),2x))') &
!                   & (real(cbuf(itht)),itht=1,4),(aimag(cbuf(itht)),itht=1,4)	
!              endif
	      CALL FFTCC(cbuf,ntht,+1.0,cbuf,work,ier)
              IF (ier .NE. 0) THEN
                write(*,'("Error calculating first FFT = ",i6)') ier
                ierr = 21
                GOTO 999
              ENDIF
!              if ((verbose == 1) .and. (inu == inu3)) then
!                write(*,'("cbuf after fft ",4(2(1pd14.4,2x),2x))') &
!                    (cbuf(itht),itht=1,4)	
!              endif

	      DO itht=2+mdata,ntht-mdata
	         cbuf(itht)=CMPLX(0.0,0.0)
	      ENDDO

	      DO itht=1,ntht
		 cbuf(itht)=cbuf(itht)/REAL(ntht)
	      ENDDO
!                write(*,'("cbuf before 2nd fft ",2(4(1pe14.4,2x),2x))') &
!                   & (real(cbuf(itht)),itht=1,4),(aimag(cbuf(itht)),itht=1,4)	
	      CALL FFTCC(cbuf,ntht,-1.0,cbuf,work,ier)
              IF (ier .NE. 0) THEN
                write(*,'("Error calculating second FFT = ",i6)') ier
                ierr = 22
                GOTO 999
              ENDIF
!              if ((verbose == 1) .and. (inu == inu3)) then
!               write(*,'("cbuf after 2nd fft",2(4(1pe14.4,2x),2x))') &
!                       (real(cbuf(itht)),itht=1,4), (aimag(cbuf(itht)),itht=1,4)	
!              endif 

!---------------------------------------------------------------
!   Subsample the data.                      	 		
!---------------------------------------------------------------
	      DO itht=1,nthts
	         powfilt(itht,ik,inu)=real(cbuf((itht-1)*j+1))
	      ENDDO
	   ENDDO
	ENDDO
!-------------------------------------------------------------------	
!   Apply thefilter to the data and write out both the filter and the	
!   filtered data if write_filter=1 and write_filt_pow=1
!--------------------------------------------------------------------

!---------------------------------------------------------------
!   Compute the filter from the theta average by selecting only those
!   Fourier components less than mfilt.						
!---------------------------------------------------------------

           DO ik = 1, nk
	      DO itht=1,ntht
	         cbuf(itht)=CMPLX(avg(itht,ik),0.0)
	      ENDDO
	      CALL FFTCC(cbuf,ntht,+1.0,cbuf,work,ier)
              IF (ier .NE. 0) THEN
                 WRITE(*,'("Error calculating the third FFT (used in &
                      computing the filter) = ", i6)') ier
                 ierr = 23
                 GOTO 999
              ENDIF

!---------------------------------------------------------------   
!   Check to see data if verbose = 1
!---------------------------------------------------------------
               
!              IF (verbose .EQ. 1) THEN
!	         DO itht=1,mfilt+1
!	          WRITE (*,15) ik,itht,REAL(cbuf(itht)),IMAG(cbuf(itht)) 
!	         ENDDO
!              ENDIF
 15           format(i5,i5,1pe13.5,2x,1pe13.5)

	      DO itht=2+mfilt,ntht-mfilt
	         cbuf(itht)=CMPLX(0.0,0.0)
	      ENDDO
	      cbuf(2)=CMPLX(0.0,0.0)
	      cbuf(ntht)=CMPLX(0.0,0.0)
	      DO itht=1,ntht
	         cbuf(itht)=cbuf(itht)/REAL(ntht)
	      ENDDO
	      CALL FFTCC(cbuf,ntht,-1.0,cbuf,work,ier)
              IF (ier .NE. 0) THEN
                 WRITE(*,'("Error calculating fourth  FFT = ", i6)') ier
                 ierr = 24
                 GOTO 999
              ENDIF
	      DO itht=1,ntht
	         buf(itht)=REAL(cbuf(itht))
	         avg(itht,ik)=buf(itht)
	      ENDDO
           ENDDO

!---------------------------------------------------------------
!   Apply the filter and renormalize to conserve power		
!---------------------------------------------------------------

	DO inu=1,nnu
          DO ik=1,nk
	     thetaavg=0.0
	     DO itht=1,nthts
                filter(itht,ik)=avg((itht-1)*j+1,ik)
	        bufs(itht)=powfilt(itht,ik,inu)/filter(itht,ik)
	        thetaavg=thetaavg+bufs(itht)/REAL(ntht)
	     ENDDO
	     DO itht=1,nthts
	        powfilt(itht,ik,inu)=bufs(itht)*nrmlz(ik,inu)/thetaavg
	     ENDDO
          ENDDO
	ENDDO
        IF (verbose .EQ. 1) THEN
          WRITE (*,20)((powfilt(1,ik,inu),ik=ik3,ik3+2),inu=inu3,inu3+2) 
          WRITE(*,'(" Finishing FourierFilter")')
        END IF
 20     format(3(3(1pe13.5,2x)))
 999    RETURN 
	END subroutine
      end module
