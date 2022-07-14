* -------------------------------------------------------------------------------------------
* localpot5new_jsoc.f <=== localpot.for <== pot4paul.ver4.for <== pot4mhd7 < == meshpotX
*                                                               Dec. 2009
*                                                               K.Hayashi
* -------------------------------------------------------------------------------------------
*
       subroutine localpot_wrapped_jsoc(imenu,jj,kk,ii,      ! imenu = 0:iteration, otherwise=green
     &                                  magmap0,bx0,by0,bz0) ! MIND ii is for z, jj for x and kk for y.
       implicit none
* number of grid for each three direction.
       integer ii, jj, kk
       integer imenu ! 1 for 
* magnetic field map at solar surface
       real    magmap0(jj,kk)
       real    bz0(jj,kk,0:ii-1)
       real    bx0(jj,kk,0:ii-1)
       real    by0(jj,kk,0:ii-1)
* local
       integer i, j, k
*
* previous common block
*
* MHD variables
       real    magmap(-2:jj+2,-2:kk+2)
       real    bz2(-2:jj+2,-2:kk+2,-2:ii+2)
       real    bx2(-2:jj+2,-2:kk+2,-2:ii+2)
       real    by2(-2:jj+2,-2:kk+2,-2:ii+2)
* mag-temp
       real    bzl(-2:jj+2,-2:kk+2,-2:ii+2)
       real    bxl(-2:jj+2,-2:kk+2,-2:ii+2)
       real    byl(-2:jj+2,-2:kk+2,-2:ii+2)
* coordinate system
       real    zz(-2:ii+2),xx(-2:jj+2),yy(-2:kk+2)
       real    zzb(-2:ii+1),xxb(-2:jj+1),yyb(-2:kk+1)
* scalar for Laplace Eq. solver
       real    lapsi1(-2:jj+2,-2:kk+2,-2:ii+2)
       real    lapsi2(-2:jj+2,-2:kk+2,-2:ii+2)

       write(*, '('' II JJ KK = '',3i5)') ii, jj, kk

       if (imenu .EQ. -100) then
         write(*,*) 'This is debug mode'
         do i = 0, ii - 1
         do k = 1, kk
         do j = 1, jj
           if (i .EQ. 0) then
             bx0(j,k,i) = float(j)
             by0(j,k,i) = float(j)
             bz0(j,k,i) = float(j)
           else
             bx0(j,k,i) = -1.0
             by0(j,k,i) = -2.0
             bz0(j,k,i) = -3.0
           endif
         enddo
         enddo
         enddo
         write(*,*) '   so local_pot calculation will be skipped..'
         return ! for debeg
       endif

* OpenMP parallel start........... any openMP directive should not apprear above.!!!

* First of all, all "large" arrays that will appear as "private"
*               should be initalized here for full subscription ranges.
       call initall(ii,jj,kk,
     &              bzl,bxl,byl,bz2,bx2,by2,
     &              zz,zzb,xx,xxb,yy,yyb,
     &              lapsi1,lapsi2)

* define grid address, position
       call setcrdnt(ii,jj,kk,zz,zzb,xx,xxb,yy,yyb)

* load file or make dummy
       do k = -2, kk + 2
       do j = -2, jj + 2
         magmap(j,k) = 0.0
       enddo
       enddo
       do k = 1, kk
       do j = 1, jj
         magmap(j,k) = magmap0(j,k)
       enddo
       enddo

* solve Lap. eq.
*       if (imenu .EQ. 0) then
       if (imenu .EQ. 1) then
         call solapeq0(ii,jj,kk,magmap,bzl,bxl,byl,
     &                 zz,zzb,xx,xxb,yy,yyb,
     &                 lapsi1,lapsi2)
       else
         call solapeq1(ii,jj,kk,magmap,bzl,bxl,byl,
     &                 zz,zzb,xx,xxb,yy,yyb,
     &                 lapsi1,lapsi2)
       endif

* calculate potential magnetic field ... B = (-) div Psi
       call get3dmag(ii,jj,kk,bz2,bx2,by2,
     &               zz,zzb,xx,xxb,yy,yyb,lapsi1,lapsi2)

* save & plot data
*       inum = 1 ! numbering...
*       call savedata(inum) ! save data
*       call plotmag(inum) ! plot scalar potential & potential field
       do i = 0, ii - 1
       do k = 1, kk
       do j = 1, jj
         bx0(j,k,i) = bx2(j,k,i)
         by0(j,k,i) = by2(j,k,i)
         bz0(j,k,i) = bz2(j,k,i)
       enddo
       enddo
       enddo


       return
       end subroutine



* -------------------------------------------------------------------------------------------
*
*                                  End of Main Program
*
* -------------------------------------------------------------------------------------------
*

* --------------------------------------------------------------------
*
* solve Laplace Eq.
*
* --------------------------------------------------------------------
*
       subroutine solapeq1(ii,jj,kk,magmap,bzl,bxl,byl,
     &                     zz,zzb,xx,xxb,yy,yyb,
     &                     lapsi1,lapsi2)
       implicit none
* interface variables
       integer ii, jj, kk
       real    magmap(-2:jj+2,-2:kk+2)
* coordinate system
       real    zz(-2:ii+2),xx(-2:jj+2),yy(-2:kk+2)
       real    zzb(-2:ii+1),xxb(-2:jj+1),yyb(-2:kk+1)
* scalar for Laplace Eq. solver
       real    lapsi1(-2:jj+2,-2:kk+2,-2:ii+2)
       real    lapsi2(-2:jj+2,-2:kk+2,-2:ii+2)
* temp-mag
       real    bzl(-2:jj+2,-2:kk+2,-2:ii+2)
       real    bxl(-2:jj+2,-2:kk+2,-2:ii+2)
       real    byl(-2:jj+2,-2:kk+2,-2:ii+2)
* local
       integer i, j, k, j2, k2
       integer idummy0
       real    aa, rr
       real    magmapl(-2:jj+2,-2:kk+2)
       real    ds, dsl

       ds = xx(1) - xx(0)
!$omp parallel do private(i,k,j,j2,k2,magmapl,aa,rr,idummy0,dsl)
       do i = 0, ii + 2
         dsl = ds
         do k =  0, kk
         do j =  0, jj
           magmapl(j,k) = magmap(j,k) ! map be at i = 0
         enddo
         enddo
*
         do k = -2, kk + 2
         do j = -2, jj + 2
           aa = 0.0
           do k2 = 0, kk
           do j2 = 0, jj
             idummy0 = (j-j2)**2 + (k-k2)**2
             rr = float(idummy0) + (float(i) + 0.5)**2
             rr = sqrt(rr)
             aa = aa + magmapl(j2,k2) / rr
           enddo
           enddo
           lapsi2(j,k,i) = aa * dsl
         enddo
         enddo
       enddo
!$omp end parallel do

** bottom
!$omp parallel do private(i,k,j)
       do i = -2, ii + 2
         if (i .LT. 0) then ! cond. for i < 0
           do k = -2, kk + 2
           do j = -2, jj + 2
             lapsi2(j,k,i) = magmap(j,k) * (zz(0) - zz(i))
     &                     + lapsi2(j,k,0)
           enddo
           enddo
         endif
       enddo
!$omp end parallel do

* copy Psi
!$omp parallel do private(i,k,j)
         do i = -2, ii + 2
           do k = -2, kk + 2
           do j = -2, jj + 2
             lapsi1(j,k,i) = lapsi2(j,k,i)
           enddo
           enddo
         enddo ! end i-loop
!$omp end parallel do

       return
       end

*
* ------------------------
*
       subroutine solapeq0(ii,jj,kk,magmap,bzl,bxl,byl,
     &                     zz,zzb,xx,xxb,yy,yyb,
     &                     lapsi1,lapsi2)
       implicit none
* interface variables
       integer ii, jj, kk
       real    magmap(-2:jj+2,-2:kk+2)
* coordinate system
       real    zz(-2:ii+2),xx(-2:jj+2),yy(-2:kk+2)
       real    zzb(-2:ii+1),xxb(-2:jj+1),yyb(-2:kk+1)
* scalar for Laplace Eq. solver
       real    lapsi1(-2:jj+2,-2:kk+2,-2:ii+2)
       real    lapsi2(-2:jj+2,-2:kk+2,-2:ii+2)
* temp-mag
       real    bzl(-2:jj+2,-2:kk+2,-2:ii+2)
       real    bxl(-2:jj+2,-2:kk+2,-2:ii+2)
       real    byl(-2:jj+2,-2:kk+2,-2:ii+2)
* local
       integer i, j, k
       integer ittime, itmax
       real    mindt
       real    maxdiff
       real    rhs, lhs
       real    dxx, dyy, dzz
       real    aa
       logical lquit
* parameter
       integer imethod
       parameter(imethod = 1) ! 1 = over-relax Jacobi; 2 & other  = normal time-relax
       real    difcrtrn
       parameter(difcrtrn = 1.0E-04)
       integer isidecnd
       parameter(isidecnd = 3) ! 1: d/ds, 2:divgrad, otherwise: d^2/ds^2 will be zero
       integer itopcond
       parameter(itopcond = 3) ! 1: Psi,  2:divgrad, otherwise: d^2/ds^2 will be zero

* give initial guess
!$omp parallel do private(i,k,j)
       do i = -2, ii + 2
         do k =  0, kk
           do j =  0, jj
             lapsi1(j,k,i) = - magmap(j,k) * (zzb(i) - zzb(ii))
           enddo
           do j = -2, -1
             lapsi1(j,k,i) = 0.0E+00
            enddo
           do j = jj + 1, jj + 2
             lapsi1(j,k,i) = 0.0E+00
            enddo
         enddo
         do j = -2, jj + 2
           lapsi1(j,k,i) = 0.0E+00
         enddo
       enddo
!$omp end parallel do

       if (imethod .EQ. 1) then
         mindt = 0.5E+00
         itmax = 1000000 ! tekitou
       else
         dxx = xx(1) - xx(0)
         dyy = yy(1) - yy(0)
         dzz = zz(1) - zz(0)
         mindt = min(dxx,dyy,dzz) ! MIND assumed here dxx etc is fixed constant ....
         mindt = mindt**2 * 0.5E+00 ! typical Courant cond.
     &         * 0.3333333E+00 ! for 3D
     &         * 0.95E+00 ! CFL number for parabolic eq.
         itmax = int(1.0E+00 / mindt)
     &         * 3 ! for safety.
     &         + 1 ! to guarantee to be positive
       endif

       write(*,*) ' Min Mesh_dt / Max Iteration count = ', mindt,itmax

* iterration process
       ittime = 0
       lquit = .false.
       do 100 while ((ittime .LE. itmax) .AND. (.NOT. lquit))
         ittime = ittime + 1

         
         if (imethod .EQ. 1) then

!$omp parallel do private(i,k,j,lhs,rhs,dxx,dyy,dzz)
           do i = 0, ii - 1
             dzz = 1.0E+00 / (zzb(i) - zzb(i-1))**2
             do k =  0, kk
               dyy = 1.0E+00 / (yy(k) - yy(k-1))**2
               do j =  0, jj
                 dxx = 1.0E+00 / (xx(j) - xx(j-1))**2

                 rhs = dzz * (lapsi1(j,k,i+1) + lapsi1(j,k,i-1)) ! Main part of  Poisson solver.
     &               + dxx * (lapsi1(j+1,k,i) + lapsi1(j-1,k,i))
     &               + dyy * (lapsi1(j,k+1,i) + lapsi1(j,k-1,i))
                 lhs = (dzz + dxx + dyy) * 2.0E+00
                 lapsi2(j,k,i) = rhs / lhs * mindt
     &                         + lapsi1(j,k,i) * (1.0E+00 - mindt)

               enddo
             enddo
           enddo
!$omp end parallel do

         else if (imethod .EQ. 2) then

!$omp parallel do private(i,k,j,lhs,rhs,dxx,dyy,dzz)
           do i = 0, ii - 1
             dzz = 1.0E+00 / (zzb(i) - zzb(i-1))**2
             do k =  0, kk
               dyy = 1.0E+00 / (yy(k) - yy(k-1))**2
               do j =  0, jj
                 dxx = 1.0E+00 / (xx(j) - xx(j-1))**2

                 rhs = dzz * (lapsi1(j,k,i+1)
     &                      - lapsi1(j,k,i  ) * 2.0E+00
     &                      + lapsi1(j,k,i-1))
     &               + dxx * (lapsi1(j+1,k,i)
     &                      - lapsi1(j  ,k,i) * 2.0E+00
     &                      + lapsi1(j-1,k,i))
     &               + dyy * (lapsi1(j,k+1,i)
     &                      - lapsi1(j,k  ,i) * 2.0E+00
     &                      + lapsi1(j,k-1,i))
                 lapsi2(j,k,i) = lapsi1(j,k,i) + rhs * mindt
               enddo
             enddo
           enddo
!$omp end parallel do

         else

!$omp parallel do private(i,k,j)
           do i = 0, ii
             dzz = 1.0E+00 / (zz(i+1) - zz(i))
             do k =  0, kk + 1
               dyy = 1.0E+00 / (yy(k+1) - yy(k))
               do j =  0, jj + 1
                 dxx = 1.0E+00 / (xx(j+1) - xx(j))

                 bzl(j,k,i) = ((lapsi1(j+1,k+1,i+1)
     &                         -lapsi1(j+1,k+1,i  ))
     &                        +(lapsi1(j  ,k+1,i+1)
     &                         -lapsi1(j  ,k+1,i  ))
     &                        +(lapsi1(j+1,k  ,i+1)
     &                         -lapsi1(j+1,k  ,i  ))
     &                        +(lapsi1(j  ,k  ,i+1)
     &                         -lapsi1(j  ,k  ,i  )))
     &                      * dzz * 0.25E+00
                 bxl(j,k,i) = ((lapsi1(j+1,k,  i+1)
     &                         -lapsi1(j  ,k,  i+1))
     &                        +(lapsi1(j+1,k  ,i+1)
     &                         -lapsi1(j  ,k  ,i+1))
     &                        +(lapsi1(j+1,k,  i  )
     &                         -lapsi1(j  ,k,  i  ))
     &                        +(lapsi1(j+1,k  ,i  )
     &                         -lapsi1(j  ,k  ,i  )))
     &                      * dxx * 0.25E+00
                 byl(j,k,i) = ((lapsi1(j+1,k+1,i+1)
     &                         -lapsi1(j+1,k  ,i+1))
     &                        +(lapsi1(j+1,k+1,i  )
     &                         -lapsi1(j+1,k  ,i  ))
     &                        +(lapsi1(j  ,k+1,i+1)
     &                         -lapsi1(j  ,k  ,i+1))
     &                        +(lapsi1(j  ,k+1,i  )
     &                         -lapsi1(j  ,k  ,i  )))
     &                      * dyy * 0.25E+00
               enddo
             enddo
           enddo
!$omp end parallel do

!$omp parallel do private(i,k,j,rhs)
           do i = 0, ii - 1
             dzz = 1.0E+00 / (zzb(i) - zzb(i-1))
             do k =  0, kk
               dyy = 1.0E+00 / (yyb(k) - yyb(k-1))
               do j =  0, jj
                 dxx = 1.0E+00 / (xxb(j) - xxb(j-1))

                 rhs =((bxl(j  ,k-1,i-1)+bxl(j  ,k  ,i-1)
     &                 +bxl(j  ,k-1,i  )+bxl(j  ,k  ,i  ))
     &                -(bxl(j-1,k-1,i-1)+bxl(j-1,k  ,i-1)
     &                 +bxl(j-1,k-1,i  )+bxl(j-1,k  ,i  )))
     &               * dxx * 0.25E+00
     &               +((byl(j-1,k  ,i-1)+byl(j-1,k  ,i  )
     &                 +byl(j  ,k  ,i-1)+byl(j  ,k  ,i  ))
     &                -(byl(j-1,k-1,i-1)+byl(j-1,k-1,i  )
     &                 +byl(j  ,k-1,i-1)+byl(j  ,k-1,i  )))
     &               * dyy * 0.25E+00
     &               +((bzl(j-1,k-1,i  )+bzl(j  ,k-1,i  )
     &                 +bzl(j-1,k  ,i  )+bzl(j  ,k  ,i  ))
     &                -(bzl(j-1,k-1,i-1)+bzl(j  ,k-1,i-1)
     &                 +bzl(j-1,k  ,i-1)+bzl(j  ,k  ,i-1)))
     &               * dzz * 0.25E+00
                 lapsi2(j,k,i) = lapsi1(j,k,i) + rhs * mindt

               enddo
             enddo
           enddo
!$omp end parallel do

        endif ! if imethod .... 

* top and bottom
!$omp parallel do private(i,k,j)
         do i = -2, ii + 2
           if (i .LT. 0) then ! cond. for i < 0
             do k =  0, kk
             do j =  0, jj
               lapsi2(j,k,i) = magmap(j,k) * (zzb(0) - zzb(i))
     &                       + lapsi2(j,k,0)
             enddo
             enddo
           endif
* 1: Psi,  2:divgrad, otherwise: d^2/ds^2 will be zero
           if (itopcond .EQ. 1) then
             if (i .GE. ii) then ! cond. for i >= ii
               do k =  0, kk
               do j =  0, jj
                 lapsi2(j,k,i) = 0.0E+00 ! source surface
               enddo
               enddo
             endif
           else if (itopcond .EQ. 2) then
             if (i .EQ. ii) then ! cond. for i >= ii
               do k =  0, kk
               do j =  0, jj
                 lapsi2(j,k,i) = lapsi2(j  ,k  ,i-1) * 6.0E+00 ! 2nd deriv = 0
     &                         - lapsi2(j  ,k  ,i-2)
     &                         - lapsi2(j-1,k  ,i-1)
     &                         - lapsi2(j+1,k  ,i-1)
     &                         - lapsi2(j  ,k-1,i-1)
     &                         - lapsi2(j  ,k+1,i-1) ! assuming dxx = dyy = dzz
               enddo
               enddo
             else if (i .GT. ii) then
               do k =  0, kk
               do j =  0, jj
                 lapsi2(j,k,i) = lapsi2(j,k,i-1) * 2.0E+00 ! 2nd deriv = 0
     &                         - lapsi2(j,k,i-2)
               enddo
               enddo
             endif
           else
             if (i .GE. ii) then ! cond. for i >= ii
               do k =  0, kk
               do j =  0, jj
                 lapsi2(j,k,i) = lapsi2(j,k,i-1) * 2.0E+00 ! 2nd deriv = 0
     &                         - lapsi2(j,k,i-2)
               enddo
               enddo
             endif
           endif
         enddo
!$omp end parallel do

* side(s) .... 1: d/ds, 2:divgrad, otherwise: d^2/ds^2 will be zero
         if (isidecnd .EQ. 1) then
!$omp parallel do private(i,k,j)
           do i = 0, ii
             do k =  0, kk
               lapsi2(  -1,k,i) = lapsi2( 0,k,i)
               lapsi2(jj+1,k,i) = lapsi2(jj,k,i)
             enddo
             do j =  0, jj
               lapsi2(j,  -1,i) = lapsi2(j, 0,i)
               lapsi2(j,kk+1,i) = lapsi2(j,kk,i)
             enddo
           enddo
!$omp end parallel do
         else if (isidecnd .EQ. 2) then
!$omp parallel do private(i,k,j)
           do i = 0, ii
             do k =  0, kk
               lapsi2(  -1,k,i) = lapsi2(   0,k  ,i  ) * 6.0E+00
     &                          - lapsi2(   1,k  ,i  )
     &                          - lapsi2(   0,k-1,i  )
     &                          - lapsi2(   0,k+1,i  )
     &                          - lapsi2(   0,k  ,i-1)
     &                          - lapsi2(   0,k  ,i+1)
               lapsi2(jj+1,k,i) = lapsi2(jj  ,k  ,i  ) * 6.0E+00
     &                          - lapsi2(jj-1,k  ,i  )
     &                          - lapsi2(jj  ,k-1,i  )
     &                          - lapsi2(jj  ,k+1,i  )
     &                          - lapsi2(jj  ,k  ,i-1)
     &                          - lapsi2(jj  ,k  ,i+1)
             enddo
             do k =  0, kk
               lapsi2(  -2,k,i) = lapsi2(  -1,k,i) * 2.0E+00
     &                          - lapsi2(   0,k,i)
               lapsi2(jj+2,k,i) = lapsi2(jj+1,k,i) * 2.0E+00
     &                          - lapsi2(jj  ,k,i)
             enddo
             do j =  0, jj
               lapsi2(j,  -1,i) = lapsi2(j  ,   0,i  ) * 6.0E+00
     &                          - lapsi2(j  ,   1,i  )
     &                          - lapsi2(j  ,   0,i-1)
     &                          - lapsi2(j  ,   0,i+1)
     &                          - lapsi2(j-1,   0,i  )
     &                          - lapsi2(j+1,   0,i  )
               lapsi2(j,kk+1,i) = lapsi2(j  ,kk  ,i  ) * 6.0E+00
     &                          - lapsi2(j  ,kk-1,i  )
     &                          - lapsi2(j  ,kk  ,i-1)
     &                          - lapsi2(j  ,kk  ,i+1)
     &                          - lapsi2(j-1,kk  ,i  )
     &                          - lapsi2(j+1,kk  ,i  )
             enddo
             do j =  0, jj
               lapsi2(j,  -2,i) = lapsi2(j,  -1,i) * 2.0E+00
     &                          - lapsi2(j,   0,i)
               lapsi2(j,kk+2,i) = lapsi2(j,kk+1,i) * 2.0E+00
     &                          - lapsi2(j,kk  ,i)
             enddo
           enddo
!$omp end parallel do
         else
!$omp parallel do private(i,k,j)
           do i = 0, ii
             do k =  0, kk
               lapsi2(  -1,k,i) = lapsi2(   0,k,i) * 2.0E+00
     &                          - lapsi2(   1,k,i)
               lapsi2(jj+1,k,i) = lapsi2(jj  ,k,i) * 2.0E+00
     &                          - lapsi2(jj-1,k,i)
             enddo
             do k =  0, kk
               lapsi2(  -2,k,i) = lapsi2(  -1,k,i) * 2.0E+00
     &                          - lapsi2(   0,k,i)
               lapsi2(jj+2,k,i) = lapsi2(jj+1,k,i) * 2.0E+00
     &                          - lapsi2(jj  ,k,i)
             enddo
             do j =  0, jj
               lapsi2(j,  -1,i) = lapsi2(j,   0,i) * 2.0E+00
     &                          - lapsi2(j,   1,i)
               lapsi2(j,kk+1,i) = lapsi2(j,kk  ,i) * 2.0E+00
     &                          - lapsi2(j,kk-1,i)
             enddo
             do j =  0, jj
               lapsi2(j,  -2,i) = lapsi2(j,  -1,i) * 2.0E+00
     &                          - lapsi2(j,   0,i)
               lapsi2(j,kk+2,i) = lapsi2(j,kk+1,i) * 2.0E+00
     &                          - lapsi2(j,kk  ,i)
             enddo
           enddo
         endif

* check convergence...
         maxdiff = -1.0E+10
!$omp parallel do private(i,k,j)
!$omp&            private(aa)
!$omp&            reduction(max:maxdiff)
         do i = 0, ii - 1
           do k =  0, kk
             do j =  0, jj
               aa = (lapsi2(j,k,i)    - lapsi1(j,k,i))**2
     &            / (lapsi2(j,k,i)**2 + lapsi1(j,k,i)**2 + 1.0E-10)
               if (aa .GT. maxdiff) maxdiff = aa
             enddo
           enddo
         enddo ! end i-loop
!$omp end parallel do

* copy Psi
!$omp parallel do private(i,k,j)
         do i = -2, ii + 2
           do k = -2, kk + 2
           do j = -2, jj + 2
             lapsi1(j,k,i) = lapsi2(j,k,i)
           enddo
           enddo
         enddo ! end i-loop
!$omp end parallel do

         if (mod(ittime,1000) .EQ. 0)
     &     write(*,'(i8,'' : '',e11.5)') ittime, maxdiff

         if (maxdiff .LT. difcrtrn) then
           write(*,*) ' converged.....'
           lquit = .true. ! frag on...
           write(*,'(i8,'' : '',e11.5)') ittime, maxdiff
         endif

 100   continue ! END of DO-WHILE-LOOP                                ----

       write(*,*) 'Iteration .. Done '
       write(*,'('' Ittime = '',i8)') ittime

       return
       end


* -----------------------------------------------------------------------
*
       subroutine get3dmag(ii,jj,kk,bz2,bx2,by2,
     &                     zz,zzb,xx,xxb,yy,yyb,lapsi1,lapsi2)
       implicit none
* number of grid for each three direction.
       integer ii, jj, kk
* variables
       real    bz2(-2:jj+2,-2:kk+2,-2:ii+2)
       real    bx2(-2:jj+2,-2:kk+2,-2:ii+2)
       real    by2(-2:jj+2,-2:kk+2,-2:ii+2)
* coordinate system
       real    zz(-2:ii+2),xx(-2:jj+2),yy(-2:kk+2)
       real    zzb(-2:ii+1),xxb(-2:jj+1),yyb(-2:kk+1)
* scalar for Laplace Eq. solver
       real    lapsi1(-2:jj+2,-2:kk+2,-2:ii+2)
       real    lapsi2(-2:jj+2,-2:kk+2,-2:ii+2)
* local
       integer i, j, k

!$omp parallel do private(i,k,j)
       do i = 0, ii
         do k = 0, kk
         do j = 0, jj
           bz2(j,k,i) =- (lapsi2(j  ,k,i  )-lapsi2(j  ,k,i-1))
     &                /  (zzb(i) - zzb(i-1))
           bx2(j,k,i) =-((lapsi2(j+1,k,i  )-lapsi2(j-1,k,i  ))
     &                  +(lapsi2(j+1,k,i-1)-lapsi2(j-1,k,i-1)))
     &                / (2.0E+00 * (xx(j+1) - xx(j-1)))
           by2(j,k,i) =-((lapsi2(j,k+1,i  )-lapsi2(j,k-1,i  ))
     &                  +(lapsi2(j,k+1,i-1)-lapsi2(j,k-1,i-1)))
     &                / (2.0E+00 * (yy(k+1) - yy(k-1)))
         enddo
         enddo
       enddo
!$omp end parallel do

       return
       end


* -----------------------------------------------------------------------
*
* Set Independent Variables, Coordinate System
*
* -----------------------------------------------------------------------
*
       subroutine setcrdnt(ii,jj,kk,zz,zzb,xx,xxb,yy,yyb)
       implicit none
*
       integer ii, jj, kk
       real    zz(-2:ii+2),xx(-2:jj+2),yy(-2:kk+2)
       real    zzb(-2:ii+1),xxb(-2:jj+1),yyb(-2:kk+1)
* local
       integer i, j, k
       real    xmax, ymax, zmax
       parameter(xmax = 1.0E+00, ymax = 1.0E+00, zmax = 1.0E+00)

!$omp parallel do private(i)
       do i = -2, ii + 2
         zz(i) = zmax / float(ii) * float(i)
       enddo
!$omp end parallel do
!$omp parallel do private(i)
       do i = -2, ii + 1
         zzb(i) = zmax / float(ii) *(float(i) + 0.5E+00)
       enddo
!$omp end parallel do
       do j = -2, jj + 2
         xx(j) = xmax / float(jj) * float(j)
       enddo
       do j = -2, jj + 1
         xxb(j) = xmax / float(jj) *(float(j) + 0.5E+00)
       enddo
       do k = -2, kk + 2
         yy(k) = ymax / float(kk) * float(k)
       enddo
       do k = -2, kk + 1
         yyb(k) = ymax / float(kk) *(float(k) + 0.5E+00)
       enddo

       return
       end


* -------------------------------------------------------------------------
*
       subroutine initall(ii,jj,kk,
     &                    bzl,bxl,byl,bz2,bx2,by2,
     &                    zz,zzb,xx,xxb,yy,yyb,
     &                    lapsi1,lapsi2)
       implicit none
* number of grid for each three direction.
       integer ii, jj, kk
* variables
       real    bz2(-2:jj+2,-2:kk+2,-2:ii+2)
       real    bx2(-2:jj+2,-2:kk+2,-2:ii+2)
       real    by2(-2:jj+2,-2:kk+2,-2:ii+2)
       real    bzl(-2:jj+2,-2:kk+2,-2:ii+2)
       real    bxl(-2:jj+2,-2:kk+2,-2:ii+2)
       real    byl(-2:jj+2,-2:kk+2,-2:ii+2)
* coordinate system
       real    zz(-2:ii+2),xx(-2:jj+2),yy(-2:kk+2)
       real    zzb(-2:ii+1),xxb(-2:jj+1),yyb(-2:kk+1)
* scalar for Laplace Eq. solver
       real    lapsi1(-2:jj+2,-2:kk+2,-2:ii+2)
       real    lapsi2(-2:jj+2,-2:kk+2,-2:ii+2)
* labor
       integer i, j, k

!$omp parallel do private(i,k,j)
       do i = -2, ii + 2
         do k = -2, kk + 2
         do j = -2, jj + 2
           lapsi1(j,k,i) = 0.0E+00
           lapsi2(j,k,i) = 0.0E+00
           bzl(j,k,i) = 0.0E+00
           bxl(j,k,i) = 0.0E+00
           byl(j,k,i) = 0.0E+00
           bz2(j,k,i) = 0.0E+00
           bx2(j,k,i) = 0.0E+00
           by2(j,k,i) = 0.0E+00
         enddo
         enddo
       enddo
!$omp end parallel do
!$omp parallel do private(i)
       do i = -2, ii + 2
         zz(i) = 0.0E+00
       enddo
!$omp end parallel do
!$omp parallel do private(i)
       do i = -2, ii + 1
         zzb(i) = 0.0E+00
       enddo
!$omp end parallel do
       do j = -2, jj + 2
         xx(j) = 0.0E+00
       enddo
       do j = -2, jj + 1
         xxb(j) = 0.0E+00
       enddo
       do k = -2, kk + 2
         yy(k) = 0.0E+00
       enddo
       do k = -2, kk + 1
         yyb(k) = 0.0E+00
       enddo

       return
       end


*--------------------------------------------------------------------
*
* END of This file
*
*--------------------------------------------------------------------
*
