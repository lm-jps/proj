!***********************************************************************
! subroutine tile                                                      *
!    Starting from a regular grid in a Mollweide equal area projection,*
! construct a set of rectangular tiles which covers the disk.          *
!***********************************************************************

subroutine tile(ntx,nty,nap)

!-----------------------------------------------------------------------
   use sizes
   use constant
   use disk_center
   use ephemeris
   use mgram_data
   use maskvec
   use pot_field
   use pad
   use pix_size
   use spherical_position
   use trnsfrm

   implicit none

   integer,intent(in) :: ntx,nty,nap
   integer :: i,j,k,l,ntile,imin,imax,jmin,jmax,nxp,nyp,ishift,jshift
   real :: theta,phi,x1,y1,y2,r1,xmin,xmax,ymin,ymax,bdotb
   real,dimension(:),allocatable :: x,y
   real,dimension(:,:),allocatable :: xs,ys,Bzwin
!-----------------------------------------------------------------------
!open(9,file='potential.dat')
! --> Allocate memory.
   allocate(x(ntx),y(nty))
   allocate(xs(ntx,nty),ys(ntx,nty))
   allocate(Bzwin(nx,ny),Bpix(nx,ny),Bpiy(nx,ny),dBpzdz(nx,ny))
   allocate(tmask(nx,ny))
   tmask=0

! --> Apply Tukey windowing function to line of sight component of the 
! --> field to reduce ringing caused by field on the limb.
   do j=1,ny
      y1=cosp(j)**2
      y2=1.-y1
      if(cosp(j).gt.-10.) then
         do i=1,nx
            if(sint(i,j).gt.-10.) then
               x1=y2*sint(i,j)**2
               r1=sqrt(y1+x1)
               if(r1.lt.0.95) then
                  Bzwin(i,j)=Bz(i,j)
               else
                  Bzwin(i,j)=0.5*Bz(i,j)*(1.-cos((1.-r1)*20.*pi))
               endif
            else
               Bzwin(i,j)=0.
            endif
         enddo
      else
         do i=1,nx
            Bzwin(i,j)=0.
         enddo
      endif
   enddo

! --> Construct regular grid in projection.
   do i=1,ntx
      x(i)=1.5*(2.*float(i-1)/float(ntx-1)-1.)
   enddo
   do j=1,nty
      y(j)=1.5*(2.*float(j-1)/float(nty-1)-1.)
   enddo

! --> Map projected points back to the sphere. Flag points which do not
! --> fall on the front side of the sphere.
   do i=1,ntx
      x1=x(i)
      do j=1,nty
         y1=y(j)
         if(x1**2+4.*y1**2.lt.8.) then
            call mollweide(x1,y1,theta,phi)
            if(abs(theta).lt.pi/2.) then
               xs(i,j)=cos(phi)*sin(theta)
               ys(i,j)=sin(phi)
            else
               xs(i,j)=-2.
               ys(i,j)=-2.
            endif
         else
            xs(i,j)=-3.
            ys(i,j)=-3.
         endif
      enddo
   enddo

! --> Pixel index of disk center, relative to 1,1 at lower left corner
! --> of FOV.
   xcen=xcen/dxi
   ycen=ycen/dyi

! --> Construct a set of rectangular boxes, each of which contains all
! --> the pixels within one of these patches. Keep track of the number
! --> of tiles.
   ntile=0
   do i=1,ntx-1
      do j=1,nty-1
         if(xs(i,j).gt.-1.and.xs(i+1,j).gt.-1.and.xs(i,j+1).gt.-1.and.xs(i+1,j+1).gt.-1.and.&
            ys(i,j).gt.-1.and.ys(i+1,j).gt.-1.and.ys(i,j+1).gt.-1.and.ys(i+1,j+1).gt.-1) then
            ntile=ntile+1
            xmin=min(xs(i,j),xs(i+1,j),xs(i,j+1),xs(i+1,j+1))
            xmax=max(xs(i,j),xs(i+1,j),xs(i,j+1),xs(i+1,j+1))
            ymin=min(ys(i,j),ys(i+1,j),ys(i,j+1),ys(i+1,j+1))
            ymax=max(ys(i,j),ys(i+1,j),ys(i,j+1),ys(i+1,j+1))

! --> "Center" point of box.
            x1=0.25*(xs(i,j)+xs(i+1,j)+xs(i,j+1)+xs(i+1,j+1))
            y1=0.5*(ys(i,j)+ys(i,j+1))

! --> Get coordinate transform information for this point.
            theta=atan(x1/sqrt(1.-x1**2-y1**2))
            phi=asin(y1)

            call transform(theta,phi)

! --> Determine which pixels are within the box.
            imin=nint(xmin*radius/dxi+xcen)
            imax=nint(xmax*radius/dxi+xcen)
            jmin=nint(ymin*radius/dyi+ycen)
            jmax=nint(ymax*radius/dyi+ycen)

! --> If the tile has some pixels in the FOV, then proceed.
            if(imax.gt.1.and.imin.lt.nx.and.jmax.gt.1.and.jmin.lt.ny) then

! --> Ensure the tile, when padded, stays within the FOV.
               if(imin.le.npad) imin=npad+1
               if(imax.ge.nx-npad) imax=nx-npad
               if(jmin.le.npad) jmin=npad+1
               if(jmax.ge.ny-npad) jmax=ny-npad

               nxp=imax-imin+2*npad+1
               nyp=jmax-jmin+2*npad+1

! --> Extract line of sight field for this tile, including applying a Tukey 
! --> windowing to limit the amount of ringing.
               allocate(blpad(nxp,nyp))

               do k=1,nap-1
                  do l=1,nap-1
                     blpad(k,l)=Bzwin(imin+k-1-npad,jmin+l-1-npad)*0.25&
                        *(1.-cos(float(k-1)/float(nap-1)*pi))*(1.-cos(float(l-1)/float(nap-1)*pi))
                  enddo
               enddo

               do k=1,nap-1
                  do l=nap,nyp-nap
                     blpad(k,l)=Bzwin(imin+k-1-npad,jmin+l-1-npad)*0.5&
                        *(1.-cos(float(k-1)/float(nap-1)*pi))
                  enddo
               enddo

               do k=1,nap-1
                  do l=nyp-nap+1,nyp
                     blpad(k,l)=Bzwin(imin+k-1-npad,jmin+l-1-npad)*0.25&
                        *(1.-cos(float(k-1)/float(nap-1)*pi))*(1.-cos(float(nyp-l)/float(nap-1)*pi))
                  enddo
               enddo

               do k=nap,nxp-nap
                  do l=1,nap-1
                     blpad(k,l)=Bzwin(imin+k-1-npad,jmin+l-1-npad)*0.5&
                        *(1.-cos(float(l-1)/float(nap-1)*pi))
                  enddo
               enddo

               do k=nap,nxp-nap
                  do l=nap,nyp-nap
                     blpad(k,l)=Bzwin(imin+k-1-npad,jmin+l-1-npad)
                  enddo
               enddo

               do k=nap,nxp-nap
                  do l=nyp-nap+1,nyp
                     blpad(k,l)=Bzwin(imin+k-1-npad,jmin+l-1-npad)*0.5&
                        *(1.-cos(float(nyp-l)/float(nap-1)*pi))
                  enddo
               enddo

               do k=nxp-nap+1,nxp
                  do l=1,nap-1
                     blpad(k,l)=Bzwin(imin+k-1-npad,jmin+l-1-npad)*0.25&
                        *(1.-cos(float(nxp-k)/float(nap-1)*pi))*(1.-cos(float(l-1)/float(nap-1)*pi))
                  enddo
               enddo

               do k=nxp-nap+1,nxp
                  do l=nap,nyp-nap
                     blpad(k,l)=Bzwin(imin+k-1-npad,jmin+l-1-npad)*0.5&
                        *(1.-cos(float(nxp-k)/float(nap-1)*pi))
                  enddo
               enddo

               do k=nxp-nap+1,nxp
                  do l=nyp-nap+1,nyp
                     blpad(k,l)=Bzwin(imin+k-1-npad,jmin+l-1-npad)*0.25&
                        *(1.-cos(float(nxp-k)/float(nap-1)*pi))*(1.-cos(float(nyp-l)/float(nap-1)*pi))
                  enddo
               enddo

! --> Calculate potential field in this tile.
               call potential(nxp,nyp)

! --> Fill in potential field (derivatives) for this tile.
               do k=1,nxp
                  do l=1,nyp
                     ishift=imin+k-1-npad
                     jshift=jmin+l-1-npad
                     x1=(float(ishift)-xcen)*dxi/radius
                     y1=(float(jshift)-ycen)*dyi/radius
                     theta=atan(x1/sqrt(1.-x1**2-y1**2))
                     phi=asin(y1)

                     call revmollweide(theta,phi,x1,y1)

                     if(x1.ge.x(i).and.x1.lt.x(i+1).and.y1.ge.y(j).and.y1.lt.y(j+1)) then
                        dBzdz(ishift,jshift)=radius*dBpzdz(k,l)
                        tmask(ishift,jshift)=1

! --> Potential field acute angle ambiguity resolution.

! --> Compute dot product of observed transverse field with potential field.
                        bdotb=Bpix(ishift,jshift)*Bx(ishift,jshift)+Bpiy(ishift,jshift)*By(ishift,jshift)

! --> Flip direction of transverse field if dot product is negative.
                        if(bdotb.lt.0.) then
                           Bx(ishift,jshift)=-Bx(ishift,jshift)
                           By(ishift,jshift)=-By(ishift,jshift)
                        endif
                     endif
                  enddo
               enddo
            endif
         endif
      enddo
   enddo
!do j=1,ny
!   write(9,*) (dBzdz(i,j),i=1,nx)
!enddo

! --> Deallocate memory.
   deallocate(x,y,xs,ys,Bzwin,dBpzdz,Bpix,Bpiy)

end subroutine tile
