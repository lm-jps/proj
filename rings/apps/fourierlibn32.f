c &&&&&&&&&&&&&&&&&&&&&
      subroutine fftrc (rldat,nrlpts,signex,hctrn,nhcpts,work,lwrk,ierr)
c
c package fft
c
c
c purpose                fast fourier transforms for data of arbitrary
c                        length.  the file fft contains three routines
c                        to handle various forms of input data as
c                        tabulated below:
c                          routine name  input form     output form
c                          ------------  ----------     -----------
c                             fftrc         real        half complex
c                             fftcr      half complex      real
c                             fftcc        complex        complex
c                        half complex refers here to the first n/2+1
c                        complex values of a conjugate symmetric array
c                        of length n.
c
c access cards           *fortran,s=ulib,n=fft
c
c
c special conditions     the efficiency of these routines is greatly
c                        affected by the number of points to be
c                        transformed.  if n is the length of the
c                        transform (nrlpts in fftrc or fftcr; ncpts in
c                        fftcc) the execution time is roughly
c                        proportional to n*sumpf where sumpf is the sum
c                        of the prime factors of n.  clearly, numbers
c                        with large prime factors should be avoided.  if
c                        n is a power of two, the package fftpow2
c                        provides even greater efficiency.
c
c     subroutine fftrc(rldat,nrlpts,signex,hctrn,nhcpts,work,lwrk,ierr)
c
c dimension of           real rldat (nrlpts)
c arguments              complex hctrn (nhcpts)
c                          where nhcpts = nrlpts/2+1
c                        real work (lwrk)
c
c latest revision        january 1978
c
c purpose                fourier transform (real to half complex).  (in
c                        case signex = +1., the computation performed by
c                        this routine is sometimes called a forward
c                        transform or fourier analysis.)  the discrete
c                        fourier transform of an array rldat, containing
c                        nrlpts real values, is a set cf of nrlpts
c                        complex values satisfying the conjugate
c                        symmetry relation cf(nrlpts+2-k) = conjg(cf(k))
c                        (k = 2,nrlpts).  due to this symmetry relation,
c                        it is only necessary to compute the first
c                        nhcpts = nrlpts/2+1 complex values, and these
c                        are returned in the complex array hctrn.
c
c usage                  call fftrc (rldat,nrlpts,signex,hctrn,nhcpts,
c                                    work,lwrk,ierr)
c                          the original values of rldat may be
c                          regenerated from hctrn by first dividing all
c                          values of hctrn by nrlpts and then calling
c                          fftcr (hctrn,nhcpts,-signex,rldat,nrlpts,
c                                 work,lwrk,ierr).
c
c note                   for these comments we assume
c                        signex = +1. or -1. and define cex(x)
c                        (for all real x) to be the complex exponential
c                        of signex*2*pi*i*x/nrlpts  where pi = 3.14...
c                        and i = sqrt(-1.).
c
c arguments
c
c on input               rldat
c                          a real array containing the nrlpts data
c                          values to be transformed.  it may be
c                          equivalenced to hctrn or work if desired.
c                          the dimension is assumed to be rldat(nrlpts).
c
c                        nrlpts
c                          the number of real data values to be
c                          transformed.  for the greatest efficiency, it
c                          should be a product of small primes.
c
c                        signex
c                          a variable whose sign determines the sign of
c                          the argument of the complex exponential used
c                          in the transform computations.  for
c                          convenience, we assume in these comments that
c                          signex is +1. or -1., but the routine in fact
c                          only uses the sign of its value.
c
c                        nhcpts
c                          the number of complex values to be returned
c                          as the half complex transform result in array
c                          hctrn.  it must be = nrlpts/2+1 or a fatal
c                          error is flagged.
c                          note:  nhcpts is not an output parameter.  it
c                                 must be set to nrlpts/2+1 by the user,
c                                 and nhcpts complex locations must be
c                                 provided for the output array hctrn.
c
c                        work
c                          a workspace of length lwrk (.ge. 4*nrlpts)
c                          for use by the routine.  either or both of
c                          the arrays rldat and hctrn may be
c                          equivalenced to work to reduce storage.  if
c                          rldat is so equivalenced, it will be
c                          destroyed.
c
c                        lwrk
c                          the length of array work.  it must be
c                          .ge. 4*nrlpts or a fatal error is flagged.
c
c on output              hctrn
c                          the complex array containing essentially the
c                          first half of the conjugate symmetric
c                          transform result.  hctrn(k)
c                          (for k = 1,nhcpts) is the complex value
c                          defined by
c                                       nrlpts
c                            hctrn(k) =  sum   rldat(j)*cex((j-1)*(k-1)
c                                       j = 1
c                          these values are also referred to as the
c                          fourier coefficients.
c                            the dimension is assumed to be
c                            complex hctrn (nhcpts) which requires
c                            2*nhcpts core locations.
c                          note:  hctrn may be equivalenced to rldat if
c                                 desired, but they are not identical in
c                                 size.  the number of core locations
c                                 required for hctrn is nrlpts+1 for odd
c                                 nrlpts and nrlpts+2 for even nrlpts.
c
c                        work
c                          the workspace containing intermediate
c                          results.
c
c                        ierr
c                          an error flag with the following meanings:
c                          =   0  no error.
c                          = 101  nrlpts is less than 1.
c                          = 102  nhcpts is not nrlpts/2+1.
c                          = 103  insufficient workspace has been
c                                 provided; lwrk is less than 4*nrlpts.
c
c common blocks          none
c
c i/o                    none
c
c precision              single
c
c required ulib          none
c routines
c
c specialist             dave fulker, ncar, boulder, colorado  80303
c
c language               fortran
c
c history                standardized march 1974 by dave fulker at ncar.
c
      real            rldat(nrlpts)          ,work(lwrk)
      complex         hctrn(nhcpts)
      dimension       nscrt(1)
c the following call is for gathering statistics on library use at ncar
      ierr = 0
      if (nrlpts .lt. 2) go to 103
      if (nhcpts .ne. nrlpts/2+1) go to 105
      if (lwrk .lt. 4*nrlpts) go to 106
      do 101 j=1,nrlpts
         work(2*j-1) = rldat(j)
         work(2*j) = 0.
  101 continue
      isign = signex
      nscrt(1) = nrlpts
      call fourt (work,nscrt,1,isign,0,work(2*nrlpts+1))
      do 102 k=1,nhcpts
         hctrn(k) = cmplx(work(2*k-1),work(2*k))
  102 continue
      return
  103 if (nrlpts .lt. 1) go to 104
      if (nhcpts .ne. 1) go to 105
      hctrn(1) = cmplx(rldat(1),0.)
      return
  104 ierr = 101
      write(2,1101)ierr
 1101 format(' ierr = ',i6,' fftrc nrlpts is .lt. 1')
      return
  105 ierr = 102
      write(2,1102)ierr
 1102 format(' ierr = ',i6,' fftrc nhcpts is not nrlpts/2+1')
      return
  106 ierr = 103
      write(2,1103)ierr
 1103 format(' ierr = ',i6,' fftrc insufficient workspace - lwrk is
     $ .lt. 4*nrlpts')
      return
      end
      subroutine fftcr (hcdat,nhcpts,signex,rltrn,nrlpts,work,lwrk,ierr)
c
c
c dimension of           complex hcdat (nhcpts)
c arguments              real rltrn (nrlpts)   where nrlpts/2+1 = nhcpts
c                        real work (lwrk)
c
c latest revision        january 1978
c
c purpose                fourier transform (half complex to real).  (in
c                        case signex = -1., the computation performed by
c                        this routine is sometimes called a backward
c                        transform or fourier synthesis.) hcdat is
c                        assumed to be the first nhcpts = nrlpts/2+1
c                        complex values of a set cf containing nrlpts
c                        complex values which satisfy the conjugate
c                        symmetry relation cf(nrlpts+2-k) = conjg(cf(k))
c                        (k = 2,nrlpts).  in addition, hctrn(1) is
c                        assumed to be real.  the discrete fourier
c                        transform of such a conjugate symmetric array
c                        is a set of nrlpts real values which are
c                        returned in the real array rltrn.
c
c usage                  call fftcr (hcdat,nhcpts,signex,rltrn,nrlpts,
c                                    work,lwrk,ierr)
c                          the original values of hcdat may be
c                          regenerated from rltrn by first dividing all
c                          values of rltrn by nrlpts and then calling
c                          fftrc (rltrn,nrlpts,-signex,hcdat,nhcpts,
c                                 work,lwrk,ierr).
c
c note                   for these comments we assume
c                        signex = +1. or -1. and define cex(x)
c                        (for all real x) to be the complex exponential
c                        of signex*2*pi*i*x/nrlpts  where pi = 3.14...
c                        and i = sqrt(-1.).
c
c arguments
c
c on input               hcdat
c                          a complex array containing nhcpts complex
c                          values which comprise essentially the first
c                          half of the conjugate symmetric data to be
c                          transformed.  it may be equivalenced to rltrn
c                          or work if desired, but it is not of
c                          identical size.
c                            the dimension is assumed to be
c                            complex hcdat(nhcpts) which requires
c                            2*hncpts core locations.
c
c                        nhcpts
c                          the number of complex values entered in
c                          hcdat.
c
c                        signex
c                          a variable whose sign determines the sign of
c                          the argument of the complex exponential used
c                          in the transform computations.  for
c                          convenience we assume in these comments that
c                          signex is +1. or -1., but the routine in fact
c                          only uses the sign of its value.
c
c                        nrlpts
c                          the number of real transform values to be
c                          returned.  nhcpts must be nrlpts/2+1 or a
c                          fatal error is flagged.
c                          note:  nrlpts is not an output parameter.  it
c                                 must satisfy nrlpts/2+1 = nhcpts, and
c                                 nrlpts real locations must be provided
c                                 for the output array rltrn.
c
c                        work
c                          a workspace of length lwrk (.ge. 4*nrlpts)
c                          for use by the routine.  either or both of
c                          the arrays hcdat or rltrn may be equivalenced
c                          to work to reduce storage.  if hcdat is so
c                          equivalenced, it will be destroyed.
c
c                        lwrk
c                          the length of array work.  it must be
c                          .ge. 4*nrlpts or a fatal error is flagged.
c
c on output              rltrn
c                          the real array containing the transform
c                          result.  rltrn(j) (for j = 1,nrlpts) is the
c                          real value defined by
c                                       nrlpts
c                            rltrn(j) =  sum   cf(k)*cex(j-1)*(k-1))
c                                       k = 1
c                          where
c                                    ( hcdat(k)/nrlpts
c                                    (for 1 .le. k .le. nhcpts
c                            cf(k) = (
c                                    (conjg(cf(nrlpts+2-k))  otherwise.
c
c                          the dimension is assumed to be rltrn(nrlpts).
c                            note:  rltrn may be equivalenced to hcdat
c                                   if desired.
c
c                        work
c                          the workspace containing intermediate
c                          results.
c
c                        ierror
c                          an error flag with the following meanings:
c                            0  no error.
c                          101  nrlpts is less than 1.
c                          102  nhcpts is not nrlpts/2+1.
c                          103  insufficient workspace has been
c                               provided:  lwrk is less than 4*nrlpts.
c
c common blocks          none
c
c i/o                    none
c
c precision              single
c
c required ulib          none
c routines
c
c specialist             dave fulker, ncar, boulder, colorado  80303
c
c language               fortran
c
c history                standardized march 1974 by dave fulker at ncar.
c
c
      complex         hcdat(nhcpts)
      real            rltrn(nrlpts)          ,work(lwrk)
      dimension       nscrt(1)
c the following call is for gathering statistics on library use at ncar
      ierr = 0
      if (nrlpts .lt. 2) go to 103
      if (nhcpts .ne. nrlpts/2+1) go to 105
      if (lwrk .lt. 4*nrlpts) go to 106
      nc = 2*nrlpts+4
      work(1) = real(hcdat(1))
      work(2) = 0.
      do 101 k=2,nhcpts
         work(2*k-1) = real(hcdat(k))
      nc2 = nc-2*k-1
      work(nc2) = real(hcdat(k))
         work(2*k) = aimag(hcdat(k))
      nc2 = nc-2*k
      work(nc2) = -aimag(hcdat(k))
  101 continue
      isign = signex
      nscrt(1) = nrlpts
      call fourt (work,nscrt,1,isign,1,work(2*nrlpts+1))
      do 102 j=1,nrlpts
         rltrn(j) = work(2*j-1)
  102 continue
      return
  103 if (nrlpts .lt. 1) go to 104
      if (nhcpts .ne. 1) go to 105
      rltrn(1) = real(hcdat(1))
      return
  104 ierr = 101
      write(2,1101)ierr
 1101 format(' ierr = ',i6,' fftcr nrlpts is .lt. 1')
      return
  105 ierr = 102
      write(2,1102)ierr
 1102 format(' ierr = ',i6,' fftcr nhcpts is not nrlpts/2+1')
      return
  106 ierr = 103
      write(2,1103)ierr
 1103 format(' ierr = ',i6,' fftcr insufficient workspace - lwrk is
     $ .lt. 4*nrlpts')
      return
      end
      subroutine fftcc (cdata,ncpts,signex,ctran,work,ierr)
c
c dimension of           complex cdata (ncpts),ctran (ncpts)
c arguments              real work(2*ncpts)
c
c latest revision        january 1978
c
c purpose                fourier transform (complex to complex).  (in
c                        case signex = +1., the computation performed by
c                        this routine is sometimes called a forward
c                        transform or fourier analysis.  for
c                        signex = -1., it is called a backward transform
c                        or fourier synthesis.)  the discrete fourier
c                        transform of an array cdata, containing ncpts
c                        complex values, is a set of ncpts complex
c                        values which are returned in the complex array
c                        ctran.
c
c usage                  call fftcc (cdata,ncpts,signex,ctran,work,ierr)
c                          the original values of cdata may be
c                          regenerated from ctran by first dividing all
c                          values of ctran by ncpts and then calling
c                          fftcc (ctran,ncpts,-signex,cdata,work,ierr).
c
c note                   for these comments we assume
c                        signex = +1. or -1. and define cex(x)
c                        (for all real x) to be the complex exponential
c                        of signex*2*pi*i*x/ncpts   where pi = 3,14...
c                        and i = sqrt(-1.).
c
c arguments
c
c on input               cdata
c                          a complex array containing the ncpts complex
c                          data values to be transformed.  it may be
c                          equivalenced to work if desired.  the
c                          dimension is assumed to be
c                          complex cdata(ncpts).
c
c                        ncpts
c                          the number of complex data values to be
c                          transformed.  it must be a positive power of
c                          2 or a fatal error is flagged.
c
c                        signex
c                          a variable whose sign determines the sign of
c                          the argument of the complex exponential used
c                          in the transform computations.  for
c                          convenience, we assume in these comments that
c                          signex is +1. or -1., but the routine in fact
c                          only uses the sign of its value.
c
c                        work
c                          a workspace of length 2*ncpts for use by the
c                          routine.  it may be equivalenced to cdata in
c                          which case cdata will be destroyed.
c
c on output              ctran
c                          the complex array in which the ncpts complex
c                          transform results are returned.  ctran(k)
c                          (for k = 1,ncpts) is the complex value
c                          defined by
c                                       ncpts
c                            ctran(k) =  sum  cdata(j)*cex(j-1)*(k-1))
c                                       j = 1
c                          these values are also referred to as the
c                          fourier coefficients.
c                            the dimension is assumed to be
c                            complex ctran(ncpts).
c
c                        work
c                          the workspace containing intermediate
c                          results.
c
c                        ierr
c                          an error flag with the following meanings:
c                            0  no error.
c                          101  ncpts is less than 1.
c
c common blocks          none
c
c i/o                    none
c
c precision              single
c
c required ulib          none
c routines
c
c specialist             dave fulker, ncar, boulder, colorado  80303
c
c language               fortran
c
c history                standardized march 1974 by dave fulker at ncar.
c
c
c
c
      complex         cdata(ncpts)           ,ctran(ncpts)           ,
     1                work(ncpts)
      dimension       nscrt(1)
c the following call is for gathering statistics on library use at ncar
      ierr = 0
      if (ncpts .lt. 1) go to 102
      do 101 j=1,ncpts
         ctran(j) = cdata(j)
  101 continue
      isign = signex
      nscrt(1) = ncpts
      call fourt (ctran,nscrt,1,isign,1,work)
      return
  102 ierr = 101
      write(2,1101)ierr
 1101 format(' ierr = ',i6,' fftcc ncpts is .lt. 1')
      return
      end
      subroutine fourt (data,nn,ndim,isign,iform,work)
c
c     the cooley-tukey fast fourier transform in usasi basic fortran
c
c     transform(j1,j2,,,,) = sum(data(i1,i2,,,,)*w1**((i2-1)*(j2-1))
c                                 *w2**((i2-1)*(j2-1))*,,,),
c     where i1 and j1 run from 1 to nn(1) and w1=exp(isign*2*pi=
c     sqrt(-1)/nn(1)), etc.  there is no limit on the dimensionality
c     (number of subscripts) of the data array.  if an inverse
c     transform (isign=+1) is performed upon an array of transformed
c     (isign=-1) data, the original data will reappear.
c     multiplied by nn(1)*nn(2)*,,,  the array of input data must be
c     in complex format.  however, if all imaginary parts are zero (i.e.
c     the data are disguised real) running time is cut up to forty per-
c     cent.  (for fastest transform of real data, nn(1) should be even.)
c     the transform values are always complex and are returned in the
c     original array of data, replacing the input data.  the length
c     of each dimension of the data array may be any integer.  the
c     program runs faster on composite integers than on primes, and is
c     particularly fast on numbers rich in factors of two.
c
c     timing is in fact given by the following formula.  let ntot be the
c     total number of points (real or complex) in the data array, that
c     is, ntot=nn(1)*nn(2)*...  decompose ntot into its prime factors,
c     such as 2**k2 * 3**k3 * 5**k5 * ...  let sum2 be the sum of all
c     the factors of two in ntot, that is, sum2 = 2*k2.  let sumf be
c     the sum of all other factors of ntot, that is, sumf = 3*k3*5*k5*..
c     the time taken by a multidimensional transform on these ntot data
c     is t = t0 + ntot*(t1+t2*sum2+t3*sumf).  on the cdc 3300 (floating
c     point add time = six microseconds), t = 3000 + ntot*(600+40*sum2+
c     175*sumf) microseconds on complex data.
c
c     implementation of the definition by summation will run in a time
c     proportional to ntot*(nn(1)+nn(2)+...).  for highly composite ntot
c     the savings offered by this program can be dramatic.  a one-dimen-
c     sional array 4000 in length will be transformed in 4000*(600+
c     40*(2+2+2+2+2)+175*(5+5+5)) = 14.5 seconds versus about 4000*
c     4000*175 = 2800 seconds for the straightforward technique.
c
c     the fast fourier transform places three restrictions upon the
c     data.
c     1.  the number of input data and the number of transform values
c     must be the same.
c     2.  both the input data and the transform values must represent
c     equispaced points in their respective domains of time and
c     frequency.  calling these spacings deltat and deltaf, it must be
c     true that deltaf=2*pi/(nn(i)*deltat).  of course, deltat need not
c     be the same for every dimension.
c     3.  conceptually at least, the input data and the transform output
c     represent single cycles of periodic functions.
c
c     the calling sequence is--
c     call fourt(data,nn,ndim,isign,iform,work)
c
c     data is the array used to hold the real and imaginary parts
c     of the data on input and the transform values on output.  it
c     is a multidimensional floating point array, with the real and
c     imaginary parts of a datum stored immediately adjacent in storage
c     (such as fortran iv places them).  normal fortran ordering is
c     expected, the first subscript changing fastest.  the dimensions
c     are given in the integer array nn, of length ndim.  isign is -1
c     to indicate a forward transform (exponential sign is -) and +1
c     for an inverse transform (sign is +).  iform is +1 if the data are
c     complex, 0 if the data are real.  if it is 0, the imaginary
c     parts of the data must be set to zero.  as explained above, the
c     transform values are always complex and are stored in array data.
c     work is an array used for working storage.  it is floating point
c     real, one dimensional of length equal to twice the largest array
c     dimension nn(i) that is not a power of two.  if all nn(i) are
c     powers of two, it is not needed and may be replaced by zero in the
c     calling sequence.  thus, for a one-dimensional array, nn(1) odd,
c     work occupies as many storage locations as data.  if supplied,
c     work must not be the same array as data.  all subscripts of all
c     arrays begin at one.
c
c     example 1.  three-dimensional forward fourier transform of a
c     complex array dimensioned 32 by 25 by 13 in fortran iv.
c     dimension data(32,25,13),work(50),nn(3)
c     complex data
c     data nn/32,25,13/
c     do 1 i=1,32
c     do 1 j=1,25
c     do 1 k=1,13
c  1  data(i,j,k)=complex value
c     call fourt(data,nn,3,-1,1,work)
c
c     example 2.  one-dimensional forward transform of a real array of
c     length 64 in fortran ii,
c     dimension data(2,64)
c     do 2 i=1,64
c     data(1,i)=real part
c  2  data(2,i)=0.
c     call fourt(data,64,1,-1,0,0)
c
c     there are no error messages or error halts in this program.  the
c     program returns immediately if ndim or any nn(i) is less than one.
c
c     program by norman brenner from the basic program by charles
c     rader,  june 1967.  the idea for the digit reversal was
c     suggested by ralph alter.
c
c     this is the fastest and most versatile version of the fft known
c     to the author.  a program called four2 is available that also
c     performs the fast fourier transform and is written in usasi basic
c     fortran.  it is about one third as long and restricts the
c     dimensions of the input array (which must be complex) to be powers
c     of two.  another program, called four1, is one tenth as long and
c     runs two thirds as fast on a one-dimensional complex array whose
c     length is a power of two.
c
c     reference--
c     ieee audio transactions (june 1967), special issue on the fft.
c
      dimension       data(1)    ,nn(1)      ,ifact(32)  ,work(1)
      data np0/0/,nprev/0/
      data twopi/6.2831853071796/,rthlf/0.70710678118655/
      if (ndim-1) 232,101,101
  101 ntot = 2
      do 103 idim=1,ndim
         if (nn(idim)) 232,232,102
  102    ntot = ntot*nn(idim)
  103 continue
c
c     main loop for each dimension
c
      np1 = 2
      do 231 idim=1,ndim
         n = nn(idim)
         np2 = np1*n
         if (n-1) 232,230,104
c
c     is n a power of two and if not, what are its factors
c
  104    m = n
         ntwo = np1
         if = 1
         idiv = 2
  105    iquot = m/idiv
         irem = m-idiv*iquot
         if (iquot-idiv) 113,106,106
  106    if (irem) 108,107,108
  107    ntwo = ntwo+ntwo
         ifact(if) = idiv
         if = if+1
         m = iquot
         go to 105
  108    idiv = 3
         inon2 = if
  109    iquot = m/idiv
         irem = m-idiv*iquot
         if (iquot-idiv) 115,110,110
  110    if (irem) 112,111,112
  111    ifact(if) = idiv
         if = if+1
         m = iquot
         go to 109
  112    idiv = idiv+2
         go to 109
  113    inon2 = if
         if (irem) 115,114,115
  114    ntwo = ntwo+ntwo
         go to 116
  115    ifact(if) = m
c
c     separate four cases--
c        1. complex transform or real transform for the 4th, 9th,etc.
c           dimensions.
c        2. real transform for the 2nd or 3rd dimension.  method--
c           transform half the data, supplying the other half by con-
c           jugate symmetry.
c        3. real transform for the 1st dimension, n odd.  method--
c           set the imaginary parts to zero.
c        4. real transform for the 1st dimension, n even.  method--
c           transform a complex array of length n/2 whose real parts
c           are the even numbered real values and whose imaginary parts
c           are the odd numbered real values.  separate and supply
c           the second half by conjugate symmetry.
c
  116    icase = 1
         ifmin = 1
         i1rng = np1
         if (idim-4) 117,122,122
  117    if (iform) 118,118,122
  118    icase = 2
         i1rng = np0*(1+nprev/2)
         if (idim-1) 119,119,122
  119    icase = 3
         i1rng = np1
         if (ntwo-np1) 122,122,120
  120    icase = 4
         ifmin = 2
         ntwo = ntwo/2
         n = n/2
         np2 = np2/2
         ntot = ntot/2
         i = 1
         do 121 j=1,ntot
            data(j) = data(i)
            i = i+2
  121    continue
c
c     shuffle data by bit reversal, since n=2**k.  as the shuffling
c     can be done by simple interchange, no working array is needed
c
  122    if (ntwo-np2) 132,123,123
  123    np2hf = np2/2
         j = 1
         do 131 i2=1,np2,np1
            if (j-i2) 124,127,127
  124       i1max = i2+np1-2
            do 126 i1=i2,i1max,2
               do 125 i3=i1,ntot,np2
                  j3 = j+i3-i2
                  tempr = data(i3)
                  tempi = data(i3+1)
                  data(i3) = data(j3)
                  data(i3+1) = data(j3+1)
                  data(j3) = tempr
                  data(j3+1) = tempi
  125          continue
  126       continue
  127       m = np2hf
  128       if (j-m) 130,130,129
  129       j = j-m
            m = m/2
            if (m-np1) 130,128,128
  130       j = j+m
  131    continue
         go to 142
c
c     shuffle data by digit reversal for general n
c
  132    nwork = 2*n
         do 141 i1=1,np1,2
            do 140 i3=i1,ntot,np2
               j = i3
               do 138 i=1,nwork,2
                  if (icase-3) 133,134,133
  133             work(i) = data(j)
                  work(i+1) = data(j+1)
                  go to 135
  134             work(i) = data(j)
                  work(i+1) = 0.
  135             ifp2 = np2
                  if = ifmin
  136             ifp1 = ifp2/ifact(if)
                  j = j+ifp1
                  if (j-i3-ifp2) 138,137,137
  137             j = j-ifp2
                  ifp2 = ifp1
                  if = if+1
                  if (ifp2-np1) 138,138,136
  138          continue
               i2max = i3+np2-np1
               i = 1
               do 139 i2=i3,i2max,np1
                  data(i2) = work(i)
                  data(i2+1) = work(i+1)
                  i = i+2
  139          continue
  140       continue
  141    continue
c
c     main loop for factors of two.  perform fourier transforms of
c     length four, with one of length two if needed.  the twiddle factor
c     w=exp(isign*2*pi*sqrt(-1)*m/(4*mmax)).  check for w=isign*sqrt(-1)
c     and repeat for w=w*(1+isign*sqrt(-1))/sqrt(2).
c
  142    if (ntwo-np1) 174,174,143
  143    np1tw = np1+np1
         ipar = ntwo/np1
  144    if (ipar-2) 149,146,145
  145    ipar = ipar/4
         go to 144
  146    do 148 i1=1,i1rng,2
            do 147 k1=i1,ntot,np1tw
               k2 = k1+np1
               tempr = data(k2)
               tempi = data(k2+1)
               data(k2) = data(k1)-tempr
               data(k2+1) = data(k1+1)-tempi
               data(k1) = data(k1)+tempr
               data(k1+1) = data(k1+1)+tempi
  147       continue
  148    continue
  149    mmax = np1
  150    if (mmax-ntwo/2) 151,174,174
  151    lmax = max0(np1tw,mmax/2)
         do 173 l=np1,lmax,np1tw
            m = l
            if (mmax-np1) 156,156,152
  152       theta = -twopi*float(l)/float(4*mmax)
            if (isign) 154,153,153
  153       theta = -theta
  154       wr = cos(theta)
            wi = sin(theta)
  155       w2r = wr*wr-wi*wi
            w2i = 2.*wr*wi
            w3r = w2r*wr-w2i*wi
            w3i = w2r*wi+w2i*wr
  156       do 169 i1=1,i1rng,2
               kmin = i1+ipar*m
               if (mmax-np1) 157,157,158
  157          kmin = i1
  158          kdif = ipar*mmax
  159          kstep = 4*kdif
               if (kstep-ntwo) 160,160,169
  160          do 168 k1=kmin,ntot,kstep
                  k2 = k1+kdif
                  k3 = k2+kdif
                  k4 = k3+kdif
                  if (mmax-np1) 161,161,164
  161             u1r = data(k1)+data(k2)
                  u1i = data(k1+1)+data(k2+1)
                  u2r = data(k3)+data(k4)
                  u2i = data(k3+1)+data(k4+1)
                  u3r = data(k1)-data(k2)
                  u3i = data(k1+1)-data(k2+1)
                  if (isign) 162,163,163
  162             u4r = data(k3+1)-data(k4+1)
                  u4i = data(k4)-data(k3)
                  go to 167
  163             u4r = data(k4+1)-data(k3+1)
                  u4i = data(k3)-data(k4)
                  go to 167
  164             t2r = w2r*data(k2)-w2i*data(k2+1)
                  t2i = w2r*data(k2+1)+w2i*data(k2)
                  t3r = wr*data(k3)-wi*data(k3+1)
                  t3i = wr*data(k3+1)+wi*data(k3)
                  t4r = w3r*data(k4)-w3i*data(k4+1)
                  t4i = w3r*data(k4+1)+w3i*data(k4)
                  u1r = data(k1)+t2r
                  u1i = data(k1+1)+t2i
                  u2r = t3r+t4r
                  u2i = t3i+t4i
                  u3r = data(k1)-t2r
                  u3i = data(k1+1)-t2i
                  if (isign) 165,166,166
  165             u4r = t3i-t4i
                  u4i = t4r-t3r
                  go to 167
  166             u4r = t4i-t3i
                  u4i = t3r-t4r
  167             data(k1) = u1r+u2r
                  data(k1+1) = u1i+u2i
                  data(k2) = u3r+u4r
                  data(k2+1) = u3i+u4i
                  data(k3) = u1r-u2r
                  data(k3+1) = u1i-u2i
                  data(k4) = u3r-u4r
                  data(k4+1) = u3i-u4i
  168          continue
               kdif = kstep
               kmin = 4*(kmin-i1)+i1
               go to 159
  169       continue
            m = m+lmax
            if (m-mmax) 170,170,173
  170       if (isign) 171,172,172
  171       tempr = wr
            wr = (wr+wi)*rthlf
            wi = (wi-tempr)*rthlf
            go to 155
  172       tempr = wr
            wr = (wr-wi)*rthlf
            wi = (tempr+wi)*rthlf
            go to 155
  173    continue
         ipar = 3-ipar
         mmax = mmax+mmax
         go to 150
c
c     main loop for factors not equal to two.  apply the twiddle factor
c     w=exp(isign*2*pi*sqrt(-1)*(j1-1)*(j2-j1)/(ifp1+ifp2)), then
c     perform a fourier transform of length ifact(if), making use of
c     conjugate symmetries.
c
  174    if (ntwo-np2) 175,201,201
  175    ifp1 = ntwo
         if = inon2
         np1hf = np1/2
  176    ifp2 = ifact(if)*ifp1
         j1min = np1+1
         if (j1min-ifp1) 177,177,184
  177    do 183 j1=j1min,ifp1,np1
            theta = -twopi*float(j1-1)/float(ifp2)
            if (isign) 179,178,178
  178       theta = -theta
  179       wstpr = cos(theta)
            wstpi = sin(theta)
            wr = wstpr
            wi = wstpi
            j2min = j1+ifp1
            j2max = j1+ifp2-ifp1
            do 182 j2=j2min,j2max,ifp1
               i1max = j2+i1rng-2
               do 181 i1=j2,i1max,2
                  do 180 j3=i1,ntot,ifp2
                     tempr = data(j3)
                     data(j3) = data(j3)*wr-data(j3+1)*wi
                     data(j3+1) = tempr*wi+data(j3+1)*wr
  180             continue
  181          continue
               tempr = wr
               wr = wr*wstpr-wi*wstpi
               wi = tempr*wstpi+wi*wstpr
  182       continue
  183    continue
  184    theta = -twopi/float(ifact(if))
         if (isign) 186,185,185
  185    theta = -theta
  186    wstpr = cos(theta)
         wstpi = sin(theta)
         j2rng = ifp1*(1+ifact(if)/2)
         do 200 i1=1,i1rng,2
            do 199 i3=i1,ntot,np2
               j2max = i3+j2rng-ifp1
               do 197 j2=i3,j2max,ifp1
                  j1max = j2+ifp1-np1
                  do 193 j1=j2,j1max,np1
                     j3max = j1+np2-ifp2
                     do 192 j3=j1,j3max,ifp2
                        jmin = j3-j2+i3
                        jmax = jmin+ifp2-ifp1
                        i = 1+(j3-i3)/np1hf
                        if (j2-i3) 187,187,189
  187                   sumr = 0.
                        sumi = 0.
                        do 188 j=jmin,jmax,ifp1
                           sumr = sumr+data(j)
                           sumi = sumi+data(j+1)
  188                   continue
                        work(i) = sumr
                        work(i+1) = sumi
                        go to 192
  189                   iconj = 1+(ifp2-2*j2+i3+j3)/np1hf
                        j = jmax
                        sumr = data(j)
                        sumi = data(j+1)
                        oldsr = 0.
                        oldsi = 0.
                        j = j-ifp1
  190                   tempr = sumr
                        tempi = sumi
                        sumr = twowr*sumr-oldsr+data(j)
                        sumi = twowr*sumi-oldsi+data(j+1)
                        oldsr = tempr
                        oldsi = tempi
                        j = j-ifp1
                        if (j-jmin) 191,191,190
  191                   tempr = wr*sumr-oldsr+data(j)
                        tempi = wi*sumi
                        work(i) = tempr-tempi
                        work(iconj) = tempr+tempi
                        tempr = wr*sumi-oldsi+data(j+1)
                        tempi = wi*sumr
                        work(i+1) = tempr+tempi
                        work(iconj+1) = tempr-tempi
  192                continue
  193             continue
                  if (j2-i3) 194,194,195
  194             wr = wstpr
                  wi = wstpi
                  go to 196
  195             tempr = wr
                  wr = wr*wstpr-wi*wstpi
                  wi = tempr*wstpi+wi*wstpr
  196             twowr = wr+wr
  197          continue
               i = 1
               i2max = i3+np2-np1
               do 198 i2=i3,i2max,np1
                  data(i2) = work(i)
                  data(i2+1) = work(i+1)
                  i = i+2
  198          continue
  199       continue
  200    continue
         if = if+1
         ifp1 = ifp2
         if (ifp1-np2) 176,201,201
c
c     complete a real transform in the 1st dimension, n even, by con-
c     jugate symmetries.
c
  201    go to (230,220,230,202),icase
  202    nhalf = n
         n = n+n
         theta = -twopi/float(n)
         if (isign) 204,203,203
  203    theta = -theta
  204    wstpr = cos(theta)
         wstpi = sin(theta)
         wr = wstpr
         wi = wstpi
         imin = 3
         jmin = 2*nhalf-1
         go to 207
  205    j = jmin
         do 206 i=imin,ntot,np2
            sumr = (data(i)+data(j))/2.
            sumi = (data(i+1)+data(j+1))/2.
            difr = (data(i)-data(j))/2.
            difi = (data(i+1)-data(j+1))/2.
            tempr = wr*sumi+wi*difr
            tempi = wi*sumi-wr*difr
            data(i) = sumr+tempr
            data(i+1) = difi+tempi
            data(j) = sumr-tempr
            data(j+1) = -difi+tempi
            j = j+np2
  206    continue
         imin = imin+2
         jmin = jmin-2
         tempr = wr
         wr = wr*wstpr-wi*wstpi
         wi = tempr*wstpi+wi*wstpr
  207    if (imin-jmin) 205,208,211
  208    if (isign) 209,211,211
  209    do 210 i=imin,ntot,np2
            data(i+1) = -data(i+1)
  210    continue
  211    np2 = np2+np2
         ntot = ntot+ntot
         j = ntot+1
         imax = ntot/2+1
  212    imin = imax-2*nhalf
         i = imin
         go to 214
  213    data(j) = data(i)
         data(j+1) = -data(i+1)
  214    i = i+2
         j = j-2
         if (i-imax) 213,215,215
  215    data(j) = data(imin)-data(imin+1)
         data(j+1) = 0.
         if (i-j) 217,219,219
  216    data(j) = data(i)
         data(j+1) = data(i+1)
  217    i = i-2
         j = j-2
         if (i-imin) 218,218,216
  218    data(j) = data(imin)+data(imin+1)
         data(j+1) = 0.
         imax = imin
         go to 212
  219    data(1) = data(1)+data(2)
         data(2) = 0.
         go to 230
c
c     complete a real transform for the 2nd or 3rd dimension by
c     conjugate symmetries.
c
  220    if (i1rng-np1) 221,230,230
  221    do 229 i3=1,ntot,np2
            i2max = i3+np2-np1
            do 228 i2=i3,i2max,np1
               imin = i2+i1rng
               imax = i2+np1-2
               jmax = 2*i3+np1-imin
               if (i2-i3) 223,223,222
  222          jmax = jmax+np2
  223          if (idim-2) 226,226,224
  224          j = jmax+np0
               do 225 i=imin,imax,2
                  data(i) = data(j)
                  data(i+1) = -data(j+1)
                  j = j-2
  225          continue
  226          j = jmax
               do 227 i=imin,imax,np0
                  data(i) = data(j)
                  data(i+1) = -data(j+1)
                  j = j-np0
  227          continue
  228       continue
  229    continue
c
c     end of loop on each dimension
c
  230    np0 = np1
         np1 = np2
         nprev = n
  231 continue
  232 return
c
c revision history---
c
c january 1978     deleted references to the  *cosy  cards and
c                  added revision history
c-----------------------------------------------------------------------
      end
c
c
c
c
c
c
c
c fftpow2    from nssl                                     08/15/79
      subroutine fft2rc (rldat,nrlpts,signex,hctrn,nhcpts,ierror)
c package fftpow2
c
c
c purpose                fast fourier transforms for data whose length
c                        is a power of two.  the file fftpow2 contains
c                        three routines to handle various forms of input
c                        data as tabulated below:
c                          routine name  input form     output form
c                          ------------  ----------     -----------
c                             fft2rc        real        half complex
c                             fft2cr     half complex      real
c                             fft2cc       complex        complex
c                        half complex refers here to the first n/2+1
c                        complex values of a conjugate symmetric array
c                        of length n.
c
c access cards           *fortran,s=ulib,n=fftpow2
c
c                        *ascent,s=ulib,n=afft2
c
c
c space required         1025 (octal)  (including both files)
c
c subroutine fft2rc (rldat,nrlpts,signex,hctrn,nhcpts,ierror)
c
c
c dimension of           real rldat (nrlpts)
c arguments              complex hctrn (nhcpts)
c                          where nhcpts = nrlpts/2+1
c
c latest revision        january 1978
c
c purpose                fourier transform for powers of 2 (real to half
c                        complex). (in case signex = +1., the
c                        computation performed by this routine is
c                        sometimes called a forward transform or fourier
c                        analysis.)  the discrete fourier transform of
c                        an array rldat, containing nrlpts (a power of
c                        2) real values, is a set cf of nrlpts complex
c                        values satisfying the conjugate symmetry
c                        relation cf(nrlpts+2-k) = conjg(cf(k)))
c                        (k = 2,nrlpts).  due to this symmetry relation
c                        it is only necessary to compute the first
c                        nhcpts = nrlpts/2+1 complex values, and these
c                        are returned in the complex array hctrn.
c
c usage                  call fft2rc (rldat,nrlpts,signex,hctrn,nhcpts,
c                                     ierror).
c                          the original values of rldat may be
c                          regenerated from hctrn by first dividing all
c                          values of hctrn by nrlpts and then calling
c                            fft2cr (hctrn,nhcpts,-signex,rldat,nrlpts,
c                                    ierror).
c
c note                   for these comments we assume
c                        signex = +1. or -1. and define cex(x)
c                        (for all real x) to be the complex exponential
c                        of signex*2*pi*i*x/nrlpts  where pi = 3.14...
c                        and i = sqrt(-1.).
c
c arguments
c
c on input               rldat
c                          a real array containing the nrlpts data
c                          values to be transformed.  it may be
c                          equivalenced to hctrn if desired.  the
c                          dimension is assumed to be rldat(nrlpts).
c
c                        nrlpts
c                          the number of real data values to be
c                          transformed.  it must be a positive power of
c                          2 or a fatal error is flagged.
c
c                        signex
c                          a variable whose sign determines the sign of
c                          the argument of the complex exponential used
c                          in the transform computations.  for
c                          convenience we assume in these comments that
c                          signex is +1. or -1., but the routine in fact
c                          only uses the sign of its value.
c
c                        nhcpts
c                          the number of complex values to be returned
c                          as the half complex transform result in array
c                          hctrn.  it must be = nrlpts/2+1 or a fatal
c                          error is flagged.
c                          note:  nhcpts is not an output parameter.  it
c                                 must be set to nrlpts/2+1 by the user,
c                                 and nhcpts complex locations must be
c                                 provided for the output array hctrn.
c
c on output              hctrn
c                          the complex array containing essentially the
c                          first half of the conjugate symmetric
c                          transform result.  hctrn(k)
c                          (for k = 1,nhcpts) is the complex value
c                          defined by
c                                       nrlpts
c                            hctrn(k) =  sum  rldat(j)*cex((j-1)*(k-1))
c                                       j = 1
c                          these values are also referred to as the
c                          fourier coefficients.
c                            the dimension is assumed to be
c                            complex hctrn (nhcpts) which requires
c                            2*nhcpts = nrlpts+2 core locations.
c                          note:  hctrn may be equivalenced to rldat if
c                                 desired, but they are not identical in
c                                 size.
c
c                        ierror
c                          an error flag with the following meanings
c                            0  no error.
c                          101  nrlpts is not a positive power of 2.
c                          102  nhcpts is not nrlpts/2+1.
c
c common blocks          fft2cm
c
c i/o                    none
c
c precision              single
c
c required ulib          fft2
c routines                 an assembly language implementation of the
c                          fast fourier transform algorithm.  the access
c                          cards are:
c                               *ascent,s=ulib,n=afft2
c
c                                 a slower fortran equivalent is
c                                 available with access cards
c                               *fortran,s=ulib,n=ffft2
c
c
c specialist             dave fulker,ncar, boulder, colorado  80303
c
c language               fortran
c
c history                developed 1971-1973 by dave fulker at ncar;
c                        standardized november 9, 1973
c
c
c
c
      real            rldat(nrlpts)
      complex         hctrn(nhcpts)
      common /fft2cm/ exp1       ,npts       ,nskip      ,mtrn       ,
     1                mskip      ,sgn        ,ierr
      complex         exp1       ,a          ,b          ,apib
      complex         args       ,argj       ,c1         ,c2
      apib(a,b) = cmplx(real(a)-aimag(b),aimag(a)+real(b))
c the following call is for gathering statistics on library use at ncar
      ierror = 0
      if (nrlpts .lt. 2) go to 103
      npts = nrlpts/2
      if (nhcpts .ne. npts+1) go to 104
      nskip = 1
      mtrn = 1
      mskip = 1
      sgn = sign(1.,signex)
      w = 2.*3.141592653589793/float(nrlpts)
      ar = cos(w)
      ai = sin(w)
      exp1 = cmplx(ar**2-ai**2,2.*ar*ai)
      args = cmplx(ar,sgn*ai)
      do 101 i=1,npts
         hctrn(i) = cmplx(.5*rldat(2*i-1),.5*rldat(2*i))
  101 continue
      call fft2 (hctrn)
      if (ierr .ne. 0) go to 103
      hctrn(npts+1) = hctrn(1)
      argj = (1.,0.)
      mc = npts+2
      indm = npts/2+1
      do 102 ind=1,indm
         indc = mc-ind
         c1 = hctrn(ind)+conjg(hctrn(indc))
         c2 = argj*(hctrn(ind)-conjg(hctrn(indc)))
         hctrn(ind) = apib(c1,-c2)
         hctrn(indc) = conjg(apib(c1,c2))
         argj = argj*args
  102 continue
      return
  103 ierror = 101
      write(2,*)' nrlpts is not a positive power of two'
      go to 105
  104 ierror = 102
      write(2,*)' nhcpts is not nrlpts/2 + 1'
  105 stop
      end
      subroutine fft2cr (hcdat,nhcpts,signex,rltrn,nrlpts,ierror)
c
c
c dimension of           complex hcdat (nhcpts)
c arguments              real rltrn (nrlpts)
c                          where nrlpts/2+1 = nhcpts
c
c latest revision        january 1978
c
c purpose                fourier transform for powers of 2 (half complex
c                        to real).  (in case signex = -1., the
c                        computation performed by this routine is
c                        sometimes called a backward transform or
c                        fourier synthesis.)  hcdat is assumed to be the
c                        first nhcpts = nrlpts/2+1 complex values of a
c                        set cf containing nrlpts (a power of 2) complex
c                        values which satisfy the conjugate symmetry
c                        relation cf(nrlpts+2-k) = conjg(cf(k))
c                        (k = 2,nrlpts).  in addition, hctrn(1) is
c                        assumed to be real.  the discrete fourier
c                        transform of such a conjugate symmetric array
c                        is a set of nrlpts real values which are
c                        returned in the real array rltrn.  for more
c                        detail, see the description of the argument
c                        rltrn.
c
c usage                  call fft2cr (hcdat,nhcpts,signex,rltrn,nrlpts,
c                                     ierror)
c                          the original values of hcdat may be
c                          regenerated from rltrn by first dividing all
c                          values of rltrn by nrlpts and then calling
c                            fft2rc (rltrn,nrlpts,-signex,hcdat,nhcpts,
c                                    ierror).
c
c note                   for these comments we assume
c                        signex = +1. or -1. and define cex(x)
c                        (for all real x) to be the complex exponential
c                        of signex*2*pi*i*x/nrlpts  where pi = 3.14...
c                        and i = sqrt(-1.).
c
c arguments
c
c on input               hcdat
c                          a complex array containing nhcpts complex
c                          values which comprise essentially the first
c                          half of the conjugate symmetric data to be
c                          transformed.  it may be equivalenced to rltrn
c                          if desired, but it is not of identical size.
c                            the dimension is assumed to be
c                            complex hcdat(nhcpts) which requires
c                            2*nhcpts = nrlpts+2 core locations.
c
c                        nhcpts
c                          the number of complex values entered in
c                          hcdat.  nrlpts = (nhcpts-1)*2 must be a
c                          positive power of 2 or a fatal error is
c                          flagged.
c
c                        signex
c                          a variable whose sign determines the sign of
c                          the argument of the complex exponential used
c                          in the transform computations.  for
c                          convenience we assume in these comments that
c                          signex is +1. or -1., but the routine in fact
c                          only uses the sign of its value.
c
c                        nrlpts
c                          the number of real transform values to be
c                          returned.  it must be (nhcpts-1)*2 or a fatal
c                          error is flagged.
c                          note:  nrlpts is not an output parameter.  it
c                                 must be set to (nhcpts-1)*2 by the
c                                 user, and nrlpts real locations must
c                                 be provided for the output array
c                                 rltrn.
c
c on output              rltrn
c                          the real array containing the transform
c                          result.  rltrn(j) (for j = 1,nrlpts) is the
c                          real value defined by
c                                       nrlpts
c                            rltran(j) =  sum  cf(k)*cex((j-1*(k-1))
c                                       k = 1
c                          where
c                                    ( hcdat(k)/nrlpts
c                                    ( for 1 .le. k .le. nhcpts
c                            cf(k) = (
c                                    ( conjg(cf(nrlpts+2-k))  otherwise.
c
c                          the dimension is assumed to be rltrn(nrlpts).
c                            note:  rltrn may be equivalenced to ncdat
c                                   if desired.
c
c                        ierror
c                          an error flag with the following meanings:
c                            0  no error.
c                          101  nrlpts is not a positive power of 2.
c                          102  nhcpts is not nrlpts/2+1.
c
c common blocks          fft2cm
c
c i/o                    none
c
c precision              single
c
c required ulib          fft2
c routines                 an assembly language implementation of the
c                          fast fourier transform algorithm.  the access
c                          cards are:
c                               *ascent,s=ulib,n=afft2
c
c                                 a slower fortran equivalent is
c                                 available with access cards.
c                               *fortran,s=ulib,n=ffft2
c
c
c specialist             dave fulker, ncar, boulder, colorado  80303
c
c language               fortran
c
c
c
c
      complex         hcdat(nhcpts)
      real            rltrn(nrlpts)
      common /fft2cm/ exp1       ,npts       ,nskip      ,mtrn       ,
     1                mskip      ,sgn        ,ierr
      complex         exp1
      complex         args       ,argj       ,c1         ,c2
c the following call is for gathering statistics on library use at ncar
      ierror = 0
      if (nrlpts .lt. 2) go to 102
      npts = nrlpts/2
      if (nhcpts .ne. npts+1) go to 103
      nskip = 1
      mtrn = 1
      mskip = 1
      sgn = sign(1.,signex)
      w = 2.*3.141592653589793/float(nrlpts)
      ar = cos(w)
      ai = sin(w)
      exp1 = cmplx(ar**2-ai**2,2.*ar*ai)
      args = cmplx(ar,sgn*ai)
      rltrn(1) = real(hcdat(1))+real(hcdat(npts+1))
      rltrn(2) = real(hcdat(1))-real(hcdat(npts+1))
      argj = args
      mc = npts+2
      indm = npts/2+1
      do 101 ind=2,indm
         indc = mc-ind
         c1 = hcdat(ind)+conjg(hcdat(indc))
         c2 = argj*(hcdat(ind)-conjg(hcdat(indc)))
         rltrn(2*ind-1) = real(c1)-aimag(c2)
         rltrn(2*ind) = aimag(c1)+real(c2)
         rltrn(2*indc-1) = real(c1)+aimag(c2)
         rltrn(2*indc) = -aimag(c1)+real(c2)
         argj = argj*args
  101 continue
      call fft2 (rltrn)
      if (ierr .ne. 0) go to 102
      return
  102 ierror = 101
      write(2,*)' nrlpts is not a positive power of two'
      go to 104
  103 ierror = 102
      write(2,*)' nhcpts is not nrlpts/2 + 1    '
  104 stop
      end
      subroutine fft2cc (cdata,ncpts,signex,ctran,ierror)
c
c
c dimension of           complex cdata (ncpts),ctran (ncpts)
c arguments
c
c latest revision        january 1978
c
c purpose                fourier transform for powers of 2 (complex to
c                        complex).  (in case signex = +1., the
c                        computation performed by this routine is
c                        sometimes called a forward transform or fourier
c                        analysis.  for signex = -1., it is called a
c                        backward transform or fourier synthesis.)  the
c                        discrete fourier transform of an array cdata,
c                        containing ncpts (a power of 2) complex values,
c                        is a set of ncpts complex values which are
c                        returned in the complex array ctran.
c
c usage                  call fft2cc (cdata,ncpts,signex,ctran,ierror)
c                          the original values of cdata may be
c                          regenerated from ctran by first dividing all
c                          values of ctran by ncpts and then calling
c                          fft2cc (ctran,ncpts,-signex,cdata,ierror).
c
c note                   for these comments we assume
c                        signex = +1. or -1. and define cex(x)
c                        (for all real x) to be the complex exponential
c                        of signex*2*pi*i*x/nrlpts  where pi = 3.14...
c                        and i = sqrt(-1.).
c
c arguments
c on input               cdata
c                          a complex array containing the ncpts complex
c                          data values to be transformed.  it may be
c                          equivalenced to ctran if desired.
c                            the dimension is assumed to be
c                            complex cdata(ncpts).
c
c                        ncpts
c                          the number of complex data values to be
c                          transformed.  it must be a positive power of
c                          2 or a fatal error is flagged.
c
c                        signex
c                          a variable whose sign determines the sign of
c                          the argument of the complex exponential used
c                          in the transform computations.  for
c                          convenience we assume in these comments that
c                          signex is +1. or -1., but the routine in fact
c                          only uses the sign of its value.
c
c on output              ctran
c                          the complex array in which the ncpts complex
c                          transform results are returned.  ctran(k)
c                          (for k = 1,ncpts) is the complex value
c                          defined by
c                                       ncpts
c                            ctran(k) =  sum  cdata(j)*cex((j-1)*(k-1))
c                                       j = 1
c                          these values are also referred to as the
c                          fourier coefficients.
c                            the dimension is assumed to be
c                            complex ctran(ncpts).
c                          note:  ctran may be equivalenced to cdata if
c                                 desired.
c
c                        ierror
c                          an error flag with the following meanings:
c                            0  no error.
c                          101  ncpts is not a positive power of 2.
c
c common blocks          fft2cm
c
c i/o                    none
c
c precision              single
c
c required ulib          ffft2
c routines                 an assembly language implementation of the
c                          fast fourier transform algorithm.  the access
c                          cards are:
c                               *ascent,s=ulib,n=afft2
c
c                                 a slower fortran equivalent is
c                                 available with access cards
c                               *fortran,s=ulib,n=ffft2
c
c
c specialist             dave fulker, ncar, boulder, colorado  80303
c
c language               fortran
c
c history                developed 1971-1973 by dave fulker at ncar,
c                        standardized october 16, 1973.
c
c
c
c
      complex         cdata(ncpts)           ,ctran(ncpts)
      common /fft2cm/ exp1       ,npts       ,nskip      ,mtrn       ,
     1                mskip      ,sgn        ,ierr
      complex         exp1
c the following call is for gathering statistics on library use at ncar
      ierror = 0
c
      do 101 j=1,ncpts
         ctran(j) = cdata(j)
  101 continue
  102 if (ncpts .lt. 2) go to 103
      npts = ncpts
      nskip = 1
      mtrn = 1
      mskip = 1
      sgn = sign(1.,signex)
      w = 2.*3.1415926535897931/float(npts)
      ar = cos(w)
      ai = sin(w)
      exp1 = cmplx(ar,ai)
      call fft2 (ctran)
      if (ierr .eq. 0) return
  103 if (ncpts .eq. 1) return
      ierror = 101
      write(2,*)' ncpts is not a positive power of two    '
      stop
c
c revision history---
c
c june 1977        rearranged common blocks and removed reference
c                  to the  loc  function to enhance portability.
c
c january 1978     deleted references to the  *cosy  cards, moved
c                  the revision histories to appear before the
c                  final end card, and moved the initial comment
c                  cards to appear after the first subroutine card
c-----------------------------------------------------------------------
      end
c
c afft2      from portlib                                  08/15/79
      subroutine fft2 (data)
c
c this is a fortran version of the assembly language subroutine fft2, to
c be used with the nssl subroutine package fftpow2.
c
      common /fft2cm/ exp1       ,npts       ,nskip      ,mtrn       ,
     1                mskip      ,sgn        ,ierr
c
c  data is a complex array of length npts which serves for both input
c         and output.
c
      complex         data(99)   ,exp1       ,exk        ,exj        ,h
      real            hh(2)
      equivalence     (hh,h)
c
c  special 1-dimensional version - nskip, mtrn, and mskip are ignored.
c
      ierr = 0
      i2k = npts
      if (i2k .eq. 1) return
      sgn1 = sign(1.,sgn)
      exk = exp1
      if (sgn1 .lt. 0.) exk = conjg(exk)
  101 continue
      i2kp = i2k
      i2k = i2k/2
      if (2*i2k .ne. i2kp) go to 117
      jli = i2k/2+1
      do 109 jl=1,i2k
         if (jl-1) 102,102,104
  102    exj = (1.,0.)
         do 103 jj=jl,npts,i2kp
            h = data(jj)-data(jj+i2k)
            data(jj) = data(jj)+data(jj+i2k)
            data(jj+i2k) = h
  103    continue
         go to 109
  104    if (jl-jli) 105,107,105
c
c  increment jl-dependent exponential factor
c
  105    exj = exj*exk
         do 106 jj=jl,npts,i2kp
            h = data(jj)-data(jj+i2k)
            data(jj) = data(jj)+data(jj+i2k)
            data(jj+i2k) = h*exj
  106    continue
         go to 109
  107    exj = cmplx(0.,sgn1)
         do 108 jj=jl,npts,i2kp
            h = data(jj)-data(jj+i2k)
            data(jj) = data(jj)+data(jj+i2k)
            data(jj+i2k) = cmplx(-sgn1*hh(2),sgn1*hh(1))
  108    continue
  109 continue
c
c  increment k-dependent exponential factor
c
      exk = exk*exk
      if (i2k-1) 110,110,101
  110 if (npts .le. 2) return
      nptsd2 = npts/2
      jmin = 0
      jmax = npts-4
      jrev = 0
      do 116 j=jmin,jmax,2
         i2k = nptsd2
         jrev2 = jrev+i2k
         if (jrev2-(j+1)) 113,113,111
  111    h = data(jrev2+1)
         data(jrev2+1) = data(j+2)
         data(j+2) = h
         if (jrev-j) 113,113,112
  112    h = data(jrev+1)
         data(jrev+1) = data(j+1)
         data(j+1) = h
  113    continue
  114    i2k = i2k/2
         jrev = jrev-i2k
         if (jrev) 115,114,114
  115    jrev = jrev+2*i2k
  116 continue
      return
  117 ierr = 101
      return
      end
c &&&&&&&&&&&&&&&&

