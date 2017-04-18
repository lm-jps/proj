c#######################################################################
      module number_types
c
      implicit none
c
c-----------------------------------------------------------------------
c ****** Basic number types.
c ****** This module is used to set the default precision for REALs.
c-----------------------------------------------------------------------
c
c ****** Set up the KIND values for the various REALs.
c
      integer, parameter :: KIND_REAL_4=kind(1.0e0)
      integer, parameter :: KIND_REAL_8=kind(1.0d0)
c
c ****** KIND values for specifying the precision of REALs.
c
      integer, private, parameter :: r4=KIND_REAL_4
      integer, private, parameter :: r8=KIND_REAL_8
c
c ****** Select the number type for REALs (one of: R4|R8).
c
      integer, parameter :: r_typ=r8
c
      end module
