c     Stand-in for hdf5_dump.c in a build made without HDF5, so that the
c     calls in output.f link.  make HDF5=1 compiles the real writer instead.
c     ss_h5_available says there is no writer, init stops a run that asks for
c     nhdf5>0, and so nothing else here is ever called.
      function ss_h5_available() bind(c,name='ss_h5_available')
      use iso_c_binding, only: c_int
      integer(c_int) :: ss_h5_available
      ss_h5_available=0
      end
************************************************************************
      function ss_h5_create() bind(c,name='ss_h5_create')
      use iso_c_binding, only: c_int
      integer(c_int) :: ss_h5_create
      ss_h5_create=1
      end
************************************************************************
      function ss_h5_put_double() bind(c,name='ss_h5_put_double')
      use iso_c_binding, only: c_int
      integer(c_int) :: ss_h5_put_double
      ss_h5_put_double=1
      end
************************************************************************
      function ss_h5_put_int() bind(c,name='ss_h5_put_int')
      use iso_c_binding, only: c_int
      integer(c_int) :: ss_h5_put_int
      ss_h5_put_int=1
      end
************************************************************************
      function ss_h5_close() bind(c,name='ss_h5_close')
      use iso_c_binding, only: c_int
      integer(c_int) :: ss_h5_close
      ss_h5_close=1
      end
