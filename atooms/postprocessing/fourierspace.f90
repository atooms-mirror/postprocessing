! Note: complex(8) in the f2py kind convention is converted to complex_double. This means both real and imaginary parts are double precision. See f2cmap_all['complexkind'] entry in numpy/f2py/capi_maps.py.
module fourierspace_module

   implicit none

contains

  subroutine sk_bare(expo,ikvec,rho)
    complex(8),intent(in)        :: expo(:,:,:)  ! (ndim, npart)
    integer, intent(in)       :: ikvec(:,:) ! (ndim, nvec), (nvec)
    integer     :: i1, i2, i3, ii
    complex(8), intent(inout) :: rho(:)
    do ii = 1,size(ikvec,2)
       i1   = ikvec(1,ii)
       i2   = ikvec(2,ii)
       i3   = ikvec(3,ii)
       rho(ii) = sum(expo(:,1,i1) * expo(:,2,i2) * expo(:,3,i3))
    end do
  end subroutine sk_bare

  ! This is the original python code
  !   acf[kk][dt] += numpy.sum(x[i0+i, :, 0, ik[0]]*x[i0, :, 0, ik[0]].conjugate() *
  !                            x[i0+i, :, 1, ik[1]]*x[i0, :, 1, ik[1]].conjugate() *
  !                            x[i0+i, :, 2, ik[2]]*x[i0, :, 2, ik[2]].conjugate()).real
  ! So we pass expo (x) and ik and do the sum
  ! We also pass i0 and i do avoid slicing in numpy
  ! This kernel vectorizes, while this is not the case for the implicit loop
  function fskt_kernel_3d(expo,t1,t2,ik) result (output)
    complex(8),intent(in)       :: expo(:,:,:,:)  ! (nsteps, npart, ndim, kvec)
    integer, intent(in)         :: t1,t2,ik(:)  ! (ndim)
    complex(8)                  :: output, tmp(size(expo,2))
    integer :: i
    do i = 1,size(expo,2)
       tmp(i) = expo(t1,i,1,ik(1)) * CONJG(expo(t2,i,1,ik(1))) * &
                expo(t1,i,2,ik(2)) * CONJG(expo(t2,i,2,ik(2))) * &
                expo(t1,i,3,ik(3)) * CONJG(expo(t2,i,3,ik(3)))
    end do
    output = SUM(tmp)
  end function fskt_kernel_3d

  function fskt_kernel_2d(expo,t1,t2,ik) result (output)
    complex(8),intent(in)       :: expo(:,:,:,:)  ! (nsteps, npart, ndim, kvec)
    integer, intent(in)         :: t1,t2,ik(:)  ! (ndim)
    complex(8)                  :: output, tmp(size(expo,2))
    integer :: i
    do i = 1,size(expo,2)
       tmp(i) = expo(t1,i,1,ik(1)) * CONJG(expo(t2,i,1,ik(1))) * &
                expo(t1,i,2,ik(2)) * CONJG(expo(t2,i,2,ik(2)))
    end do
    output = SUM(tmp)
  end function fskt_kernel_2d

  function fskt_kernel_nd(expo,t1,t2,ik) result (output)
    complex(8),intent(in)       :: expo(:,:,:,:)  ! (nsteps, npart, ndim, kvec)
    integer, intent(in)         :: t1,t2,ik(:)  ! (ndim)
    complex(8)                  :: output, tmp(size(expo,2))
    integer :: i, j
    tmp(:) = expo(t1,:,1,ik(1)) * CONJG(expo(t2,:,1,ik(1)))
    do j = 2,SIZE(expo,3)
       do i = 1,size(expo,2)
          tmp(i) = tmp(i) * expo(t1,i,j,ik(j)) * CONJG(expo(t2,i,j,ik(j)))
       end do
    end do
    output = SUM(tmp)
  end function fskt_kernel_nd

end module fourierspace_module
