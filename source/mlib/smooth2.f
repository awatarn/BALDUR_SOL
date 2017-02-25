c@smooth2.f
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine smooth2 ( array, idim, ztemp1, ztemp2
     &  , index1, ilower, iupper, iorder, cmix )
c
c    sbrtn smooth2 applies a smoothing operator to an array
c
c  array  = 2-D array to be smoothed
c  idim   = first dimension of array
c  ztemp1, ztemp2 = temporary arrays for internal use
c           these must have at least iupper elements
c  index1 = index of first dimension of array, held fixed
c             Must be less than idim
c  ilower = lower bound of the part of the array to be smoothed
c  iupper = upper bound of the part of the array to be smoothed
c           Note:  array(ilower) and array(iupper) will be held fixed
c                  these refer to the second dimension of array
c  iorder = number of smoothings to be applied
c  cmix   = fraction of the original array
c             to be added to the smoothed array
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      dimension array(idim,*), ztemp1(*), ztemp2(*)
c
      do j=ilower,iupper
        ztemp1(j) = array(index1,j)
      enddo
c
      do js=1,iorder
c
        do j=ilower,iupper
          ztemp2(j) = ztemp1(j)
        enddo
c
        do j=ilower+1,iupper-1
          ztemp1(j) = 0.25*( ztemp2(j-1) + 2.0*ztemp2(j) + ztemp2(j+1) )
        enddo
c
      enddo
c
      do j=ilower+1,iupper-1
        array(index1,j) = cmix * array(index1,j)
     &               + ( 1.0 - cmix ) * ztemp1(j)
      enddo
c
      return
      end
