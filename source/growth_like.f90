!Routine to compute likelihood from growth of structure data 
!measured by Redshift space distortions in the WiggleZ survey

!November 2010 (David Parkinson)

module growth_module
  use cmbtypes
  use Precision
  use MatrixUtils
  implicit none

  logical :: Use_growth = .false.
  integer, parameter :: WigglezbinN = 4
  real(dl), parameter :: Pi_num = 3.14159265359D0 
  real(dl) :: growth_z(14), growth_fsig8(WigglezbinN), growth_diagerr(WigglezbinN)
  real(dl) :: growth_AP_F(WigglezbinN), growth_AP_F_diagerr(WigglezbinN), fsig8_AP_corr(WigglezbinN)
  real(dl) :: growth_Ninv(WigglezbinN,WigglezbinN), growth_Ninvmarge(WigglezbinN,WigglezbinN)
  real(dl) :: growth_z_eps
  logical, save :: do_growth_init = .true.
  real(dl) :: omegam_0, omegak_0, omegav_0, w0
contains


  subroutine growth_SetTransferRedshifts(redshifts)
    real, intent(inout) :: redshifts(*)
    integer i,j
   !input is default log z spacing; can change here; check for consistency with other (e.g. lya)
         
   !Note internal ordering in CAMB is the opposite to that used in cosmomc transfer arrays (as here)
   !first index here must be redshift zero
   
    if(Use_growth.and.matter_power_lnzsteps<13) & 
         call MpiStop('For growth matter_power_lnzsteps should be set to at least 13 (hardcoded in cmbtypes)')
    growth_z(1) = 0.d0
    if(matter_power_lnzsteps==1 .or. .not. use_growth) return
    do i=1,4
       j = 3*i
       if(i.eq.1) then
          growth_z(j) = 0.22
       else if(i.eq.2) then
          growth_z(j) = 0.41
       else if(i.eq.3) then
          growth_z(j) = 0.60
       else if(i.eq.4) then
          growth_z(j) = 0.78
       endif
       growth_z_eps = 1.d-2*growth_z(j)
       growth_z(j-1) = growth_z(j)-growth_z_eps
       growth_z(j+1) = growth_z(j)+growth_z_eps
       if(i.eq.4) growth_z(j+2) = growth_z(j+1)+growth_z_eps
    enddo
    do i=1,14
       print*, i, growth_z(i)
    enddo
    If (Feedback > 0) print*, 'Seting Transfer Redshifts'
    redshifts(1:matter_power_lnzsteps) = growth_z(1:matter_power_lnzsteps)
    return

  end subroutine growth_SetTransferRedshifts

 subroutine growth_init()
   character (LEN=1200) :: InLine
   character (LEN=20) :: names(186)
   real(dl) :: input (186,3)
   real dum
   character(len=32) dum_char
   integer i, num_redshifts_initial
!   real(dl) growth_ADD_func, growth_H_func
!   real(dl) growth_ADD_func, growth_H_func

   if (Feedback > 0) write (*,*) 'reading: growth data'
   call OpenTxtFile('data/cb_growth.dat',tmp_file_unit)
   read(tmp_file_unit,*) dum_char
   do i=1,WigglezbinN
      read(tmp_file_unit, *) dum, growth_fsig8(i), growth_diagerr(i), growth_AP_F(i),&
          growth_AP_F_diagerr(i), fsig8_AP_corr(i)

      growth_diagerr(i)=growth_diagerr(i)**2
      growth_AP_F_diagerr(i)=growth_AP_F_diagerr(i)**2
   end do
   close(tmp_file_unit)
   
   do_growth_init = .false.

 end subroutine growth_init

 function Growth_LnLike(CMB,T)
  !Assume this is called just after CAMB with the correct model

  implicit none
  type(CMBParams) CMB
  type(CosmoTheory) T
  real Growth_LnLike
  integer i,j
  real(dl) z
  real z_tmp(3), a, growth_rate(WigglezbinN),sigma8_theory(WigglezbinN)
  real chisq_bin(WigglezbinN), chisq, k_in
  real diffs(2), step1(2), covmatrix(2,2)
  real z_out, a_tmp(3)
  real delta(3), dd, da, fsig8_theory, delta_0
  real(dl) angdiam, hz
  real AP_F


  if(feedback.gt.0) print*, "calling growth like"
  omegam_0 = CMB%omc+CMB%omb
  omegav_0 = CMB%omv
  omegak_0 = CMB%omk
  w0       = CMB%w
  print*, "do growth init", do_growth_init
  if (do_growth_init) call growth_init
!  delta_0 = sqrt(MatterPowerAt_Z(T,0.1,0.))

     do i=1, WigglezbinN
        !Obviously this is not v efficient...
        z= growth_z(3*i)
        a = 1.d0/(1.d0+z)
        k_in = 0.1d0/(CMB%H0/100.d0)
        do j=1,3
           z_tmp(j) = growth_z(1+((i-1)*3)+j)
           a_tmp(j) = 1./(1.+z_tmp(j))
!           delta(j) = sqrt(k_in**2*MatterPowerAt_Z(T,k_in,z_tmp(j)))
           delta(j) = T%growth_function(1+((i-1)*3)+j)
!           print*, j, z_tmp(j), z_tmp(j), delta(j)
        enddo
        dd = delta(1)-delta(3)
!        da = (1.d0/(1.d0+(z_tmp(1))))-(1.d0/(1.d0+(z_tmp(3))))
        da = a_tmp(1) - a_tmp(3)
        growth_rate(i) = (dd/da)*(a/delta(2))
        sigma8_theory(i) = T%sigma_8(3*i)
        fsig8_theory = growth_rate(i)*sigma8_theory(i)
        diffs(1) = fsig8_theory-growth_fsig8(i)
        covmatrix(1,1) = growth_diagerr(i)
        angdiam = growth_ADD_func(z)
        hz = growth_H_func(z)
!        AP_F = angdiam*hz/(growth_D_A_fid(i)*growth_H_z_fid(i))
        AP_F = angdiam*hz*(1.+z)
        diffs(2) = AP_F-growth_AP_F(i)
        covmatrix(2,2) = growth_AP_F_diagerr(i)
        covmatrix(1,2) = fsig8_AP_corr(i)*sqrt(covmatrix(1,1)*covmatrix(2,2))
        covmatrix(2,1) = covmatrix(1,2)
        call Matrix_Inverse(covmatrix)
        step1 = MATMUL(covmatrix,diffs)
        chisq_bin(i) = dot_product(diffs,step1)
     end do


     chisq = sum(chisq_bin)



     Growth_LnLike = chisq/2
     if (Feedback > 1) write (*,*) 'Growth chisq: ', chisq

!     stop
 end function Growth_LnLike

 function growth_ADD_func(z)
! Angular Diameter Distance function for growth of structure module
! Can't use inbuilt CAMB one as will need to call for fiducial model
  implicit none
  real(dl) growth_ADD_func,z, ADD, r_flat
  real(dl), external :: rombint

  r_flat = rombint(growth_ADD_int_func,0.d0,z,1.d-5)
  if(omegak_0.lt.0.d0) then
    ! elliptical universe
    ADD = (1.d0/sqrt(-omegak_0))*sin(sqrt(-omegak_0)*r_flat)
  else if(omegak_0.gt.0.d0) then
    ! hyperbolic universe
    ADD = (1.d0/sqrt(omegak_0))*sinh(sqrt(omegak_0)*r_flat)
  else
    ADD = r_flat
  endif
  growth_ADD_func = ADD/(1.d0+z)
  return
 end function growth_ADD_func


 function growth_H_func(z)
  ! Normalised Hubble rate as a function of redshift
  ! negelcting contributions from radiation density
  implicit none
  real(dl) growth_H_func, z

  growth_H_func = sqrt(omegam_0*(1.d0+z)**3+omegak_0*(1.d0+z)**2+omegav_0*(1.d0+z)**(3.d0*(1+w0))) 

  return
 end function growth_H_func


 function growth_ADD_int_func(z)
! Function to be interated to get comoving distance in a flat universe
  implicit none
  real(dl) growth_ADD_int_func, z
  

  growth_ADD_int_func = 1.d0/growth_H_func(z)

  return
 end function growth_ADD_int_func

 end module growth_module
