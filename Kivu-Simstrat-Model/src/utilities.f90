! ---------------------------------------------------------------------------------
!     Kivu-Simstrat a physical 1D model for simulating Lake Kivu (Rwanda/Dep. Rep. Congo)
!
!     Developed by:  Group of Applied System Analysis
!                    Dept. of Surface Waters - Research and Management
!                    Eawag - Swiss Federal institute of Aquatic Science and Technology
!
!     Copyright (C) 2025, Eawag
!
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>. 
! ---------------------------------------------------------------------------------
!<    +---------------------------------------------------------------+
!     | Generic utilities that are used throughout the code
!<    +---------------------------------------------------------------+


module utilities
   use strat_kinds
   use strat_consts
   implicit none

   type, public :: string
      character(len=:), allocatable :: str
   end type

   interface toStr
      module procedure str_int, str_real
   end interface

contains

   !> Interpolation of yi on grid zi (based on the given y on grid z)
   subroutine Interp(z, y, num_z, zi, yi, num_zi)
      implicit none
      real(RK), dimension(:), intent(in) :: z, y, zi
      real(RK), dimension(:), intent(out) :: yi
      integer, intent(in) :: num_z, num_zi
      integer :: posk1, posk2, posi, i

      if (num_z == 1) then ! only one value
         yi(1:num_zi) = y(1)
         return
      end if

      ! Assign closest value if out of given grid
      posk1 = 1
      do while (zi(posk1) <= z(1))
         yi(posk1) = y(1)
         posk1 = posk1 + 1
      end do
      posk2 = num_zi
      do while (zi(posk2) >= z(num_z))
         yi(posk2) = y(num_z)
         posk2 = posk2 - 1
      end do
      ! Linear interpolation
      posi = 1
      do i = posk1, posk2
         do while (zi(i) > z(posi + 1))
            posi = posi + 1
         end do
         yi(i) = y(posi) + ((zi(i) - z(posi)) / (z(posi + 1) - z(posi))) * (y(posi + 1) - y(posi))
      end do

      return
   end

   ! Interpolation of yi on grid zi (based on the given y on grid z)
   !####################################################################
   pure subroutine Interp_nan(z, y, num_z, zi, yi, num_zi)
      !####################################################################
      !use, intrinsic :: iso_fortran_env
      use, intrinsic :: ieee_arithmetic
      implicit none

      real(RK), dimension(:), intent(in) :: z, y, zi
      real(RK), dimension(:), intent(out) :: yi
      integer, intent(in) :: num_z, num_zi

      integer posk1, posk2, posi, i

      ! Assign NaN if out of given grid
      posk1 = 1
      do while (zi(posk1) < z(1))
         yi(posk1) = 0.0_RK
         yi(posk1) = ieee_value(yi(posk1), ieee_quiet_nan) ! NaN
         posk1 = posk1 + 1
      end do
      posk2 = num_zi
      do while (zi(posk2) > z(num_z))
         yi(posk2) = 0.0_RK
         yi(posk2) = ieee_value(yi(posk2), ieee_quiet_nan) ! NaN
         posk2 = posk2 - 1
      end do

      ! Linear interpolation
      posi = 1
      do i = posk1, posk2
         do while (zi(i) > z(posi + 1))
            posi = posi + 1
         end do
         yi(i) = y(posi) + ((zi(i) - z(posi)) / (z(posi + 1) - z(posi))) * (y(posi + 1) - y(posi))
      end do

      return
   end subroutine Interp_nan

   !! Integrate discrete function y[x] using the trapezoidal rule
   !!####################################################################
   subroutine Integrate(x, y, inty, num)
      !!####################################################################
      implicit none

      integer :: num, i
      real(RK), dimension(:), intent(in) :: x, y
      real(RK), dimension(:), intent(inout) :: inty

      inty(1) = 0
      do i = 2, num
         inty(i) = inty(i - 1) + 0.5_RK*(x(i) - x(i - 1))*(y(i) + y(i - 1))
      end do
      return
   end

      ! Assign nan to values out of current grid
   !####################################################################
   pure subroutine Assign_nan(y, ubnd, ubnd_grid)
      !####################################################################
      !use, intrinsic :: iso_fortran_env
      use, intrinsic :: ieee_arithmetic
      implicit none

      real(RK), dimension(:), intent(out) :: y
      integer, intent(in) :: ubnd, ubnd_grid

      integer :: i

      ! Assign NaN if out of given grid
      i = ubnd_grid
      do while (i > ubnd)
         y(i) = 0.0_RK
         y(i) = ieee_value(y(i), ieee_quiet_nan) ! NaN
         i = i - 1
      end do

      return
   end subroutine Assign_nan


   pure function linspace(x0, xend, n, endpoint) result(x)
      implicit none

      integer, intent(in) :: n
      real(RK), intent(in) :: x0, xend
      logical, optional, intent(in) :: endpoint
      real(RK), dimension(n) :: x

      real(RK) :: dx
      real(RK) :: denom
      integer :: i

      denom = real(n - 1, RK)
      if (present(endpoint) .and. (.not. endpoint)) then
         denom = real(n, RK)
      end if

      dx = (xend - x0)/denom

      x = [(real(i, RK)*dx + x0, i=0, n - 1)]

      return
   end function linspace

   pure subroutine diff(d, a, N)

      implicit none

      integer, intent(in) :: N
      real(RK), dimension(N - 1), intent(out) :: d
      real(RK), dimension(N), intent(in) :: a

      d = a(2:N) - a(1:N - 1)

      return
   end subroutine diff

   subroutine check_file_exists(fname)
      implicit none
      character(len=*), intent(in) :: fname

      logical :: file_exists
      if (fname == '') then
         call error('Filename is empty')
      else
         inquire (file=fname, exist=file_exists)
         if (.not. file_exists) then
            call error('File '//fname//' does not exist')
         end if
      end if
   end subroutine check_file_exists

   subroutine ok(message)
      implicit none
      character(len=*), intent(in) :: message
      write(*, *) '[OK] '//message
   end subroutine ok

   subroutine error(message)
      implicit none
      character(len=*), intent(in) :: message
      write(6, *) '[ERROR] '//message
      stop
   end subroutine error

   subroutine warn(message)
      implicit none
      character(len=*), intent(in) :: message
      write(*, *) '[WARNING] '//message
   end subroutine warn

   pure function find_index_ordered(array, target_value) result(idx)
      implicit none
      real(RK), dimension(:), intent(in) :: array
      real(RK), intent(in) :: target_value

      integer :: idx

      do idx = 1, size(array)
         if (array(idx) > target_value) exit
      end do
   end function

   pure function linear_interpolate(t_start, t_end, v_start, v_end, t) result(v)
      implicit none
      real(RK), intent(in) :: t_start, t_end, v_start, v_end, t
      real(RK) :: v

      v = v_start + t*(v_end - v_start)/(t_end - t_start)
   end function

   pure function convert2height_above_sed(z, z_zero) result(h)
      implicit none
      real(RK), dimension(:), intent(in) :: z
      real(RK), intent(in) :: z_zero

      real(RK), dimension(size(z)) :: h
      integer :: n

      n = size(z)
      h = -z_zero + z(n:1:-1)
   end function


   ! Reverse an array without allocating a second array
   subroutine reverse_in_place(in_arr)
      implicit none
      real(RK), intent(inout) :: in_arr(:)
      real(RK) :: temp

      integer :: first, last, i, len

      first = lbound(in_arr, dim=1)
      last = ubound(in_arr, dim=1)
      len = size(in_arr)

      ! Works for even and odd sized arrays
      !(as len/2 is always integer and not rounded, but cutoff)
      do i = last, first + int(len/2), -1
         temp = in_arr(i)
         in_arr(i) = in_arr(len + 1 - i)
         in_arr(len + 1 - i) = temp
      end do

   end subroutine

   character(len=20) function str_int(k)
      implicit none
      ! "Convert an integer to string."
      integer, intent(in) :: k
      write (str_int, '(a)') k
      str_int = adjustl(str_int)
   end function str_int

   character(len=20) function str_real(k)
      implicit none
      ! "Convert an integer to string."
      real(RK), intent(in) :: k
      write (str_real, '(a)') k
      str_real = adjustl(str_real)
   end function str_real

   character(len=20) function real_to_str(k, fmt)
      implicit none
      real(RK), intent(in) :: k
      character(len=*), intent(in) :: fmt
      write (real_to_str, fmt) k
      real_to_str = adjustl(real_to_str)
   end function

   pure logical function is_leap_year(year)
      integer, intent(in) :: year
      is_leap_year = mod(year, 4) == 0 .and. (.not. mod(year, 100) == 0 .or. mod(year, 400) == 0)
   end function

   pure function calc_days_per_month(year) result(days_per_month)
      integer, intent(in) :: year
      integer, dimension(12) :: days_per_month

      if (is_leap_year(year)) then
         days_per_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      else
         days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      end if
   end function

   ! Initialize the calender (day, month and year) based on the given starting year and and start datum
   subroutine init_calendar(start_year, datum, current_year, current_month, current_day)
      implicit none
      integer, intent(in) :: start_year
      real(RK), intent(in) :: datum
      integer, intent(out) :: current_year, current_month
      real(RK), intent(out) :: current_day

      ! Local variables
      real(RK) :: elapsed_days, days_left, days_per_year
      integer :: i
      integer, dimension(12) :: days_per_month

      ! Determine current year
      current_year = start_year
      elapsed_days = 0
      do
         if (is_leap_year(current_year)) then
            days_per_year = 366
         else
            days_per_year = 365
         end if

         if ((elapsed_days + days_per_year) < datum) then
            elapsed_days = elapsed_days + days_per_year
            current_year = current_year + 1
         else
            exit
         end if
      end do
      days_left = datum - elapsed_days

      ! Determine current month and day
      days_per_month = calc_days_per_month(current_year)

      do i=1,12
         if (days_left > days_per_month(i)) then
            days_left = days_left - days_per_month(i)
         else
            current_month = i
            current_day = days_left
            exit
         end if
      end do
   end subroutine


   ! Update calendar month and day (used for albedo assignment)
   subroutine update_calendar(current_year, current_month, current_day, dt)
      implicit none
      integer, intent(inout) :: current_year, current_month
      real(RK), intent(inout) :: current_day
      real(RK), intent(in) :: dt

      ! Local variables
      integer, dimension(12) :: days_per_month
      real(RK) :: current_day_new


      ! Prepare day per month arrays
      days_per_month = calc_days_per_month(current_year)

      ! Update current day
      current_day_new = current_day + dt/24/60/60

      ! If new month is reached
      if (ceiling(current_day_new) > days_per_month(current_month)) then
         current_month = current_month + 1
         current_day = current_day_new - floor(current_day_new)

         ! If new year is reached
         if (current_month > 12) then
            current_month = 1
            current_year = current_year + 1
         end if
      else
         ! If not a new month, just go on counting
         current_day = current_day_new
      end if

   end subroutine

   ! Calculate density as a function of T and S
   subroutine calc_density(rho, T, S)
      implicit none

      ! Arguments
      real(RK), intent(in) :: T, S
      real(RK), intent(out) :: rho

      ! Local variables
      real(RK) :: rho0t, rho0st

      ! According to Chen Millero, changed according "Double diffusion in Lake Kivu" from Schmid et al., 2010
      rho0t = 0.9998395_RK + T*(6.7914e-5_RK + T*(-9.0894e-6_RK + T*(1.0171e-7_RK + T*(-1.2846e-9_RK + T*(1.1592e-11_RK + T*(-5.0125e-14_RK))))))
      rho0st = (7.5e-4_RK + T*(-3.85e-6_RK + T*(4.96e-8_RK)))*S

      rho = rho_0*(rho0t + rho0st)
   end subroutine

   ! Transform CO2 to DIC (analogous to AED2 carbon module)
   subroutine calc_co2_from_dic(CO2, T, S, DIC, pH)
      implicit none

      ! Arguments
      real(RK), intent(in) :: T, S, DIC, pH
      real(RK), intent(out) :: CO2

      ! Local variables
      real(RK) :: T_K, pK1, pK2, K1, K2, H

      T_K = T + 273.15

      ! Calculation taken from AED2 carbon module
      pK1 = 3670.7/T_K - 62.008+9.7944*log(T_K) - 0.0118*S + 0.000116*(S**2)
      K1  = 10.**(-pK1)  ! this is on the SWS pH scale in mol/kg-SW
      pK2 = 1394.7/T_K + 4.777 - 0.0184*S + 0.000118*(S**2)
      K2  = 10.**(-pK2)

      H = 10.**(-pH)
      CO2 = DIC*H*H/(H*H + K1*H + K1*K2)
      CO2 = CO2/1e6 ! Transform to mol/L

   end subroutine

   ! Compute density as a function of T, S, CO2 and CH4
   subroutine calc_density_aed2(rho, T, S, co2, ch4)
      implicit none

      ! Arguments
      real(RK), intent(in) :: T, S, co2, ch4
      real(RK), intent(out) :: rho

      ! Local variables
      real(RK) :: rho0t, rho0st, rho0_co2, rho0_ch4

      ! According to Chen Millero, changed according "Double diffusion in Lake Kivu" from Schmid et al., 2010
      rho0t = 0.9998395_RK + T*(6.7914e-5_RK + T*(-9.0894e-6_RK + T*(1.0171e-7_RK + T*(-1.2846e-9_RK + T*(1.1592e-11_RK + T*(-5.0125e-14_RK))))))
      rho0st = 7.5e-4_RK*S

      rho0_co2 = 0.0125*co2   ! Schmid et al., 2010
      rho0_ch4 = -0.02*ch4    ! Schmid et al., 2010

      rho = rho_0*(rho0t + rho0st + rho0_co2 + rho0_ch4)
   end subroutine

   subroutine calc_rho_ratio(r_rho, T1, T2, S1, S2, co2_1, co2_2, ch4_1, ch4_2, depth)
      implicit none

      ! Arguments
      real(RK), intent(in) :: T1, T2, S1, S2, co2_1, co2_2, ch4_1, ch4_2, depth
      real(RK), intent(out) :: r_rho

      ! Local variables
      real(RK) :: alpha, beta_S, beta_co2, beta_ch4, T, S, P

      T = (T1 + T2)/2
      S = (S1 + S2)/2
      P = depth/9.81

      ! According to Chen Millero, changed according "Double diffusion in Lake Kivu" from Schmid et al., 2004 Acta
      alpha = -68.00_RK + 18.2091_RK*T - 0.30866_RK*T**2 + 5.3445e-3_RK*T**3 - 6.0721e-5_RK*T**4 + 3.1441e-7_RK*T**5 + &
       (4.599_RK - 0.1999_RK*T + 2.790e-3_RK*T**2)*S + (0.3682_RK - 1.52e-2_RK*T+ 1.91e-4_RK*T**2 - 4.613e-3_RK*S)*P
      alpha = alpha*1e-6
      beta_S = 7.5e-4_RK

      beta_co2 = 0.0125   ! Schmid et al., 2010
      beta_ch4 = -0.02    ! Schmid et al., 2010

      r_rho = (beta_S*(S1-S2) + beta_co2*(co2_1-co2_2) + beta_ch4*(ch4_1-ch4_2))/(alpha*(T1-T2))
   end subroutine

   subroutine calculate_pH(pH,TEM,SAL,DIC)      ! Added by Modeste, 2025
      implicit none

      real(RK),    intent(in)   :: TEM, SAL, DIC
      real(RK),    intent(out)  :: pH              
      ! LOCAL
      real(RK)                  :: TA0, PRE, K0, KS, KF, fH, KB, KW, KP1, KP2, KP3, KSi = 0., K1, K2, TB, TP, TS, TF, TSi, TC, TA

      !===========Initialize the conditions =========================!

      ! compute total alkalinity from salinity
      !IF (data%kivu_mode) THEN
      TA0 = 11960 * SAL       ! Relationship from Wüest et al., 2009, added by FB, 2020
      TA0 = TA0 / 1.0D6      ! change unit to mol/kgSW
      !ENDIF

      TB  = 0.
      TP  = 0./1.e6
      TS  = 0.
      TF  = 0.
      TSi = 0./1.e6

      TC  = DIC /1.e6 ! convert from mmol/m3 to mol/L
      TA  = TA0 

      PRE = 0.

      Call Cal_constants(TEM, SAL, PRE, K0, KS, KF, fH, KB, KW, KP1, KP2, KP3, &
                              & KSi,  K1, K2, TB, TP, TS, TF)

      Call Cal_pHfromTATC(TA, TC, pH, K1, K2, TB, KB, KW, KP1, KP2, KP3,&
                                 & TP, TSi, TS, KS, KSi, TF, KF)

   end subroutine calculate_pH
   
   !---- Calculate pH from total alkalinity and carbon- aed_carbon module
   subroutine Cal_pHfromTATC(TAx, TCx, pHx, K1F, K2F, TBF, KBF, KWF, KP1F, KP2F, KP3F,&
                           & TPF, TSiF, TSF, KSF, KSiF, TFF, KFF)
      implicit none

      real(RK),     intent(in) :: TAx, TCx
      real(RK),     intent(in) :: K1F, K2F, TBF, KBF, KWF, KP1F, KP2F, KP3F, TPF, TSiF, TSF, KSF, TFF, KFF, KSiF
      real(RK),     intent(out):: pHx

      !local
      real(RK)                 :: H, Denom, CAlk, BAlk, OH, PhosTop, PhosBot, PAlk, SiAlk, &
                                 & FREEtoTOT, Hfree, HSO4, HF, Residual, Slope, deltapH, pHTol, ln10

      integer :: count
      integer :: max_count

      
      pHx      = 8.        !This is the first guess
      pHTol    = 0.0001    !tolerance for iterations end
      ln10     = log(10.)
      deltapH  = pHTol + 1

      count = 0
      max_count = 100

      DO WHILE (abs(deltapH) > pHTol .AND. count < max_count)
         H         = 10.**(-pHx)
         Denom     = (H*H + K1F*H + K1F*K2F)
         CAlk      = TCx*K1F*(H + 2.*K2F)/Denom
         BAlk      = TBF*KBF/(KBF + H)
         OH        = KWF/H
         PhosTop   = KP1F*KP2F*H + 2.*KP1F*KP2F*KP3F - H*H*H
         PhosBot   = H*H*H + KP1F*H*H + KP1F*KP2F*H + KP1F*KP2F*KP3F
         PAlk      = TPF*PhosTop/PhosBot
         SiAlk     = TSiF*KSiF/(KSiF + H)
         FREEtoTOT = (1 + TSF/KSF) !% pH scale conversion factor
         Hfree     = H/FREEtoTOT !% for H on the total scale
         HSO4      = TSF/(1 + KSF/Hfree) !% since KS is on the free scale
         HF        = TFF/(1 + KFF/Hfree) !% since KF is on the free scale
         Residual  = TAx - CAlk - BAlk - OH - PAlk - SiAlk + Hfree + HSO4 + HF
         Slope     = ln10*(TCx*K1F*H*(H*H + K1F*K2F + 4.*H*K2F)/Denom/Denom &
                  & + BAlk*H/(KBF + H) + OH + H)
         deltapH   = Residual/Slope !% this is Newton's method
         ! to keep the jump from being too big;
         DO WHILE (abs(deltapH) > 1)
            deltapH = deltapH/2
         END DO

         pHx       = pHx + deltapH  !% Is on the same scale as K1 and K2 were calculated...
         count = count + 1
      END DO

   end subroutine Cal_pHfromTATC

   !--------Calculate constants to use for pH----------------------------------------
   subroutine Cal_constants(TempC, Sal, PRE, K0, KS, kF, fH, KB, KW, KP1, KP2, KP3, &
                         & KSi,  K1, K2, TB, TP, TS, TF)
      implicit none

      real(RK),   intent(in)      :: TempC, Sal, PRE
      real(RK),   intent(out)     :: K0, KS, kF, fH, KB, KW, KP1, KP2, KP3, KSi,  K1, K2
      real(RK),   intent(inout)   :: TB, TP, TS, TF

      ! local
      real(RK)                    :: TempK, RT, logTempK, sqrSal, TempK100, Pbar
      real(RK)                    :: lnK0, IonS, lnKS, lnKF, lnKBtop, lnKB, lnKW, lnKP1, &
                                    & lnKP2, lnKP3, lnKSi, lnK1, lnK2, pK1, pK2
      real(RK)                    :: SWStoTOT, FREEtoTOT
      real(RK)                    ::   deltaV, kappa, lnK1fac, lnK2fac, lnKWfac, lnKFfac,&
                                    & lnKSfac, lnKP1fac, lnKP2fac, lnKP3fac, lnKSifac,    &
                                    & K1fac, K2fac, KWfac, KFfac, KSfac, KP1fac, KP2fac,  &
                                    & KP3fac, KSifac, pHfactor, Delta, b, P1atm, FugFac,  &
                                    & VPWP, VPCorrWP, VPSWWP, VPFac, lnKBfac, KBfac

      TempK = TempC + 273.15
      RT    = 83.1451*TempK
      logTempK = log(TempK)
      sqrSal   = sqrt(Sal)
      Pbar     = PRE/10.

      TB    = (0.000232/10.811)*(Sal/1.80655) ! Total Borate, mol/kg-sw
      TF    = (0.000067/18.998)*(Sal/1.80655) ! Total Fluoride, mol/kg-sw
      TS    = (0.14/96.062)*(Sal/1.80655)     ! Total Sulfate, mol/kg-sw

      !---------- CO2 solubility -------------------!
      TempK100 = TempK/100.
      lnK0 = -60.2409 + 93.4517 / TempK100 + 23.3585 * log(TempK100) + Sal* &
         & (0.023517 - 0.023656*TempK100 + 0.0047036 * (TempK100**2.))
      K0   = exp(lnK0)   ! mol/kg-sw/atm

      !----------- Ion Strength --------------------!
      IonS     = 19.924 * Sal / (1000. - 1.005 * Sal)

      !-------- KS for the bisulfate ion -----------!
      lnKS = -4276.1/TempK + 141.328 - 23.093*logTempK + &
         & (-13856./TempK + 324.57 - 47.986*logTempK)*sqrt(IonS) + &
         & (35474./TempK - 771.54 + 114.723*logTempK)*IonS + &
         & (-2698./TempK)*sqrt(IonS)*IonS + (1776./TempK)*(IonS**2)
      KS   = exp(lnKS) * (1. - 0.001005*Sal)  ! mol/kg-sw

      !-------- KF for hydrogen fluoride -----------!
      lnKF = 1590.2/TempK - 12.641 + 1.525*sqrt(IonS)
      KF   = exp(lnKS) * (1. - 0.001005*Sal)  ! mol/kg-sw
      SWStoTOT  = (1 + TS/KS)/(1 + TS/KS + TF/KF)
      FREEtoTOT =  1 + TS/KS

      !---fH, activity coefficient of the H+ ion ---!
      fH   = 1.2948 - 0.002036*TempK + (0.0004607 - 0.000001475*TempK)*(Sal**2)

      !--------KB for boric acid --------------------!
      lnKBtop = -8966.9 - 2890.53*sqrSal - 77.942*Sal + 1.728*sqrSal*Sal - 0.0996*(Sal**2)
      lnKB    = lnKBtop/TempK + 148.0248 + 137.1942*sqrSal + 1.62142*Sal + (-24.4344 - 25.085*sqrSal &
            & - 0.2474*Sal)*logTempK + 0.053105*sqrSal*TempK
      KB      = exp(lnKB)/SWStoTOT          ! convert to SWS pH scale


      !--------- KW for water -----------------------!
      lnKW = 148.9802 - 13847.26/TempK - 23.6521*logTempK + (-5.977 + &
            & 118.67/TempK + 1.0495*logTempK)*sqrSal - 0.01615*Sal
      KW   = exp(lnKW)

      !------KP1, KP2, KP3 for phosphoric acid-------!
      lnKP1 = -4576.752/TempK + 115.54 - 18.453*logTempK + (-106.736/TempK +  &
            & 0.69171)*sqrSal + (-0.65643/TempK - 0.01844)*Sal
      KP1   = exp(lnKP1)
      lnKP2 = -8814.715/TempK + 172.1033 - 27.927*logTempK + (-160.34/TempK + &
            & 1.3566)*sqrSal + (0.37335/TempK - 0.05778)*Sal
      KP2   = exp(lnKP2)
      lnKP3 = -3070.75/TempK - 18.126 + (17.27039/TempK + 2.81197)*sqrSal + &
            & (-44.99486/TempK - 0.09984)*Sal
      KP3   = exp(lnKP3)

      !-------KSi for silicic acid--------------------!
      lnKSi = -8904.2/TempK + 117.4 - 19.334*logTempK + (-458.79/TempK + 3.5913)*sqrt(IonS) + &
            & (188.74/TempK - 1.5998)*IonS + (-12.1652/TempK + 0.07871)*(IonS**2)

      !--------K1 and K2 for carbonic acid------------!
      pK1 = 3670.7/TempK-62.008+9.7944*logTempK-0.0118*Sal+0.000116*(Sal**2)
      K1  = 10.**(-pK1)  ! this is on the SWS pH scale in mol/kg-SW
      pK2 = 1394.7/TempK + 4.777 - 0.0184*Sal + 0.000118*(Sal**2)
      K2  = 10.**(-pK2)

      !============correct constants for pressure=================!
      !------correct K1, k2, kB for pressure----------!
      deltaV  = -25.5 + 0.1271*TempC
      Kappa   = (-3.08 + 0.0877*TempC)/1000.
      lnK1fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

      deltaV  = -15.82 - 0.0219*TempC
      Kappa   = (1.13 - 0.1475*TempC)/1000
      lnK2fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

      deltaV  = -29.48 + 0.1622*TempC - 0.002608*(TempC**2)
      Kappa   = -2.84/1000.
      lnKBfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

      !----correct other constants for pressure--------!

      deltaV  = -20.02 + 0.1119*TempC - 0.001409*(TempC**2)
      Kappa   = (-5.13 + 0.0794*TempC)/1000
      lnKWfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

      deltaV  = -9.78 - 0.009*TempC - 0.000942*(TempC**2)
      Kappa   = (-3.91 + 0.054*TempC)/1000
      lnKFfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

      deltaV  = -18.03 + 0.0466*TempC + 0.000316*(TempC**2)
      Kappa   = (-4.53 + 0.09*TempC)/1000
      lnKSfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

      deltaV  = -14.51 + 0.1211*TempC - 0.000321*(TempC**2)
      Kappa   = (-2.67 + 0.0427*TempC)/1000
      lnKP1fac= (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

      deltaV  = -23.12 + 0.1758*TempC - 0.002647*(TempC**2)
      Kappa   = (-5.15 + 0.09  *TempC)/1000
      lnKP2fac= (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

      deltaV  = -26.57 + 0.202 *TempC - 0.003042*(TempC**2)
      Kappa   = (-4.08 + 0.0714*TempC)/1000
      lnKP3fac= (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

      deltaV  = -29.48 + 0.1622*TempC - 0.002608*(TempC**2)
      Kappa   = -2.84/1000
      lnKSifac= (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT

      K1fac  = exp(lnK1fac);  K1  = K1 *K1fac
      K2fac  = exp(lnK2fac);  K2  = K2 *K2fac
      KWfac  = exp(lnKWfac);  KW  = KW *KWfac
      KBfac  = exp(lnKBfac);  KB  = KB *KBfac
      KFfac  = exp(lnKFfac);  KF  = KF *KFfac
      KSfac  = exp(lnKSfac);  KS  = KS *KSfac
      KP1fac = exp(lnKP1fac); KP1 = KP1*KP1fac
      KP2fac = exp(lnKP2fac); KP2 = KP2*KP2fac
      KP3fac = exp(lnKP3fac); KP3 = KP3*KP3fac
      KSifac = exp(lnKSifac); KSi = KSi*KSifac

      SWStoTOT  = (1 + TS/KS)/(1 + TS/KS + TF/KF)
      FREEtoTOT =  1 + TS/KS

      pHfactor = SWStoTOT

      !----convert from SWS pH scale to chosen scale------!
      K1  = K1* pHfactor; K2  = K2* pHfactor;
      KW  = KW* pHfactor; KB  = KB* pHfactor;
      KP1 = KP1*pHfactor; KP2 = KP2*pHfactor;
      KP3 = KP3*pHfactor; KSi = KSi*pHfactor;

      !-----Calculate Fugacity constants------------------!
      Delta = (57.7 - 0.118*TempK)
      b = -1636.75 + 12.0408*TempK - 0.0327957*(TempK**2) + 3.16528*0.00001*(TempK**3)

      P1atm = 1.01325 ! in bar
      FugFac = exp((b + 2.*Delta)*P1atm/RT)

      !------Calculate VP faction
      VPWP = exp(24.4543 - 67.4509*(100./TempK) - 4.8489*log(TempK/100.))
      VPCorrWP = exp(-0.000544*Sal)
      VPSWWP = VPWP*VPCorrWP
      VPFac = 1. - VPSWWP

   end subroutine Cal_constants

   subroutine save_array(output_unit, array)
      implicit none
      integer, intent(in) :: output_unit
      real(RK), dimension(:), allocatable, intent(in) :: array

      write(output_unit) lbound(array), ubound(array)
      write(output_unit) array
   end subroutine

   subroutine save_integer_array(output_unit, array)
      implicit none
      integer, intent(in) :: output_unit
      integer, dimension(:), allocatable, intent(in) :: array

      write(output_unit) lbound(array), ubound(array)
      write(output_unit) array
   end subroutine

   subroutine save_logical_array(output_unit, array)
      implicit none
      integer, intent(in) :: output_unit
      logical, dimension(:), allocatable, intent(in) :: array

      write(output_unit) lbound(array), ubound(array)
      write(output_unit) array
   end subroutine

   subroutine save_array_pointer(output_unit, array)
      implicit none
      integer, intent(in) :: output_unit
      real(RK), dimension(:), pointer, intent(in) :: array

      write(output_unit) lbound(array), ubound(array)
      write(output_unit) array
   end subroutine

   subroutine read_array(input_unit, array)
      implicit none
      integer, intent(in) :: input_unit
      real(RK), dimension(:), allocatable, intent(inout) :: array
      integer :: array_lbound, array_ubound

      read(input_unit) array_lbound, array_ubound
      read(input_unit) array(array_lbound:array_ubound)
   end subroutine

   subroutine read_integer_array(input_unit, array)
      implicit none
      integer, intent(in) :: input_unit
      integer, dimension(:), allocatable, intent(inout) :: array
      integer :: array_lbound, array_ubound

      read(input_unit) array_lbound, array_ubound
      read(input_unit) array(array_lbound:array_ubound)
   end subroutine

   subroutine read_logical_array(input_unit, array)
      implicit none
      integer, intent(in) :: input_unit
      logical, dimension(:), allocatable, intent(inout) :: array
      integer :: array_lbound, array_ubound

      read(input_unit) array_lbound, array_ubound
      read(input_unit) array(array_lbound:array_ubound)
   end subroutine

   subroutine read_array_pointer(input_unit, array)
      implicit none
      integer, intent(in) :: input_unit
      real(RK), dimension(:), pointer, intent(inout) :: array
      integer :: array_lbound, array_ubound

      read(input_unit) array_lbound, array_ubound
      read(input_unit) array(array_lbound:array_ubound)
   end subroutine


   subroutine save_matrix(output_unit, matrix)
      implicit none
      integer, intent(in) :: output_unit
      real(RK), dimension(:, :), allocatable, intent(in) :: matrix

      write(output_unit) lbound(matrix, 1), ubound(matrix, 1), lbound(matrix, 2), ubound(matrix, 2)
      write(output_unit) matrix
   end subroutine

   subroutine save_matrix_pointer(output_unit, matrix)
      implicit none
      integer, intent(in) :: output_unit
      real(RK), dimension(:, :), pointer, intent(in) :: matrix

      write(output_unit) lbound(matrix, 1), ubound(matrix, 1), lbound(matrix, 2), ubound(matrix, 2)
      write(output_unit) matrix
   end subroutine

   subroutine read_matrix(input_unit, matrix)
      implicit none
      integer, intent(in) :: input_unit
      real(RK), dimension(:, :), allocatable, intent(inout) :: matrix
      integer :: matrix_lbound_1, matrix_ubound_1, matrix_lbound_2, matrix_ubound_2

      read(input_unit) matrix_lbound_1, matrix_ubound_1, matrix_lbound_2, matrix_ubound_2
      if (.not. allocated(matrix)) then
         allocate (matrix(matrix_lbound_1:matrix_ubound_1, matrix_lbound_2:matrix_ubound_2))
      end if
      read(input_unit) matrix(matrix_lbound_1:matrix_ubound_1, matrix_lbound_2:matrix_ubound_2)
   end subroutine

   subroutine read_matrix_pointer(input_unit, matrix)
      implicit none
      integer, intent(in) :: input_unit
      real(RK), dimension(:, :), pointer, intent(inout) :: matrix
      integer :: matrix_lbound_1, matrix_ubound_1, matrix_lbound_2, matrix_ubound_2

      read(input_unit) matrix_lbound_1, matrix_ubound_1, matrix_lbound_2, matrix_ubound_2
      read(input_unit) matrix(matrix_lbound_1:matrix_ubound_1, matrix_lbound_2:matrix_ubound_2)
   end subroutine

   pure real(RK) function datum(start_datum, simulation_time)
      implicit none
      real(RK), intent(in) :: start_datum
      integer(8), dimension(2), intent(in) :: simulation_time

      ! Simulation time is a tuple of 2 integers (days, seconds)
      datum = start_datum + real(simulation_time(1), RK) + real(simulation_time(2), RK) / SECONDS_PER_DAY
   end function

end module utilities
