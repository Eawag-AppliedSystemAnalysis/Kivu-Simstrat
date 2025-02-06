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
!     |  Advection module
!     |  - Based on already read in/outflows, calculates Advection
!     |  - Might grow or shrink grid (methods merge/add_box)
!<    +---------------------------------------------------------------+

module strat_advection
   use strat_kinds
   use strat_simdata
   use strat_consts
   use strat_grid
   use simstrat_aed2
   use utilities
   implicit none
   private

   type, public :: AdvectionModule
      class(ModelConfig), pointer :: cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param

   contains
      procedure, pass :: init => advection_init
      procedure, pass :: update => advection_update
      procedure, pass :: merge_box => advection_merge_box
      procedure, pass :: add_box => advection_add_box
   end type

contains

   subroutine advection_init(self, state, model_config, model_param, grid)
      implicit none
      class(AdvectionModule) :: self
      class(ModelState) :: state
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(ModelParam), target :: model_param

      ! Local variables
      integer :: i

      self%cfg => model_config
      self%param => model_param
      self%grid => grid

      if (self%cfg%couple_aed2) then
         do i = 1, state%n_AED2
            select case(trim(state%AED2_names(i)))
            case('CAR_pH')
               state%n_pH = i
            end select
         end do
      end if

   end subroutine

   ! A lot of code that is hard to test - might be refactored in the future
   subroutine advection_update(self, state)
      implicit none
      class(AdvectionModule) :: self
      class(ModelState) :: state

      real(RK) :: top_z, top_h
      real(RK) :: top
      real(RK) :: dh, dh_i(1:2), h_div_2, h_mult_2 ! depth differences
      real(RK) :: dU(self%grid%nz_grid), dV(self%grid%nz_grid), dTemp(self%grid%nz_grid), dS(self%grid%nz_grid)
      real(RK) :: dTr(self%grid%nz_grid), dHO(self%grid%nz_grid), dD(self%grid%nz_grid), dLA(self%grid%nz_grid)
      real(RK) :: dHe(self%grid%nz_grid), dNe(self%grid%nz_grid), dAr(self%grid%nz_grid), dKr(self%grid%nz_grid)
      real(RK) :: dt_i(1:2) ! first and second time step
      real(RK) :: AreaFactor_adv(1:self%grid%nz_grid)
      integer :: i, t_i

      associate (grid=>self%grid, &
                 nz_occupied=>self%grid%nz_occupied, &
                 dt=>state%dt, &
                 h=>self%grid%h, &
                 Q_vert=>state%Q_vert, &
                 AED2_state=>state%AED2_state, &
                 ubnd_vol=>self%grid%ubnd_vol, &
                 ubnd_fce=>self%grid%ubnd_fce)

         !Depth difference compared to previous timestep
         top_z = grid%z_face(ubnd_fce)
         top_h = grid%h(ubnd_vol)

         dh = state%Q_vert(ubnd_fce)/grid%Az(ubnd_fce)*state%dt
         h_div_2 = 0.5_RK*h(nz_occupied - 1) ! Take second highest box since the top box might not be at the full height
         h_mult_2 = 2_RK*h(nz_occupied - 1)

         ! Calculate timestep splitting
         !Split timestep depending on situation
         if (dh == 0.) then ! If volume does not change, take one normal time step
            dt_i(1) = dt
         else if (top_z == grid%max_depth) then ! If we are already at the maximum lake level
            dt_i(1) = dt
         else if ((dh + top_z) >= grid%max_depth) then ! If the full timestep would lead to a lake level higher than maximum allowed lake level, split the timestep.
            dt_i(1) = (grid%max_depth - top_z)/dh*dt
         else if (((dh + top_h) > h_div_2) .and. & ! If top box>0.5*lower box and <2*lower box, take one time step
                  ((dh + top_h) < h_mult_2)) then
            dt_i(1) = dt
         else if ((dh + top_h) <= h_div_2) then ! If top box<=0.5*lower box, first step until top box=0.5*lower box
            dt_i(1) = abs((top_h - h_div_2)/dh)*dt
         else ! If top box>=2*lower box, first step until top box = 2*lower box
            dt_i(1) = abs((2*h(nz_occupied - 1) - top_h)/dh)*dt
         end if
         dt_i(2) = dt - dt_i(1) ! Rest of timestep

         ! FB 2016/2019: Revisions
         do t_i = 1, 2 !First and (if needed) second timestep
            AreaFactor_adv(1:nz_occupied) = dt_i(t_i)/((grid%Az(1:nz_occupied) + grid%Az(2:nz_occupied+1))/2*grid%h(1:nz_occupied)) ! Area factor for dt(t_i)
            dh_i(t_i) = dh*dt_i(t_i)/dt ! Depth difference for dt(t_i)

            ! Update Simstrat variables U, V, T and S
            call do_update_statvars(self, state, AreaFactor_adv(1:nz_occupied), dh_i(t_i))

            ! Update AED2 variables
            if(self%cfg%couple_aed2) call do_update_statvars_AED2(self, state, AreaFactor_adv(1:nz_occupied), dh_i(t_i))

            ! Adjust boxes (Horrible if/else construction - replace!)
            if (t_i == 1) then
               if (dh == 0) then ! If volume does not change, return
                  return
               else if ((dh + top_z) >= grid%max_depth) then ! If surface level reached
                  call grid%modify_top_box(grid%max_depth - top_z)
                  return
               else if (((dh_i(t_i) + top_h) > h_div_2) .and. &
                        ((dh_i(t_i) + top_h) < (h_mult_2))) then ! and top box<2*lower box
                  call grid%modify_top_box(dh_i(t_i))
                  return
               else if (t_i == 1 .and. (dh + top_h) <= h_div_2) then ! If top box<=0.5*lower box, merge 2 boxes
                  call self%merge_box(state, dh_i(t_i))
               else if (t_i == 1 .and. (dh + top_h) >= h_mult_2) then ! If top box>=2*lower box, add one box
                  call self%add_box(state, dh_i(t_i))
               end if ! dh==0
            end if

         end do !end do t_i=1,2
      end associate
   end subroutine

   subroutine do_update_statvars(self, state, AreaFactor_adv, dh)
      ! Arguments
      class(AdvectionModule) :: self
      class(ModelState) :: state
      real(RK), dimension(:) :: AreaFactor_adv
      real(RK) :: dh

      ! Local variables
      integer :: i, top
      real(RK) :: dU(self%grid%nz_grid), dV(self%grid%nz_grid), dTemp(self%grid%nz_grid), dS(self%grid%nz_grid)
      real(RK) :: dTr(self%grid%nz_grid), dHO(self%grid%nz_grid), dD(self%grid%nz_grid), dLA(self%grid%nz_grid)
      real(RK) :: dHe(self%grid%nz_grid), dNe(self%grid%nz_grid), dAr(self%grid%nz_grid), dKr(self%grid%nz_grid)
      real(RK) :: Tr_evap, Tr_rain, Tr_p, Tr_lake, alpha_Tr, f, T_abs, hum, E
      real(RK) :: HO_evap, HO_rain, HO_p, HO_lake, alpha_18O, delta_eps_18O
      real(RK) :: D_evap, D_rain, D_p, D_lake, alpha_D, delta_eps_D
      real(RK) :: Vap_atm, Vap_lake
      integer :: outflow_above, outflow_below

      associate(ubnd_vol => self%grid%ubnd_vol, &
         Q_vert => state%Q_vert, &
         h => self%grid%h)

            ! Calculate changes
            do i = 1, ubnd_vol
               ! For the top-most cell, if Q_vert at the upper face is positive, there is still no outflow (the cell is simply growing, but this is done elsewhere)
               if ((i == ubnd_vol) .and. Q_vert(i + 1) > 0) then
                  top = 0
               else
                  top = 1
               end if

               ! If Q_vert at the upper face of cell i is positive, then there is outflow to the cell above
               if (Q_vert(i + 1) > 0) then
                     outflow_above = 1
               else 
                     outflow_above = 0
               end if

               ! If Q_vert at the lower face of cell i is negative, then there is outflow to the cell below
               if (Q_vert(i) < 0) then
                     outflow_below = 1
               else
                     outflow_below = 0
               end if

               ! Calculate advective flow out of cell (thus negative sign in the front) i to the cells above and below
               dU(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%U(i)
               dV(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%V(i)
               dTemp(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%T(i)
               dS(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%S(i)

               dTr(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%Tr(i)
               dHO(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%heavy_oxygen(i)
               dD(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%deuterium(i)
               dLA(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%light_ar(i)

               dHe(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%He(i)
               dNe(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%Ne(i)
               dAr(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%Ar(i)
               dKr(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%Kr(i)

               ! Calculate the advective flow into cell i from below
               if (i > 1 .and. Q_vert(i ) > 0) then
                  dU(i) = dU(i) + Q_vert(i)*state%U(i - 1)
                  dV(i) = dV(i) + Q_vert(i)*state%V(i - 1)
                  dTemp(i) = dTemp(i) + Q_vert(i)*state%T(i - 1)
                  dS(i) = dS(i) + Q_vert(i)*state%S(i - 1)

                  dTr(i) = dTr(i) + Q_vert(i)*state%Tr(i - 1)
                  dHO(i) = dHO(i) + Q_vert(i)*state%heavy_oxygen(i - 1)
                  dD(i) = dD(i) + Q_vert(i)*state%deuterium(i - 1)
                  dLA(i) = dLA(i) + Q_vert(i)*state%light_ar(i - 1)

                  dHe(i) = dHe(i) + Q_vert(i)*state%He(i - 1)
                  dNe(i) = dNe(i) + Q_vert(i)*state%Ne(i - 1)
                  dAr(i) = dAr(i) + Q_vert(i)*state%Ar(i - 1)
                  dKr(i) = dKr(i) + Q_vert(i)*state%Kr(i - 1)
               end if

               ! Calculate the advective flow into cell i from above (- sign in front because Q_vert is negative if there is inflow)
               if (i < ubnd_vol .and. Q_vert(i + 1) < 0) then
                  dU(i) = dU(i) - Q_vert(i + 1)*state%U(i + 1)
                  dV(i) = dV(i) - Q_vert(i + 1)*state%V(i + 1)
                  dTemp(i) = dTemp(i) - Q_vert(i + 1)*state%T(i + 1)
                  dS(i) = dS(i) - Q_vert(i + 1)*state%S(i + 1)

                  dTr(i) = dTr(i) - Q_vert(i + 1)*state%Tr(i + 1)
                  dHO(i) = dHO(i) - Q_vert(i + 1)*state%heavy_oxygen(i + 1)
                  dD(i) = dD(i) - Q_vert(i + 1)*state%deuterium(i + 1)
                  dLA(i) = dLA(i) - Q_vert(i + 1)*state%light_ar(i + 1)

                  dHe(i) = dHe(i) - Q_vert(i + 1)*state%He(i + 1)
                  dNe(i) = dNe(i) - Q_vert(i + 1)*state%Ne(i + 1)
                  dAr(i) = dAr(i) - Q_vert(i + 1)*state%Ar(i + 1)
                  dKr(i) = dKr(i) - Q_vert(i + 1)*state%Kr(i + 1)
               end if
            end do

            ! Compute change
            ! dT = dT(vertical advection) + dT(inflow) + dT(outflow), units: °C*m^3/s
            dTemp(1:ubnd_vol) = dTemp(1:ubnd_vol) + state%Q_inp(3, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%T(1:ubnd_vol)
            ! dS = dS(vertical advection) + dS(inflow) + dS(outflow), units: ‰*m^3/s
            dS(1:ubnd_vol) = dS(1:ubnd_vol) + state%Q_inp(4, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%S(1:ubnd_vol)
            ! dTr = dTr(vertical advection) + dTr(inflow) + dTr(outflow), units: TU*m^3/s
            dTr(1:ubnd_vol) = dTr(1:ubnd_vol) + state%Q_inp(5, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%Tr(1:ubnd_vol)
            ! dHO = dHO(vertical advection) + dHO(inflow) + dHO(outflow), units: [-]*m^3/s
            dHO(1:ubnd_vol) = dHO(1:ubnd_vol) + state%Q_inp(6, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%heavy_oxygen(1:ubnd_vol)
            ! dD = dD(vertical advection) + dD(inflow) + dD(outflow), units: [-]*m^3/s
            dD(1:ubnd_vol) = dD(1:ubnd_vol) + state%Q_inp(7, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%deuterium(1:ubnd_vol)
            ! dLA = dLA(vertical advection) + dLA(inflow) + dLA(outflow), units: [-]*m^3/s
            dLA(1:ubnd_vol) = dLA(1:ubnd_vol) + state%Q_inp(8, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%light_ar(1:ubnd_vol)

            ! dNG = dNG(vertical advection) + dNG(inflow) + dNG(outflow), units: [ccSTP/g]*m^3/s
            dHe(1:ubnd_vol) = dHe(1:ubnd_vol) + state%Q_inp(9, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%He(1:ubnd_vol)
            dNe(1:ubnd_vol) = dNe(1:ubnd_vol) + state%Q_inp(10, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%Ne(1:ubnd_vol)
            dAr(1:ubnd_vol) = dAr(1:ubnd_vol) + state%Q_inp(11, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%Ar(1:ubnd_vol)
            dKr(1:ubnd_vol) = dKr(1:ubnd_vol) + state%Q_inp(12, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%Kr(1:ubnd_vol)


            ! Add change to the state variable
            state%U(1:ubnd_vol) = state%U(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dU(1:ubnd_vol)
            state%V(1:ubnd_vol) = state%V(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dV(1:ubnd_vol)
            state%T(1:ubnd_vol) = state%T(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dTemp(1:ubnd_vol)
            state%S(1:ubnd_vol) = state%S(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dS(1:ubnd_vol)

            state%Tr(1:ubnd_vol) = state%Tr(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dTr(1:ubnd_vol)
            state%heavy_oxygen(1:ubnd_vol) = state%heavy_oxygen(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dHO(1:ubnd_vol)
            state%deuterium(1:ubnd_vol) = state%deuterium(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dD(1:ubnd_vol)
            state%light_ar(1:ubnd_vol) = state%light_ar(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dLA(1:ubnd_vol)

            state%He(1:ubnd_vol) = state%He(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dHe(1:ubnd_vol)
            state%Ne(1:ubnd_vol) = state%Ne(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dNe(1:ubnd_vol)
            state%Ar(1:ubnd_vol) = state%Ar(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dAr(1:ubnd_vol)
            state%Kr(1:ubnd_vol) = state%Kr(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dKr(1:ubnd_vol)

            ! Preparation
            T_abs = state%T(ubnd_vol) + 273.15
            f = state%Vap_atm/1013/exp(24.4543 - 67.4509*100/T_abs - 4.8489*log(T_abs/100)) ! Weiss, 1980

            Vap_atm = 10**((0.7859_RK + 0.03477_RK*state%T_atm)/(1 + 0.00412_RK*state%T_atm))
            Vap_atm = Vap_atm*(1 + 1e-6_RK*self%param%p_air*(4.5_RK + 0.00006_RK*state%T_atm**2))

            Vap_lake = 10**((0.7859_RK + 0.03477_RK*state%T(ubnd_vol))/(1 + 0.00412_RK*state%T(ubnd_vol)))
            Vap_lake = Vap_lake*(1 + 1e-6_RK*self%param%p_air*(4.5_RK + 0.00006_RK*state%T(ubnd_vol)**2))

            hum = f*Vap_atm/Vap_lake

            E = 115_RK  ! Computed from water balance in Lake Kivu (Muvundja et al, 2014)

            ! Change due to rain/evaporation for Tritium
            alpha_Tr = 0.89_RK ! according to W. Aeschbach (kinetic contribution is negligible)

            Tr_p = state%Q_inp(5,ubnd_vol)/state%Q_inp(1,ubnd_vol) ! Assuming surface inflow as same Tr as rain
            Tr_lake = state%Tr(ubnd_vol)

            Tr_evap = E*alpha_Tr*(hum*Tr_p - Tr_lake)/(1.0_RK - hum)
            Tr_rain = E*Tr_p ! E=R for Lake Kivu (Muvundja et al., 2014)

            state%Tr(ubnd_vol) = state%Tr(ubnd_vol) + AreaFactor_adv(ubnd_vol)*(Tr_evap + Tr_rain)

            ! Change due to rain/evaporation for 18O/D
            HO_p = state%Q_inp(6,ubnd_vol)/state%Q_inp(1,ubnd_vol)   ! Precipitation is assumed to be the same as surface inflow
            HO_lake = state%heavy_oxygen(ubnd_vol)

            D_p = state%Q_inp(7,ubnd_vol)/state%Q_inp(1,ubnd_vol)
            D_lake = state%deuterium(ubnd_vol)

            ! Contributions of equilibrium and kinetic fractionation
            alpha_18O = 1/exp(-0.00207_RK - 0.4156_RK/T_abs + 1137.0_RK/T_abs**2) ! Majoube 1971 in Gonfiantini 1986
            alpha_D = 1/exp(0.05261_RK - 76.248_RK/T_abs + 24844.0_RK/T_abs**2) ! Majoube 1971 in Gonfiantini 1986
            delta_eps_18O = 14.2*(1.0_RK - hum)/1000 ! Kinetic contribution according to Gonfiantini, 1986
            delta_eps_D = 12.5*(1.0_RK - hum)/1000 ! Kinetic contribution according to Gonfiantini, 1986

            HO_evap = E*alpha_18O*(hum*HO_p - HO_lake)/(1.0_RK - hum + delta_eps_18O)
            HO_rain = E*HO_p

            D_evap = E*alpha_D*(hum*D_p - D_lake)/(1.0_RK - hum + delta_eps_D)
            D_rain = E*D_p

            state%heavy_oxygen(ubnd_vol) = state%heavy_oxygen(ubnd_vol) + AreaFactor_adv(ubnd_vol)*(HO_evap + HO_rain)
            state%deuterium(ubnd_vol) = state%deuterium(ubnd_vol) + AreaFactor_adv(ubnd_vol)*(D_evap + D_rain)

            ! Variation of variables due to change in volume
            state%U(ubnd_vol) = state%U(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%V(ubnd_vol) = state%V(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%T(ubnd_vol) = state%T(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%S(ubnd_vol) = state%S(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%Tr(ubnd_vol) = state%Tr(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%heavy_oxygen(ubnd_vol) = state%heavy_oxygen(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%deuterium(ubnd_vol) = state%deuterium(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%light_ar(ubnd_vol) = 1 ! Always equal to atmosphere (actually a boundary condition)

            state%He(ubnd_vol) = 3.83e-8!state%He(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%Ne(ubnd_vol) = 1.55e-7!state%Ne(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%Ar(ubnd_vol) = 8.33e-7!state%Ar(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
            state%Kr(ubnd_vol) = 5.37e-8!state%Kr(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh)
      end associate
   end subroutine

   subroutine do_update_statvars_AED2(self, state, AreaFactor_adv, dh)
      ! Arguments
      class(AdvectionModule) :: self
      class(ModelState) :: state
      real(RK), dimension(:) :: AreaFactor_adv(1:self%grid%ubnd_vol)
      real(RK) :: dh

      ! Local variables
      integer :: i, top, outflow_above, outflow_below
      real(RK) :: dAED2(self%grid%nz_grid, state%n_AED2)

      associate(ubnd_vol => self%grid%ubnd_vol, &
         Q_vert => state%Q_vert, &
         AED2_state => state%AED2_state, &
         h => self%grid%h)

         ! Calculate changes
         do i = 1, ubnd_vol
            ! For the top-most cell, if Q_vert at the upper face is positive, there is still no outflow (the cell is simply growing, but this is done elsewhere)
            if ((i == ubnd_vol) .and. Q_vert(i + 1) > 0) then
               top = 0
            else
               top = 1
            end if

            ! If Q_vert at the upper face of cell i is positive, then there is outflow to the cell above
            if (Q_vert(i + 1) > 0) then
                  outflow_above = 1
            else 
                  outflow_above = 0
            end if

            ! If Q_vert at the lower face of cell i is negative, then there is outflow to the cell below
            if (Q_vert(i) < 0) then
                  outflow_below = 1
            else
                  outflow_below = 0
            end if

            ! Calculate advective flow out of cell (thus negative sign in the front) i to the cells above and below
            dAED2(i,:) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*AED2_state(i,:)

            ! Calculate the advective flow into cell i from below
            if (i > 1 .and. Q_vert(i) > 0) then
               dAED2(i,:) = dAED2(i,:) + Q_vert(i)*AED2_state(i - 1,:)
            end if
            ! Calculate the advective flow into cell i from above (- sign in front because Q_vert is negative if there is inflow)
            if (i < ubnd_vol .and. Q_vert(i + 1) < 0) then
               dAED2(i,:) = dAED2(i,:) - Q_vert(i + 1)*AED2_state(i + 1,:)
            end if
         end do

         ! Add change to state variables
         ! dAED2 = dAED2(vertical advection) + dAED2(inflow) + Outflow(negative)*AED2, units: C*m^3/s

         ! Add change to the state variable
         do i=1,state%n_AED2
            dAED2(1:ubnd_vol,i) = dAED2(1:ubnd_vol,i) + state%Q_inp(n_simstrat + i, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*AED2_state(1:ubnd_vol,i)
            AED2_state(1:ubnd_vol,i) = AED2_state(1:ubnd_vol,i) + AreaFactor_adv(1:ubnd_vol)*dAED2(1:ubnd_vol,i)
         end do

         ! Variation of variables due to change in volume
         AED2_state(ubnd_vol,:) = AED2_state(ubnd_vol,:)*h(ubnd_vol)/(h(ubnd_vol) + dh)

         ! Transform [H] back to pH
         if(self%cfg%couple_aed2) then
            AED2_state(:,state%n_pH) = -log10(AED2_state(:,state%n_pH))
         end if

      end associate
   end subroutine

   ! Merges two boxes
   ! - Takes care of calculating the new state variable for this box
   ! - Calls grid methods to modify grid spacing etc
   subroutine advection_merge_box(self, state, dh)
      implicit none
      class(AdvectionModule) :: self
      class(ModelState) :: state
      real(RK) :: dh
      real(RK) :: w_a, w_b
      associate (ubnd_fce=>self%grid%ubnd_fce, ubnd_vol=>self%grid%ubnd_vol, AED2_state=>state%AED2_state)

         ! New values of the state variables are weighted averages
         !determine weighting an normalization connstant
         w_a = 0.5_RK*self%grid%Az(ubnd_fce)
         w_b = self%grid%Az(ubnd_fce - 1)

         ! shrink grid by one (this also updates ubnd_fce/vol)
         call self%grid%shrink(dh)

         ! update quantities in new top box (based on former top box and current value)
         state%U(ubnd_vol) = (w_a*state%U(ubnd_vol + 1) + w_b*state%U(ubnd_vol))/(w_a + w_b)
         state%V(ubnd_vol) = (w_a*state%V(ubnd_vol + 1) + w_b*state%V(ubnd_vol))/(w_a + w_b)
         state%T(ubnd_vol) = (w_a*state%T(ubnd_vol + 1) + w_b*state%T(ubnd_vol))/(w_a + w_b)
         state%S(ubnd_vol) = (w_a*state%S(ubnd_vol + 1) + w_b*state%S(ubnd_vol))/(w_a + w_b)
         state%Tr(ubnd_vol) = (w_a*state%Tr(ubnd_vol + 1) + w_b*state%Tr(ubnd_vol))/(w_a + w_b)
         state%heavy_oxygen(ubnd_vol) = (w_a*state%heavy_oxygen(ubnd_vol + 1) + w_b*state%heavy_oxygen(ubnd_vol))/(w_a + w_b)
         state%deuterium(ubnd_vol) = (w_a*state%deuterium(ubnd_vol + 1) + w_b*state%deuterium(ubnd_vol))/(w_a + w_b)
         state%light_ar(ubnd_vol) = (w_a*state%light_ar(ubnd_vol + 1) + w_b*state%light_ar(ubnd_vol))/(w_a + w_b)

         state%He(ubnd_vol) = (w_a*state%He(ubnd_vol + 1) + w_b*state%He(ubnd_vol))/(w_a + w_b)
         state%Ne(ubnd_vol) = (w_a*state%Ne(ubnd_vol + 1) + w_b*state%Ne(ubnd_vol))/(w_a + w_b)
         state%Ar(ubnd_vol) = (w_a*state%Ar(ubnd_vol + 1) + w_b*state%Ar(ubnd_vol))/(w_a + w_b)
         state%Kr(ubnd_vol) = (w_a*state%Kr(ubnd_vol + 1) + w_b*state%Kr(ubnd_vol))/(w_a + w_b)

         state%k(ubnd_fce) = (w_a*state%k(ubnd_fce + 1) + w_b*state%k(ubnd_fce))/(w_a + w_b)
         state%eps(ubnd_fce) = (w_a*state%eps(ubnd_fce + 1) + w_b*state%eps(ubnd_fce))/(w_a + w_b)
         state%Q_vert(ubnd_fce) = (w_a*state%Q_vert(ubnd_fce + 1) + w_b*state%Q_vert(ubnd_fce))/(w_a + w_b)

         ! AED2
         if (self%cfg%couple_AED2) AED2_state(ubnd_vol,:) = (w_a*AED2_state(ubnd_vol + 1,:) + w_b*AED2_state(ubnd_vol,:))/(w_a + w_b)

         ! update area factors
         call self%grid%update_area_factors()

      end associate
   end subroutine

   ! Adds a new box
   ! - Takes care of calculating the new state variable for this box
   ! - Calls grid methods to modify grid spacing etc
   subroutine advection_add_box(self, state, dh)
      implicit none
      class(AdvectionModule) :: self
      class(ModelState) :: state
      real(RK) :: dh
      associate (ubnd_fce=>self%grid%ubnd_fce, ubnd_vol=>self%grid%ubnd_vol)

         ! extend grid by one (also updates ubnd_vol etc)
         call self%grid%grow(dh)

         ! Update quantities in new grid element
         state%U(ubnd_vol) = state%U(ubnd_vol - 1)
         state%V(ubnd_vol) = state%V(ubnd_vol - 1)
         state%T(ubnd_vol) = state%T(ubnd_vol - 1)
         state%S(ubnd_vol) = state%S(ubnd_vol - 1)
         state%Tr(ubnd_vol) = state%Tr(ubnd_vol - 1)
         state%heavy_oxygen(ubnd_vol) = state%heavy_oxygen(ubnd_vol - 1)
         state%deuterium(ubnd_vol) = state%deuterium(ubnd_vol - 1)
         state%light_ar(ubnd_vol) = state%light_ar(ubnd_vol - 1)
         state%Q_vert(ubnd_fce) = state%Q_vert(ubnd_fce - 1) ! Vertical discharge of new box

         state%He(ubnd_vol) = state%He(ubnd_vol - 1)
         state%Ne(ubnd_vol) = state%Ne(ubnd_vol - 1)
         state%Ar(ubnd_vol) = state%Ar(ubnd_vol - 1)
         state%Kr(ubnd_vol) = state%Kr(ubnd_vol - 1)

         state%k(ubnd_fce) = state%k(ubnd_fce - 1)
         state%eps(ubnd_fce) = state%eps(ubnd_fce - 1)

         if (self%cfg%couple_AED2) state%AED2_state(ubnd_vol,:) = state%AED2_state(ubnd_vol - 1,:)

         call self%grid%update_area_factors()

      end associate
   end subroutine

end module
