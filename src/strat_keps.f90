   !<    +---------------------------------------------------------------+
!     |  ModelVariable implementation for K and EPS variable
!<    +---------------------------------------------------------------+

module strat_keps
   use strat_kinds
   use strat_consts
   use strat_simdata
   use strat_statevar
   use strat_grid
   use strat_solver
   implicit none
   private

   ! Class for K variable
   type, extends(ModelVariable), public :: KModelVar
   contains
      procedure, pass(self), public :: calc_terms => k_var_calc_terms
      procedure, pass(self), public :: post_solve => k_var_post_solve
   end type

   ! Class for Eps variable
   type, extends(ModelVariable), public :: EpsModelVar
   contains
      procedure, pass(self), public :: calc_terms => eps_var_calc_terms
      procedure, pass(self), public :: post_solve => eps_var_post_solve
   end type

contains

  ! Calculates terms for linear system for a K variable
   subroutine k_var_calc_terms(self, state, param, sources, boundaries)
      class(KModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), dimension(:) ::  sources, boundaries
      real(RK) :: rhs_0, rhs_ubnd
      real(RK) :: pminus(2:self%grid%ubnd_fce - 1), pplus(2:self%grid%ubnd_fce - 1), Prod, Buoy, Diss
      integer :: i
      associate (grid=>self%grid, &
                 ubnd_fce=>self%grid%ubnd_fce, &
                 ubnd_vol=>self%grid%ubnd_vol)

         !!!!!!!! Precalculations !!!!!!!!
         sources = 0.0_RK
         boundaries = 0.0_RK
         state%ko(1:ubnd_fce) = state%k(1:ubnd_fce) ! ko = TKE at old time step

         ! Diffusivity for k is located on volume grid
         state%avh(2:ubnd_vol - 1) = 0.5_RK/sig_k*(state%num(2:ubnd_fce - 2) + state%num(3:ubnd_fce - 1))

         if (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) then
            state%avh(1) = 0.0_RK
            state%avh(ubnd_vol) = 0.0_RK
         else
            state%avh(1) = 2*state%u_taub**4/(state%eps(1) + state%eps(2)) ! = 0 for no shear stress
            state%avh(ubnd_vol) = 2*state%u_taus**4/(state%eps(ubnd_fce) + state%eps(ubnd_fce - 1)) ! = 0 for no shear stress
         end if

         do i = 2, ubnd_fce - 1
            Prod = state%P(i) + state%P_Seiche(i) ! Add seiche energy
            Buoy = state%B(i)
            Diss = state%eps(i)
            if (Prod + Buoy > 0) then
               pplus(i) = Prod + Buoy
               pminus(i) = Diss
            else
               pplus(i) = Prod
               pminus(i) = Diss - Buoy
            end if
         end do

         !!!!!!!! Define sources !!!!!!!!
         sources(2:ubnd_fce - 1) = pplus(2:ubnd_fce - 1)

         !!!!!! Define boundary conditions !!!!
         boundaries(2:ubnd_fce - 1) = pminus(2:ubnd_fce - 1)/state%k(2:ubnd_fce - 1)

         if (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) then
            ! K(0) and K(ubnd_fce) are assigned in the post processing function
         else ! no fluxes, unity A-matrix + condition on RHS
            rhs_0 = state%u_taub**2/sqrt(state%cm0*state%cde)
            rhs_ubnd = state%u_taus**2/sqrt(state%cm0*state%cde)

            ! Trick to have rhs(1) = rhs_0 and rhs(ubnd_fce) = rhs_ubnd
            ! in discretization, rhs is calculated as rhs(1) = var(1) + sources(1)*dt
            ! Given the equation below, the following results:
            ! rhs(1) = var(1) + (-var(1)/dt + rhs_0/dt)*dt = var(1) -var(1) + rhs_0 = rhs_0
            sources(1) = -(self%var(1)/state%dt) + (rhs_0/state%dt)
            sources(ubnd_fce) = -(self%var(ubnd_fce)/state%dt) + (rhs_ubnd/state%dt)
         end if
      end associate
   end subroutine

   ! Updates variable boundaries after solve has been executed
   subroutine k_var_post_solve(self, state)
      class(KModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state

      integer :: i
      associate (grid=>self%grid, &
                 ubnd_fce=>self%grid%ubnd_fce, &
                 ubnd_vol=>self%grid%ubnd_vol)

         if (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) then
            ! Define TKE at boundary (no flux)
            self%var(1) = self%var(2)
            self%var(self%grid%ubnd_fce) = self%var(self%grid%ubnd_fce - 1)
         end if

         ! check lower limit of k
         do i = 1, ubnd_fce
            if (self%var(i) < k_min) self%var(i) = k_min ! Lower limit of TKE
         end do

      end associate
   end subroutine

   ! Calculates terms for linear system for a EPS variable
   subroutine eps_var_calc_terms(self, state, param, sources, boundaries)
      class(EpsModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), dimension(:) ::  sources, boundaries
      real(RK) :: cee3, rhs_0, rhs_ubnd
      real(RK) :: flux(1:self%grid%ubnd_fce)
      real(RK) :: pminus(2:self%grid%ubnd_fce - 1), pplus(2:self%grid%ubnd_fce - 1), Prod, Buoy, Diss

      integer :: i
      associate (grid=>self%grid, &
                 ubnd_fce=>self%grid%ubnd_fce, &
                 ubnd_vol=>self%grid%ubnd_vol)

         !!!!!!!! Precalculations !!!!!!!!
         sources = 0.0_RK
         boundaries = 0.0_RK

         ! Diffusivity for eps is located on volume grid
         state%avh(1:ubnd_vol) = 0.5_RK/sig_e*(state%num(1:ubnd_fce - 1) + state%num(2:ubnd_fce))

         if (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) then
            flux(1) = state%avh(1)*(state%cde*((state%ko(2))**1.5_RK))/(kappa*(K_s + 0.5_RK*grid%h(1)))**2
            flux(ubnd_fce) = state%avh(ubnd_vol)*(state%cde*((state%ko(ubnd_fce - 1))**1.5_RK))/(kappa*(z0 + 0.5_RK*grid%h(ubnd_vol)))**2

            do i = 2, ubnd_fce - 1
              flux(i) = state%num(i)/sig_e*(state%cde*((state%ko(i))**1.5_RK)/(kappa*(z0 + 0.25_RK*(grid%h(i - 1) + grid%h(i))))**2)
            end do

            state%avh(1) = 0.0_RK
            state%avh(ubnd_vol) = 0.0_RK
         else
            state%avh(1) = 2*state%u_taub**4/sig_e/(state%eps(1) + state%eps(2)) ! = 0 for no shear stress
            state%avh(ubnd_vol) = 2*state%u_taus**4/sig_e/(state%eps(ubnd_fce) + state%eps(ubnd_fce - 1)) ! = 0 for no shear stress
         end if

         do i = 2, ubnd_fce - 1
            if (state%B(i) > 0) then
               cee3 = 1.
            else
               cee3 = ce3
            end if
            if(grid%z_volume(i) < 365) then 
               Prod = ce1*state%eps(i)/state%ko(i)*state%P(i)
            else
               Prod = ce1*state%eps(i)/state%ko(i)*(state%P(i) + state%P_Seiche(i)) ! New code plus seiche
            end if
            Buoy = cee3*state%eps(i)/state%ko(i)*state%B(i)
            Diss = ce2*state%eps(i)*state%eps(i)/state%ko(i)

            if (Prod + Buoy > 0) then
               pplus(i) = Prod + Buoy
               pminus(i) = Diss
            else
               pplus(i) = Prod
               pminus(i) = Diss - Buoy
            end if
         end do

         !!!!!!!! Define sources !!!!!!!!
         sources(2:ubnd_fce - 1) = pplus(2:ubnd_fce - 1)

         !!!!!! Define boundary conditions !!!!
         boundaries(2:ubnd_fce - 1) = pminus(2:ubnd_fce - 1)/state%eps(2:ubnd_fce - 1)
         if (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) then
            sources(2:ubnd_fce - 1) = sources(2:ubnd_fce - 1) + flux(2:ubnd_fce - 1)*grid%AreaFactor_eps(1:ubnd_fce - 2) ! AreaFactor_eps = 1/A * dA/dz (at epsilon posions)
            if (grid%Az(1) /= 0) then ! Flux from bottom only!
               sources(2) = sources(2) + flux(1)*(grid%Az(1) + grid%Az(2))/(grid%Az(2)*(grid%h(1) + grid%h(2)))
            end if
            sources(ubnd_fce-1)= sources(ubnd_fce-1)+flux(ubnd_fce)*(grid%Az(ubnd_fce)+grid%Az(ubnd_fce-1))/(grid%Az(ubnd_fce-1)*(grid%h(ubnd_vol)+grid%h(ubnd_vol-1)))
            ! eps(0) and eps(ubnd_fce) are assigned in the post processing function
         else ! no fluxes, unity A-matrix + condition on RHS
            rhs_0 = state%u_taub**2/sqrt(state%cm0*state%cde)
            rhs_ubnd = state%u_taus**2/sqrt(state%cm0*state%cde)

            ! Trick to have rhs(1) = rhs_0 and rhs(ubnd_fce) = rhs_ubnd
            ! in discretization, rhs is calculated as rhs(1) = var(1) + sources(1)*dt
            ! Given the equation below, the following results:
            ! rhs(1) = var(1) + (-var(1)/dt + rhs_0/dt)*dt = var(1) -var(1) + rhs_0 = rhs_0
            sources(1) = -(self%var(1)/state%dt) + (rhs_0/state%dt)
            sources(ubnd_fce) = -(self%var(ubnd_fce)/state%dt) + (rhs_ubnd/state%dt)
         end if

      end associate
   end subroutine


   subroutine eps_var_post_solve(self, state)
      class(EpsModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      integer :: i,j
      real(RK) ::epslim
      real(RK) :: log_NN
      associate (grid=>self%grid, &
                 ubnd_fce=>self%grid%ubnd_fce, &
                 ubnd_vol=>self%grid%ubnd_vol)

         if (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) then
            ! Define eps at boundary (no flux)
            self%var(1) = self%var(2) + (state%cde*((state%ko(2))**1.5_RK)/(kappa*(K_s + grid%h(1)))**2)*grid%h(1)
            self%var(ubnd_fce)=   self%var(ubnd_fce-1)+(state%cde*((state%ko(ubnd_fce-1))**1.5_RK)/(kappa*(z0 +grid%h(ubnd_vol)))**2)*grid%h(ubnd_vol)
         end if

         state%keps_counter = state%keps_counter + 1
         ! check lower limit of eps / update num
         do i = 1, ubnd_fce
            if (state%NN(i) > 0) then
               epslim = 0.212_RK*state%k(i)*sqrt(state%NN(i))
            else
               epslim = eps_min
            end if
            if (state%eps(i) < epslim) then
               state%eps(i) = epslim
            end if
            if (state%eps(i) < 0) then
               write(6, *) 'Dissipation negative'
            end if


            if (self%cfg%apparent_diffusivity) then
               ! Below 120 m (where double diffusive staircases exist)
               if (grid%z_volume(i) < 365) then

                  ! If stratification is negative, then use k-eps-model to mix
                  if ((state%NN(i) < 0.0_RK)) then
                     state%num(i) = state%cmue1(i)*state%k(i)*state%k(i)/state%eps(i) + num_mol
                     state%nuh(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nuh_mol
                     state%nus(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nus_mol
                     state%nug(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nug_mol
                     state%nut(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nut_mol

                     state%nu_he(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_he_mol
                     state%nu_ne(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_ne_mol
                     state%nu_ar(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_ar_mol
                     state%nu_kr(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_kr_mol
                  else
                     ! Otherwise use the apparent diffusivity parameterization of double diffusive staircases

                     ! Use molecular diffusion of momentum (no real reason for this...)
                     state%num(i) = state%cmue1(i)*state%k(i)*state%k(i)/state%eps(i) + num_mol

                     if (state%NN(i) <= 0) state%NN(i) = 1e-15_RK
                     log_NN = log(state%NN(i))


                     ! Between 120 and 170 m
                     if (grid%z_volume(i) > 3316) then
                        ! From Kapp-N2 loglog regression new
                        !state%nuh(i) = exp(-0.5573*log_NN - 18.86)
                        !state%nus(i) = exp(-0.6511*log_NN - 23.24)

                        ! From Kapp-N2 loglog regression new new (29.4.2020)
                        state%nuh(i) = exp(-0.7184*log_NN - 21.22)
                        state%nus(i) = exp(-0.8202*log_NN - 25.79)

                        !if (state%nuh(i) < 1e-6+nuh_mol) state%nuh(i) = 1e-6+nuh_mol
                        !if (state%nus(i) < 1e-8+nus_mol) state%nus(i) = 1e-8+nus_mol

                     else
                        ! From Kapp-N2 loglog regression new new (8.4.2020)
                        !state%nuh(i) = 0.4e-9*(log_NN + 1)**4
                        !state%nus(i) = exp(-0.6587*log_NN - 22.95)

                        ! From Kapp-N2 loglog regression new new new (29.4.2020)
                        !state%nuh(i) = 0.4e-9*(log_NN + 1)**4
                        state%nuh(i) = exp(-0.7184*log_NN - 19.90)
                        state%nus(i) = exp(-0.8202*log_NN - 24.41)

                     end if

                  end if
               else
                  state%num(i) = state%cmue1(i)*state%k(i)*state%k(i)/state%eps(i) + num_mol
                  state%nuh(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nuh_mol
                  state%nus(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nus_mol
                  state%nug(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nug_mol
                  state%nut(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nut_mol

                  state%nu_he(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_he_mol
                  state%nu_ne(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_ne_mol
                  state%nu_ar(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_ar_mol
                  state%nu_kr(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_kr_mol      
               end if

            else
               state%num(i) = state%cmue1(i)*state%k(i)*state%k(i)/state%eps(i) + num_mol
               state%nuh(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nuh_mol
               state%nus(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nus_mol
               state%nug(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nug_mol
               state%nut(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nut_mol

               state%nu_he(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_he_mol
               state%nu_ne(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_ne_mol
               state%nu_ar(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_ar_mol
               state%nu_kr(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + nu_kr_mol
            end if

            ! Diffusion cannot be lower than molecular
            if (state%nuh(i) < nuh_mol) state%nuh(i) = nuh_mol
            if (state%nus(i) < nus_mol) state%nus(i) = nus_mol

            ! We assume that diffusivity of gases is equal to diffusivity of salts
            state%nug(i) = state%nus(i)!/nus_mol*nug_mol
            state%nut(i) = state%nug(i)

            ! Calculate heat flux
            state%heat_flux(i) = -cp*rho_0*state%nuh(i)*(state%T(i) - state%T(i - 1))/(grid%h(i) + grid%h(i - 1))*2           
            state%salt_flux(i) = -state%nus(i)*(state%S(i) - state%S(i - 1))/(grid%h(i) + grid%h(i - 1))*2

         end do

         state%num(1) = kappa*state%u_taub*K_s + avh_min
         state%num(ubnd_fce) = kappa*state%u_taus*z0 + avh_min
      end associate
   end subroutine

end module strat_keps
