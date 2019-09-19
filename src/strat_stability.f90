!<    +---------------------------------------------------------------+
!     |  Stability Module
!     |  - Contains methods to update cmue_cn /qe and NN
!<    +---------------------------------------------------------------+

module strat_stability
   use strat_kinds
   use strat_consts
   use strat_grid
   use strat_simdata
   implicit none
   private

   ! Common Types
   type, public :: StabilityModule
      class(StaggeredGrid), pointer :: grid
      class(ModelConfig), pointer :: model_cfg
      class(ModelParam), pointer :: model_param
      logical couple_aed2

   contains
      procedure, pass :: init => stability_module_init
      procedure, pass :: update => stability_module_update
      procedure, pass :: update_cmue_cn => stability_module_update_cmue_cn
      procedure, pass :: update_cmue_qe => stability_module_update_cmue_qe
      procedure, pass :: update_NN => stability_module_update_NN

   end type
contains

   subroutine stability_module_init(self, grid, model_cfg, model_param)
      implicit none
      class(StabilityModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_cfg
      class(ModelParam), target :: model_param

      self%grid => grid
      self%model_cfg => model_cfg
      self%model_param => model_param

      self%couple_aed2 = model_cfg%couple_aed2
   end subroutine

   ! Update state variables
   subroutine stability_module_update(self, state)
      implicit none
      class(StabilityModule) :: self
      class(ModelState) :: state
      real(RK), dimension(self%grid%ubnd_fce) :: beta

      ! Do buoyancy update (update NN)
      call self%update_NN(state)

      ! Update cmue depending on selected stabilty function
      if (self%model_cfg%stability_func == 1) then
         call self%update_cmue_cn(state%cmue1, state%cmue2)

      else if (self%model_cfg%stability_func == 2) then
         beta(1:self%grid%ubnd_fce) = state%NN(1:self%grid%ubnd_fce)*(state%k(1:self%grid%ubnd_fce)/state%eps(1:self%grid%ubnd_fce))**2

         beta(1) = 0
         beta(self%grid%ubnd_fce) = 0

         call self%update_cmue_qe(beta, state%cmue1, state%cmue2, state%cde)

      end if

   end subroutine

   ! Compute NN from T and salinity
   subroutine stability_module_update_NN(self, state)
      implicit none
      class(StabilityModule) :: self
      class(ModelState) :: state

      ! Local variables
      real(RK) :: buoy(self%grid%length_vol)
      real(RK) :: rho0t(self%grid%length_vol), rho0st(self%grid%length_vol), rho0_co2(self%grid%length_vol), rho0_ch4(self%grid%length_vol)
      real(RK) :: ch4(self%grid%length_vol), co2(self%grid%length_vol), dic(self%grid%length_vol), pH(self%grid%length_vol)
      real(RK) :: k1(self%grid%length_vol)
      integer :: i

      associate (grd=>self%grid, T=>state%T, S=>state%S, NN=>state%NN, rho=>state%rho)

         ! Get gas concentrations from AED2
         if (self%couple_aed2) then
            do i = 1, state%n_AED2
               select case(trim(state%AED2_names(i)))
               case('CAR_ch4')
                  ch4 = state%AED2_state(:,i)/1000
               end select
            end do
            do i = 1, state%n_AED2_diag
               select case(trim(state%AED2_diagnames(i)))
               case('CAR_CO2')
                  co2 = state%AED2_diagstate(:,i)/1000
               end select
            end do
         end if

         do i = 1, grd%ubnd_vol
            ! According to Chen Millero, changed according "Double diffusion in Lake Kivu" from Schmid et al., 2010
            rho0t(i) = 0.9998395_RK + T(i)*(6.7914e-5_RK + T(i)*(-9.0894e-6_RK + T(i)* &
                                              (1.0171e-7_RK + T(i)*(-1.2846e-9_RK + T(i)*(1.1592e-11_RK + T(i)*(-5.0125e-14_RK))))))
            rho0st(i) = (7.5e-4_RK + T(i)*(-3.85e-6_RK + T(i)*(4.96e-8_RK)))*S(i)

            if (self%couple_aed2) then
               rho0_co2(i) = 0.0125*co2(i)
               rho0_ch4(i) = -0.02*ch4(i)
               rho0t(i) = rho0t(i) + rho0_co2(i) + rho0_ch4(i)
            end if

            buoy(i) = -g*(rho(i) - rho_0)/rho_0
         end do

         NN(2:grd%ubnd_fce - 1) = grd%meanint(1:grd%ubnd_vol - 1)*(buoy(2:grd%ubnd_vol) - buoy(1:grd%ubnd_vol - 1))
         NN(1) = NN(2)
         NN(grd%ubnd_fce) = NN(grd%ubnd_fce - 1)

      end associate
   end subroutine

   subroutine stability_module_update_cmue_cn(self, cmue1, cmue2)
      implicit none

      ! Global variables
      class(StabilityModule) :: self
      real(RK), dimension(:), intent(inout) :: cmue1, cmue2

      ! Standard version of k-eps model
      cmue1 = cmue
      cmue2 = cmue/Prndtl

      ! Set boundaries
      cmue1(1) = cmue1(2)
      cmue2(1) = cmue2(2)
      cmue1(self%grid%ubnd_fce) = cmue1(self%grid%ubnd_fce - 1)
      cmue2(self%grid%ubnd_fce) = cmue2(self%grid%ubnd_fce - 1)
   end subroutine

   subroutine stability_module_update_cmue_qe(self, beta, cmue1, cmue2, cde)
      implicit none
      class(StabilityModule) :: self
      real(RK), dimension(:), intent(in) :: beta
      real(RK), dimension(:), intent(inout) ::cmue1, cmue2
      real(RK) :: cde
      real(RK) :: gh, sm, sh
      integer :: i

      do i = 2, self%grid%ubnd_fce - 1
         gh = -cde**2*0.5_RK*beta(i)
         if (gh > 0.02) gh = gh - (gh - 0.02_RK)**2/(gh + 0.0233_RK - 2*0.02_RK)
         if (gh < -0.28) gh = -0.28_RK

         sm = 1.0_RK - 3*c1 - 6*a1/b1 - 3*a2*gh*((b2 - 3*a2)*(1.0_RK - 6*a1/b1) - 3*c1*(b2 + 6*a1))
         sm = a1*sm/((1.0_RK - 3*a2*gh*(6*a1 + b2))*(1.0_RK - 9*a1*a2*gh))
         sh = a2*(1.0_RK - 6*a1/b1)/(1.0_RK - 3*a2*gh*(6*a1 + b2))

         cmue1(i) = sqrt(2.0_RK)*cde*sm
         cmue2(i) = sqrt(2.0_RK)*cde*sh
      end do

      ! Set boundaries
      cmue1(1) = cmue1(2)
      cmue2(1) = cmue2(2)
      cmue1(self%grid%ubnd_fce) = cmue1(self%grid%ubnd_fce - 1)
      cmue2(self%grid%ubnd_fce) = cmue2(self%grid%ubnd_fce - 1)
   end subroutine

end module
