!<    +---------------------------------------------------------------+
!     | Implementation of a statevar for a generic transport variable
!     | At the moment used for S, but could be used for any other biogeochemical type
!<    +---------------------------------------------------------------+


module strat_transport_reaction
   use strat_kinds
   use strat_consts
   use strat_simdata
   use strat_statevar
   use strat_grid
   use strat_solver
   implicit none
   private

   type, extends(ModelVariable), public :: TransportReactionModVar
      real(RK), dimension(:), pointer :: dVar
      real(RK), pointer               :: decay
   contains
      procedure, pass(self), public :: assign_external_source => transreact_assign_external_source
      procedure, pass(self), public :: assign_decay_constant => transreact_assign_decay_constant
      procedure, pass(self), public :: calc_terms => transreact_var_calc_terms
   end type

contains

  ! Method to assign external source
   subroutine transreact_assign_external_source(self, dVar)
      class(TransportReactionModVar), intent(inout) :: self
      real(RK), dimension(:), target :: dVar
      self%dVar => dVar
   end subroutine

   ! Method to assign decay constant
   subroutine transreact_assign_decay_constant(self, decay)
      class(TransportReactionModVar), intent(inout) :: self
      real(RK), target :: decay
      self%decay => decay
   end subroutine

   ! Calculate source terms according to external source
   subroutine transreact_var_calc_terms(self, state, param, sources, boundaries)
      class(TransportReactionModVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), dimension(:) ::  sources, boundaries
      associate (grid=>self%grid, &
                 ubnd_fce=>self%grid%ubnd_fce, &
                 ubnd_vol=>self%grid%ubnd_vol)

         !!!!!!!! Define sources !!!!!!!!
         sources = self%dVar ! We only have a generic, external source!

         ! Decay (if applicable)
         self%Var = self%Var*exp(-self%decay*state%dt)

         ! No explicit boundary conditions
         boundaries(1:ubnd_vol) = 0

      end associate
   end subroutine

end module strat_transport_reaction
