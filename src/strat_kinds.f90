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
!     |  Type definitions for some global interfaces
!<    +---------------------------------------------------------------+

module strat_kinds
   implicit none
   private

   ! Common constants
   integer, parameter, public :: RK = kind(0.d0) !Real kind
   character(len=3), parameter, public :: version = 'Kivu'
   integer, parameter, public :: n_simstrat = 12

   ! Common Types
   type, abstract, public :: LinSysSolver
   contains
      procedure(generic_linsyssolver_solve), deferred, nopass :: solve
   end type

   type, abstract, public :: SimModule
   contains
      procedure(generic_simmodule_update), deferred, nopass :: update
   end type

   type, abstract, extends(SimModule), public :: StateVariable
contains
   procedure(generic_statevariable_solve), deferred, nopass :: solve
end type

contains

subroutine generic_linsyssolver_solve(ld, md, ud, rhs, x, ubnd)
   implicit none

   ! Arguments
   real(RK), intent(in) :: ld(:), md(:), ud(:) ! Diagonals (A)
   real(RK), intent(in) :: rhs(:) ! right-hand side (b)
   real(RK), intent(out) :: x(:) ! solution (x)
   integer, intent(in) :: ubnd

end subroutine

subroutine generic_statevariable_solve()
   implicit none
end subroutine

subroutine generic_simmodule_update()
   implicit none
end subroutine
end module strat_kinds
