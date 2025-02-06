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
!     |  Lateral module
!     |     Reads and processes inflows/outflows such that
!     |     Advection can be calculated in the next step
!     |     There are at least two different implementations possible:
!     |         - LateralRhoModule : Inflow plunges according to density
!     |         - LateralModule:     Inflow affects layer as configured in file
!<    +---------------------------------------------------------------+


module strat_lateral
   use strat_kinds
   use strat_simdata
   use strat_consts
   use strat_grid
   use utilities
   implicit none
   private

   ! Generic base class for both modules (LateralRho and Normal)
   type, abstract, public :: GenericLateralModule
      class(ModelConfig), pointer :: cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param

      ! Variables that where either marked with "save" before, or that have been
      ! global, but only used in the lateral environment:
      real(RK), dimension(:, :), allocatable   :: z_Inp, Q_start, Qs_start, Q_end, Qs_end, Q_read_start, Q_read_end
      real(RK), dimension(:, :), allocatable   :: Inp_read_start, Inp_read_end, Qs_read_start, Qs_read_end
      real(RK), dimension(:), allocatable  :: tb_start, tb_end ! Input depths, start time, end time
      integer, dimension(:), allocatable  :: eof, nval, nval_deep, nval_surface, fnum
      logical, dimension(:), allocatable :: has_surface_input, has_deep_input
      integer :: n_vars, n_ch4, n_dic, max_n_inflows
      logical :: couple_aed2
      character(len=100) :: simstrat_path(n_simstrat), aed2_path

   contains
      procedure, pass :: init => lateral_generic_init
      procedure(lateral_generic_update), deferred, pass :: update
   end type

   ! Subclasses
   type, extends(GenericLateralModule), public :: LateralRhoModule
   contains
      procedure, pass, public :: update => lateral_rho_update
   end type

   type, extends(GenericLateralModule), public:: LateralModule
   contains
      procedure, pass, public :: update => lateral_update
   end type

contains
   subroutine lateral_generic_update(self, state)
      implicit none
      class(GenericLateralModule) :: self
      class(ModelState) :: state
   end subroutine

   subroutine lateral_generic_init(self, state, model_config, input_config, aed2_config, model_param, grid)
      implicit none
      class(GenericLateralModule) :: self
      class(ModelState) :: state
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(InputConfig), target :: input_config
      class(AED2Config), target :: aed2_config
      class(ModelParam), target :: model_param

      ! Locals
      integer :: i

      self%cfg => model_config
      self%param => model_param
      self%grid => grid

      self%n_vars = n_simstrat
      self%simstrat_path(1) = input_config%QinpName
      self%simstrat_path(2) = input_config%QoutName
      self%simstrat_path(3) = input_config%TinpName
      self%simstrat_path(4) = input_config%SinpName
      self%simstrat_path(5) = input_config%TrinpName
      self%simstrat_path(6) = input_config%HOinpName
      self%simstrat_path(7) = input_config%DinpName
      self%simstrat_path(8) = input_config%LAinpName
      self%simstrat_path(9) = input_config%HeinpName
      self%simstrat_path(10) = input_config%NeinpName
      self%simstrat_path(11) = input_config%ArinpName
      self%simstrat_path(12) = input_config%KrinpName

      self%couple_aed2 = model_config%couple_aed2
      if (self%couple_aed2) then
         self%n_vars = self%n_vars + state%n_AED2_state
         self%aed2_path = aed2_config%path_aed2_inflow
      end if

      self%max_n_inflows = model_config%max_length_input_data

      allocate(self%eof(self%n_vars))
      allocate(self%fnum(self%n_vars))
      allocate(self%nval(self%n_vars))
      allocate(self%nval_deep(self%n_vars))
      allocate(self%nval_surface(self%n_vars))
      allocate(self%tb_start(self%n_vars))
      allocate(self%tb_end(self%n_vars))

      allocate (self%z_Inp(1:self%n_vars,1:self%max_n_inflows)) ! Input depths read from file
      allocate (self%Inp_read_start(1:self%n_vars,1:self%max_n_inflows))  ! Input read from file
      allocate (self%Inp_read_end(1:self%n_vars,1:self%max_n_inflows))  ! Input read from file
      allocate (self%Q_read_start(1:self%n_vars, 1:self%max_n_inflows)) ! Integrated input
      allocate (self%Q_read_end(1:self%n_vars, 1:self%max_n_inflows)) ! Integrated input
      allocate (self%Qs_read_start(1:self%n_vars,1:self%max_n_inflows))  ! Integrated surface input
      allocate (self%Qs_read_end(1:self%n_vars,1:self%max_n_inflows))  ! Integrated surface input
      allocate (self%Q_start(1:self%n_vars,1:grid%nz_grid + 1)) ! Input interpolated on grid
      allocate (self%Q_end(1:self%n_vars,1:grid%nz_grid + 1)) ! Input interpolated on grid
      allocate (self%Qs_start(1:self%n_vars,1:grid%nz_grid + 1)) ! Surface input interpolated on grid
      allocate (self%Qs_end(1:self%n_vars,1:grid%nz_grid + 1)) ! Surface input interpolated on grid

      allocate(state%Q_inp(1:self%n_vars,1:grid%nz_grid + 1))
      allocate(self%has_surface_input(1:self%n_vars))
      allocate(self%has_deep_input(1:self%n_vars))

      self%n_ch4 = 0
      self%n_dic = 0

      ! Get gas concentrations from AED2
      if (self%couple_aed2) then
         do i = 1, state%n_AED2_state
            select case(trim(state%AED2_state_names(i)))
            case('CAR_ch4')
               self%n_ch4 = i
            end select
         end do
         do i = 1, state%n_AED2_state
            select case(trim(state%AED2_state_names(i)))
            case('CAR_dic')
               self%n_dic = i
            end select
         end do
         do i = 1, state%n_AED2_state
            select case(trim(state%AED2_state_names(i)))
            case('CAR_pH')
               state%n_pH = i
            end select
         end do
      end if
   end subroutine
      
      
   ! Implementation for lateral rho
   subroutine lateral_rho_update(self, state)
      implicit none
      class(LateralRhoModule) :: self
      class(ModelState) :: state

      ! Local Declarations
      real(RK) :: Inp(1:self%n_vars,1:self%max_n_inflows)
      real(RK) :: dummy
      real(RK) :: Q_in(1:self%grid%ubnd_vol), h_in(1:self%grid%ubnd_vol)
      real(RK) :: T_in, S_in, Tr_in, HO_in, D_in, LA_in, co2_in, ch4_in, rho_in, CD_in, g_red, slope, Ri, E, Q_inp_inc
      real(RK) :: He_in, Ne_in, Ar_in, Kr_in
      real(RK) :: AED2_in(state%n_AED2_state)

      integer :: i, j, k, i1, i2, status
      character(len=100) :: fname

      associate (datum=>state%datum, &
                 idx=>state%first_timestep, &
                 Q_inp=>state%Q_inp, & ! Q_inp is the input at each depth for each time step
                 Q_vert=>state%Q_vert, & ! Q_vert is the integrated net water input at each depth (integrated inflow - outflow)
                 grid=>self%grid, &
                 ubnd_vol=>self%grid%ubnd_vol, &
                 ubnd_fce=>self%grid%ubnd_fce)


         do i=1, self%n_vars
            if (idx) then
               self%max_n_inflows = 0

               if (i > n_simstrat) then
                  fname = trim(self%aed2_path)//trim(state%AED2_state_names(i - n_simstrat))//'_inflow.dat'
               else
                  fname = trim(self%simstrat_path(i))
               end if
               self%fnum(i) = i + 60  ! Should find a better way to manage unit numbers
               open(self%fnum(i), action='read', status='old', file=fname)
               
               if (status .ne. 0) then
                  call error('File '//fname//' not found.')
               else
                  write(6,*) 'Reading ', fname
               end if

               ! Default values
               self%Q_start(i,:) = 0.0_RK
               self%Q_end(i,:) = 0.0_RK
               self%Qs_start(i,:) = 0.0_RK
               self%Qs_end(i,:) = 0.0_RK

               ! End of file is not reached
               self%eof(i) = 0

               ! Read input depths
               read(self%fnum(i),*,end=9)

               ! Read number of deep and surface columns
               read(self%fnum(i), *, end=9) self%nval_deep(i), self%nval_surface(i)
               if (self%nval_deep(i) > 0) then
                  self%has_deep_input(i) = .true.
               else
                  self%has_deep_input(i) = .false.
               end if
               if (self%nval_surface(i) > 0) then
                  self%has_surface_input(i) = .true.
               else
                  self%has_surface_input(i) = .false.
               end if

               ! Total number of values to read
               self%nval(i) = self%nval_deep(i) + self%nval_surface(i)

               if (self%nval(i) > self%max_n_inflows) self%max_n_inflows = self%nval(i)

               ! Read input depths
               read(self%fnum(i),*,end=9) dummy, (self%z_Inp(i,j),j=1,self%nval(i))
               ! Convert deep input depths
               self%z_Inp(i,1:self%nval_deep(i)) = grid%z_zero + self%z_Inp(i,1:self%nval_deep(i))
               ! Convert surface input depths
               self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = grid%lake_level + self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i))

               !Read first input values
               read(self%fnum(i),*,end=9) self%tb_start(i),(self%Inp_read_start(i,j),j=1,self%nval(i))

               ! If there is deep outflow (i==2)
               if (i==2 .and. self%has_deep_input(i)) then
                  call Integrate(self%z_Inp(i,1:self%nval_deep(i)),self%Inp_read_start(i,1:self%nval_deep(i)),self%Q_read_start(i,:),self%nval_deep(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i,1:self%nval_deep(i)),self%Q_read_start(i,:),self%nval_deep(i),self%Q_start(i,:))
               end if
               ! If there is any surface inflow
               if (self%has_surface_input(i)) then
                  call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_start(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
               end if


               ! Read next line
               read(self%fnum(i),*,end=7) self%tb_end(i),(self%Inp_read_end(i,j),j=1,self%nval(i))

               ! If there is deep outflow (i==2)
               if (i==2 .and. self%has_deep_input(i)) then
                  call Integrate(self%z_Inp(i,1:self%nval_deep(i)),self%Inp_read_end(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i),self%Q_end(i,:))
               end if

               ! If there is any surface inflow
               if (self%has_surface_input(i)) then
                  call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
               end if
               call ok('Input file successfully read: '//fname)
            end if ! if idx

            ! If lake level changes and if there is surface inflow, adjust inflow depth to keep relative inflow depth constant
            if ((.not. grid%lake_level == grid%lake_level_old) .and. self%has_surface_input(i)) then

               ! Readjust surface input depths
               self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i)) - grid%lake_level_old + grid%lake_level

               ! Adjust surface inflow to new lake level
               if (self%has_surface_input(i)) then
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
               end if
            end if ! end if not lake_level...


            if ((datum<=self%tb_start(i)).or.(self%eof(i)==1)) then    ! if datum before first date or end of file reached
               goto 8
            else
               do while (.not.((datum>=self%tb_start(i)).and.(datum<=self%tb_end(i)))) ! do until datum between dates
                  ! Move one step in time
                  self%tb_start(i) = self%tb_end(i)
                  self%Qs_start(i, :) = self%Qs_end(i, :)
                  self%Qs_read_start(i, :) = self%Qs_read_end(i, :)

                  ! For outflow, take Q_start; for inflow, temperature and salinity, take Inp for the plunging algorithm
                  if (i==2) then
                     self%Q_start(i, :) = self%Q_end(i, :)
                  else
                     self%Inp_read_start(i,1:self%nval_deep(i)) = self%Inp_read_end(i,1:self%nval_deep(i))
                  end if

                  ! Read next line
                  read(self%fnum(i),*,end=7) self%tb_end(i),(self%Inp_read_end(i,j),j=1,self%nval(i))

                  ! If there is deep outflow (i==2)
                  if (i==2 .and. self%has_deep_input(i)) then
                     call Integrate(self%z_Inp(i,1:self%nval_deep(i)),self%Inp_read_end(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i),self%Q_end(i,:))
                  end if

                  ! If there is any surface inflow
                  if (self%has_surface_input(i)) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  end if
               end do

               if(self%tb_end(i)<=self%tb_start(i)) then
                  call error('Dates in '//trim(fname)//' file must always be increasing.')
               end if

               ! Linearly interpolate value at correct datum
               if (i/=2) then
                  ! For plunging input, Inp will be needed later
                  Inp(i,1:self%nval_deep(i)) = self%Inp_read_start(i,1:self%nval_deep(i)) + (datum-self%tb_start(i)) &
                  * (self%Inp_read_end(i,1:self%nval_deep(i)) - self%Inp_read_start(i,1:self%nval_deep(i)))/(self%tb_end(i) - self%tb_start(i))

                  ! Surface input is already added to Q_inp; the plunging algorithm will add the deep input further below
                  do j=1,ubnd_fce
                     Q_inp(i,j) = (self%Qs_start(i,j)) + (datum - self%tb_start(i)) &
                     * (self%Qs_end(i,j) - self%Qs_start(i,j))/(self%tb_end(i) - self%tb_start(i))
                  end do
               else
                  ! For outflow (i==2), both surface and deep inputs are added to Q_inp
                  do j=1,ubnd_fce
                     Q_inp(i,j) = (self%Q_start(i,j) + self%Qs_start(i,j)) + (datum-self%tb_start(i)) &
                     * (self%Q_end(i,j) + self%Qs_end(i,j) - self%Q_start(i,j) - self%Qs_start(i,j))/(self%tb_end(i)-self%tb_start(i))
                  end do
               end if
            end if
            goto 11

7           self%eof(i) = 1
8           if(i/=2) Inp(i,1:self%nval_deep(i)) = self%Inp_read_start(i,1:self%nval_deep(i)) ! Set to closest available value
            Q_inp(i,1:ubnd_vol) = self%Q_start(i,1:ubnd_vol) + self%Qs_start(i,1:ubnd_vol) ! Set to closest available value
            goto 11

9           write(6,*) '[WARNING] ','No data found in ',trim(fname),' file. Check number of depths. Values set to zero.'
            self%eof(i) = 1
            if(i/=2) Inp(i,1:self%nval_deep(i)) = 0.0_RK
            if(i/=2) self%Inp_read_start(i,1) = 0.0_RK
            if(i==2) Q_inp(i,1:ubnd_vol) = 0.0_RK
            if(i==2) self%Q_start(i,1:ubnd_fce) = 0.0_RK
            self%Qs_start(i,1:ubnd_fce) = 0.0_RK

11          continue
         end do      ! end do i=1,n_vars

         !Set Q_inp to the differences (from the integrals)
         do i = 1,self%n_vars
            do j = 1, ubnd_vol
               Q_inp(i, j) = Q_inp(i, j + 1) - Q_inp(i, j)
            end do
            Q_inp(i,ubnd_vol + 1) = 0
         end do

         ! Only if biochemistry enabled: Transform pH to [H] for physical mixing processes
         if (self%couple_aed2 .and.state%n_pH > 0) then
            ! current pH profile
            state%AED2_state(:,state%n_pH) = 10.**(-state%AED2_state(:,state%n_pH))
            do i=1,ubnd_vol
               if (Q_inp(n_simstrat + state%n_pH,i) > 0) then
                  ! Surface inflows: pH is given as pH*m2/s, so before transforming to [H], we need to get rid of the m2/s temporarily
                  Q_inp(n_simstrat + state%n_pH,i) = Q_inp(n_simstrat + state%n_pH,i)/Q_inp(1,i)
                  Q_inp(n_simstrat + state%n_pH,i) = 10.**(-Q_inp(n_simstrat + state%n_pH,i))
                  Q_inp(n_simstrat + state%n_pH,i) = Q_inp(n_simstrat + state%n_pH,i)*Q_inp(1,i)
               end if
            end do
         end if

         ! Plunging algorithm
         do j = 1,self%nval_deep(1)  ! nval_deep needs to be the same for all i
            if (Inp(1,j) > 1E-15) then
               k = ubnd_vol
               do while (grid%z_volume(k) > self%z_Inp(1,j)) ! Find the place where the plunging inflow enters the lake (defined in file)
                  k = k - 1
               end do

               ! Get initial Q, T and S for the plunging inflow (before entrainment of ambient water)
               Q_in(k) = Inp(1,j) !Inflow flow rate [m3/s]
               T_in = Inp(3,j) !Inflow temperature [°C]
               S_in = Inp(4,j) !Inflow salinity [‰]
               Tr_in = Inp(5,j) ! Inflow tritium [TU]
               HO_in = Inp(6,j) ! Inflow heavy oxygen [‰]
               D_in = Inp(7,j) ! Inflow deuterium [‰]
               LA_in = Inp(8,j) ! Inflow 39Ar [‰]

               He_in = Inp(9,j) ! Inflow He [ccSTP]
               Ne_in = Inp(10,j) ! Inflow Ne [ccSTP]
               Ar_in = Inp(11,j) ! Inflow Ar [ccSTP]
               Kr_in = Inp(12,j) ! Inflow Kr [ccSTP]

               ! Only if biochemistry enabled
               if (self%couple_aed2) then
                  ! Get AED2 values for the plunging inflow (before entrainment of ambient water)
                  AED2_in = Inp(n_simstrat + 1 : self%n_vars,j)
                  ! Calculate CO2 [mol/L] from DIC and pH (analogous to AED2 carbon module) for density calculation
                  call calc_co2_from_dic(co2_in, T_in, S_in, AED2_in(self%n_dic), AED2_in(state%n_pH))
                  ! Transform pH to [H] for physical mixing processes
                  if (state%n_pH > 0) AED2_in(state%n_pH) = 10.**(-AED2_in(state%n_pH))
                  ! Get CH4 [mol/L] for density calculation
                  ch4_in = AED2_in(self%n_ch4)/1e6 ! Inflow CH4 [mol/L]
                  ! Compute density as a function of T, S, CO2 and CH4
                  call calc_density_aed2(rho_in, T_in, S_in, co2_in, ch4_in) !Inflow density [kg/m3]
               else
                  ! Compute density as a function of T and S
                  call calc_density(rho_in, T_in, S_in)
               end if
               g_red = g*(rho_in - state%rho(k))/rho_in !Reduced gravity [m/s2]

               slope = pi/72 !Slope of inflow
               !hang = pi/3 !Stream half-angle
               CD_in = self%param%CD*10 !Inflow drag coefficient
               !Ri = CD_in*(1+0.21*CD_in**0.5*sin(hang))/(sin(hang)*tan(slope)) !Richardson number
               !Ri = CD_in/tan(slope)*(1/sin(hang)+0.21*CD_in**0.5) !Richardson number
               Ri = CD_in/tan(slope)*(1.15 + 0.21*CD_in**0.5) !Richardson number (assuming an inflow half-angle of pi/3)
               E = 1.6*CD_in**1.5/Ri !Entrainment coefficient
               h_in(k) = (2*Q_in(k)**2*Ri*tan(slope)**2/abs(g_red))**0.2 !Inflow thickness [m]

               if (g_red > 0) then !Inflow plunges
                  do while ((rho_in > state%rho(k)).and.(k > 1))
                     h_in(k - 1) = 1.2*E*(grid%z_volume(k) - grid%z_volume(k-1))/sin(slope) + h_in(k)
                     Q_in(k - 1) = Q_in(k)*(h_in(k - 1)/h_in(k))**(5./3.)
                     Q_inp(2,k) = Q_inp(2,k) - (Q_in(k-1) - Q_in(k))
                     T_in = (T_in*Q_in(k) + state%T(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     S_in = (S_in*Q_in(k) + state%S(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     Tr_in = (Tr_in*Q_in(k) + state%Tr(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     HO_in = (HO_in*Q_in(k) + state%heavy_oxygen(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     D_in = (D_in*Q_in(k) + state%deuterium(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     LA_in = (LA_in*Q_in(k) + state%light_ar(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)

                     He_in = (He_in*Q_in(k) + state%He(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     Ne_in = (Ne_in*Q_in(k) + state%Ne(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     Ar_in = (Ar_in*Q_in(k) + state%Ar(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     Kr_in = (Kr_in*Q_in(k) + state%Kr(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     if (self%couple_aed2) then
                        AED2_in = (AED2_in*Q_in(k) + state%AED2_state(k,:)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     end if
                     rho_in = (rho_in*Q_in(k) + state%rho(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     k = k - 1
                  end do
                  i2 = k
                  do i1 = k,ubnd_vol !extend upwards
                     if(i1 == ubnd_vol) exit
                     if(grid%z_volume(i1 + 1) > (grid%z_volume(k) + h_in(k))) exit
                  end do
               else if (g_red < 0) then !Inflow rises
                  do while ((rho_in < state%rho(k)) .and. (k < ubnd_vol))
                     h_in(k + 1) = 1.2*E*(grid%z_volume(k + 1) - grid%z_volume(k))/sin(slope) + h_in(k)
                     Q_in(k + 1) = Q_in(k)*(h_in(k + 1)/h_in(k))**(5./3.)
                     Q_inp(2,k) = Q_inp(2,k) - (Q_in(k + 1) - Q_in(k))
                     T_in = (T_in*Q_in(k) + state%T(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     S_in = (S_in*Q_in(k) + state%S(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     Tr_in = (Tr_in*Q_in(k) + state%Tr(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     HO_in = (HO_in*Q_in(k) + state%heavy_oxygen(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     D_in = (D_in*Q_in(k) + state%deuterium(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     LA_in = (LA_in*Q_in(k) + state%light_ar(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)

                     He_in = (He_in*Q_in(k) + state%He(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     Ne_in = (Ne_in*Q_in(k) + state%Ne(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     Ar_in = (Ar_in*Q_in(k) + state%Ar(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     Kr_in = (Kr_in*Q_in(k) + state%Kr(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     if (self%couple_aed2) then
                        AED2_in = (AED2_in*Q_in(k) + state%AED2_state(k,:)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     end if
                     rho_in = (rho_in*Q_in(k) + state%rho(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     k = k + 1
                  end do
                  i1 = k
                  do i2 = k,1,-1 !extend downwards
                     if(i2 == 1) exit
                     if(grid%z_volume(i2 - 1) < (grid%z_volume(k) - h_in(k))) exit
                  end do
               end if

               ! Deep plunging input is added to Q_inp for i=1,3,4 (inflow, temperature, salinity)
               do i = i2,i1
                  Q_inp_inc = Q_in(k)/(grid%z_face(i1 + 1) - grid%z_face(i2))*grid%h(i)
                  Q_inp(1,i) = Q_inp(1,i) + Q_inp_inc
                  Q_inp(3,i) = Q_inp(3,i) + T_in*Q_inp_inc
                  Q_inp(4,i) = Q_inp(4,i) + S_in*Q_inp_inc
                  Q_inp(5,i) = Q_inp(5,i) + Tr_in*Q_inp_inc
                  Q_inp(6,i) = Q_inp(6,i) + HO_in*Q_inp_inc
                  Q_inp(7,i) = Q_inp(7,i) + D_in*Q_inp_inc
                  Q_inp(8,i) = Q_inp(8,i) + LA_in*Q_inp_inc
                  Q_inp(9,i) = Q_inp(9,i) + He_in*Q_inp_inc
                  Q_inp(10,i) = Q_inp(10,i) + Ne_in*Q_inp_inc
                  Q_inp(11,i) = Q_inp(11,i) + Ar_in*Q_inp_inc
                  Q_inp(12,i) = Q_inp(12,i) + Kr_in*Q_inp_inc
                  if (self%couple_aed2) Q_inp(n_simstrat + 1 : self%n_vars,i) = Q_inp(n_simstrat + 1 : self%n_vars,i) + AED2_in*Q_inp_inc
               end do
            end if
         end do
         ! Q_vert is the integrated difference between in- and outflow (starting at the lake bottom)
         ! Q_vert is located on the face grid, m^3/s
         Q_vert(1) = 0
         do i = 2,ubnd_fce
            Q_vert(i) = Q_vert(i - 1) + Q_inp(1,i - 1) + Q_inp(2,i - 1)
            state%w(i) = Q_vert(i)/grid%Az(i)
            state%lateral_input(i) = Q_inp(1,i - 1)
         end do
      end associate
   end subroutine

   ! "Normal" Implementation
   subroutine lateral_update(self, state)
      implicit none
      class(LateralModule) :: self
      class(ModelState) :: state

      ! Local Declarations
      real(RK) :: dummy
      integer :: i, j, status
      character(len=100) :: fname

      associate (datum=>state%datum, &
                 idx=>state%first_timestep, &
                 Q_inp=>state%Q_inp, & ! Q_inp is the input at each depth for each time step
                 Q_vert=>state%Q_vert, & ! Q_vert is the integrated net water input
                 grid=>self%grid, &
                 ubnd_vol=>self%grid%ubnd_vol, &
                 ubnd_Fce=>self%grid%ubnd_fce)

         do i=1, self%n_vars
            if (idx) then
               if (i > n_simstrat) then
                  fname = trim(self%aed2_path)//trim(state%AED2_state_names(i - n_simstrat))//'_inflow.dat'
               else
                  fname = trim(self%simstrat_path(i))
               end if
               self%fnum(i) = i + 60  ! Should find a better way to manage unit numbers
               open(self%fnum(i), action='read', status='old', file=fname)
               
               if (status .ne. 0) then
                  call error('File '//fname//' not found.')
               else
                  write(6,*) 'Reading ', fname
               end if

               ! Default values
               self%Q_start(i,:) = 0.0_RK
               self%Q_end(i,:) = 0.0_RK
               self%Qs_start(i, :) = 0.0_RK
               self%Qs_end(i, :) = 0.0_RK

               ! Open file and start to read
               self%eof(i) = 0
               read (self%fnum(i), *, end=9) ! Skip first row: description of columns

               ! Number of deep (fixed) and surface inputs to read
               read (self%fnum(i), *, end=9) self%nval_deep(i), self%nval_surface(i)
               if (self%nval_deep(i) > 0) then
                  self%has_deep_input(i) = .true.
               else
                  self%has_deep_input(i) = .false.
               end if
               if (self%nval_surface(i) > 0) then
                  self%has_surface_input(i) = .true.
               else
                  self%has_surface_input(i) = .false.
               end if

               ! Total number of values to read
               self%nval(i) = self%nval_deep(i) + self%nval_surface(i)

               ! Read input depths
               read (self%fnum(i), *, end=9) dummy, (self%z_Inp(i, j), j=1, self%nval(i))

               ! Convert input depths
               self%z_Inp(i, 1:self%nval_deep(i)) = grid%z_zero + self%z_Inp(i, 1:self%nval_deep(i))

               if (self%has_surface_input(i)) then
                  ! Convert surface input depths
                  self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = grid%lake_level + self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i))
               end if

               ! Read first input line
               read (self%fnum(i), *, end=9) self%tb_start(i), (self%Inp_read_start(i, j), j=1, self%nval(i))

               if (self%has_deep_input(i)) then
                  ! Cumulative integration of input
                  call Integrate(self%z_Inp(i, :), self%Inp_read_start(i, :), self%Q_read_start(i, :), self%nval_deep(i))
                  ! Interpolation on face grid
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_start(i, :), self%nval_deep(i), self%Q_start(i, :))
               end if

               ! If there is surface input, integrate and interpolate
               if (self%has_surface_input(i)) then
                  call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_start(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
               end if


               ! Read second line and treatment of deep inflow
               read (self%fnum(i), *, end=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))
               if (self%has_deep_input(i)) then
                  call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))
               end if
               ! If there is surface input, integrate and interpolate
               if (self%has_surface_input(i)) then
                  call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                  call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
               end if

               call ok('Input file successfully read: '//fname)
            end if ! idx==1



            ! If lake level changes and if there is surface inflow, adjust inflow depth to keep them at the surface
            if ((.not. grid%lake_level == grid%lake_level_old) .and. (self%has_surface_input(i))) then

               ! Readjust surface input depths
               self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i)) - grid%lake_level_old + grid%lake_level

               ! Adjust surface inflow to new lake level
               call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
               call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))

            end if ! end if not lake_level...



            ! Temporal treatment of inflow
            if ((datum <= self%tb_start(i)) .or. (self%eof(i) == 1)) then ! if datum before first date or end of file reached
               goto 8
            else
               do while (.not. ((datum >= self%tb_start(i)) .and. (datum <= self%tb_end(i)))) ! Do until datum between dates
                  self%tb_start(i) = self%tb_end(i) ! Move one step in time
                  self%Q_start(i, :) = self%Q_end(i, :)
                  self%Qs_start(i, :) = self%Qs_end(i, :)
                  self%Qs_read_start(i, :) = self%Qs_read_end(i, :)

                  read (self%fnum(i), *, end=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))

                  if (self%has_deep_input(i)) then
                    call Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                    call grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))
                  end if

                  if (self%has_surface_input(i)) then
                     call Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     call grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  end if
               end do ! end do while
            end if

            ! Linearly interpolate value at correct datum (Q_inp is on face grid)
            do j = 1, ubnd_fce
               Q_inp(i,j) = (self%Q_start(i,j) + self%Qs_start(i,j)) + (datum-self%tb_start(i)) &
               * (self%Q_end(i,j) + self%Qs_end(i,j) - self%Q_start(i,j) - self%Qs_start(i,j))/(self%tb_end(i)-self%tb_start(i))
            end do
            goto 11

            ! If end of file reached, set to closest available value
 7          self%eof(i) = 1
 8          Q_inp(i,1:ubnd_fce) = self%Q_start(i,1:ubnd_fce) + self%Qs_start(i,1:ubnd_fce)
            goto 11

            ! If no data available
 9          write(6,*) '[WARNING] ','No data found in ',trim(fname),' file. Check number of depths. Values set to zero.'
            self%eof(i) = 1
            Q_inp(i, 1:ubnd_fce) = 0.0_RK
            self%Q_start(i, 1:ubnd_fce) = 0.0_RK
            self%Qs_start(i, 1:ubnd_fce) = 0.0_RK
            11        continue

         end do ! end do i=1,self%n_vars
         ! Q_vert is the integrated difference between in- and outflow (starting at the lake bottom)
         ! Q_vert is located on the face grid, m^3/s
         Q_vert(1)=0
         Q_vert(2:ubnd_fce) = Q_inp(1, 2:ubnd_fce) + Q_inp(2, 2:ubnd_fce)

         ! The final Q_inp is located on the volume grid
         do i = 1, self%n_vars
            do j = 1, ubnd_vol
               Q_inp(i, j) = Q_inp(i, j + 1) - Q_inp(i, j)
            end do
            Q_inp(i,ubnd_vol + 1) = 0
         end do

         ! Only if biochemistry enabled: Transform pH to [H] for physical mixing processes
         if (self%couple_aed2 .and.state%n_pH > 0) then
            ! current pH profile
            state%AED2_state(:,state%n_pH) = 10.**(-state%AED2_state(:,state%n_pH))
            do i=1,ubnd_vol
               if (Q_inp(n_simstrat + state%n_pH,i) > 0) then
                  ! Surface inflows: pH is given as pH*m2/s, so before transforming to [H], we need to get rid of the m2/s temporarily
                  Q_inp(n_simstrat + state%n_pH,i) = Q_inp(n_simstrat + state%n_pH,i)/Q_inp(1,i)
                  Q_inp(n_simstrat + state%n_pH,i) = 10.**(-Q_inp(n_simstrat + state%n_pH,i))
                  Q_inp(n_simstrat + state%n_pH,i) = Q_inp(n_simstrat + state%n_pH,i)*Q_inp(1,i)
               end if
            end do
         end if

      end associate
   end subroutine

end module
