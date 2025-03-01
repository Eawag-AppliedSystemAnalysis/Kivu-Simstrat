!     +---------------------------------------------------------------+
!     | Inputfile  module
!     |  - Reads configuration and initial conditions
!     |  - Sets up simulation data structure!
!     +---------------------------------------------------------------+

module strat_inputfile
   use strat_kinds
   use strat_simdata
   use strat_grid
   use strat_consts
   use utilities
   use json_kinds, only: CK
   use json_module
   use csv_module
   implicit none

   private

   !##################################################
   !# Inputfile
   !##################################################
   type, public :: SimstratSimulationFactory
      private
      class(SimulationData), pointer :: simdata

   contains
      procedure, pass(self), public :: initialize_model
      procedure, pass(self), public :: read_json_par_file
      procedure, pass(self), public :: read_initial_data
      procedure, pass(self), public :: read_grid_config
      procedure, pass(self), public :: setup_model
      procedure, pass(self), public :: setup_output_conf
   end type SimstratSimulationFactory

contains

   subroutine initialize_model(self, fname, simdata)
      implicit none
      class(SimstratSimulationFactory) :: self
      class(SimulationData), pointer, intent(out) :: simdata
      character(len=*) :: fname

      ! Allocate model
      ! If(associated(simdata)) deallocate(simdata)

      allocate (SimulationData :: self%simdata)
      simdata => self%simdata

      ! Parse inputfile
      call self%read_json_par_file(fname)

      ! Set up grid
      call self%read_grid_config

      call self%simdata%model%init(self%simdata%grid%nz_grid)

      ! Read initial data
      call self%read_initial_data

      ! Update area factors
      call self%simdata%grid%update_area_factors()

      ! Init rest of model
      call self%setup_model()

      ! Set output configuration
      call self%setup_output_conf()

   end subroutine initialize_model


   ! Set up logger configuration
   subroutine setup_output_conf(self)
      implicit none
      class(SimstratSimulationFactory) :: self
      type(csv_file) :: f
      logical :: status_ok
      integer :: i

      associate (model=>self%simdata%model, &
                 output_cfg=>self%simdata%output_cfg)

         ! Read output depths from file if path is specified in parfile
         if (output_cfg%output_depth_type == 7) then
            call check_file_exists(output_cfg%zoutName)

            call f%read (output_cfg%zoutName, header_row=1, status_ok=status_ok)
            if (.not. status_ok) then
               call f%destroy()
               call error('Unable to read output depths: '//output_cfg%zoutName)
            end if
            ! Values are read into array zout_read instead of zout to keep zout for the final output depths. For example,
            ! zout_read will contain only one value here if the output frequency is given, but zout will contain the output
            ! depths produced from this frequency.
            call f%get(1, output_cfg%zout_read, status_ok)
            call f%destroy()

            ! Decide whether an interval or a list of values is given
            if(size(output_cfg%zout_read)==1) then
               output_cfg%depth_interval = output_cfg%zout_read(1)
            else
               output_cfg%depth_interval = 0
               if (output_cfg%output_depth_reference == 'surface') then
                  call reverse_in_place(output_cfg%zout_read)
               end if

               ! Check if depths are monotonously in- or decreasing
               do i = 1,size(output_cfg%zout_read)
                  if (i>1) then
                     if (output_cfg%zout_read(i) < output_cfg%zout_read(i-1)) then
                        call error('Output depths are not monotonous. Maybe you chose the wrong output depth reference?')
                     end if
                  end if
               end do

            end if
         end if

         ! Read output times if path is specified in parfile
         if (output_cfg%output_time_type == 7) then
            call check_file_exists(output_cfg%toutName)

            call f%read (output_cfg%toutName, header_row=1, status_ok=status_ok)
            if (.not. status_ok) then
               call f%destroy()
               call error('Unable to read output depths: '//output_cfg%toutName)
            end if
            call f%get(1, output_cfg%tout, status_ok)
            call f%destroy()

            ! Determine output mode: If only one element in tout, then this number is the thinning interval.
            ! If there are more than one number in tout, these are the output times and thinning_interval is set to 0.
            if (size(output_cfg%tout)==1) then
               output_cfg%thinning_interval = int(output_cfg%tout(1))
            else
               output_cfg%thinning_interval = 0
            endif
         end if

         ! Define variables that should be written
         if (output_cfg%output_all) then
            output_cfg%number_output_vars = 44
            output_cfg%output_var_names = [character(len=12) :: 'V','U','T','S','co2','ch4','num','nuh','nus','nug','NN','k','eps','P','B','Ps','HA','HW','HK','HV','Rad0','TotalIceH','BlackIceH','WhiteIceH','SnowH','WaterH','Qvert','rho','DiffHeatFlux','BuoyHeatFlux','AdvHeatFlux','DiffSaltFlux','AdvSaltFlux','VertVel','DensityRatio','LateralInput','Tr','HeavyOxygen','Deuterium','LightAr','He','Ne','Ar','Kr']
         else
            output_cfg%number_output_vars = size(output_cfg%output_var_names)
         end if
         allocate (self%simdata%output_cfg%output_vars(output_cfg%number_output_vars))

         do i=1,output_cfg%number_output_vars
            select case(trim(output_cfg%output_var_names(i)))
               case('V')
                  ! Horizontal velocity in y direction [m s-1]
                  self%simdata%output_cfg%output_vars(i)%name = "V"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%V
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('U')
                  ! Horizontal velocity in x direction [m s-1]
                  self%simdata%output_cfg%output_vars(i)%name = "U"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%U
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('co2')
                  ! CO2 concentration in [mmol/m3] (from AED2)
                  self%simdata%output_cfg%output_vars(i)%name = "co2"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%co2
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('ch4')
                  ! CH4 concentration in [mmol/m3] (from AED2)
                  self%simdata%output_cfg%output_vars(i)%name = "ch4"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%ch4
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('T')
                  ! Temperature [°C]
                  self%simdata%output_cfg%output_vars(i)%name = "T"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%T
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('S')
                  ! Salinity [‰]
                  self%simdata%output_cfg%output_vars(i)%name = "S"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%S
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('num')
                  ! Turbulent diffusivity for momentum  (viscosity) [m2 s]
                  self%simdata%output_cfg%output_vars(i)%name = "num"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%num
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('nuh')
                  ! Turbulent diffusivity for temperature [m2 s]
                  self%simdata%output_cfg%output_vars(i)%name = "nuh"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%nuh
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('nus')
                  ! Turbulent diffusivity for salinity [m2 s]
                  self%simdata%output_cfg%output_vars(i)%name = "nus"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%nus
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('nug')
                  ! Turbulent diffusivity for gases [m2 s]
                  self%simdata%output_cfg%output_vars(i)%name = "nug"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%nug
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('NN')
                  ! Brunt-Väisälä frequency [s-1]
                  self%simdata%output_cfg%output_vars(i)%name = "NN"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%NN
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('k')
                  ! Turbulent kinetic energy (TKE) [J kg-1]
                  self%simdata%output_cfg%output_vars(i)%name = "k"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%k
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('eps')
                  ! TKE dissipation rate [W kg-1]
                  self%simdata%output_cfg%output_vars(i)%name = "eps"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%eps
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('P')
                  ! Shear stress production P [W kg-1]
                  self%simdata%output_cfg%output_vars(i)%name = "P"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%P
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('B')
                  ! Buoyancy production [W kg-1]
                  self%simdata%output_cfg%output_vars(i)%name = "B"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%B
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('Ps')
                  ! Production of TKE due to internal seiching [W kg-1]
                  self%simdata%output_cfg%output_vars(i)%name = "Ps"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%P_seiche
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('HA')
                  ! Infrared radiation from sky [W m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "HA"
                  self%simdata%output_cfg%output_vars(i)%values_surf => self%simdata%model%ha
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('HW')
                  ! Infrared ratidation from water [W m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "HW"
                  self%simdata%output_cfg%output_vars(i)%values_surf => self%simdata%model%hw
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('HK')
                  ! Sensible heat flux from water [W m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "HK"
                  self%simdata%output_cfg%output_vars(i)%values_surf => self%simdata%model%hk
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('HV')
                  ! Latent heat flux from water [W m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "HV"
                  self%simdata%output_cfg%output_vars(i)%values_surf => self%simdata%model%hv
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('Rad0')
                  ! Solar radiation at surface [W m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "Rad0"
                  self%simdata%output_cfg%output_vars(i)%values_surf => self%simdata%model%rad0
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('TotalIceH')
                  ! Total ice thickness [m]
                  self%simdata%output_cfg%output_vars(i)%name = "TotalIceH"
                  self%simdata%output_cfg%output_vars(i)%values_surf => self%simdata%model%total_ice_h
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('BlackIceH')
                  ! Black ice thickness [m]
                  self%simdata%output_cfg%output_vars(i)%name = "BlackIceH"
                  self%simdata%output_cfg%output_vars(i)%values_surf => self%simdata%model%black_ice_h
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('WhiteIceH')
                  ! White ice (snow ice) thickness [m]
                  self%simdata%output_cfg%output_vars(i)%name = "WhiteIceH"
                  self%simdata%output_cfg%output_vars(i)%values_surf => self%simdata%model%white_ice_h
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('SnowH')
                  ! Snow layer thickness [m]
                  self%simdata%output_cfg%output_vars(i)%name = "SnowH"
                  self%simdata%output_cfg%output_vars(i)%values_surf => self%simdata%model%snow_h
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('WaterH')
                  ! Water level [m]
                  self%simdata%output_cfg%output_vars(i)%name = "WaterH"
                  self%simdata%output_cfg%output_vars(i)%values_surf => self%simdata%grid%lake_level
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('Qvert')
                  ! Vertical advection [m3 s-1]
                  self%simdata%output_cfg%output_vars(i)%name = "Qvert"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%Q_vert
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('rho')
                  ! Vertical advection [m3 s-1]
                  self%simdata%output_cfg%output_vars(i)%name = "rho"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%rho
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('DiffHeatFlux')
                  ! Vertical heat flux [W m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "DiffHeatFlux"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%diff_heat_flux
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('BuoyHeatFlux')
                  ! Vertical heat flux [W m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "BuoyHeatFlux"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%buoy_heat_flux
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('AdvHeatFlux')
                  ! Vertical advective heat flux [W m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "AdvHeatFlux"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%adv_heat_flux
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('DiffSaltFlux')
                  ! Vertical diffusive salt flux [‰ m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "DiffSaltFlux"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%diff_salt_flux
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('AdvSaltFlux')
                  ! Vertical advective salt flux [‰ m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "AdvSaltFlux"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%adv_salt_flux
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('VertVel')
                  ! Vertical velocity [m s-1]
                  self%simdata%output_cfg%output_vars(i)%name = "VertVel"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%w
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('DensityRatio')
                  ! Vertical heat flux [W m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "DensityRatio"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%R_rho
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('LateralInput')
                  ! Vertical heat flux [W m-2]
                  self%simdata%output_cfg%output_vars(i)%name = "LateralInput"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%lateral_input
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .false.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .true.

               case('Tr')
                  ! Tritium concentration [TU]
                  self%simdata%output_cfg%output_vars(i)%name = "Tr"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%Tr
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('HeavyOxygen')
                  ! 18O isotope concentration ratio [18_O/16_O]
                  self%simdata%output_cfg%output_vars(i)%name = "HeavyOxygen"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%heavy_oxygen
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('Deuterium')
                  ! 2H isotope concentration ratio [2_H/1_H]
                  self%simdata%output_cfg%output_vars(i)%name = "Deuterium"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%deuterium
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('LightAr')
                  ! 39Ar concentration as percent modern
                  self%simdata%output_cfg%output_vars(i)%name = "LightArgon"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%light_ar
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('He')
                  ! Tritium concentration [TU]
                  self%simdata%output_cfg%output_vars(i)%name = "He"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%He
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('Ne')
                  ! Tritium concentration [TU]
                  self%simdata%output_cfg%output_vars(i)%name = "Ne"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%Ne
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('Ar')
                  ! Tritium concentration [TU]
                  self%simdata%output_cfg%output_vars(i)%name = "Ar"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%Ar
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case('Kr')
                  ! Tritium concentration [TU]
                  self%simdata%output_cfg%output_vars(i)%name = "Kr"
                  self%simdata%output_cfg%output_vars(i)%values => self%simdata%model%Kr
                  self%simdata%output_cfg%output_vars(i)%volume_grid = .true.
                  self%simdata%output_cfg%output_vars(i)%face_grid = .false.

               case default
                  call warn('Output variable specified in config file not found: ' // trim(output_cfg%output_var_names(i)))
            end select
         end do

      end associate
   end subroutine

   ! Setup model configuration and state vars
   subroutine setup_model(self)
      implicit none
      class(SimstratSimulationFactory) :: self
      ! Integer :: i
      associate (simdata=>self%simdata, &
                 model_cfg=>self%simdata%model_cfg, &
                 model_param=>self%simdata%model_param, &
                 model=>self%simdata%model, &
                 grid=>self%simdata%grid)

         ! Initialize some more values
         if (model_cfg%stability_func == 1) model%cm0 = 0.5625_RK
         if (model_cfg%stability_func == 2) model%cm0 = 0.556171_RK
         model%cde = model%cm0**3
         sig_e = (kappa/model%cm0)**2/(ce2 - ce1)

         model%num(1:grid%nz_grid + 1) = 0.0_RK
         model%nuh(1:grid%nz_grid + 1) = 0.0_RK
         model%nus(1:grid%nz_grid + 1) = 0.0_RK
         model%nug(1:grid%nz_grid + 1) = 0.0_RK

         model%tx = 0.0_RK
         model%ty = 0.0_RK

         model%drag = (kappa/log(1.0_RK + 30/K_s*grid%h(1)/2))**2

         model%gamma = grid%Az(grid%ubnd_fce)/(grid%volume**1.5_RK)/sqrt(rho_0)*model_param%CD

         ! Geothermal heat flux
         if (model_param%fgeo /= 0) then
            allocate(model%fgeo_add(grid%nz_grid))

            model%fgeo_add(1:grid%nz_grid) = model_param%fgeo/rho_0/cp*grid%dAz(1:grid%nz_grid)/grid%Az(2:grid%nz_grid + 1) ! calculation per kg
            if (grid%Az(1) /= 0) then
               model%fgeo_add(1) = model%fgeo_add(1) + 2*model_param%fgeo/rho_0/cp*grid%Az(1)/((grid%Az(1) + grid%Az(2))*grid%h(1))
            end if
         end if

         ! Set up timing
         model%datum = self%simdata%sim_cfg%start_datum
         model%dt = self%simdata%sim_cfg%timestep
         model%model_step_counter = 0
         model%output_counter = 1
      end associate
   end subroutine

   ! Read config of grid and init grid
   subroutine read_grid_config(self)
      implicit none
      class(SimstratSimulationFactory) :: self
      type(GridConfig) :: grid_config
      real(RK), dimension(:) :: z_tmp(self%simdata%model_cfg%max_length_input_data)
      real(RK), dimension(:) :: A_tmp(self%simdata%model_cfg%max_length_input_data)
      integer :: num_read, i
      associate (simdata=>self%simdata, &
                 max_length_input_data=>grid_config%max_length_input_data)

         grid_config%max_length_input_data = self%simdata%model_cfg%max_length_input_data

         if (simdata%input_cfg%grid_input_type == 7) then
            allocate (grid_config%grid_read(grid_config%max_length_input_data))
            ! Read grid
            open (12, status='old', file=simdata%input_cfg%GridName)
            read (12, *)
            do i = 1, max_length_input_data
               read (12, *, end=69) grid_config%grid_read(i)
            end do

69          if(i==max_length_input_data) then
               write(6,*) '[WARNING] ','Only first ',max_length_input_data,' values of file read.'
            else
               call ok('Grid file successfully read')
            end if
            close (12)

            if (i == 2) then ! Constant spacing
               grid_config%nz_grid = int(grid_config%grid_read(1))
               grid_config%equidistant_grid = .TRUE.
            else ! Variable spacing
               grid_config%nz_grid = i - 2
               grid_config%equidistant_grid = .FALSE.
            end if

         ! If grid was read from json
         ! If an array of grid values is given
         else if (simdata%input_cfg%grid_input_type == 3) then
            ! Determine nz_grid from grid
            grid_config%nz_grid = size(simdata%input_cfg%read_grid_array_from_json) - 1
            ! Set flag
            grid_config%equidistant_grid = .FALSE.
            ! Store json-read values in grid_read array to be compatible with rest of the code
            allocate (grid_config%grid_read(grid_config%max_length_input_data))
            grid_config%grid_read(1:grid_config%nz_grid + 1) = simdata%input_cfg%read_grid_array_from_json(1:grid_config%nz_grid + 1)
            call ok('Grid file successfully read')
         ! If the spacing is given
         else
            grid_config%nz_grid = int(simdata%input_cfg%read_grid_value_from_json)
            grid_config%equidistant_grid = .TRUE.
            call ok('Grid file successfully read')
         end if

         ! Read Morphology
         open (11, status='old', file=simdata%input_cfg%MorphName)
         read (11, *) ! Skip header
         do i = 1, max_length_input_data ! Read depth and area
            read (11, *, end=86) z_tmp(i), A_tmp(i)

            ! Check that the uppermost depth is not negative
            if (i==1) then
               if(z_tmp(i)<0) then
                  call error('The uppermost depth in the morphology file is negative, it should be at least 0!')
               end if
            else
            ! Check that depth and area are monotonous
               if (z_tmp(i)>z_tmp(i-1)) then
                  call error('The depths of the morphology input file are not decreasing monotonously.')
               else if(A_tmp(i)>A_tmp(i-1)) then
                  call warn('The lake area increases with depth somewhere. Do you really want this?')
               end if
            end if
         end do

86       if(i==max_length_input_data) then
            write(6,*) '[WARNING] ','Only first ',max_length_input_data,' values of file read.'
         else
            call ok('Morphology file successfully read')
         end if
         close (11)

         num_read = i - 1 ! Number of area values

         allocate (grid_config%z_A_read(num_read), grid_config%A_read(num_read))

         ! Reverse order of values
         do i = 1, num_read
            grid_config%z_A_read(i) = -z_tmp(num_read - i + 1)
            grid_config%A_read(i) = A_tmp(num_read - i + 1)
         end do

         grid_config%max_depth = grid_config%z_A_read(1) - grid_config%z_A_read(num_read) ! depth = max - min depth

         ! Initialize Grid of simdata
         call simdata%grid%init(grid_config)
      end associate
   end subroutine

   ! Read Par file and setup rest of config
   !#######################################################################
   subroutine read_json_par_file(self, ParName)
      !#######################################################################
      implicit none
      class(SimstratSimulationFactory) :: self
      character(kind=CK, len=*), intent(in) :: ParName

      type(json_file) :: par_file
      logical :: found
      integer :: n_children_dummy, index_bs

      ! gfortran cannot handle type bound allocatable character that are passed to subroutine as intent(out)
      ! as a workaround we have to store the values in a local scope allocatable character
      character(kind=CK, len=:), allocatable          :: MorphName, InitName, ForcingName, AbsorpName
      character(kind=CK, len=:), allocatable          :: GridName, zoutName, toutName, PathOut
      character(kind=CK, len=:), allocatable          :: QinpName, QoutName, TinpName, SinpName, TrinpName, HOinpName, DinpName, LAinpName
      character(kind=CK, len=:), allocatable          :: HeinpName, NeinpName, ArinpName, KrinpName
      character(len=20), dimension(:), allocatable :: output_var_names

      associate (input_cfg=>self%simdata%input_cfg, &
                 output_cfg=>self%simdata%output_cfg, &
                 model_cfg=>self%simdata%model_cfg, &
                 aed2_cfg=>self%simdata%aed2_cfg, &
                 sim_cfg=>self%simdata%sim_cfg, &
                 model_param=>self%simdata%model_param)

         ! model%ParName = ParName
         ! Check if inputfile SimstratModelexists
         call check_file_exists(ParName)

         call par_file%initialize(comment_char='!')

         ! Load file or stop if fail
         call par_file%load_file(filename=ParName)
         if (par_file%failed()) then
            call error('Could not read inputfile '//ParName)
         end if

         ! Names of Inputfile
         call par_file%get('Input.Morphology', MorphName, found); input_cfg%MorphName = MorphName; call check_field(found, 'Input.Morphology', ParName)
         call par_file%get('Input.Initial conditions', InitName, found); input_cfg%InitName = InitName; call check_field(found, 'Input.Initial conditions', ParName)
         call par_file%get('Input.Forcing', ForcingName, found); input_cfg%ForcingName = ForcingName; call check_field(found, 'Input.Forcing', ParName)
         call par_file%get('Input.Absorption', AbsorpName, found); input_cfg%AbsorpName = AbsorpName; call check_field(found, 'Input.Absorption', ParName)

         ! Grid information can also be included in par-file
         ! Check type of json entry
         call par_file%info('Input.Grid',found,input_cfg%grid_input_type,n_children_dummy)

         if (input_cfg%grid_input_type == 7) then ! Path name
            call par_file%get('Input.Grid', GridName, found); input_cfg%GridName = GridName; call check_field(found, 'Input.Grid', ParName)
         else if (input_cfg%grid_input_type == 3) then ! Grid depths are given
            call par_file%get('Input.Grid', input_cfg%read_grid_array_from_json, found); call check_field(found, 'Input.Grid', ParName)
         else if (input_cfg%grid_input_type == 5 .or. input_cfg%grid_input_type == 6) then
            call par_file%get('Input.Grid', input_cfg%read_grid_value_from_json, found); call check_field(found, 'Input.Grid', ParName)
         end if

         call par_file%get('Input.Inflow', QinpName, found); input_cfg%QinpName = QinpName; call check_field(found, 'Input.Inflow', ParName)
         call par_file%get('Input.Outflow', QoutName, found); input_cfg%QoutName = QoutName; call check_field(found, 'Input.Outflow', ParName)
         call par_file%get('Input.Inflow temperature', TinpName, found); input_cfg%TinpName = TinpName; call check_field(found, 'Input.Inflow temperature', ParName)
         call par_file%get('Input.Inflow salinity', SinpName, found); input_cfg%SinpName = SinpName; call check_field(found, 'Input.Inflow salinity', ParName)
         call par_file%get('Input.Inflow tritium', TrinpName, found); input_cfg%TrinpName = TrinpName; call check_field(found, 'Input.Inflow tritium', ParName)
         call par_file%get('Input.Inflow heavy oxygen', HOinpName, found); input_cfg%HOinpName = HOinpName; call check_field(found, 'Input.Inflow heavy oxygen', ParName)
         call par_file%get('Input.Inflow deuterium', DinpName, found); input_cfg%DinpName = DinpName; call check_field(found, 'Input.Inflow deuterium', ParName)
         call par_file%get('Input.Inflow light ar', LAinpName, found); input_cfg%LAinpName = LAinpName; call check_field(found, 'Input.Inflow light ar', ParName)
         
         call par_file%get('Input.Inflow helium', HeinpName, found); input_cfg%HeinpName = HeinpName; call check_field(found, 'Input.Inflow helium', ParName)
         call par_file%get('Input.Inflow neon', NeinpName, found); input_cfg%NeinpName = NeinpName; call check_field(found, 'Input.Inflow neon', ParName)
         call par_file%get('Input.Inflow argon', ArinpName, found); input_cfg%ArinpName = ArinpName; call check_field(found, 'Input.Inflow argon', ParName)
         call par_file%get('Input.Inflow krypton', KrinpName, found); input_cfg%KrinpName = KrinpName; call check_field(found, 'Input.Inflow krypton', ParName)

         ! Path to output folder
         call par_file%get('Output.Path', PathOut, found); call check_field(found, 'Output.Path', ParName)

         ! Transform backslashes to slash
         do while(scan(PathOut,'\\')>0)
            index_bs = scan(PathOut,'\\')
            PathOut(index_bs:index_bs) = '/'
         end do

         ! Remove trailing slashes at the end
         if (len(PathOut) == scan(trim(PathOut),"/", BACK= .true.)) then
            output_cfg%PathOut = PathOut(1:len(PathOut) - 1)
         else
            output_cfg%PathOut = trim(PathOut)
         end if

         ! Output depth reference
         call par_file%get("Output.OutputDepthReference", output_cfg%output_depth_reference, found); call check_field(found, 'Output.OutputDepthReference', ParName)
         if (.not.(output_cfg%output_depth_reference == 'surface' .or. output_cfg%output_depth_reference == 'bottom')) then
            call error('Invalid field "Output.OutputDepthReference" in par-file.')
         end if

         ! Output depths
         ! Check type of input: string for path, array for depth list and integer for interval
         call par_file%info('Output.Depths',found,output_cfg%output_depth_type,n_children_dummy)
         ! Treat different input possibilities for output depths
         if (output_cfg%output_depth_type == 7) then ! Path name
            call par_file%get('Output.Depths', zoutName, found); output_cfg%zoutName = zoutName; call check_field(found, 'Output.Depths', ParName)
         else if (output_cfg%output_depth_type == 3) then ! Output depths are given
            call par_file%get('Output.Depths', output_cfg%zout_read, found); call check_field(found, 'Output.Depths', ParName)
            if (output_cfg%output_depth_reference == 'surface') then
               call reverse_in_place(output_cfg%zout_read)
            end if
            output_cfg%depth_interval = 0
         else if (output_cfg%output_depth_type == 5 .or. output_cfg%output_depth_type == 6) then ! Output interval is given
            call par_file%get('Output.Depths', output_cfg%depth_interval, found); call check_field(found, 'Output.Depths', ParName)

         else
            call error('Invalid field "Output.Depths" in par-file.')
         end if

         ! Output times
         ! Check type of input: string for path, array for depth list and integer for interval
         call par_file%info('Output.Times',found,output_cfg%output_time_type,n_children_dummy)
         ! Treat different input possibilities for output times
         if (output_cfg%output_time_type == 7) then ! Path name
            call par_file%get('Output.Times', toutName, found); output_cfg%toutName = toutName; call check_field(found, 'Output.Times', ParName)
         else if (output_cfg%output_time_type == 3) then ! Output depths are given
            call par_file%get('Output.Times', output_cfg%tout, found); call check_field(found, 'Output.Times', ParName)
            output_cfg%thinning_interval = 0

         else if (output_cfg%output_time_type == 5 .or. output_cfg%output_time_type == 6) then ! Output interval is given
            call par_file%get('Output.Times', output_cfg%thinning_interval_read, found); call check_field(found, 'Output.Times', ParName)
            output_cfg%thinning_interval = int(output_cfg%thinning_interval_read)

            ! This code is needed for line 141 in simstrat.f90. Not very elegant.. might change in the future.
            allocate(output_cfg%tout(1))
            output_cfg%tout(1) = output_cfg%thinning_interval

         else
            call error('Invalid field "Output.Times" in par-file.')
         end if
         call par_file%get('Output.All', output_cfg%output_all, found); call check_field(found, 'Output.All', ParName)
         if (.not. output_cfg%output_all) then
            call par_file%get('Output.Variables', output_var_names, found); output_cfg%output_var_names = output_var_names; call check_field(found, 'Output.Variables', ParName)
         end if

         ! Model configuration
         call par_file%get("ModelConfig.MaxLengthInputData", model_cfg%max_length_input_data, found);
         if (.not. found) then
            model_cfg%max_length_input_data = 1000
            call warn('Variable "ModelConfig.MaxLengthInputData" is not set. Assume a value of 1000')
         end if
         call par_file%get("ModelConfig.CoupleAED2", model_cfg%couple_aed2,found);
         if (.not. found) then
            model_cfg%couple_aed2 = .false.
            call warn('Variable "ModelConfig.CoupleAED2" is not set. Assume you do not want to couple simstrat with aed2.')
         end if
         call par_file%get("ModelConfig.TurbulenceModel", model_cfg%turbulence_model, found); call check_field(found, 'ModelConfig.TurbulenceModel', ParName)
         call par_file%get("ModelConfig.ApparentDiffusivity", model_cfg%apparent_diffusivity, found); call check_field(found, 'ModelConfig.ApparentDiffusivity', ParName)
         call par_file%get("ModelConfig.SplitSeicheParameter", model_cfg%split_a_seiche, found); call check_field(found, 'ModelConfig.SplitSeicheParameter', ParName)
         call par_file%get("ModelConfig.StabilityFunction", model_cfg%stability_func, found); call check_field(found, 'ModelConfig.StabilityFunction', ParName)
         call par_file%get("ModelConfig.FluxCondition", model_cfg%flux_condition, found); call check_field(found, 'ModelConfig.FluxCondition', ParName)
         call par_file%get("ModelConfig.Forcing", model_cfg%forcing_mode, found); call check_field(found, 'ModelConfig.Forcing', ParName)
         call par_file%get("ModelConfig.UserDefinedWaterAlbedo", model_cfg%user_defined_water_albedo, found); call check_field(found, 'ModelConfig.UserDefinedWaterAlbedo', ParName)
         call par_file%get("ModelConfig.UseFilteredWind", model_cfg%use_filtered_wind, found); call check_field(found, 'ModelConfig.UseFilteredWind', ParName)
         call par_file%get("ModelConfig.SeicheNormalization", model_cfg%seiche_normalization, found); call check_field(found, 'ModelConfig.SeicheNormalization', ParName)
         call par_file%get("ModelConfig.WindDragModel", model_cfg%wind_drag_model, found); call check_field(found, 'ModelConfig.WindDragModel', ParName)
         call par_file%get("ModelConfig.InflowMode", model_cfg%inflow_mode, found); call check_field(found, 'ModelConfig.InflowMode', ParName)
         call par_file%get("ModelConfig.PressureGradients", model_cfg%pressure_gradients, found); call check_field(found, 'ModelConfig.PressureGradients', ParName)
         call par_file%get("ModelConfig.IceModel", model_cfg%ice_model, found); call check_field(found, 'ModelConfig.IceModel', ParName)
         call par_file%get("ModelConfig.SnowModel", model_cfg%snow_model, found); call check_field(found, 'ModelConfig.SnowModel', ParName)

         ! AED2 configuration (or another biogeochemical model if implemented)
         if (model_cfg%couple_aed2) then
            call par_file%get("AED2Config.AED2ConfigFile", aed2_cfg%aed2_config_file,found); call check_field(found, 'AED2Config.AED2ConfigFile', ParName)
            call par_file%get("AED2Config.PathAED2initial", aed2_cfg%path_aed2_initial,found); call check_field(found, 'AED2Config.PathAED2initial', ParName)
            call par_file%get("AED2Config.PathAED2inflow", aed2_cfg%path_aed2_inflow,found); call check_field(found, 'AED2Config.PathAED2inflow',ParName)
            call par_file%get("AED2Config.PathAED2output", aed2_cfg%path_aed2_output,found); call check_field(found, 'AED2Config.PathAED2output',ParName)
            call par_file%get("AED2Config.ParticleMobility", aed2_cfg%particle_mobility,found); call check_field(found, 'AED2Config.ParticleMobility', ParName)
            call par_file%get("AED2Config.BioshadeFeedback", aed2_cfg%bioshade_feedback,found); call check_field(found, 'AED2Config.BioshadeFeedback', ParName)
            call par_file%get("AED2Config.BackgroundExtinction", aed2_cfg%background_extinction,found); call check_field(found, 'AED2Config.BackgroundExtinction', ParName)
            call par_file%get("AED2Config.BenthicMode", aed2_cfg%benthic_mode,found); call check_field(found, 'AED2Config.BenthicMode', ParName)
            call par_file%get("AED2Config.NZones", aed2_cfg%n_zones,found); call check_field(found, 'AED2Config.NZones', ParName)
            call par_file%get("AED2Config.ZoneHeights", aed2_cfg%zone_heights,found); call check_field(found, 'AED2Config.ZoneHeights', ParName)
         end if

         !Model Parameter
         call par_file%get("ModelParameters.lat", model_param%Lat, found); call check_field(found, 'ModelParameters.lat', ParName)
         call par_file%get("ModelParameters.p_air", model_param%p_air, found); call check_field(found, 'ModelParameters.p_air', ParName)
         call par_file%get("ModelParameters.a_seiche", model_param%a_seiche, found); call check_field(found, 'ModelParameters.a_seiche', ParName)
         if (model_cfg%split_a_seiche) then
            call par_file%get("ModelParameters.a_seiche_w", model_param%a_seiche_w, found); call check_field(found, 'ModelParameters.a_seiche_w', ParName)
            call par_file%get("ModelParameters.strat_sumr", model_param%strat_sumr, found); call check_field(found, 'ModelParameters.strat_sumr', ParName)
         end if 
         call par_file%get("ModelParameters.q_nn", model_param%q_NN, found); call check_field(found, 'ModelParameters.q_nn', ParName)
         call par_file%get("ModelParameters.f_wind", model_param%f_wind, found); call check_field(found, 'ModelParameters.f_wind', ParName)

         ! C10 is a physical parameter on the order of e-3 if wind drag model is 1. Conversely, it is a calibration parameter
         ! with a value around 1 for wind drag models 2 and 3. This fact can lead to confusion and thus we check for realistic
         ! input values depending on the wind drag model chosen.
         call par_file%get("ModelParameters.c10", model_param%C10_constant, found); call check_field(found, 'ModelParameters.c10', ParName)
         if (model_cfg%wind_drag_model == 1 .and. model_param%C10_constant > 0.1) then
            call error('The input value of C10 is too high to be physically possible. Choose a lower value or change the wind drag model to 2 or 3 if you meant to use C10 as a calibration parameter.')
         else if (model_cfg%wind_drag_model > 1 .and. model_param%C10_constant < 0.1) then
            call error('The input value of C10 is too low to serve as calibration parameter. Maybe you intended it as physical parameter but then you need to use wind drag model = 1.')
         end if

         call par_file%get("ModelParameters.cd", model_param%CD, found); call check_field(found, 'ModelParameters.cd', ParName)
         call par_file%get("ModelParameters.hgeo", model_param%fgeo, found); call check_field(found, 'ModelParameters.hgeo', ParName)
         call par_file%get("ModelParameters.p_sw", model_param%p_sw, found); call check_field(found, 'ModelParameters.p_sw', ParName)
         call par_file%get("ModelParameters.p_lw", model_param%p_lw, found); call check_field(found, 'ModelParameters.p_lw', ParName)
         call par_file%get("ModelParameters.p_windf", model_param%p_windf, found); call check_field(found, 'ModelParameters.p_windf', ParName)
         call par_file%get("ModelParameters.beta_sol", model_param%beta_sol, found); call check_field(found, 'ModelParameters.beta_sol', ParName)

         ! Get water albedo value if switch is on
         if (model_cfg%user_defined_water_albedo) then
            call par_file%get("ModelParameters.wat_albedo", model_param%wat_albedo, found); call check_field(found, 'ModelParameters.wat_albedo', ParName)
         end if
         if (model_cfg%ice_model == 1) then
           call par_file%get("ModelParameters.p_albedo", model_param%p_albedo, found); call check_field(found, 'ModelParameters.p_albedo', ParName)
           call par_file%get("ModelParameters.freez_temp", model_param%freez_temp, found); call check_field(found, 'ModelParameters.freez_temp', ParName)
           call par_file%get("ModelParameters.snow_temp", model_param%snow_temp, found); call check_field(found, 'ModelParameters.snow_temp', ParName)
         end if

         ! Simulation Parameter
         call par_file%get("Simulation.Timestep s", sim_cfg%timestep, found); call check_field(found, 'Simulation.Timestep s', ParName)
         call par_file%get("Simulation.Start year", sim_cfg%start_year, found); call check_field(found, 'Simulation.Start year', ParName)
         call par_file%get("Simulation.Start d", sim_cfg%start_datum, found); call check_field(found, 'Simulation.Start d', ParName)
         call par_file%get("Simulation.End d", sim_cfg%end_datum, found); call check_field(found, 'Simulation.End d', ParName)
         call par_file%get("Simulation.DisplaySimulation", sim_cfg%disp_simulation, found); call check_field(found, 'Simulation.DisplaySimulation', ParName)

         call par_file%destroy()

         call ok('Configuration: '//trim(ParName))
         !  end if
      end associate

   end subroutine read_json_par_file

   ! Read initial data and set in state
   subroutine read_initial_data(self)
      implicit none
      class(SimstratSimulationFactory) :: self

      ! Local variables
      real(RK) :: z_read(self%simdata%model_cfg%max_length_input_data), U_read(self%simdata%model_cfg%max_length_input_data), V_read(self%simdata%model_cfg%max_length_input_data)
      real(RK) :: T_read(self%simdata%model_cfg%max_length_input_data), S_read(self%simdata%model_cfg%max_length_input_data), Tr_read(self%simdata%model_cfg%max_length_input_data)
      real(RK) :: HO_read(self%simdata%model_cfg%max_length_input_data), D_read(self%simdata%model_cfg%max_length_input_data), LA_read(self%simdata%model_cfg%max_length_input_data), k_read(self%simdata%model_cfg%max_length_input_data), eps_read(self%simdata%model_cfg%max_length_input_data)
      real(RK) :: He_read(self%simdata%model_cfg%max_length_input_data), Ne_read(self%simdata%model_cfg%max_length_input_data), Ar_read(self%simdata%model_cfg%max_length_input_data), Kr_read(self%simdata%model_cfg%max_length_input_data)
      real(RK) :: z_ini_depth
      integer :: i, num_read

      associate (grid=>self%simdata%grid, &
                 model=>self%simdata%model, &
                 nz_occupied=>self%simdata%grid%nz_occupied, &
                 max_length_input_data=>self%simdata%model_cfg%max_length_input_data)

         ! Read file
         open (13, status='old', file=self%simdata%input_cfg%InitName) ! Opens initial conditions file
         read (13, *) ! Skip header
         do i = 1, max_length_input_data ! Read initial u,v,T, etc
            read (13, *, end=99) z_read(i), U_read(i), V_read(i), T_read(i), S_read(i), Tr_read(i), HO_read(i), D_read(i), LA_read(i), He_read(i), Ne_read(i), Ar_read(i), Kr_read(i), k_read(i), eps_read(i)
            if (z_read(i)>0) then
               call error('One or several input depths of initial conditions are positive.')
            end if
         end do
99       num_read = i-1                               ! Number of valuInitNamees
         if (num_read < 1) then
            call error('Unable to read initial conditions files (no data found).')
         else if(num_read==max_length_input_data) then
            write(6,*) '[ERROR] ','Only first ',max_length_input_data,' values of initial data file read.'
         end if

         close (13)
         do i = 1, num_read
            z_read(i) = abs(z_read(i)) ! Make depths positive
         end do
         z_ini_depth = z_read(1) ! Initial depth (top-most)

         ! Update actual filled z in grid
         call grid%update_depth(z_ini_depth)

         ! Set initial lake level
         allocate(grid%lake_level)  ! Allocation is needed because lake_level is defined as pointer for use with AED2
         grid%lake_level = grid%z_face(grid%ubnd_fce)

         grid%lake_level_old = grid%lake_level

         ! Reverse arrays
         call reverse_in_place(z_read(1:num_read))
         z_read(1:num_read) = grid%z_zero - z_read(1:num_read)
         call reverse_in_place(U_read(1:num_read))
         call reverse_in_place(V_read(1:num_read))
         call reverse_in_place(T_read(1:num_read))
         call reverse_in_place(S_read(1:num_read))
         call reverse_in_place(Tr_read(1:num_read))
         call reverse_in_place(HO_read(1:num_read))
         call reverse_in_place(D_read(1:num_read))
         call reverse_in_place(LA_read(1:num_read))
         call reverse_in_place(k_read(1:num_read))
         call reverse_in_place(eps_read(1:num_read))

         call reverse_in_place(He_read(1:num_read))
         call reverse_in_place(Ne_read(1:num_read))
         call reverse_in_place(Ar_read(1:num_read))
         call reverse_in_place(Kr_read(1:num_read))

         if (num_read == 1) then
            call warn('Only one row! Water column will be initially homogeneous.')
            model%U = U_read(1)
            model%V = V_read(1)
            model%T = T_read(1)
            model%S = S_read(1)
            model%Tr = Tr_read(1)
            model%heavy_oxygen = HO_read(1)
            model%deuterium = D_read(1)
            model%light_ar = LA_read(1)

            model%He = He_read(1)
            model%Ne = Ne_read(1)
            model%Ar = Ar_read(1)
            model%Kr = Kr_read(1)

            model%k = k_read(1)
            model%eps = eps_read(1)
         else
            ! Interpolate variables UVTS on central grid and store
            call grid%interpolate_to_vol(z_read, U_read, num_read, model%U)
            call grid%interpolate_to_vol(z_read, V_read, num_read, model%V)
            call grid%interpolate_to_vol(z_read, T_read, num_read, model%T)
            call grid%interpolate_to_vol(z_read, S_read, num_read, model%S)
            call grid%interpolate_to_vol(z_read, Tr_read, num_read, model%Tr)
            call grid%interpolate_to_vol(z_read, HO_read, num_read, model%heavy_oxygen)
            call grid%interpolate_to_vol(z_read, D_read, num_read, model%deuterium)
            call grid%interpolate_to_vol(z_read, LA_read, num_read, model%light_ar)

            call grid%interpolate_to_vol(z_read, He_read, num_read, model%He)
            call grid%interpolate_to_vol(z_read, Ne_read, num_read, model%Ne)
            call grid%interpolate_to_vol(z_read, Ar_read, num_read, model%Ar)
            call grid%interpolate_to_vol(z_read, Kr_read, num_read, model%Kr)

            ! Interpolate k/eps on upper grid and store
            call grid%interpolate_to_face(z_read, k_read, num_read, model%k)
            call grid%interpolate_to_face(z_read, eps_read, num_read, model%eps)
         end if

         call ok('Initial data file successfully read')

      end associate
   end subroutine

   subroutine check_field(found, field_name, file_name)
      implicit none

      logical, intent(in) :: found
      character(len=*), intent(in) :: field_name, file_name

      if (.not. found) then
         call error('Field '//field_name//' not found in '//file_name)
      end if
   end subroutine check_field

end module strat_inputfile
