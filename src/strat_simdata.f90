!     +---------------------------------------------------------------+
!     |  Data structure definitions for simulation data
!     +---------------------------------------------------------------+

module strat_simdata
   use strat_kinds
   use strat_grid
   use strat_consts
   implicit none
   private

   ! All Input files
   type, public :: InputConfig
      character(len=:), allocatable          :: MorphName
      character(len=:), allocatable          :: InitName
      character(len=:), allocatable          :: ForcingName
      character(len=:), allocatable          :: AbsorpName
      character(len=:), allocatable          :: GridName
      character(len=:), allocatable          :: QinpName
      character(len=:), allocatable          :: QoutName
      character(len=:), allocatable          :: TinpName
      character(len=:), allocatable          :: SinpName
      character(len=:), allocatable          :: TrinpName
      character(len=:), allocatable          :: HOinpName
      character(len=:), allocatable          :: DinpName
      character(len=:), allocatable          :: LAinpName
      character(len=:), allocatable          :: HeinpName
      character(len=:), allocatable          :: NeinpName
      character(len=:), allocatable          :: ArinpName
      character(len=:), allocatable          :: KrinpName
      real(RK), dimension(:), allocatable    :: read_grid_array_from_json
      real(RK) :: read_grid_value_from_json
      integer :: grid_input_type
   end type

   ! Definition of a variable to log
   type, public :: LogVariable
      character(len=:), allocatable :: name
      real(RK), dimension(:), pointer :: values
      real(RK), pointer :: values_surf 
      logical :: volume_grid, face_grid
   end type

   ! Definition of a AED2 variable to log
   type, public :: LogVariableAED2
      character(len=48), pointer, dimension(:) :: names
      real(RK), dimension(:,:), pointer :: values
   end type

   ! Logging configuration
   type, public :: OutputConfig
      character(len=:), allocatable :: PathOut
      character(len=:), allocatable :: zoutName
      character(len=:), allocatable :: toutName
      character(len=:), allocatable :: output_depth_reference
      real(RK), dimension(:), allocatable :: zout, zout_read
      real(RK), dimension(:), allocatable :: tout
      real(RK), dimension(:), allocatable :: n_timesteps_between_tout
      real(RK), dimension(:), allocatable :: adjusted_timestep
      logical :: write_to_file, output_all
      integer :: number_output_vars
      character(len=20), dimension(:), allocatable :: output_var_names ! Names of output variables
      class(LogVariable), dimension(:), allocatable :: output_vars
      class(LogVariableAED2), allocatable :: output_vars_aed2

      integer :: output_time_type, output_depth_type, thinning_interval
      real(RK) :: depth_interval, thinning_interval_read ! thinning_interval_read is a real to make sure that also values
      ! like 72.0 can be read (which are interpreted as a double)
   end type

   ! Simulation configuration
   type, public :: SimConfig
      integer :: timestep
      integer :: start_year
      real(RK) :: start_datum
      real(RK) :: end_datum
      integer :: disp_simulation
   end type

   ! Model configuration (read from file)
   type, public :: ModelConfig
      integer :: max_length_input_data
      logical :: couple_aed2
      integer :: turbulence_model
      logical :: apparent_diffusivity
      logical :: split_a_seiche
      integer :: stability_func
      integer :: flux_condition
      integer :: forcing_mode
      logical :: user_defined_water_albedo
      logical :: use_filtered_wind
      integer :: seiche_normalization
      integer :: wind_drag_model
      integer :: inflow_mode
      integer :: pressure_gradients
      logical :: salinity_transport
      integer :: ice_model
      integer :: snow_model
   end type

   ! AED2 configuration (read from file)
   type, public :: AED2Config
      character(len=:), allocatable :: aed2_config_file
      character(len=:), allocatable :: path_aed2_initial
      character(len=:), allocatable :: path_aed2_inflow
      character(len=:), allocatable :: path_aed2_output
      logical :: particle_mobility
      logical :: bioshade_feedback
      real(RK) :: background_extinction
      integer :: benthic_mode
      integer :: n_zones
      real(RK), dimension(:), allocatable :: zone_heights
   end type

   ! Model params (read from file)
   type, public :: ModelParam
      real(RK) :: Lat
      real(RK) :: p_air
      real(RK) :: a_seiche
      real(RK) :: a_seiche_w
      real(RK) :: strat_sumr
      real(RK) :: q_NN
      real(RK) :: f_wind
      real(RK) :: C10_constant
      real(RK) :: CD
      real(RK) :: fgeo
      real(RK) :: p_sw
      real(RK) :: p_lw
      real(RK) :: p_windf
      real(RK) :: beta_sol
      real(RK) :: wat_albedo
      real(RK) :: p_albedo
      real(RK) :: freez_temp
      real(RK) :: snow_temp
   end type

   ! Model state (self is actually the simulation data!!!)
   type, public :: ModelState
      ! Iteration variables
      integer :: i, j, output_counter, model_step_counter
      real(RK) :: datum, dt
      logical :: first_timestep = .true.

      ! Variables located on z_cent grid
      ! Note that for these variables the value at 0 z.b. U(0) is not used
      real(RK), dimension(:), allocatable :: U, V, co2, ch4 ! Water velocities [m/s]
      real(RK), dimension(:), pointer :: T, S, Tr, heavy_oxygen, deuterium, light_AR, R_rho ! Temperature [°C], Salinity [‰], 18O [-], D [-], 39Ar [-], Density ratio [-]
      real(RK), dimension(:), pointer :: He, Ne, Ar, Kr ! Noble gas concentrations [ccSTP/g]
      real(RK), dimension(:), allocatable :: dS, dTr, dHO, dD, dLA ! Source/sink for salinity, tritium, heavy oxygen, deuterium, 39Ar
      real(RK), dimension(:), pointer :: dHe, dNe, dAr, dKr ! Source/sink for noble gas concentrations [ccSTP/g]
      real(RK), dimension(:, :), allocatable :: Q_inp ! Horizontal inflow [m^3/s]
      real(RK), dimension(:), pointer :: rho ! Water density [kg/m^3]
      real(RK), dimension(:), allocatable :: diff_heat_flux, buoy_heat_flux, adv_heat_flux ! Heat flux [W m-2]
      real(RK), dimension(:), allocatable :: diff_salt_flux, adv_salt_flux ! Salt_flux [permil m-2]
      real(RK), dimension(:,:), pointer :: AED2_state ! State matrix of AED2 variables
      real(RK), dimension(:,:), pointer :: AED2_diagstate ! State matrix of AED2 variables
      character(len=48), dimension(:), pointer :: AED2_names ! Names of AED2 state variables used in the simulation
      character(len=48), dimension(:), pointer :: AED2_diagnames ! Names of AED2 state variables used in the simulation
      integer :: n_pH
      integer :: keps_counter
   
      ! Variables located on z_upp grid
      real(RK), dimension(:), allocatable :: k, ko ! Turbulent kinetic energy (TKE) [J/kg]
      real(RK), dimension(:), allocatable :: avh
      real(RK), dimension(:), allocatable :: eps ! TKE dissipation rate [W/kg]
      real(RK), dimension(:), allocatable :: num, nuh, nus, nug, nut ! Turbulent viscosity (momentum) and diffusivity (temperature, salinity, gases, water)
      real(RK), dimension(:), allocatable :: nu_he, nu_ne, nu_ar, nu_kr ! Noble gas diffusivities
      real(RK), dimension(:), allocatable :: P, B ! Shear stress production [W/kg], buoyancy production [W/kg]
      real(RK), dimension(:), allocatable :: NN ! Brunt-Väisälä frequency [s-2]
      real(RK), dimension(:), allocatable :: cmue1, cmue2 ! Model constants
      real(RK), dimension(:), allocatable :: P_Seiche ! Production of TKE [W/kg] and seiche energy [J]
      real(RK) :: E_Seiche
      real(RK) :: gamma ! Proportionality constant for loss of seiche energy

      real(RK), dimension(:), allocatable :: absorb ! Absorption coeff [m-1]
      real(RK), dimension(:), pointer :: absorb_vol ! Absorption coeff on vol grid [m-1]
      real(RK) :: u10, v10, Wf, Vap_atm ! Wind speeds, wind factor
      real(RK), pointer :: uv10 ! pointer attribute needed for AED2
      real(RK), pointer :: rain ! pointer attribute needed for AED2, rain is not calculated in Simstrat for the moment, but required by AED2
      real(RK) :: drag, u_taus ! Drag
      real(RK), pointer :: u_taub ! pointer attribute needed for AED2
      real(RK) :: tx, ty ! Shear stress
      real(RK) :: C10 ! Wind drag coefficient
      real(RK) :: SST, heat, heat_snow, heat_ice, heat_snowice! Sea surface temperature and heat flux

      real(RK) :: T_atm ! Air temp at surface
      real(RK), dimension(:), allocatable :: rad, rad_vol ! Solar radiation (in water)
      real(RK), dimension(:), allocatable :: Q_vert, lateral_input, w ! Vertical exchange between boxes (integrated and not-integrated over depth)
      real(RK), dimension(9,12) :: albedo_data  ! Experimental monthly albedo data for determination of current water albedo
      real(RK) :: albedo_water   ! Current water albedo
      integer :: lat_number ! Latitude band (used for determination of albedo)

      ! Snow and Ice
      real(RK), allocatable :: snow_h ! Snow layer height [m]
      real(RK), allocatable :: total_ice_h ! Total ice layer height [m]
      real(RK), allocatable :: black_ice_h ! Black ice layer height [m]
      real(RK), allocatable :: white_ice_h ! Snowice layer height [m]
      real(RK) :: snow_dens ! Snow density [kg m-3]
      real(RK) :: ice_temp ! Ice temperature [°C]
      real(RK) :: precip ! Precipiation in water eqvivalent hight [m]

      !For saving heatflux
      real(RK), allocatable :: ha ! Incoming long wave [W m-2]
      real(RK), allocatable :: hw ! Outgoing long wave [W m-2]
      real(RK), allocatable :: hk ! Sensible flux [W m-2]
      real(RK), allocatable :: hv ! Latent heat [W m-2]
      real(RK), pointer :: rad0 !  Solar radiation at surface  [W m-2]
   
      real(RK) :: cde, cm0
      real(RK) ::  fsed
      real(RK), dimension(:), allocatable     :: fgeo_add
      integer :: n_AED2, n_AED2_diag


   contains
      procedure, pass :: init => model_state_init
   end type

   ! Structure that encapsulates a full program state
   type, public :: SimulationData
      type(InputConfig), public   :: input_cfg
      type(OutputConfig), public  :: output_cfg
      type(SimConfig), public     :: sim_cfg
      type(ModelConfig), public   :: model_cfg
      type(AED2Config), public    :: aed2_cfg
      type(ModelParam), public    :: model_param
      type(ModelState), public    :: model
      type(StaggeredGrid), public :: grid
   contains
      procedure, pass :: init => simulation_data_init
   end type

contains
   subroutine simulation_data_init(self, state_size)
      class(SimulationData), intent(inout) :: self
      integer, intent(in) :: state_size
      ! Init model data structures
      call self%model%init(state_size)

   end subroutine

   ! Allocates all arrays of the model state in the correct size
   subroutine model_state_init(self, state_size)
      class(ModelState), intent(inout) :: self
      integer, intent(in) :: state_size

      ! Values on volume grid
      ! Important: Size is smaller than vars on upper grid.
      !            https://en.wikipedia.org/wiki/Off-by-one_error#Fencepost_error ;-)
      allocate (self%U(state_size))
      allocate (self%V(state_size))
      allocate (self%co2(state_size))
      allocate (self%ch4(state_size))
      allocate (self%T(state_size))
      allocate (self%S(state_size))
      allocate (self%Tr(state_size))
      allocate (self%heavy_oxygen(state_size))
      allocate (self%deuterium(state_size))
      allocate (self%light_ar(state_size))
      allocate (self%He(state_size))
      allocate (self%Ne(state_size))
      allocate (self%Ar(state_size))
      allocate (self%Kr(state_size))
      allocate (self%R_rho(state_size))
      allocate (self%dS(state_size))
      allocate (self%dTr(state_size))
      allocate (self%dHO(state_size))
      allocate (self%dD(state_size))
      allocate (self%dLA(state_size))
      allocate (self%dHe(state_size))
      allocate (self%dNe(state_size))
      allocate (self%dAr(state_size))
      allocate (self%dKr(state_size))
      allocate (self%rho(state_size))
      allocate (self%avh(state_size))
      allocate (self%diff_heat_flux(state_size))
      allocate (self%buoy_heat_flux(state_size))
      allocate (self%adv_heat_flux(state_size))
      allocate (self%diff_salt_flux(state_size))
      allocate (self%adv_salt_flux(state_size))

      ! Values on z_upp grid
      allocate (self%k(state_size + 1))
      allocate (self%ko(state_size + 1))
      allocate (self%eps(state_size + 1))
      allocate (self%num(state_size + 1))
      allocate (self%nuh(state_size + 1))
      allocate (self%nus(state_size + 1))
      allocate (self%nug(state_size + 1))
      allocate (self%nut(state_size + 1))
      allocate (self%nu_he(state_size + 1))
      allocate (self%nu_ne(state_size + 1))
      allocate (self%nu_ar(state_size + 1))
      allocate (self%nu_kr(state_size + 1))
      allocate (self%P(state_size + 1))
      allocate (self%B(state_size + 1))
      allocate (self%NN(state_size + 1))
      allocate (self%cmue1(state_size + 1))
      allocate (self%cmue2(state_size + 1))
      allocate (self%P_Seiche(state_size + 1))

      allocate (self%absorb(state_size + 1))
      allocate (self%absorb_vol(state_size))
      allocate (self%rad(state_size + 1))
      allocate (self%rad_vol(state_size))
      allocate (self%Q_vert(state_size + 1))
      allocate (self%lateral_input(state_size + 1))
      allocate (self%w(state_size + 1))

      allocate (self%snow_h)
      allocate (self%total_ice_h)
      allocate (self%black_ice_h)
      allocate (self%white_ice_h)

      allocate (self%ha)
      allocate (self%hw)
      allocate (self%hk)
      allocate (self%hv)
      allocate (self%rad0)

      ! Init to zero
      self%U = 0.0_RK
      self%V = 0.0_RK
      self%T = 0.0_RK
      self%S = 0.0_RK
      self%Tr = 0.0_RK
      self%heavy_oxygen = 0.0_RK
      self%deuterium = 0.0_RK
      self%light_ar = 0.0_RK
      self%He = 0.0_RK
      self%Ne = 0.0_RK
      self%Ar = 0.0_RK
      self%Kr = 0.0_RK
      self%R_rho = 0.0_RK
      self%dS = 0.0_RK
      self%dTr = 0.0_RK
      self%dHO = 0.0_RK
      self%dD = 0.0_RK
      self%dLA = 0.0_RK
      self%dHe = 0.0_RK
      self%dNe = 0.0_RK
      self%dAr = 0.0_RK
      self%dKr = 0.0_RK
      self%rho = 0.0_RK
      self%co2 = 0.0_RK
      self%ch4 = 0.0_RK
      self%buoy_heat_flux = 0.0_RK
      self%adv_heat_flux = 0.0_RK
      self%diff_salt_flux = 0.0_RK
      self%adv_salt_flux = 0.0_RK

      self%k = 0.0_RK
      self%ko = 0.0_RK
      self%eps = 0.0_RK
      self%num = 0.0_RK
      self%nuh = 0.0_RK
      self%nus = 0.0_RK
      self%nug = 0.0_RK
      self%nut = 0.0_RK
      self%nu_he = 0.0_RK
      self%nu_ne = 0.0_RK
      self%nu_ar = 0.0_RK
      self%nu_kr = 0.0_RK
      self%P = 0.0_RK
      self%B = 0.0_RK
      self%NN = 0.0_RK
      self%cmue1 = 0.0_RK
      self%cmue2 = 0.0_RK
      self%P_Seiche = 0.0_RK
      self%E_Seiche = 0.0_RK

      self%absorb = 0.0_RK
      self%absorb_vol = 0.0_RK
      self%rad = 0.0_RK
      self%rad_vol = 0.0_RK
      self%Q_vert = 0.0_RK
      self%lateral_input = 0.0_RK
      self%w = 0.0_RK

      self%snow_h = 0.0_RK
      self%total_ice_h = 0.0_RK
      self%black_ice_h = 0.0_RK
      self%white_ice_h = 0.0_RK
      self%ice_temp = 0.0_RK
      self%snow_dens = rho_s_0
      self%precip = 0.0_RK 
   
      self%ha = 0.0_RK
      self%hw = 0.0_RK
      self%hk = 0.0_RK 
      self%hv = 0.0_RK 
      self%rad0 = 0.0_RK
      self%n_pH = 0
      self%keps_counter = 0

      ! init pointers
      allocate(self%uv10)
      self%uv10 = 0.0_RK
      allocate(self%rain)
      self%rain = 0.0_RK
      allocate(self%u_taub)
      self%u_taub = 0.0_RK

   end subroutine

end module strat_simdata
