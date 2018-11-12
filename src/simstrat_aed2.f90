!     +---------------------------------------------------------------+
!     |  Simstrat - AED2 interface
!     +---------------------------------------------------------------+

module simstrat_aed2
   use strat_simdata
   use strat_grid
   use utilities
   use aed2_common
   use aed2_core

   implicit none
   private

   type, public :: SimstratAED2
      class(AED2Config), pointer :: aed2_cfg
      class(StaggeredGrid), pointer :: grid

      real(RK),allocatable,dimension(:) :: lKw    !# background light attenuation (m**-1)

      !# Arrays for state and diagnostic variables
      real(RK),allocatable,dimension(:,:) :: cc !# water quality array: nlayers, nvars
      real(RK),allocatable,dimension(:,:) :: cc_diag
      real(RK),allocatable,dimension(:) :: cc_diag_hz
      real(RK),allocatable,dimension(:) :: tss
      real(RK),allocatable,dimension(:) :: sed_zones

      !# Arrays for work, vertical movement, and cross-boundary fluxes
      real(RK),allocatable,dimension(:,:) :: ws
      real(RK),allocatable,dimension(:)   :: total
      real(RK),allocatable,dimension(:)   :: local
      real(RK),allocatable,dimension(:) :: dz

      !# Arrays for environmental variables not supplied externally.
      real(RK),allocatable,dimension(:) :: par, pres
      real(RK),allocatable,dimension(:) :: uva, uvb, nir

      character(len=48),allocatable :: names(:)
      character(len=48),allocatable :: bennames(:)

      integer,allocatable,dimension(:) :: externalid

      real(RK),allocatable,dimension(:) :: min_, max_

      integer :: n_aed2_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet

   contains
      procedure, pass(self), public :: init
      procedure, pass(self), public :: update
   end type SimstratAED2

contains

   subroutine init(self, aed2_cfg, grid)
      implicit none
      class(SimstratAED2) :: self
      class(AED2Config), target :: aed2_cfg
      class(StaggeredGrid), target :: grid

      ! Local variables
      character(len=80) :: fname
      type(aed2_variable_t),pointer :: tvar

      character(len=64) :: models(64)
      namelist /aed2_models/ models
      integer i, j, status, rc

      self%grid => grid
      self%aed2_cfg => aed2_cfg

      associate (n_aed2_vars => self%n_aed2_vars, &
                 n_vars => self%n_vars, &
                 n_vars_ben => self%n_vars_ben, &
                 n_vars_diag => self%n_vars_diag, &
                 n_vars_diag_sheet => self%n_vars_diag_sheet)

         fname = 'aed2.nml'

         if ( aed2_init_core('.') /= 0 ) call error("Initialisation of aed2_core failed")
         call aed2_print_version

         ! Create model tree
         write (6,*) "     Processing aed2_models config from ", trim(fname)
         open(50,file=fname,action='read',status='old',iostat=status)
         if ( status /= 0 ) then
            call error("Cannot open file " // trim(fname))
            stop
         end if

         models = ''
         read(50, nml=aed2_models, iostat=status)
         if ( status /= 0 ) then
            call error("Cannot read namelist entry aed2_models")
            stop
         end if

         do i=1,size(models)
            if (models(i)=='') exit
            call aed2_define_model(models(i), 50)
         end do

         !# should be finished with this file
         close(50)
         write (6,*) "      AED2 file parsing completed."

         n_aed2_vars = aed2_core_status(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)

         print "(/,5X,'AED2 : n_aed2_vars = ',I3,' ; MaxLayers         = ',I4)",n_aed2_vars,grid%nz_grid
         print "(  5X,'AED2 : n_vars      = ',I3,' ; n_vars_ben        = ',I3)",n_vars,n_vars_ben
         print "(  5X,'AED2 : n_vars_diag = ',I3,' ; n_vars_diag_sheet = ',I3,/)",n_vars_diag,n_vars_diag_sheet

         call check_data(self)

         !# names = grab the names from info
         allocate(self%names(n_vars),stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (names)'
         allocate(self%bennames(n_vars_ben),stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (bennames)'

         !# Now that we know how many vars we need, we can allocate space for them
         allocate(self%cc(grid%nz_grid, (n_vars + n_vars_ben)),stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (CC)'
         self%cc = 0.         !# initialise to zeroFarc

         allocate(self%min_((n_vars + n_vars_ben)))
         allocate(self%max_((n_vars + n_vars_ben)))
         print "(5X,'Configured variables to simulate:')"

         j = 0
         do i=1,self%n_aed2_vars
            if ( aed2_get_var(i, tvar) ) then
               if ( .not. (tvar%sheet .or. tvar%diag .or. tvar%extern) ) then
                  j = j + 1
                  self%names(j) = trim(tvar%name)
                  self%min_(j) = tvar%minimum
                  self%max_(j) = tvar%maximum
                  print *,"     S(",j,") AED2 pelagic(3D) variable: ", trim(self%names(j))
            end if
         end if
      end do

      j = 0
      do i=1,n_aed2_vars
         if ( aed2_get_var(i, tvar) ) then
               if ( tvar%sheet .and. .not. (tvar%diag .or. tvar%extern) ) then
                  j = j + 1
                  self%bennames(j) = trim(tvar%name)
                  self%min_(n_vars+j) = tvar%minimum
                  self%max_(n_vars+j) = tvar%maximum
                  print *,"     B(",j,") AED2 benthic(2D) variable: ", trim(self%bennames(j))
               end if
            end if
         end do

         j = 0
         do i=1,n_aed2_vars
            if ( aed2_get_var(i, tvar) ) then
               if ( tvar%diag ) then
                  if ( .not.  tvar%sheet ) then
                     j = j + 1
                     print *,"     D(",j,") AED2 diagnostic 3Dvariable: ", trim(tvar%name)
                  end if
               end if
            end if
         end do

         j = 0
         do i=1,n_aed2_vars
            if ( aed2_get_var(i, tvar) ) then
               if ( tvar%diag ) then
                  if (tvar%sheet ) then
                     j = j + 1
                     print *,"     D(",j,") AED2 diagnostic 2Dvariable: ", trim(tvar%name)
                  end if
               end if
            end if
         enddo

         allocate(self%externalid(n_aed2_vars))

         !# Allocate diagnostic variable array and set all values to zero.
         !# (needed because time-integrated/averaged variables will increment rather than set the array)
         allocate(self%cc_diag(grid%nz_grid, n_vars_diag),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (cc_diag)'
         self%cc_diag = zero_

         !# Allocate diagnostic variable array and set all values to zero.
         !# (needed because time-integrated/averaged variables will increment rather than set the array)
         allocate(self%cc_diag_hz(n_vars_diag_sheet),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (cc_diag_hz)'
         self%cc_diag_hz = zero_

         !# Allocate array with vertical movement rates (m/s, positive for upwards),
         !# and set these to the values provided by the model.
         allocate(self%ws(grid%nz_grid, n_vars),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (ws)'
         self%ws = zero_

         !# Allocate array for photosynthetically active radiation (PAR).
         !# This will be calculated internally during each time step.
         allocate(self%par(grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (par)'
         self%par = zero_

         allocate(self%nir(grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (nir)'
         self%nir = zero_
         allocate(self%uva(grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (uva)'
         self%uva = zero_
         allocate(self%uvb(grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (uvb)'
         self%uvb = zero_

         allocate(self%dz(grid%nz_grid),stat=rc)
         self%dz = zero_

         !# Allocate array for local pressure.
         !# This will be calculated [approximated] from layer depths internally
         !# during each time step.
         allocate(self%pres(grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (pres)'
         self%pres = zero_

         allocate(self%tss(grid%nz_grid),stat=rc)
         if (rc /= 0) stop 'allocate_memory(): Error allocating (tss)'
         self%tss = zero_

         write(*,"(/,5X,'----------  AED2 config : end  ----------',/)")
      end associate
   end subroutine


   subroutine update(self, aed2_cfg)
      implicit none
      class(SimstratAED2) :: self
      class(AED2Config), target :: aed2_cfg


   end subroutine

   subroutine check_data(self)
   !-------------------------------------------------------------------------------
   ! Check that all variable dependencies have been met
   !-------------------------------------------------------------------------------
   !ARGUMENTS
   class(SimstratAED2) :: self
   !LOCALS
      integer :: av !, i
      integer :: v, d, sv, sd, ev, err_count
      type(aed2_variable_t),pointer :: tvar
   !-------------------------------------------------------------------------------
   !BEGIN
      v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
      err_count = 0

      do av=1,self%n_aed2_vars
      if ( .not.  aed2_get_var(av, tvar) ) then
         call error("Error getting variable info")
         stop
      end if

      if ( tvar%extern ) then !# global variable
         ev = ev + 1
         select case (tvar%name)
            case ( 'temperature' ) ; tvar%found = .true.
            case ( 'salinity' )    ; tvar%found = .true.
            case ( 'density' )     ; tvar%found = .true.
            case ( 'layer_ht' )    ; tvar%found = .true.
            case ( 'extc_coef' )   ; tvar%found = .true.
            case ( 'tss' )         ; tvar%found = .true.
            case ( 'par' )         ; tvar%found = .true.
            case ( 'nir' )         ; tvar%found = .true.
            case ( 'uva' )         ; tvar%found = .true.
            case ( 'uvb' )         ; tvar%found = .true.
            case ( 'pressure' )    ; tvar%found = .true.
            case ( 'depth' )       ; tvar%found = .true.
            case ( 'sed_zone' )    ; tvar%found = .true.
            case ( 'wind_speed' )  ; tvar%found = .true.
            case ( 'par_sf' )      ; tvar%found = .true.
            case ( 'taub' )        ; tvar%found = .true.
            case ( 'lake_depth' )  ; tvar%found = .true.
            case ( 'layer_area' )  ; tvar%found = .true.
            case default ; call error("ERROR: external variable "//trim(tvar%name)//" not found.")
         end select
      elseif ( tvar%diag ) then  !# Diagnostic variable
         if ( tvar%sheet ) then
            sd = sd + 1
         else
            d = d + 1
         end if
      else    !# state variable
         if ( tvar%sheet ) then
            sv = sv + 1
         else
            v = v + 1
         end if
      end if
      if ( .not. tvar%found ) then
         call error("Undefined variable " //trim(tvar%name))
         err_count = err_count + 1
      end if
   enddo

   if ( self%n_vars < v ) print *,"More vars than expected",v,self%n_vars
   if ( self%n_vars_ben < sv ) print *,"More sheet vars than expected"
   if ( self%n_vars_diag < d ) print *,"More diag vars than expected"
   if ( self%n_vars_diag_sheet < sd ) print *,"More sheet diag vars than expected"

   if ( err_count > 0 ) then
      call error("In AED2 configuration")
      stop
   end if
END SUBROUTINE check_data

end module simstrat_aed2