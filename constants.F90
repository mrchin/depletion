module constants

  implicit none

  ! ============================================================================
  ! VERSIONING NUMBERS

  ! OpenMC major, minor, and release numbers
  integer, parameter :: VERSION_MAJOR   = 0
  integer, parameter :: VERSION_MINOR   = 5
  integer, parameter :: VERSION_RELEASE = 1

  ! Revision numbers for binary files
  integer, parameter :: REVISION_SOURCE     = 1
  integer, parameter :: REVISION_STATEPOINT = 8

  ! ============================================================================
  ! ADJUSTABLE PARAMETERS 

  ! NOTE: This is the only section of the constants module that should ever be
  ! adjusted. Modifying constants in other sections may cause the code to fail.

  ! Monoatomic ideal-gas scattering treatment threshold
  real(8), parameter :: FREE_GAS_THRESHOLD = 400.0

  ! Significance level for confidence intervals
  real(8), parameter :: CONFIDENCE_LEVEL = 0.95_8

  ! Used for surface current tallies
  real(8), parameter :: TINY_BIT = 1e-8_8

  ! User for precision in geometry
  real(8), parameter :: FP_PRECISION = 1e-14_8
  real(8), parameter :: FP_REL_PRECISION = 1e-5_8
  real(8), parameter :: FP_COINCIDENT = 1e-12_8

  ! Maximum number of collisions/crossings
  integer, parameter :: MAX_EVENTS = 10000
  integer, parameter :: MAX_SAMPLE = 100000

  ! Maximum number of words in a single line, length of line, and length of
  ! single word
  integer, parameter :: MAX_WORDS    = 500
  integer, parameter :: MAX_LINE_LEN = 250
  integer, parameter :: MAX_WORD_LEN = 150
  integer, parameter :: MAX_FILE_LEN = 255

  ! ============================================================================
  ! PHYSICAL CONSTANTS

  real(8), parameter ::            &
       PI           = 3.1415926535898_8, & ! pi
       MASS_NEUTRON = 1.0086649156,      & ! mass of a neutron
       MASS_PROTON  = 1.00727646677,     & ! mass of a proton
       AMU          = 1.66053873e-27,    & ! 1 amu in kg
       N_AVOGADRO   = 0.602214179,       & ! Avogadro's number in 10^24/mol
       K_BOLTZMANN  = 8.617342e-11,      & ! Boltzmann constant in MeV/K
       INFINITY     = huge(0.0_8),       & ! positive infinity
       ZERO         = 0.0_8,             &
       ONE          = 1.0_8,             &
       TWO          = 2.0_8

  ! ============================================================================
  ! GEOMETRY-RELATED CONSTANTS

  ! Boundary conditions
  integer, parameter ::  &
       BC_TRANSMIT = 0,  & ! Transmission boundary condition (default)
       BC_VACUUM   = 1,  & ! Vacuum boundary condition
       BC_REFLECT  = 2,  & ! Reflecting boundary condition
       BC_PERIODIC = 3     ! Periodic boundary condition

  ! Logical operators for cell definitions
  integer, parameter ::              &
       OP_LEFT_PAREN  = huge(0),     & ! Left parentheses
       OP_RIGHT_PAREN = huge(0) - 1, & ! Right parentheses
       OP_UNION       = huge(0) - 2, & ! Union operator
       OP_DIFFERENCE  = huge(0) - 3    ! Difference operator

  ! Cell types
  integer, parameter ::  &
       CELL_NORMAL  = 1, & ! Cell with a specified material
       CELL_FILL    = 2, & ! Cell filled by a separate universe
       CELL_LATTICE = 3    ! Cell filled with a lattice

  ! Void material
  integer, parameter :: MATERIAL_VOID = -1

  ! Lattice types
  integer, parameter ::  &
       LATTICE_RECT = 1, & ! Rectangular lattice
       LATTICE_HEX  = 2    ! Hexagonal lattice

  ! Lattice boundary crossings
  integer, parameter ::    &
       LATTICE_LEFT   = 1, & ! Flag for crossing left (x) lattice boundary
       LATTICE_RIGHT  = 2, & ! Flag for crossing right (x) lattice boundary
       LATTICE_BACK   = 3, & ! Flag for crossing back (y) lattice boundary
       LATTICE_FRONT  = 4, & ! Flag for crossing front (y) lattice boundary
       LATTICE_BOTTOM = 5, & ! Flag for crossing bottom (z) lattice boundary
       LATTICE_TOP    = 6    ! Flag for crossing top (z) lattice boundary

  ! Surface types
  integer, parameter ::  &
       SURF_PX     =  1, & ! Plane parallel to x-plane 
       SURF_PY     =  2, & ! Plane parallel to y-plane 
       SURF_PZ     =  3, & ! Plane parallel to z-plane 
       SURF_PLANE  =  4, & ! Arbitrary plane
       SURF_CYL_X  =  5, & ! Cylinder along x-axis
       SURF_CYL_Y  =  6, & ! Cylinder along y-axis
       SURF_CYL_Z  =  7, & ! Cylinder along z-axis
       SURF_SPHERE =  8, & ! Sphere
       SURF_CONE_X =  9, & ! Cone parallel to x-axis
       SURF_CONE_Y = 10, & ! Cone parallel to y-axis
       SURF_CONE_Z = 11    ! Cone parallel to z-axis

  ! ============================================================================
  ! CROSS SECTION RELATED CONSTANTS

  ! Interpolation flag
  integer, parameter ::   &
       HISTOGRAM     = 1, & ! y is constant in x
       LINEAR_LINEAR = 2, & ! y is linear in x
       LINEAR_LOG    = 3, & ! y is linear in ln(x)
       LOG_LINEAR    = 4, & ! ln(y) is linear in x
       LOG_LOG       = 5    ! ln(y) is linear in ln(x)

  ! Particle type
  integer, parameter :: &
       NEUTRON  = 1, &
       PHOTON   = 2, &
       ELECTRON = 3

  ! Angular distribution type
  integer, parameter :: & 
       ANGLE_ISOTROPIC = 1, & ! Isotropic angular distribution
       ANGLE_32_EQUI   = 2, & ! 32 equiprobable bins
       ANGLE_TABULAR   = 3    ! Tabular angular distribution

  ! Secondary energy mode for S(a,b) inelastic scattering
  integer, parameter :: &
       SAB_SECONDARY_EQUAL  = 0, & ! Equally-likely outgoing energy bins
       SAB_SECONDARY_SKEWED = 1    ! Skewed outgoing energy bins

  ! Elastic mode for S(a,b) elastic scattering
  integer, parameter :: &
       SAB_ELASTIC_DISCRETE = 3, & ! Sample from discrete cosines
       SAB_ELASTIC_EXACT    = 4    ! Exact treatment for coherent elastic

  ! Reaction types
  integer, parameter :: &
       TOTAL_XS    = 1, &
       ELASTIC     = 2, &
       N_LEVEL     = 4, &
       MISC        = 5, &
       N_2ND       = 11, &
       N_2N        = 16, &
       N_3N        = 17, &
       N_FISSION   = 18, &
       N_F         = 19, &
       N_NF        = 20, &
       N_2NF       = 21, &
       N_NA        = 22, &
       N_N3A       = 23, &
       N_2NA       = 24, &
       N_3NA       = 25, &
       N_NP        = 28, &
       N_N2A       = 29, &
       N_2N2A      = 30, &
       N_ND        = 32, &
       N_NT        = 33, &
       N_N3HE      = 34, &
       N_ND2A      = 35, &
       N_NT2A      = 36, &
       N_4N        = 37, &
       N_3NF       = 38, &
       N_2NP       = 41, &
       N_3NP       = 42, &
       N_N2P       = 44, &
       N_NPA       = 45, &
       N_N1        = 51, &
       N_N40       = 90, &
       N_NC        = 91, &
       N_DISAPPEAR = 101, &
       N_GAMMA     = 102, &
       N_P         = 103, &
       N_D         = 104, &
       N_T         = 105, &
       N_3HE       = 106, &
       N_A         = 107, &
       N_2A        = 108, &
       N_3A        = 109, &
       N_2P        = 111, &
       N_PA        = 112, &
       N_T2A       = 113, &
       N_D2A       = 114, &
       N_PD        = 115, &
       N_PT        = 116, &
       N_DA        = 117

  ! ACE table types
  integer, parameter :: &
       ACE_NEUTRON   = 1, & ! continuous-energy neutron
       ACE_THERMAL   = 2, & ! thermal S(a,b) scattering data
       ACE_DOSIMETRY = 3    ! dosimetry cross sections

  ! Fission neutron emission (nu) type
  integer, parameter ::   &
       NU_NONE       = 0, & ! No nu values (non-fissionable)
       NU_POLYNOMIAL = 1, & ! Nu values given by polynomial
       NU_TABULAR    = 2    ! Nu values given by tabular distribution

  ! Cross section filetypes
  integer, parameter :: &
       ASCII  = 1, & ! ASCII cross section file
       BINARY = 2    ! Binary cross section file

  ! Probability table parameters
  integer, parameter :: &
       URR_CUM_PROB = 1, &
       URR_TOTAL    = 2, &
       URR_ELASTIC  = 3, &
       URR_FISSION  = 4, &
       URR_N_GAMMA  = 5, &
       URR_HEATING  = 6

  ! Maximum number of partial fission reactions
  integer, parameter :: PARTIAL_FISSION_MAX = 4

  ! ENDF decay modes
  integer, parameter :: &
       DECAY_GAMMA    = 0, & ! Gamma decay
       DECAY_BETA     = 1, & ! Beta decay
       DECAY_EC       = 2, & ! Electron capture and/or positron emission
       DECAY_IT       = 3, & ! Isomeric transition
       DECAY_ALPHA    = 4, & ! Alpha decay
       DECAY_NEUTRON  = 5, & ! Neutron emission
       DECAY_SF       = 6, & ! Spontaneous fission
       DECAY_PROTON   = 7, & ! Proton emission
       DECAY_ELECTRON = 8, & ! "Discrete electrons"
       DECAY_XRAY     = 9, & ! X-rays and annihilation radiation
       DECAY_UNKNOWN  = 10

  ! ============================================================================
  ! TALLY-RELATED CONSTANTS

  ! Tally type
  integer, parameter :: &
       TALLY_VOLUME          = 1, &
       TALLY_SURFACE_CURRENT = 2

  ! Tally estimator types
  integer, parameter :: &
       ESTIMATOR_ANALOG      = 1, &
       ESTIMATOR_TRACKLENGTH = 2

  ! Event types for tallies
  integer, parameter :: &
       EVENT_SURFACE = -2, &
       EVENT_LATTICE = -1, &
       EVENT_SCATTER =  1, &
       EVENT_ABSORB  =  2, &
       EVENT_FISSION =  3 

  ! Tally score type
  integer, parameter :: N_SCORE_TYPES = 15
  integer, parameter :: &
       SCORE_FLUX          = -1,  & ! flux
       SCORE_TOTAL         = -2,  & ! total reaction rate
       SCORE_SCATTER       = -3,  & ! scattering rate
       SCORE_NU_SCATTER    = -4,  & ! scattering production rate
       SCORE_SCATTER_N     = -5,  & ! arbitrary scattering moment
       SCORE_SCATTER_PN    = -6,  & ! system for scoring 0th through nth moment
       SCORE_TRANSPORT     = -7,  & ! transport reaction rate
       SCORE_DIFFUSION     = -8,  & ! diffusion coefficient
       SCORE_N_1N          = -9,  & ! (n,1n) rate
       SCORE_ABSORPTION    = -10, & ! absorption rate
       SCORE_FISSION       = -11, & ! fission rate
       SCORE_NU_FISSION    = -12, & ! neutron production rate
       SCORE_KAPPA_FISSION = -13, & ! fission energy production rate
       SCORE_CURRENT       = -14, & ! partial current
       SCORE_EVENTS        = -15    ! number of events
       
  ! Maximum scattering order supported
  integer, parameter :: SCATT_ORDER_MAX = 10
  character(len=*), parameter :: SCATT_ORDER_MAX_PNSTR = "scatter-p10"

  ! Tally map bin finding
  integer, parameter :: NO_BIN_FOUND = -1

  ! Tally filter and map types
  integer, parameter :: N_FILTER_TYPES = 8
  integer, parameter :: &
       FILTER_UNIVERSE  = 1, &
       FILTER_MATERIAL  = 2, &
       FILTER_CELL      = 3, &
       FILTER_CELLBORN  = 4, &
       FILTER_SURFACE   = 5, &
       FILTER_MESH      = 6, &
       FILTER_ENERGYIN  = 7, &
       FILTER_ENERGYOUT = 8

  ! Tally surface current directions
  integer, parameter :: &
       IN_RIGHT  = 1,   &
       OUT_RIGHT = 2,   &
       IN_FRONT  = 3,   &
       OUT_FRONT = 4,   &
       IN_TOP    = 5,   &
       OUT_TOP   = 6

  ! Global tallY parameters
  integer, parameter :: N_GLOBAL_TALLIES = 4
  integer, parameter :: &
       K_COLLISION   = 1, &
       K_ABSORPTION  = 2, &
       K_TRACKLENGTH = 3, &
       LEAKAGE       = 4

  ! ============================================================================
  ! EXTERNAL SOURCE PARAMETERS

  ! Source spatial distribution types
  integer, parameter :: &
       SRC_SPACE_BOX   = 1, & ! Source in a rectangular prism
       SRC_SPACE_POINT = 2    ! Source at a single point

  ! Source angular distribution types
  integer, parameter :: &
       SRC_ANGLE_ISOTROPIC = 1, & ! Isotropic angular 
       SRC_ANGLE_MONO      = 2, & ! Monodirectional source
       SRC_ANGLE_TABULAR   = 3    ! Tabular distribution

  ! Source energy distribution types
  integer, parameter :: &
       SRC_ENERGY_MONO    = 1, & ! Monoenergetic source
       SRC_ENERGY_MAXWELL = 2, & ! Maxwell fission spectrum
       SRC_ENERGY_WATT    = 3, & ! Watt fission spectrum
       SRC_ENERGY_TABULAR = 4    ! Tabular distribution
       
  ! ============================================================================
  ! MISCELLANEOUS CONSTANTS

  ! indicates that an array index hasn't been set
  integer, parameter :: NONE = 0

  ! Codes for read errors -- better hope these numbers are never used in an
  ! input file!
  integer, parameter :: ERROR_INT  = -huge(0)
  real(8), parameter :: ERROR_REAL = -huge(0.0_8) * 0.917826354_8

  ! Energy grid methods
  integer, parameter :: &
       GRID_NUCLIDE  = 1, & ! non-unionized energy grid
       GRID_UNION    = 2, & ! union grid with pointers
       GRID_LETHARGY = 3    ! lethargy mapping

  ! Running modes
  integer, parameter ::        &
       MODE_FIXEDSOURCE = 1, & ! Fixed source mode
       MODE_EIGENVALUE  = 2, & ! K eigenvalue mode
       MODE_PLOTTING    = 3, & ! Plotting mode
       MODE_TALLIES     = 4, & ! Tally results mode
       MODE_PARTICLE    = 5, & ! Particle restart mode
       MODE_DEPLETION   = 6    ! Depletion mode

  ! Unit numbers
  integer, parameter :: UNIT_SUMMARY  = 11 ! unit # for writing summary file
  integer, parameter :: UNIT_TALLY    = 12 ! unit # for writing tally file
  integer, parameter :: UNIT_PLOT     = 13 ! unit # for writing plot file
  integer, parameter :: UNIT_XS       = 14 ! unit # for writing xs summary file
  integer, parameter :: UNIT_SOURCE   = 15 ! unit # for writing source file
  integer, parameter :: UNIT_STATE    = 16 ! unit # for writing state point
  integer, parameter :: CMFD_BALANCE  = 17 ! unit # for writing cmfd balance file
  integer, parameter :: UNIT_PARTICLE = 18 ! unit # for writing particle restart
  integer, parameter :: UNIT_ENDF     = 19 ! unit # for reading ENDF files

  !=============================================================================
  ! CMFD CONSTANTS

  ! for non-accelerated regions on coarse mesh overlay
  integer, parameter :: CMFD_NOACCEL = 99999

  ! constant to represent a zero flux "albedo"
  real(8), parameter :: ZERO_FLUX = 999.0_8

  ! constant to represent albedo rejection
  real(8), parameter :: ALBEDO_REJECT = 999.0_8

  ! constant for writing out no residual
  real(8), parameter :: CMFD_NORES = 99999.0_8

  !=============================================================================
  ! DEPLETION CONSTANTS

  ! Chebyshev coefficients from M. Pusa, "Corrections to partial fraction
  ! decomposition coefficients for Chebyshev rational approximation on the
  ! negative real axis," arXiv:1206.2880v1 [math.NA] 13 Jun 2012.

  integer, parameter :: CRAM_ORDER = 14
  complex(8), parameter :: CRAM_ALPHA0 = (1.8321743782540412751e-14_8, 0.0_8)
  complex(8), parameter :: CRAM_ALPHA(7) = (/ &
       (-7.1542880635890672853e-5_8, +1.4361043349541300111e-4_8), &
       (+9.4390253107361688779e-3_8, -1.7184791958483017511e-2_8), &
       (-3.7636003878226968717e-1_8, +3.3518347029450104214e-1_8), &
       (-2.3498232091082701191e+1_8, -5.8083591297142074004e+0_8), &
       (+4.6933274488831293047e+1_8, +4.5643649768827760791e+1_8), &
       (-2.7875161940145646468e+1_8, -1.0214733999056451434e+2_8), &
       (+4.8071120988325088907e+0_8, -1.3209793837428723881e+0_8) /)
  complex(8), parameter :: CRAM_THETA(7) = (/ &
       (-8.8977731864688888199_8, +1.6630982619902085304e+1_8), &
       (-3.7032750494234480603_8, +1.3656371871483268171e+1_8), &
       (-0.2087586382501301251_8, +1.0991260561901260913e+1_8), &
       (+3.9933697105785685194_8, +6.0048316422350373178_8), &
       (+5.0893450605806245066_8, +3.5888240290270065102_8), &
       (+5.6231425727459771248_8, +1.1940690463439669766_8), &
       (+2.2697838292311127097_8, +8.4617379730402214019_8) /)

  ! Indicate a blank column in sparse matrix
  integer, parameter :: NULL_COLUMN = -1

  ! Depletion step units
  integer, parameter :: &
       DEPLETION_STEP_DAY   = 1, & ! Step in days
       DEPLETION_STEP_MWDKG = 2    ! Step in MWd/kgU

  ! Indicate that depletion step is decay only
  real(8), parameter :: DEPLETION_DECAY = -1.0_8

end module constants
