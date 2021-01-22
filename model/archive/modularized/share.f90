module share

    IMPLICIT NONE

    ! Basic model parameters
    REAL(8), parameter :: delta=.022 !depreciation rate of durable  was .03
    REAL(8), parameter :: elasticity= 2  !elasticity of substitution
    REAL(8), parameter :: F = 0.05 ! fixed cost
    REAL(8), parameter :: r=.024 ! rental rate?
    REAL(8), parameter :: rborrow=.024
    REAL(8), parameter :: psi=(r/(1+r))**(-elasticity)
    REAL(8), parameter :: rent=1-(1-delta)/(1+r)
    REAL(8), parameter :: thetamatlab=.15
    REAL(8), parameter :: theta=1-(1-thetamatlab)*(1-delta)/(1+r) !down payment

    ! Parameters formerly defined in durables. Assigning initial values here
    REAL(8), parameter :: r_rentalO = 0.05 + 0.0105
    REAL(8), parameter :: r_rental_retireO= r_rentalO
    REAL(8), parameter :: elasticity2O= .83 + .03 ! 1 - (share of expenditure in housing)
    REAL(8), parameter :: beta2O = 0.92 + 0.002   !Quarterly discount factor
    REAL(8), parameter :: beta2retireO = beta2O
    REAL(8), parameter :: psi0 = 1.0
    REAL(8), parameter :: ConsElas = 2.5
    REAL(8) :: r_rental, r_rental_retire, elasticity2, beta2, beta2retire

    ! The asset and durable grids
    REAL(8), parameter :: Dmin=0
    REAL(8), parameter :: Dmax=12

    !REAL(8), parameter :: amin=(exp(hpmin)-exp(hpmax))*(1-theta)*Dmax    !minimum asset value (actually voluntary equity if theta!=1).  Note that if theta<1, we can have assets <0 even with amin=0 since a is vol. equity.  But with theta=1 & amin=0, no borrowing.
    !REAL(8), parameter :: amin=-.89
    REAL(8), parameter :: amin=-0.40/(1+r)  !.40 is roughly min y, so this is basically a no default condition
    !REAL(8), parameter :: amax2=8   !max asset
    REAL(8), parameter :: amax=40  !max asset

    integer, parameter :: agridsize = 60
    integer, parameter :: Dgridsize = 75

    !idiosyncratic earnings shocks
    REAL(8), parameter :: sigma_z= .21
    REAL(8), parameter :: rho_z=.91
    REAL(8), parameter :: sigma_eta_init=(sigma_z**2.0/(1-rho_z**2.0))**.5
    REAL(8), parameter :: zmin=-2.5*sigma_eta_init
    REAL(8), parameter :: zmax=2.5*sigma_eta_init

    integer, parameter :: zgridsize=13  !size of ar(1) income shock grid
    REAL(8), parameter :: borrowconstraint=0.0

    !House price subsidy shocks and policy settings
    REAL(8), parameter :: hpmin=-0.05
    REAL(8), parameter :: hpmax=0.0
    integer, parameter :: hpgridsize=3 !size of house price grid
    integer, parameter :: PolYrs=1 ! Number of years subsidy is active
    integer, parameter :: PolYrDiff=0 ! Years btwn Repeat and Homebuyer policies
    integer, parameter :: EligYrsF=3 ! Years under which agent is not owning to be first-time
    integer, parameter :: EligYrsR=5 ! Years under which agent owns the same durable to be repeat
    REAL(8), DIMENSION(hpgridsize) :: hpnodes, thetanodes, hpdelta
    REAL(8), DIMENSION(hpgridsize, hpgridsize) :: Probhp
    REAL(8), DIMENSION(hpgridsize-1) :: eta ! 0-1: Smoothing how much subsidy applies to down

    integer, parameter :: Tretire=38 !40 !35
    integer, parameter :: Tdie=Tretire+25


    ! Below are the state and policy holders


    REAL(8), DIMENSION(agridsize) :: anodes  !Nodes for idiosyncratic assets (voluntary equity)
    REAL(8), DIMENSION(Dgridsize) :: Dnodes  !Nodes for idiosyncratic durable

    REAL(8), DIMENSION(zgridsize) :: Probinit, Probinitcum



    REAL(8), DIMENSION(zgridsize) :: znodes  !Nodes for idiosyncratic productivity
    REAL(8), DIMENSION(zgridsize, zgridsize) :: Probz, Probzcum  !while probability of being at any z will vary with time, transition probability will not, so max we need to loop over wont
    INTEGER, DIMENSION(zgridsize) :: maxexpectationz, minexpectationz


    REAL(8), DIMENSION(Tretire) :: ageearnings
    REAL(8), DIMENSION(Tdie-Tretire) :: deathprob
    REAL(8), DIMENSION(zgridsize) :: predictedlifetimeearnings, retirementincome
    REAL(8), DIMENSION(zgridsize, Tdie) :: income
    REAL(8) :: averagelifetimeearnings

    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Vnoadjust,  Vadjust, Vrent  !Conjectured value of not adjusting and adjusting, and new values of not adjusting and adjusting
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: achoiceadjust,  achoicenoadjust, achoicerent,  Dchoiceadjust, Dchoicenoadjust, Dchoicerent, achoice, Dchoice, cchoice, cchoiceadjust, cchoicenoadjust, cchoicerent, cchoiceshock, Dchoiceshock  !Policy functions when solving problem
    REAL(8) :: priceshock

    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: rentalindicator, choiceindicator !choice indicator = 1 if adjust, 2 if noadjust, 3 if rent
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: EV


    INTEGER :: i, j, k, l, m, n, t, iter
    INTEGER, parameter :: numhouseholds=8000   !Size of household panel
    REAL(8), dimension(numhouseholds, Tdie, 13) :: householdresultsholder
    !integer, parameter :: burnin=1000    !Initial periods dropped in simulation


    REAL(8) :: diffmoments, diffmoments1, diffmoments2, diffmoments3
    REAL(8), parameter :: difftol=.01
    REAL(8) :: numiter

    REAL(8), parameter :: costlyequity=1

    REAL(8), parameter :: rentelasticity=1

    ! GE transition variables
    REAL(8), DIMENSION(Tdie+1) :: Ps,gs
    REAL(8) :: g

    integer :: totalTranslength=40

    REAL(8) :: P1
    REAL(8) :: P2

    INTEGER :: TransTime ! To solveDP?

    REAL(8) :: age

    integer :: ageatbubblestart ! To transition?



    CONTAINS
    
    subroutine share_arrays
        IMPLICIT NONE

        eta(1) = 1
        eta(2) = 1

        ! House price grid - three for different policy states, fill manually
        write(*,*) "House price grid:"
        hpnodes(1)=hpmin
        hpnodes(hpgridsize-1)=hpmin
        hpnodes(hpgridsize)=hpmax

        !Housing shocks
        Probhp=0
        DO i=1, hpgridsize-1
            Probhp(i, i)=0.0
            Probhp(i, hpgridsize)=1.0
        END DO
        Probhp(hpgridsize, hpgridsize)=1.0

    end subroutine share_arrays

end module share
