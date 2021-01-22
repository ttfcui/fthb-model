module share

    IMPLICIT NONE

    ! Basic model parameters
    REAL(8), parameter :: elasticity2= .759 ! 1 - (share of expenditure in housing)
    REAL(8), parameter :: beta2= 0.94! 0.9175 !Quarterly discount factor
    REAL(8), parameter :: e_of_s= 2 ! elasticity of substitution (unused)
    REAL(8), parameter :: r=.010 ! return on safe asset
    REAL(8), parameter :: rborrow=r+.008 ! Plus 0.8% mortgage premium
    REAL(8), parameter :: delta=.022 !depreciation rate of durable  was .03
    REAL(8), parameter :: dtau=.000 !Property tax on housing
    REAL(8), parameter :: tax1=0.175 !Income tax on labour
    REAL(8), parameter :: tax2=.15 ! Progressivity parameter on income tax
    REAL(8), parameter :: F= 0.06 ! fixed cost (e.g. realtor fees)
    REAL(8), parameter :: F2= 0.00 ! quadratic trans cost term (like furniture?)
    REAL(8) :: Dmin= 1.075 ! Minimum owned house size (not parameter to fit to D grid)
    REAL(8), parameter :: thetamatlab=.20 ! Down payment proportion
    REAL(8), parameter :: ConsElas= 2.5 ! Price elasticity of agg. construction

    ! Current parameters to be calibrated
    REAL(8), parameter :: elasticity= 2.5 ! intertemporal elasticity of substitution
    REAL(8), parameter :: rentPrem= 7.7e-2 ! Rental operation costs on top of user cost
    REAL(8), parameter :: rentPremRetire= rentPrem ! Different premium for retirement
    REAL(8), parameter :: ret_wealth=0.03 ! Lump-sum transfer at retirement, relative to labour income
    REAL(8), parameter :: psi= 5 ! End-of-life bequest parameter
    REAL(8), parameter :: rentUtil= .93 ! Disutility of renting relative to owning
    REAL(8), parameter :: beq_base= 0.00 ! affects MU of a unit increase in bequests

    ! Derived parameters.
    REAL(8), parameter :: beta2retire=beta2  ! Different discount rate in retirement (deprecated)
    REAL(8), parameter :: usercost = (1.0-delta-dtau)/(1.0+r)  ! Traditional user cost of housing
    REAL(8), parameter :: rent=1.0-usercost  ! rent equals user cost in frictionless market
    ! Downpayment incorporating user cost - from 1-theta = (1-thetamatlab)*usercost
    REAL(8), parameter :: theta=1-(1-thetamatlab)*usercost
    REAL(8), parameter :: psi2 = ConsElas/(1+ConsElas) ! Parameter for housing supply fn
    INTEGER, parameter :: transgridsize= 1

    ! The asset and durable grids
    REAL(8), parameter :: Dmax=7.5  ! Maximum housing level (measured in units of $67,250)
    REAL(8) :: Dmax_rent=Dmax ! Maximum rental housing available
    !REAL(8), parameter :: amin=(exp(hpmin)-exp(hpmax))*(1-theta)*Dmax    !minimum asset value (actually voluntary equity if theta!=1).  Note that if theta<1, we can have assets <0 even with amin=0 since a is vol. equity.  But with theta=1 & amin=0, no borrowing.
    REAL(8), parameter :: amin=-0.40/(1+r)  !.40 is roughly min y, so this is basically a no default condition
    REAL(8), parameter :: amax=9 !max asset (measured in units of $67,250)
    REAL(8), parameter :: borrowconstraint=0.0

    integer, parameter :: agridsize= 90.000000 
    integer, parameter :: Dgridsize= 55.000000 

    ! Permanent income process
    REAL(8), parameter :: sigma_z= 0.21  ! Standard deviation of productivity shocks
    REAL(8), parameter :: rho_z=.91  ! Autocorrelation of AR(1) process
    integer, parameter :: zgridsize= 19 !size of ar(1) income shock grid
    REAL(8), parameter :: poismean=0.1  ! Mean of initial poisson dist (deprecated)

    ! Job separation / transitory income process
    REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: empprob
    REAL(8), parameter :: unempprob= 0.118  ! Prob. of going into unemp. (deprecated)
    REAL(8), parameter :: movprob= 0.015  ! Prob. of exogenous move for owners
    REAL(8), parameter :: movprobR= 0.000 ! Prob. of exogenous home purchase for renters (deprecated)
    REAL(8), parameter :: unemptrans= -0.05  ! Unemploymeny benefits (deprecated)

    !House price subsidy shocks and policy settings
    integer, parameter :: hptransLength=1 ! Total number of years in transition
    integer, parameter :: PolStart=0 ! Period when policy starts (but policy is announced at period 0)
    integer, parameter :: PolEnd=1 ! Period when policy stops
    integer, parameter :: EligYrsF= 40.000000 ! Years under which agent is not owning to be first-time
    integer, parameter :: EligYrsR= 99.000000 ! Years under which agent owns the same durable to be repeat
    REAL(8), DIMENSION(:), ALLOCATABLE :: hpnodes, r_rental, r_rental_retire
    REAL(8), DIMENSION(2,hptransLength) :: startprice
    INTEGER, dimension(9) :: polParam ! Whether or not other params change in policy period
    REAL(8), dimension(9), TARGET :: polRead, polLevel ! Level of *difference* in the parameters
    INTEGER :: TransTime

    ! Life duration (assume simulation starts at age 22 for HH)
    integer, parameter :: Tretire=38 !40 !35
    integer, parameter :: Tdie=Tretire+15

    ! Flat subsidy/tax to every living agent to ensure government budget balance
    REAL(8), DIMENSION(hptransLength+1) :: balancer
    REAL(8), DIMENSION(zgridsize, Tdie) :: adjtransfer  ! Main adjustment transfer matrix
    REAL(8) :: trans_phase_int, trans_phase_out ! Parameters for controlling phase-out of policy
    REAL(8) :: trans_phasein_int, trans_phasein_out ! Parameters for controlling phase-in of policy
    REAL(8) :: renttransfer, noadjtransfer ! Proportionate or level change in prices when adjustment occurs
    REAL(8) :: renttransfer_internal, noadjtransfer_internal
    LOGICAL, DIMENSION(zgridsize, Tdie) :: pctageflag ! Controls if transfer params are level change or proportion changes
    REAL(8), parameter :: Eta=1.0 ! 0-1: Smoothing how much subsidy applies to down
    REAL(8), parameter :: Eta_transfer=1.0 ! Eta but for monetary transfers
    REAL(8) :: hpdelta ! (risky) scalar rewritten during transition period loops
                       ! to record one-period changes in house prices
    REAL(8) :: balancer_internal
    REAL(8), DIMENSION(hptransLength, hptransLength) :: Probhp


    ! Below are the state and policy holders


    REAL(8), DIMENSION(agridsize), TARGET :: anodes  !Nodes for idiosyncratic assets (voluntary equity)
    REAL(8), DIMENSION(Dgridsize), TARGET :: Dnodes  !Nodes for idiosyncratic durable

    REAL(8), DIMENSION(zgridsize) :: znodes  !Nodes for idiosyncratic productivity
    REAL(8), DIMENSION(transgridsize), TARGET :: transcst, transcstProbs !Nodes for idiosyncratic transaction costs
    REAL(8), DIMENSION(zgridsize, zgridsize) :: Probz, Probzcum  !while probability of being at any z will vary with time, transition probability will not, so max we need to loop over wont
    REAL(8), DIMENSION(zgridsize) :: Probinit, Probinitcum
    INTEGER, DIMENSION(zgridsize), TARGET :: maxexpectationz, minexpectationz

    REAL(8), DIMENSION(Tdie) :: utilscale  ! Utility scale by age
    REAL(8), DIMENSION(Tdie-Tretire), TARGET :: deathprob  ! Probability of death in retirements
    REAL(8), DIMENSION(Tretire) :: ageearnings  ! Deterministic age component of income
    REAL(8), DIMENSION(zgridsize, Tdie) :: income, inctax  ! Levels (not logs) of income and taxes


    ! Algorithm parameters.
    INTEGER, PARAMETER :: vars=2
    REAL(8), DIMENSION(vars,vars+1,2) :: amoebaGrid
    REAL(8) :: steady_conv_thres, transition_conv_thres, pthres, ftol

    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: achoice, Dchoice, cchoice !Policy functions when solving problem
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: achoiceMov, DchoiceMov, cchoiceMov
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: achoiceMovR, DchoiceMovR, cchoiceMovR
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: choiceindicator, choiceindicatorMov !choice indicator= 1 if adjust, 2 if noadjust, 3 if rent
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: achoiceshock, Dchoiceshock, &
        cchoiceshock, choiceindicatorshock
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: achoiceshockMov, DchoiceshockMov, cchoiceshockMov, choiceindicatorshockMov
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: achoiceshockMovR, DchoiceshockMovR, cchoiceshockMovR
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: achoiceexpect, Dchoiceexpect, &
        cchoiceexpect, choiceindicatorexpect
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: achoiceexpectMov, DchoiceexpectMov, cchoiceexpectMov, choiceindicatorexpectMov
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: achoiceexpectMovR, DchoiceexpectMovR, cchoiceexpectMovR
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: chindMR, chindshockMR, chindexpectMR
    REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, TARGET :: EV, EVMov, EVMovR, EVpol, EVpolMov, EVpolMovR

    ! TODO: get rid of this, let each module have their own tickers.
    !INTEGER :: i, j, k, l, m, n, t, iter

    INTEGER, parameter :: numhouseholds=30000   !Size of household panel
    REAL(8), dimension(:,:), ALLOCATABLE :: shock_calibrate
    REAL(8), dimension(:,:,:), ALLOCATABLE :: shock
    REAL(8), dimension(:,:), ALLOCATABLE :: incshocks
    REAL(8), dimension(:,:,:), ALLOCATABLE :: householdresultsholder
    REAL(8), dimension(:,:,:,:), ALLOCATABLE :: householdtransitionholder
    REAL(8), dimension(:,:,:), ALLOCATABLE :: housingflow, housingstock, rentalflow,&
        bequestflow, newownerflow, resaleflow
    INTEGER, dimension(:,:,:), ALLOCATABLE :: numadjust
    integer, dimension(33) :: seedvalue


    REAL(8) :: diffmoments, diffmoments1, diffmoments2, diffmoments3
    REAL(8), parameter :: difftol=.01
    INTEGER :: numiter

    LOGICAL, parameter :: costlyequity=.FALSE.  ! Something about ability to borrow?
    REAL(8), parameter :: maint= 1.0  ! Exogeneous ability to maintain durable
    REAL(8), parameter :: scrapped=0 ! Nothing received for sale of durable
    REAL(8), parameter :: scrapvalue=0.02  ! Fixed value when depreciated durable is sold
    REAL(8), parameter :: downflag=1  ! Subsidy applies only to down or not
    REAL(8), parameter :: rebatedflag=0  ! Transaction costs rebated
    REAL(8), parameter :: discountflag=0  ! Subsidy applies to face value but not down
    REAL(8), parameter :: himdflag=0  ! Home Interest Mortgage Deduction in effect

    ! If TRUE, solves the model in GE for hptransLength periods
    LOGICAL, parameter :: ge_start=.FALSE.  
    ! If TRUE, solves the model in PE given GE steady-state
    LOGICAL, parameter :: pe_start=.TRUE.  
    ! If TRUE, solves only GE steady-state and not transition (takes too long)
    LOGICAL, parameter :: ss_only=.FALSE.  
    ! If TRUE, loads in the price law of motion rather than solving for it
    LOGICAL, parameter :: loadpricepath=.FALSE.  
    ! If TRUE, activates parameter changes taking place in policy period
    LOGICAL :: agg_policies
    ! If TRUE, draw from custom starting values for asset/housing states
    LOGICAL, parameter :: custom_start=.FALSE.  
    ! If TRUE, draw from custom starting values for income states
    LOGICAL, parameter :: custom_incs=.FALSE.  
    ! If TRUE, no unemployment (and one less state dimension)
    LOGICAL, parameter :: no_unemp=.TRUE.  
    ! If TRUE, no deterministic age component is read
    LOGICAL, parameter :: no_lifeinc=.FALSE.  
    ! If TRUE, no one dies after retirement
    LOGICAL, parameter :: no_death=.FALSE.  
    ! If TRUE, consumption in certain periods not weighed more than others
    LOGICAL, parameter :: no_utilscale=.FALSE.  
    ! If TRUE, displays a set of policy functions when solving DP problem
    LOGICAL, parameter :: show_pol=.FALSE.  
    ! If TRUE, exports the largest microdata files from simulations to txt files
    LOGICAL, parameter :: print_micro=.TRUE.  

    REAL(8), parameter :: rentelasticity=1


    CONTAINS
    
    subroutine open_log
    OPEN (UNIT=0, FILE="model_log.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    ! TODO: Cram all the output files in this subroutine?

    end subroutine
    
    subroutine share_arrays
        IMPLICIT NONE

        ! House price grid - initial guesses for stationary equilibrium prices
        ALLOCATE(hpnodes(hptransLength+1))
        hpnodes(1)=0.0
        ALLOCATE(r_rental(hptransLength), r_rental_retire(hptransLength))
        write(0,*) "// Equilibrium price guess: ", exp(hpnodes(1))

        startprice(1,1)= 0.08
        startprice(2,1)= -0.08
        startprice(:,1)= hpnodes(1) + startprice(:,1)

        startprice(1,2:hptransLength) = 0.05
        startprice(2,2:hptransLength) = -0.05

        ! Different policy regimes (prices, subsidies)
        adjTransfer= 0
        rentTransfer= 0
        noadjTransfer= 0
        polLevel = 0
        polRead = 0

        ! Switches entire policy regime btwn level and proportion terms
        ! (Later programs can vary regime)
        agg_policies= .FALSE.
        pctageflag= .FALSE.

        trans_phase_int= 0! 1.33
        trans_phase_out= 0 ! 1.6
        trans_phasein_int= 0 ! 0.03
        trans_phasein_out= 0 ! 1.2

        ! Initial guesses for budget balance transfers
        balancer= 0.0000

        seedvalue(:)=195423
        CALL random_seed(put=seedvalue(:))

        ALLOCATE(shock(5, numhouseholds, TDie),shock_calibrate(4, numhouseholds),&
                 householdresultsholder(15, numhouseholds, TDie))

        CALL random_number(shock_calibrate)
        CALL random_number(shock)

        steady_conv_thres= 5e-2
        transition_conv_thres= 3e-2
        pthres= 5e-4
        ftol= 5e-9


        !Policy functions when solving problem
        ALLOCATE(achoice(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1),&
             Dchoice(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1),&
             cchoice(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1))
        ALLOCATE(choiceindicator(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1))
        ALLOCATE(EV(zgridsize, 2, Dgridsize, agridsize, Tdie+1, hptransLength+1))

        !Policy functions for storing subsidy policy response
        ALLOCATE(achoiceshock(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             Dchoiceshock(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             cchoiceshock(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd))
        ALLOCATE(choiceindicatorshock(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd))

        !Policy functions for storing response given policy in the future
        ALLOCATE(achoiceexpect(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             Dchoiceexpect(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             cchoiceexpect(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd))
        ALLOCATE(choiceindicatorexpect(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd))
        ALLOCATE(EVpol(zgridsize, 2, Dgridsize, agridsize, Tdie+1, hptransLength+1))

        ! Moving shock arrays (max over fewer things?)
        ALLOCATE(achoiceMov(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1),&
             achoiceMovR(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1),&
             achoiceshockMov(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),& 
             achoiceshockMovR(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),& 
             achoiceexpectMov(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             achoiceexpectMovR(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             DchoiceMov(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1),&
             DchoiceMovR(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1),&
             DchoiceshockMov(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),& 
             DchoiceshockMovR(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),& 
             DchoiceexpectMov(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             DchoiceexpectMovR(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             cchoiceMov(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1),&
             cchoiceMovR(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1),&
             cchoiceshockMov(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),& 
             cchoiceshockMovR(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),& 
             cchoiceexpectMov(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             cchoiceexpectMovR(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             choiceindicatorMov(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1),&
             choiceindicatorshockMov(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             choiceindicatorexpectMov(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             chindMR(zgridsize, 2, Dgridsize, agridsize, Tdie, hptransLength+1),&
             chindshockMR(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             chindexpectMR(zgridsize, 2, Dgridsize, agridsize, Tdie, PolEnd),&
             EVMov(zgridsize, 2, Dgridsize, agridsize, Tdie+1, hptransLength+1),&
             EVMovR(zgridsize, 2, Dgridsize, agridsize, Tdie+1, hptransLength+1),&
             EVpolMov(zgridsize, 2, Dgridsize, agridsize, Tdie+1, hptransLength+1),&
             EVpolMovR(zgridsize, 2, Dgridsize, agridsize, Tdie+1, hptransLength+1))


        if (scrapped == 1) then
        amoebaGrid(:,:,1)= RESHAPE( (/ .01, .8, .25, .25, .5, .2 /), &
                                    shape(amoebaGrid(:,:,1) ))
        amoebaGrid(:,:,2)= RESHAPE( (/ .0, .01, .05, .1, .4, .21 /), &
                                    shape(amoebaGrid(:,:,2) ))
        else
        !Work in progress
        amoebaGrid(:,:,1)= RESHAPE( (/ 0.99, 0.45, .25, 0.01, 0.10, 1.30 /), &
                                    shape(amoebaGrid(:,:,1) ))
        amoebaGrid(:,:,2)= RESHAPE( (/ 0.1, 0.01, 0.05, 0.25, 0.01, 0.05 /), &
                                    shape(amoebaGrid(:,:,2) ))
        end if

    end subroutine share_arrays

end module share
