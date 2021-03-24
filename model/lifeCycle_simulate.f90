module lifecycle_simulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE SIMULATE contains the steady-state simulation, outputs calibration
! statistics and also contains routines used to find stationary equilibrium.
! problem, as well as certain debug/test subroutines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    USE lifecycle_calibrate
    USE lifecycle_solveDP
    USE share
    USE OMP_LIB
    IMPLICIT NONE
    REAL(8), dimension(:,:), ALLOCATABLE :: alive, consumption, incomeholder, mortint, mpc
    REAL(8), dimension(:,:), ALLOCATABLE :: numrent, numowntotal, numresale, numbequest
    REAL(8), dimension(:,:), ALLOCATABLE :: posttaxholder, posttaxincome
    REAL(8) :: conditionalindicator, constrained
    REAL(8) :: exanteEVborn, exanteEVoverall, numobswelfare
    REAL(8) :: ratiocrent, ratiocown, numrent_t, numown_t, numbuy_t, numsell_t, fracown
    REAL(8) :: leverage_t, a_t, leverageown_t, aown_t, leveragerent_t, arent_t

    CONTAINS

    subroutine simulate(price, numHH, household_nextperiod, shock, write_files)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    The core of the module: this subroutine executes the simulation
        !    of agents in steady-state, subject to shocks, and writes microdata.
        !
        !    Modified: 08/09/2019
        !
        !    PARAMETERS
        !
        !    price: Real: records the price per unit of housing services.
        !
        !    numHH: Marks the number of agents simulated. A total of
        !    numHH*TDie cycle are done, and the number of simulations is lesser
        !    than or equal that due to death.
        !
        !    household_nextperiod: Array for storing the next-period policy
        !    functions of all simulated agents, as well as other statistics
        !    and policy-eligibility status. Used in transition subroutine.
        !
        !    shock: The array of shocks generated before execution to be used
        !    in the simulation. Ensures the same shocks are used every time
        !    to minimize noise.
        !
        !    write_files: If TRUE, overwrite existing data files outputted
        !    from the simulation. E.g. this is turned off as the stationary
        !    equilibrium is found.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: price
        INTEGER, INTENT(IN) :: numHH
        INTEGER :: i, j, t
        LOGICAL, INTENT(IN) :: write_files
        REAL(8), DIMENSION(:,:,:), INTENT(INOUT) :: household_nextperiod, shock
        REAL(8), dimension(:,:), ALLOCATABLE :: currenthouseholdstate, newhouseholdstate
        REAL(8), dimension(:,:), ALLOCATABLE :: consumptionmpc, rental_mpc, consumptionhpshock, durable_mpc, currenttemp, currentperm, rentalind, choiceind, movedind
        REAL(8), dimension(:,:), ALLOCATABLE :: durableinvestment, actualwealth, financialwealth, housingnet, housinggross, taxrate, welfare, diffc, ratioc, qchg
        REAL(8), dimension(:,:), ALLOCATABLE :: adjustmean
        REAL(8), dimension(numHH) :: insurancetemp, insuranceperm, cov1, cov2, cov3, cov4, insurancedtemp, insurancedperm, insurancedtempconditional, insurancedpermconditional, numtrans
        real(8), dimension(numHH, 3) :: fthb_flag
        REAL(8) :: mpc_shifter
        ! 1: age when purchase fthb, 2: last period when house sold,
        ! 3: last period when house bought

        ! %< Setting out output arrays
        if (write_files) then

        OPEN (UNIT=41, FILE="fthb.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        OPEN (UNIT=42, FILE="dist_fthb.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        OPEN (UNIT=43, FILE="lifecycle_adjust.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        OPEN (UNIT=44, FILE="lifecycleprofiles.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        OPEN (UNIT=45, FILE="housingstock.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        end if
        ! %>

        ! Allocation of automatic arrays for OpenMP (skip this) %<
        ALLOCATE(mpc(numHH, Tdie), consumption(numHH, Tdie), consumptionmpc(numHH, Tdie), &
        consumptionhpshock(numHH, Tdie), durable_mpc(numHH, Tdie), rental_mpc(numHH, Tdie), mortint(numHH, Tdie),&
        currenttemp(numHH, Tdie), currentperm(numHH, Tdie), incomeholder(numHH, Tdie), posttaxincome(numHH, Tdie),& 
        rentalind(numHH, Tdie), choiceind(numHH, Tdie), movedind(numHH, Tdie), durableinvestment(numHH, Tdie), actualwealth(numHH, Tdie),&
        financialwealth(numHH, Tdie), housingnet(numHH, Tdie), housinggross(numHH, Tdie), taxrate(numHH, Tdie), &
        welfare(numHH, Tdie), diffc(numHH, Tdie), ratioc(numHH, Tdie), qchg(numHH, Tdie))

        ALLOCATE(adjustmean(numHH, Tdie))

        ALLOCATE(alive(numHH, Tdie), numowntotal(numHH, Tdie),&
                 numrent(numHH, Tdie), numresale(numHH, Tdie),&
                 numbequest(numHH, Tdie))

        ALLOCATE(currenthouseholdstate(6, numHH),newhouseholdstate(6, numHH))
        ! %>

        ! STATE DEFINITION
        ! state 1: idiosyncratic income shock: z
        ! state 2: unemployment status 'ue'
        ! state 3: durable assets d
        ! state 4: liquid assets a
        ! state 5: Transition period (taken as subroutine argument)
        ! state 6: age: t

        ! %< Initializion and allocation for simulation
        housingstock=0
        housingflow=0
        rentalflow=0
        resaleflow=0
        newownerflow=0
        bequestflow=0
        numadjust=0
        alive=1

        ! flag for when purchasing a house (Tdie+1 = not buying a house in lifetime)
        fthb_flag = Tdie + 1 
        rentalind = 0
        movedind = 0

        call cohort_calibrate(numHH, price, currenthouseholdstate(1:4,:),&
                              fthb_flag(:,1), rentalind(:,1), shock_calibrate)
        fthb_flag(:,3) = fthb_flag(:,1) ! Durable owned at start

        ! Placeholder data if we need state vars for first age cohort
        if (write_files) then
        do i=1, numHH
        write(41, '(2I6.2, 22F16.6)') i, 0, currenthouseholdstate(3, i),&
            currenthouseholdstate(4, i)-(1-theta)*exp(price)*currenthouseholdstate(3, i),&
            0d0, currenthouseholdstate(4, i)+theta*exp(price)*currenthouseholdstate(3, i),&
            rentalind(i,1),2d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,currenthouseholdstate(1, i),0d0,0d0,0d0,0d0 
        end do
        end if

        consumption=0
        incomeholder=0
        posttaxincome=0
        mortint=0
        consumptionmpc=0
        durable_mpc=0
        rental_mpc=0
        mpc=0
        mpc_shifter=7.5e-2

        ! 5th state is age (or equiv) time
        currenthouseholdstate(5, :)=1
        ! 1 indicates starting/ending steady-state values
        currenthouseholdstate(6, :)=1
        newhouseholdstate(6, :)=1

        numowntotal=0
        numrent=0
        numtrans=0
        numresale=0
        numbequest=0
        constrained=0
        ! %>

        do t=1, Tdie

        ! %< More legacy allocation
        ratiocrent=0
        ratiocown=0
        numrent_t=0
        numown_t=0
        numbuy_t=0
        numsell_t=0

        leverage_t=0
        a_t=0
        leverageown_t=0
        aown_t=0
        leveragerent_t=0
        leverageown_t=0
        ! %>

!$OMP PARALLEL DO
            do i=1, numHH

                ! %< first thing is that if you are old, do you die or not
                if ((t>Tretire) .AND. alive(i, t) == 1) then
                   if (shock(1,i,t) < deathprob(t-Tretire)) then
                       alive(i, t:Tdie)=0
                       household_nextperiod(5,i,t)=1
                       household_nextperiod(5,i,t+1:Tdie)=0
                       bequestflow(i, t, 1)=currenthouseholdstate(3, i)
                       numbequest(i, t) = 1
                       cycle
                   end if
                end if
                ! %>

                if (alive(i, t)==1) then

                actualwealth(i, t)=currenthouseholdstate(4, i)+theta*currenthouseholdstate(3, i)*exp(price)  !total wealth is voluntary equity plus equity in durable
                financialwealth(i, t)=currenthouseholdstate(4, i)-(1-theta)*currenthouseholdstate(3, i)*exp(price)
                if (financialwealth(i, t) < 0) mortint(i, t) = -(rborrow-r)*financialwealth(i, t)
                housingnet(i, t)=theta*currenthouseholdstate(3, i)*exp(price)
                housinggross(i, t)=currenthouseholdstate(3, i)*exp(price)

                !write(*,*) "here"

                ! %< Interpolation of next period policy. Everything else in
                ! the simulation subroutine is just bean counting.
                currenthouseholdstate(4, i) = currenthouseholdstate(4, i)+mpc_shifter
                ! Exogenous shock if HH currently renter
                if (shock(2,i,t) < movProbR .AND. currenthouseholdstate(3, i) == 0) then
                    call pol_linworking(currenthouseholdstate(:,i),achoiceMovR,DchoiceMovR,cchoiceMovR,&
                        EVmovR,chindMR,newhouseholdstate(4, i), newhouseholdstate(3, i),&
                        consumptionmpc(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    currenthouseholdstate(4, i) = currenthouseholdstate(4, i)-mpc_shifter
                    call pol_linworking(currenthouseholdstate(:,i),achoiceMovR,DchoiceMovR,cchoiceMovR,&
                        EVmovR,chindMR,newhouseholdstate(4, i), newhouseholdstate(3, i),&
                        consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    movedind(i, t) = 1
                ! Exogeneous shock if HH currently owner
                else if (shock(2,i,t) < movprob .AND. currenthouseholdstate(3, i) > 0) then
                    call pol_linworking(currenthouseholdstate(:,i),achoiceMov,DchoiceMov,cchoiceMov,EVmov,&
                        choiceindicatorMov,newhouseholdstate(4, i), newhouseholdstate(3, i),&
                        consumptionmpc(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    durable_mpc(i, t) = newhouseholdstate(3, i)
                    rental_mpc(i, t) = rentalind(i, t)
                    currenthouseholdstate(4, i) = currenthouseholdstate(4, i)-mpc_shifter
                    call pol_linworking(currenthouseholdstate(:,i),achoiceMov,DchoiceMov,cchoiceMov,EVmov,&
                        choiceindicatorMov,newhouseholdstate(4, i), newhouseholdstate(3, i),&
                        consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    movedind(i, t) = 2
                else 
                ! Derive MPC (dC/dW)
                call pol_linworking(currenthouseholdstate(:,i),achoice,Dchoice,cchoice,EV,&
                    choiceindicator,newhouseholdstate(4, i), newhouseholdstate(3, i),&
                    consumptionmpc(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                durable_mpc(i, t) = newhouseholdstate(3, i)
                rental_mpc(i, t) = rentalind(i, t)
                ! Actual a/D policy
                    currenthouseholdstate(4, i) = currenthouseholdstate(4, i)-mpc_shifter
                call pol_linworking(currenthouseholdstate(:,i),achoice,Dchoice,cchoice,EV,&
                    choiceindicator,newhouseholdstate(4, i), newhouseholdstate(3, i),&
                    consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    movedind(i, t) = 0
                end if

                mpc(i, t) = (consumptionmpc(i, t) - consumption(i, t))/mpc_shifter

                ! To prevent any case of bad interpolation (but wouldn't
                ! pol_linworking catch this?
                if (newhouseholdstate(3, i)>Dmax) newhouseholdstate(3, i)=(1.0-1e-7)*Dmax
                if (newhouseholdstate(4, i)>amax) newhouseholdstate(4, i)=(1.0-1e-7)*amax
                ! %>

                ! %< Impute pure transitory shock going into next period
                if (shock(5,i,t) < 0.5 ) then
                    newhouseholdstate(4, i) = newhouseholdstate(4, i) - sigma_temp
                else
                    newhouseholdstate(4, i) = newhouseholdstate(4, i) + sigma_temp
                end if

                ! %>

                ! %< Recover next period income shock
                incomeholder(i, t)=(currenthouseholdstate(2, i) - 1.0)*&
                                   income(currenthouseholdstate(1, i), t)
                ! Note the post-tax income does not include transfers for stability?
                posttaxincome(i, t)=(currenthouseholdstate(2, i) - 1.0)*&
                                    (income(currenthouseholdstate(1, i), t)-&
                                    inctax(currenthouseholdstate(1, i), t))-&
                                    dtau*exp(price)*currenthouseholdstate(3, i)
                                    
                if (t<=Tretire) then
                    if (currenthouseholdstate(2, i) == 2 .AND. &
                        shock(3,i,t) < empprob(t+1)) then
                        ! Remaining employed
                        newhouseholdstate(2, i) = currenthouseholdstate(2, i)
                        ! Permanent income state
                        newhouseholdstate(1, i) = incshocks(i, t+1)
                    else
                        newhouseholdstate(1, i) = currenthouseholdstate(1, i)
                        if (currenthouseholdstate(2, i) == 2) then
                            ! Entering unemployment
                            newhouseholdstate(2, i) = 1
                        else if (currenthouseholdstate(2, i) == 1 .AND. &
                                 shock(3,i,t) >= unempprob) then
                            ! Leaving unemployment
                            newhouseholdstate(2, i) = 2
                        else
                            ! Remaining in unemployment
                            newhouseholdstate(2, i) = currenthouseholdstate(2, i)
                        end if
                        if (t == Tretire) newhouseholdstate(2, i) = 2
                    end if

                else
                    newhouseholdstate(1:2, i)=currenthouseholdstate(1:2, i)
                end if
                ! %>

                ! Generate fixed cost shock
                ! Notice we allow arbitrary shifts for assets in the future
                if (choiceind(i,t) == 1.0) newhouseholdstate(4, i) = newhouseholdstate(4,i) - &
                   F*shock(4,i,t)*newhouseholdstate(3, i)*exp(price)


                ! Track housing decisions going into period
                householdresultsholder(12, i, t)=fthb_flag(i, 1)
                householdresultsholder(13, i, t)=fthb_flag(i, 2)
                householdresultsholder(14, i, t)=fthb_flag(i, 3)
                householdresultsholder(15, i, t)=rentalflow(i, max(t-1,1), 1)
                if (newhouseholdstate(4, i)<1e-4 .and. t==1) constrained=constrained+1

                if (t > 1) then
                    ! Record lags of consumption, durables, assets
                    household_nextperiod(11,i,t)=consumption(i, t-1)
                    household_nextperiod(12,i,t)=financialwealth(i, t-1)
                    household_nextperiod(13,i,t)=welfare(i, t-1)
                end if

                durableinvestment(i, t)=newhouseholdstate(3, i)-(1-delta-dtau)*currenthouseholdstate(3, i)

                ! %< Legacy code tracking renters/owners at start of each period,
                ! plus buyers/sellers...  could be useful again.
                !if (currenthouseholdstate(3, i)==0) then
                !    numrent_t=numrent_t+1
                !    ratiocrent=ratiocrent+(ratioc(i, t)-1)/(exp(hpmin)-1)
                !    arent_t=arent_t+currenthouseholdstate(4, i)-(1-theta)*currenthouseholdstate(3, i)
                    !if (currenthouseholdstate(3, i)/(currenthouseholdstate(4, i)+theta*currenthouseholdstate(3, i)) .ne. 0/0) then
                !        leveragerent_t=leveragerent_t+currenthouseholdstate(3, i)/(currenthouseholdstate(4, i)+theta*currenthouseholdstate(3, i))
                    !end if
                !else
                !    ratiocown=ratiocown+(ratioc(i, t)-1)/(exp(hpmin)-1)
                    !aown_t=aown_t+currenthouseholdstate(4, i)-(1-theta)*currenthouseholdstate(3, i)
                    !if (currenthouseholdstate(3, i)/(currenthouseholdstate(4, i)+theta*currenthouseholdstate(3, i)) .ne. 0/0) then
                    !    leverageown_t=leverageown_t+currenthouseholdstate(3, i)/(currenthouseholdstate(4, i)+theta*currenthouseholdstate(3, i))
                    !end if
                !end if

                !if (newhouseholdstate(3, i)-currenthouseholdstate(3, i)>.001 .and. rentalind(i, t)<=.995) then
                !numbuy_t=numbuy_t+1
                !elseif (newhouseholdstate(3, i)-currenthouseholdstate(3, i)<-.001) then
                !numsell_t=numsell_t+1
                !end if
                !a_t=a_t+currenthouseholdstate(4, i)-(1-theta)*currenthouseholdstate(3, i)
                !if (currenthouseholdstate(3, i)/(currenthouseholdstate(4, i)+theta*currenthouseholdstate(3, i)) .ne. 0/0) then
                !leverage_t=leverage_t+currenthouseholdstate(3, i)/(currenthouseholdstate(4, i)+theta*currenthouseholdstate(3, i))
                !end if

                ! %>

                ! %< write(*,*) "here 2"
                householdresultsholder(1:4, i, t)=currenthouseholdstate(1:4, i)
                householdresultsholder(5, i, t)=alive(i, t)
                householdresultsholder(6, i, t)=consumption(i, t)
                householdresultsholder(7, i, t)=newhouseholdstate(4, i)
                householdresultsholder(8, i, t)=newhouseholdstate(3, i)
                householdresultsholder(9, i, t)=rentalind(i, t)
                householdresultsholder(10, i, t)=incomeholder(i, t)
                householdresultsholder(11, i, t)=welfare(i, t)
                ! %>

                ! Decision tree begins with unambiguous renting
                if ((rentalind(i, t) > .995) .or.  (rentalind(i, t) /= rentalind(i, t))) then
                    ! %<
                    rentalflow(i, t, 1)=newhouseholdstate(3, i)
                    newhouseholdstate(3, i)=0
                    fthb_flag(i, 3) = Tdie+1
                    numrent(i, t)=1
                    if (currenthouseholdstate(3, i) > 0) then
                        resaleflow(i, t, 1) = currenthouseholdstate(3, i)
                        newownerflow(i, t, 1) = newhouseholdstate(3, i)
                        numresale(i, t) = 1
                    end if 
                ! Starts the FTHB policy timer.
                    if ((fthb_flag(i, 1) /= Tdie+1) .and. (fthb_flag(i, 2) == Tdie+1)) then
                        fthb_flag(i, 2) = t
                    end if
                ! Track agents again as FTHBs after three consecutive years of renting
                    if (t >= fthb_flag(i, 2) + (EligYrsF - 1)) then
                        fthb_flag(i, 1) = Tdie+1
                        fthb_flag(i, 2) = Tdie+1
                    end if
                    ! %>
                else
                    ! %<
                    ! Housingstock measures demand at end of period. This line
                    ! is for the nonadjusters
                    housingstock(i, t, 1)=currenthouseholdstate(3, i)
                    numowntotal(i, t)=1
                    ! Any house purchase resets timer
                    fthb_flag(i, 2) = Tdie+1
                ! Track agents as eligible repeat buyer if their last purchase old enough
                        if (t >= fthb_flag(i, 3) + EligYrsR) fthb_flag(i, 3) = -1 
                        if (choiceind(i, t) == 1.0) then
                ! %< Start the repeat buyer policy countdown whenever a durables purchase is made
                            numtrans(i) = 1
                            numadjust(i, t, 1) = 1
                            housingflow(i, t, 1)=newhouseholdstate(3, i)
                            housingstock(i, t, 1)=newhouseholdstate(3, i)
                            newownerflow(i, t, 1) = newhouseholdstate(3, i)
                            if (currenthouseholdstate(3, i) > 0) then
                                resaleflow(i, t, 1) = currenthouseholdstate(3, i)
                                numresale(i, t) = 1
                            end if
                ! Record stats for agent if they are eligible for policy
                ! (The final negative indicators are factors for which policy
                !  was taken since Matlab doesn't handle mixed types well)
                            if (fthb_flag(i, 1) == Tdie+1) then
                                if (write_files) then
                                write(42, '(2I6.2, 14F19.6, I6.2)') i, t, newhouseholdstate(3, i),&
                                newhouseholdstate(4, i) - (1-theta)*exp(price)*newhouseholdstate(3, i),&
                                rentalflow(i, max(t-1,1), 1), 0.0, movedind(i,t), &
                                financialwealth(i,t), consumption(i,t), welfare(i,t),&
                                householdresultsholder(1, i, t), householdresultsholder(10, i, t),&
                                householdresultsholder(1, i, max(t-1,1)), householdresultsholder(10, i, max(t-1,1)),&
                                household_nextperiod(11:12,i,t),-1
                                end if
                            else if (fthb_flag(i, 3) == -1) then
                                if (write_files) then
                                write(42, '(2I6.2, 14F19.6, I6.2)') i, t, newhouseholdstate(3, i),&
                                newhouseholdstate(4, i) - (1-theta)*exp(price)*newhouseholdstate(3, i),&
                                currenthouseholdstate(3, i), fthb_flag(i, 1),  movedind(i,t), &
                                financialwealth(i,t), consumption(i,t), welfare(i,t),&
                                householdresultsholder(1, i, t), householdresultsholder(10, i, t),&
                                householdresultsholder(1, i, max(t-1,1)), householdresultsholder(10, i, max(t-1,1)),&
                                household_nextperiod(11:12,i,t),-2
                                end if
                            end if
                            fthb_flag(i, 1) = t
                            fthb_flag(i, 3) = t
                        ! %>
                        end if
                    ! %>
                end if

                if (write_files .AND. print_micro) then
                write(41, '(2I6.2, 22F19.6)') i, t, max(newhouseholdstate(3, i),rentalflow(i,t,1)),&
                    newhouseholdstate(4, i)-(1-theta)*exp(price)*newhouseholdstate(3, i),&
                    max(currenthouseholdstate(3, i),rentalflow(i, max(t-1,1), 1)), &
                    newhouseholdstate(4, i)+theta*exp(price)*newhouseholdstate(3, i),&
                    rentalind(i,t),3-choiceind(i,t), movedind(i, t), financialwealth(i,t), consumption(i,t), &
                    consumptionmpc(i, t), durable_mpc(i, t), rental_mpc(i, t), welfare(i,t),&
                    householdresultsholder(12, i, t), householdresultsholder(1, i, t), householdresultsholder(10, i, t),&
                    householdresultsholder(1, i, max(t-1,1)), householdresultsholder(10, i, max(t-1,1)),&
                    householdresultsholder(1, i, max(t-2,1)), householdresultsholder(10, i, max(t-2,1)),&
                    household_nextperiod(11:12,i,t)
                end if

                ! %<
                household_nextperiod(1:5,i,t)=householdresultsholder(1:5,i,t)
                household_nextperiod(6:9,i,t)=householdresultsholder(12:15,i,t)
                household_nextperiod(10,i,t)=householdresultsholder(10,i,t)

                ! End of life means liquidation of housing
                if (t==TDie) bequestflow(i, t, 1)=newhouseholdstate(3, i)
                if (t==TDie) numbequest(i, t) = 1
                ! Degree to which durable can be maintained
                newhouseholdstate(3, i) = (1-delta+maint*delta)*newhouseholdstate(3, i)
                ! %>

                newhouseholdstate(5, i)= t+1
                end if ! Checking if person is alive
            end do
!$OMP END PARALLEL DO

            currenthouseholdstate = newhouseholdstate
            if (write_files) then
            write(44, '(17F16.6)') t*1.0, sum(alive(:, t)), &
            sum(numrent(:, t)), sum(alive(:, t)*(1-rentalind(:, t))), &
            sum(alive(:,t)*consumption(:, t))/sum(alive(:, t)), sum(alive(:,t)*mpc(:, t))/sum(alive(:, t)),&
            sum(alive(:,t)*consumption(:, t)*(1-rentalind(:, t)))/sum(alive(:, t)*(1-rentalind(:, t))), &
            sum(alive(:,t)*exp(price)*newhouseholdstate(3, :)*(1-rentalind(:, t)))/sum(alive(:, t)*(1-rentalind(:, t))), &
            sum(alive(:,t)*r_rental(1)*rentalflow(1:numHH,t,1)*rentalind(:, t))/sum(numrent(:, t)), &
            sum(alive(:,t)*newhouseholdstate(4, :))/sum(alive(:, t)), &
            sum(alive(:,t)*(newhouseholdstate(4, :)*(1-rentalind(:, t))))/sum(alive(:, t)*(1-rentalind(:, t))), &
            sum(alive(:,t)*(newhouseholdstate(4, :)*rentalind(:,t)))/sum(numrent(:, t)), &
            sum(alive(:,t)*(newhouseholdstate(4, :)-(1-theta)*exp(price)*newhouseholdstate(3, :))*(1-rentalind(:, t)))/sum(alive(:, t)*(1-rentalind(:, t))), &
            sum(alive(:,t)*(exp(price)*newhouseholdstate(3, :)/(newhouseholdstate(4, :)+theta*exp(price)*newhouseholdstate(3, :)))*(1-rentalind(:, t)))/sum(alive(:, t)*(1-rentalind(:, t))), &
            sum(alive(:,t)*(newhouseholdstate(4, :)+theta*exp(price)*newhouseholdstate(3, :))*(1-rentalind(:, t))), &
            sum(alive(:, t)*(incomeholder(:, t))), 1-(sum(numrent(:, t))/sum(alive(:, t)))
            end if

            if (write_files) then
            write(45, '(I6.2, 3F16.6)') t, sum(housingstock(:, t, 1)), sum(housingflow(:, t, 1)), sum(numtrans(:))
            end if
            numtrans = 0

        end do

!$OMP PARALLEL DO
        do t=1, Tdie
        do i=1, numHH
            if (write_files) write(43, '(3I6.2)') i, t, sum(numadjust(i,t:min(t+Tretire-1,Tdie), 1))
        end do
        end do
!$OMP END PARALLEL DO

        DEALLOCATE(consumptionmpc, consumptionhpshock, durable_mpc, &
        currenttemp, currentperm, rentalind, choiceind, movedind, &
        durableinvestment, actualwealth, financialwealth, &
        housingnet, housinggross, taxrate, welfare, &
        diffc, ratioc, qchg)

        DEALLOCATE(numowntotal,numrent,numresale,numbequest)

        DEALLOCATE(currenthouseholdstate,newhouseholdstate)

        close(41)
        close(42)
        close(43)
        close(44)
        close(45)

    end subroutine simulate ! %>

    REAL(8) FUNCTION SimSteady(price, numHH, iter)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    This function is called in every iteration of the algorithm used
        !    to find the stationary equilibrium. It solves the DP problem,
        !    runs the simulation and checks market clearing.
        !
        !    Modified: 06/22/2017
        !
        !    PARAMETERS
        !
        !    price: Real: records the price per unit of housing services.
        !
        !    numHH: Marks the number of agents simulated. A total of
        !
        !    iter: Total iterations of the algorithm so far. May be necessary
        !    to track in subroutine for efficiency.
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE lifecycle_solveDP
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: price
        INTEGER, INTENT(IN) :: numHH, iter
        REAL(8) :: avghousing, construction, totrev

        write(*,*) "// Price: ", EXP(price)
        write(0,*) NEW_LINE('A') // "// Price: ", EXP(price)
        ! Agents face the competitive price for rental housing
        r_rental = rentPrem + exp(price)*rent
        r_rental_retire = rentPremRetire + exp(price)*rent
        call bequests((/price/), 1)
        write(0,*) NEW_LINE('A') // "Begin value function iteration:"
        call solveworkingproblem((/price/), (/balancer_internal/),&
            achoice(:,:,:,:,:,1:1), Dchoice(:,:,:,:,:,1:1), rentchoice(:,:,:,:,:,1:1),&
            cchoice(:,:,:,:,:,1:1), choiceindicator(:,:,:,:,:,1:1),&
            achoiceMov(:,:,:,:,:,1:1), DchoiceMov(:,:,:,:,:,1:1),&
            cchoiceMov(:,:,:,:,:,1:1), choiceindicatorMov(:,:,:,:,:,1:1),&
            achoiceMovR(:,:,:,:,:,1:1), DchoiceMovR(:,:,:,:,:,1:1),&
            cchoiceMovR(:,:,:,:,:,1:1),chindMR(:,:,:,:,:,1:1),.TRUE., .TRUE., 'None',&
            EV(:,:,:,:,:,1:1), EVMov(:,:,:,:,:,1:1), EVMovR(:,:,:,:,:,1:1))
        write(0,*) "Value function iteration complete"

        write(0,*) "Simulation starting..."
        call simulate(price, numHH, householdtransitionholder(:,1:numHH,:,1),&
                      shock(:,1:numHH,:),.FALSE.)
        call get_mktclear(price, 1, avghousing, construction, totrev)

        ! Balancer calibrated in every iteration, so no need to throw into
        ! objective fn
        ! This is an ad-hoc "smoothing method": to make sure algorithm isn't
        ! all dependent on balanced budget instead of prices, decisions are made
        ! by putting more weight on the last iterated budget transfer value
        if (iter > 1) then
            balancer_internal = (0.25*balancer(1)*SUM(alive)+0.75*(totrev))/SUM(alive)
            balancer(1) = (0.25*balancer(1)*SUM(alive)+0.75*(totrev))/SUM(alive)
        end if
        write(*,*) balancer(1) - balancer_internal

        DEALLOCATE(alive, consumption, incomeholder, mortint, posttaxincome, mpc)
        SimSteady = (avghousing - EXP(construction))**2.0
        
    END FUNCTION SimSteady ! %>

    subroutine get_mktclear(price, period, avghousing, construction, totrev)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    This function computes and approximates market clearing conditions
        !    in the stationary equilibrium using some algorithms.
        !
        !    Modified: 08/09/2019
        !
        !    PARAMETERS
        !
        !    price: Real: records the price per unit of housing services.
        !
        !    period: Integer: current period of the transition
        !    (base period is 1 for steady state)
        !
        !    avghousing: Output. A real number indicating the number of units
        !    of excess demand for housing services.
        !
        !    construction: Output. Real indicates the number of units of
        !    housing construction supplied at the market price.
        !
        !    totrev: Output. The approximate level of subsidy/lump-sum
        !    taxes needed per capita to balance the government's budget.    
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE lifecycle_solveDP
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: price
        INTEGER, INTENT(IN) :: period
        REAL(8) :: avghousing, construction, totrev
        REAL(8) :: laborshare, taxshare, netrevenue,&
                   bequestres, lifeinc, privrevenue, psi0, coef, pindex

        !coef = ((balancer(1) - balancer_internal)/balancer_internal)*&
        !        &SUM(mpc)/SUM(alive)
        !coef = SUM(mpc)/SUM(alive)
        pindex = EXP((1.0/(1.0-psi2)*(price+psi2*LOG(psi2))))

        ! Government budget calibration
        taxshare = (sum(incomeholder(:,1:Tretire) - posttaxincome(:,1:Tretire)) - balancer(1)*sum(alive))&
                   /sum(incomeholder(:,1:Tretire))
        ! Impose market clearing in consumption goods market to back out
        ! the share of labour going into construction industry (IN PROGRESS)
        lifeinc = sum(posttaxincome(:,1:Tretire)) !+ sum(incomeholder(:,Tretire+1:Tdie))
        write(0,*) lifeinc, balancer(1)*SUM(alive)
        privrevenue = F*SUM(housingflow(:,1:Tretire,period)) + SUM(mortint(:,1:Tretire))
        write(0,*) sum(consumption(:,1:Tretire)), privrevenue
        laborshare = (lifeinc + balancer(1)*SUM(alive)- (sum(consumption(:,1:Tretire)) + privrevenue))&
                     &/(lifeinc + balancer(1)*SUM(alive))
        ! The following equation must be from when trying to solve laborshare
        ! and netrevenue simultaneously - but subbing in post tax income may approximate it enough
        !laborshare = (lifeinc - (coef)**(-1.0)*(sum(consumption) + privrevenue) - pindex)/&
        !             &(SUM(incomeholder(:,1:Tretire))*(1.0-coef**(-1.0)) - pindex)
        !laborshare = 1.0 - laborshare
        if (laborshare <= 0.0) then
            write(*,*) "Market clearing condition violated, continuing..."
            laborshare = 1e-3
        end if
        psi0 = (laborshare*SUM(incomeholder(:,1:Tretire)))**(1.0-psi2)
        write(0,'(2A30)') "// Share of construction labor", "  |  Factor value"
        write(0,'(2F30.5)') laborshare, psi0
        construction = LOG(psi0) + psi2/(1.0-psi2)*LOG((psi2*psi0*exp(price)))

        ! First budget balance formula: government's budget is financed by
        ! an income tax and captured revenue from land, and its outlays
        ! are pension payments plus purchases of consumption good
        netrevenue = (taxshare*sum(incomeholder(:,1:Tretire))) + &
                     (dtau*EXP(price)*SUM(housingstock(:,:,period))) + &
                     (EXP(price)*EXP(construction) - &
                      laborshare*SUM(incomeholder(:,1:Tretire)))- &
                     (sum(incomeholder(:,Tretire+1:Tdie)) + ret_wealth*SUM(incomeholder(:,Tretire)))
        write(0,'(A30, F30.5)') "// Taxes as share of income:", taxshare
        write(0,'(A30, F30.5)') "// Transfer to agents:", balancer(1)
        ! Note netrevenue is actually not a residual - it is the total amount
        ! of expenses transferred toward agents or taken away. Hence taking
        ! the difference below should mean the value -> 0
        write(0,'(A40, F30.5)') "// Net government revenue per capita:", (netrevenue +&
            privrevenue - balancer(1)*SUM(alive))/SUM(alive)
        totrev = netrevenue + privrevenue

        ! Calculate new housing construction - does it
        ! 1) make up for depreciation in both the homeowner and rental
        ! industry housing stock;
        ! 2) Balances the amount of housing entering the market from homeowners
        ! entering rental and dead homeowners?
        bequestres = SUM(bequestflow(:,:,period)) - SUM(householdtransitionholder(3, :, 1, period))&
                     - SUM(householdtransitionholder(4, :, 1, period))/EXP(price)
        avghousing =  (SUM(newownerflow(:,:,period)) - SUM(resaleflow(:,:,period))) - bequestres
        avghousing = avghousing + maint*delta*(SUM(housingstock(:,:,period)) + SUM(rentalflow(:,:,period)))
        if (avghousing < 0) avghousing = 1e-6
        if (construction < 0) construction = LOG(1e-6)
        write(0,'(3A30)') "// Housing in SS", "  |  New homeowners in SS", "  |  Resold housing in SS"
        write(0,'(3F30.4)') SUM(housingstock(:,:,period))/SUM(alive), SUM(newownerflow(:,:,period))/SUM(alive), SUM(resaleflow(:,:,period))/SUM(alive)
        write(0,'(2A30)') "// Final rental homes in SS", "  |  Bequested housing in SS"
        write(0,'(2F30.4)') SUM(rentalflow(:,:,period))/SUM(alive), bequestres/SUM(alive)
        write(0,'(A40,A21)') "// Demand for housing less resale supply", "  |  New construction"
        write(0,'(F40.10,F20.10)') avghousing/SUM(alive), EXP(construction)/SUM(alive)
        
    end subroutine

end module lifecycle_simulate
