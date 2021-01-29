module lifecycle_transition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE TRANSITION is similar to the simulate module, but runs code pertaining
! to the policy transition every period. Subroutines used to find the
! Rational Expectations transition path are included.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    USE lifecycle_calibrate
    USE lifecycle_solveDP
    USE lifecycle_simulate
    USE share
    USE OMP_LIB
    IMPLICIT NONE

    REAL(8) :: m_eps2, m_eta2, m_cFy, m_yFy, m_epsc, m_etac, m_cIVy, m_yIVy, m_eps, &
    m_eta, m_c, m_y, m_IVy, m_Fy, m_d, m_dFy, m_dIVy, m_epsd, m_etad
    REAL(8) :: actualwealthtotal
    REAL(8), dimension(hptransLength, 1) :: exanteEVbornstate, numobsstate, &
    overallwelfarebystate, overallobsbystate, consumptionbystate, &
    overallwelfarebystateyoung, overallwelfarebystatemiddle, overallwelfarebystateold, overallobsbystateyoung, overallobsbystatemiddle, overallobsbystateold, &
    consumptionbystateyoung, consumptionbystatemiddle, consumptionbystateold
    REAL(8), dimension(hptransLength, 350) :: overallwelfarebystateCE, &
    overallwelfarebystateyoungCE, overallwelfarebystatemiddleCE, overallwelfarebystateoldCE
    REAL(8), dimension(hptransLength) :: cutoff
    REAL(8), dimension(350, 1) :: CE

    CONTAINS

    subroutine transition(price, numHH, period, household_byperiod,&
                          household_nextperiod, shock, partial, write_files)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    The core of the module: this subroutine executes the policy transition
        !    simulation, period by period. Additional output is written.
        !
        !    Modified: 08/06/2019
        !
        !    PARAMETERS
        !
        !    price: Real: records the price per unit of housing services.
        !
        !    numHH: Marks the number of agents simulated. A total of
        !    numHH*TDie cycle are done, and the number of simulations is lesser
        !    than or equal that due to death.
        !
        !    Period: the number of periods since the simulation began.
        !    The first period is denoted as 0, the second "1" and so forth.
        !
        !    household_byperiod: Array for reading the states, incomes
        !    and policy eligibility statuses of agents in the current policy
        !    period.
        !
        !    household_nextperiod: Array for storing the next-period policy
        !    functions of all simulated agents, as well as other statistics
        !    and policy-eligibility status. Used in transition subroutine.
        !
        !    shock: The array of shocks generated before execution to be used
        !    in the simulation. Ensures the same shocks are used every time
        !    to minimize noise.
        !
        !    partial: Boolean. If TRUE, the policy transition is run in
        !    partial equilibrium, i.e. prices remain constant during the
        !    whole transition path.
        !
        !    write_files: If TRUE, overwrite existing data files outputted
        !    from the simulation. E.g. this is turned off as the RE transition
        !    path is found.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: price
        INTEGER, INTENT(IN) :: numHH, period
        INTEGER :: i, t
        LOGICAL, INTENT(IN) :: partial, write_files
        REAL(8), DIMENSION(:,:,:), INTENT(INOUT) :: household_byperiod, household_nextperiod
        REAL(8), dimension(:,:,:), INTENT(INOUT) :: shock
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: currenthouseholdstate, transhouseholdstate, newhouseholdstate
        !REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: newteststate
        REAL(8), dimension(:,:), ALLOCATABLE :: aggregatecontrib, mpc, consumptionmpc, consumptionhpshock, durableconsumption, currenttemp, currentperm, rentalind, choiceind, subval
        REAL(8), dimension(:,:), ALLOCATABLE :: durableinvestment, actualwealth, financialwealth, housingnet, housinggross, actualwealthdividedbyincome, totalnetworthdividedbyincome, taxrate, welfare, diffc, ratioc
        REAL(8), dimension(:,:), ALLOCATABLE :: qchg, numtrans
        REAL(8), dimension(numHH) :: insurancetemp, insuranceperm, cov1, cov2, cov3, cov4, &
        insurancedtemp, insurancedperm, insurancedtempconditional, insurancedpermconditional
        real(8), dimension(numHH, 3) :: fthb_flag



        ! %< Setting out (additional) output arrays
        ! TODO: Add a flag that writes files only if prices have converged.
        if ((write_files) .AND. (period == 0)) then
        OPEN (UNIT=71, FILE="housing_transit.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        OPEN (UNIT=72, FILE="transition_adjust.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        OPEN (UNIT=73, FILE="transition_fthb.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        OPEN (UNIT=74, FILE="transition_debug.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        end if
        ! %>

        ! Allocation of automatic arrays for OpenMP (skip this) %<
        ALLOCATE(aggregatecontrib(numHH, Tdie), mpc(numHH, Tdie), &
        consumptionmpc(numHH, Tdie), consumptionhpshock(numHH, Tdie), durableconsumption(numHH, Tdie),&
        currenttemp(numHH, Tdie), currentperm(numHH, Tdie), posttaxincome(numHH, Tdie),&
        rentalind(numHH, Tdie), choiceind(numHH, Tdie), durableinvestment(numHH, Tdie), actualwealth(numHH, Tdie),&
        financialwealth(numHH, Tdie), housingnet(numHH, Tdie), housinggross(numHH, Tdie), taxrate(numHH, Tdie), &
        welfare(numHH, Tdie), diffc(numHH, Tdie), ratioc(numHH, Tdie))

        ALLOCATE(qchg(numHH, Tdie),numtrans(numHH,TDie))

        ALLOCATE(alive(numHH, Tdie), numowntotal(numHH, Tdie),&
                 numrent(numHH, Tdie), numresale(numHH, Tdie),&
                 numbequest(numHH, Tdie), subval(numHH, Tdie))

        ALLOCATE(currenthouseholdstate(7, numHH),transhouseholdstate(6, numHH),&
                 newhouseholdstate(6, numHH))
        ! %>

        consumption=0
        incomeholder=0
        posttaxincome=0
        subval=0
        mortint=0

        numowntotal=0
        numrent = 0
        numresale=0
        numbequest=0
        numtrans = 0  ! Recording transactions between cohorts for one period


        ! The loop here is solving for the entire population in each transition period
        ! Data from the Tdie-1 cohorts surviving from the preceding period is read in,
        ! And a new cohort of age 1 is initialized.
        do t =Tdie,1,-1

        ! STATE DEFINITION
        ! state 1: idiosyncratic income shock: z
        ! state 2: unemployment shock
        ! state 3: durable assets d
        ! state 4: liquid assets a
        ! state 5: age cohort in current period
        ! state 6: Transition period (taken as subroutine argument)
        ! STATE 7: Policy eligibility; controls which policy functions read

        ! Notice we load in stats from earlier simulation, where they exist
        if ((t == 1) .AND. (period > 0)) THEN
            fthb_flag(:,1) = Tdie+1
            rentalind(:,t) = 0
            call cohort_calibrate(numHH, price, currenthouseholdstate(1:4,:),&
                                  fthb_flag(:,1),rentalind(:,t), shock_calibrate)
            fthb_flag(:,3) = fthb_flag(:,1) ! Durable owned at start
            alive(:,t)=1
        else
        ! TODO: Change order of fthb_flag to get rid of loop
            currenthouseholdstate(1:4, :) = household_byperiod(1:4, :, t)
            !fthb_flag(1:3, :) = household_byperiod(5:7, :, t)
            alive(:,t)=household_byperiod(5, :, t)
            do i=1,numHH
            fthb_flag(i, 1:3) = household_byperiod(6:8, i, t)
            ! UNEMPLOYMENT SHOCK SHOULD BE AT JUST 0.25
            ! if (period == 0 .AND. t <= Tretire .AND. shock(2,i,t) < 0.25) currenthouseholdstate(2, i) = 1.0
            ! if (period == 0 .AND. currenthouseholdstate(4, i) >= 0.0 .AND. fthb_flag(i, 1) == Tdie+1) currenthouseholdstate(4, i) = 1e-4
            end do
        end if
        currenthouseholdstate(5, :) = t
        ! +1 for not being steady state, +1 for indexing shock as period 0
        ! In PE case (one transition period), stick with one state space
        if (hptransLength == 1) THEN
            currenthouseholdstate(6, :) = 2
        else
            currenthouseholdstate(6, :) = period+2
        end if
        do i = 1, numHH
            ! First-time eligibility: never bought a durable or has
            ! been renting durable for EligYrsF years. The EligYrsF
            ! <= condition ensures policy turned off if that parameter
            ! is large enough.
            if ((fthb_flag(i, 1)>=Tdie+1) &
                .and. (currenthouseholdstate(1, i) > 0) &
                .and. (EligYrsF <= Tdie) .and. (period < PolEnd) &
                .and. (period >= PolStart)) then
                currenthouseholdstate(7, i) = 1
                if (period > 0) rentalflow(i, max(t-1,1), period+1) &
                                = household_byperiod(9, i, t)
           ! Repeat eligibility: has owned the same durable level
           ! for EligYrsR years.
            else if ((t >= fthb_flag(i, 3) + EligYrsR) &
                     .and. (period < PolEnd) .and. (period >= PolStart)) then
                currenthouseholdstate(7, i) = 2
            else
                currenthouseholdstate(7, i) = 3
            end if

        end do

        !$OMP PARALLEL DO
            do i=1, numHH

                ! first thing is that if you are old, do you die or not
                if ((t>Tretire) .AND. (alive(i,t) == 1)) then
                   if (shock(1,i,t)<deathprob(t-Tretire)) then
                       alive(i, t)=0
                       household_nextperiod(5, i, t+1:Tdie)=0
                       bequestflow(i, t, period+2)=currenthouseholdstate(3, i)
                       numbequest(i, t) = 1
                       cycle
                   end if
                end if

                if (alive(i, t) == 1) then

                !total wealth is voluntary equity plus equity in durable
                actualwealth(i, t)=currenthouseholdstate(4, i)+theta*currenthouseholdstate(3, i)*&
                    &exp(price)
                financialwealth(i, t)=currenthouseholdstate(4, i)-(1-theta)*currenthouseholdstate(3, i)*&
                    &exp(price)
                if (financialwealth(i, t) < 0) mortint(i, t) = -(rborrow-r)*financialwealth(i, t)
                housingnet(i, t)=theta*currenthouseholdstate(3, i)*exp(price)
                housinggross(i, t)=currenthouseholdstate(3, i)*exp(price)


                ! %< Interpolation of next period policy. Everything else in
                ! the simulation subroutine is just bean counting.
                ! TODO: shock and expect arrays need to be indexed
                qchg(i, t) = currenthouseholdstate(4, i) + &
                    (1-(1-polLevel(3))*theta)*currenthouseholdstate(3, i)*hpdelta
                transhouseholdstate(1:6, i) = currenthouseholdstate(1:6, i)
                transhouseholdstate(4, i) = qchg(i, t)
                ! Minimizes population by changing the HP state when policy is
                ! active (so the index is 1 not 2 like other matrices)
                if (period < PolEnd) then
                    transhouseholdstate(6, i) = currenthouseholdstate(6, i) - 1
                end if
                ! The control flow here dictates if agent responds to policy shock
                if (currenthouseholdstate(7, i) /= 3 .AND. period == PolEnd - 1) then ! %<
                    if (shock(2,i,t) < movProbR .AND. currenthouseholdstate(3, i) == 0) then
                    ! These are "inattentive renters" so they do not take up policy, period
                    call pol_linworking(transhouseholdstate(:,i),&
                            achoiceMovR,DchoiceMovR,cchoiceMovR,EVpolMovR,chindMR,&
                        newhouseholdstate(4, i), newhouseholdstate(3, i),&
                        consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    else if (shock(2,i,t) < movprob .AND. currenthouseholdstate(3, i) > 0) then
                        call pol_linworking(transhouseholdstate(:,i),&
                            achoiceshockMov,DchoiceshockMov,cchoiceshockMov,EVpolMov,choiceindicatorshockMov,&
                            newhouseholdstate(4, i), newhouseholdstate(3, i),&
                            consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    else
                        call pol_linworking(transhouseholdstate(:,i),&
                        achoiceshock,Dchoiceshock,cchoiceshock,EVpol,choiceindicatorshock,&
                            newhouseholdstate(4, i), newhouseholdstate(3, i),&
                            consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    end if

                    ! Fixed durable level adjustment if turned on
                    if (agg_policies .AND. polParam(2) == 1.0 .AND. choiceind(i, t) == 1.0) then
                       newhouseholdstate(3,i) = polLevel(2)
                    end if

                    ! compare utility of scrapping to not adjusting?
                    if (scrapped > 0) then
                        if (shock(2,i,t) < movProbR .AND. currenthouseholdstate(3, i) == 0) then
                        call adj_func_comp(newhouseholdstate(3:4, i), &
                                transhouseholdstate(:,i), achoiceMovR, DchoiceMovR, cchoiceMovR, chindMR, EVpolmovR,&
                            (/.FALSE., .FALSE./), price, consumption(i,t), &
                            rentalind(i,t), choiceind(i,t), welfare(i,t))
                        else if (shock(2,i,t) < movprob .AND. currenthouseholdstate(3, i) > 0) then
                            call adj_func_comp(newhouseholdstate(3:4, i), &
                                transhouseholdstate(:,i), achoiceMov, DchoiceMov, cchoiceMov, choiceindicatorMov, EVpolmov,&
                                (/.FALSE., .FALSE./), price, consumption(i,t), &
                                rentalind(i,t), choiceind(i,t), welfare(i,t))
                        else
                            call adj_func_comp(newhouseholdstate(3:4, i), &
                            transhouseholdstate(:,i), achoice, Dchoice, cchoice, choiceindicator, EVpol,&
                                (/.FALSE., .FALSE./), price, consumption(i,t), &
                                rentalind(i,t), choiceind(i,t), welfare(i,t))
                        end if
                    end if 
                ! %>
                ! Case with two (maybe more?) policy periods: compare utility of
                ! considering policy now or considering policy later %<
                else if (currenthouseholdstate(7, i) /= 3 .AND. period >= PolStart &
                         .AND. period < PolEnd - 1) then
                    if (shock(2,i,t) < movProbR .AND. currenthouseholdstate(3, i) == 0) then
                    call pol_linworking(transhouseholdstate(:,i),&
                           achoiceMovR,DchoiceMovR,cchoiceMovR,EVpolMovR,chindMR,&
                        newhouseholdstate(4,i), newhouseholdstate(3,i),&
                        consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    else if (shock(2,i,t) < movprob .AND. currenthouseholdstate(3, i) > 0) then
                    call pol_linworking(transhouseholdstate(:,i),&
                           achoiceshockMov,DchoiceshockMov,cchoiceshockMov,EVpolMov,choiceindicatorshockMov,&
                        newhouseholdstate(4,i), newhouseholdstate(3,i),&
                        consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    else
                       call pol_linworking(transhouseholdstate(:,i),&
                        achoiceshock,Dchoiceshock,cchoiceshock,EVpol,choiceindicatorshock,&
                           newhouseholdstate(4,i), newhouseholdstate(3,i),&
                           consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    end if

                    ! Fixed durable level adjustment if turned on
                    if (agg_policies .AND. polParam(2) == 1.0 .AND. choiceind(i, t) == 1.0) then
                       newhouseholdstate(3,i) = polLevel(2)
                    end if
                    ! But comparison only necessary if policy actually taken up
                    if (choiceind(i, t) == 1) then
                        if (shock(2,i,t) < movProbR .AND. currenthouseholdstate(3, i) == 0) then
                            call adj_func_comp(newhouseholdstate(3:4, i), &
                                transhouseholdstate(:,i), achoiceMovR, DchoiceMovR,&
                                cchoiceMovR, chindMR, EVpolMovR,&
                            (/.FALSE., .TRUE./), price, consumption(i,t), &
                            rentalind(i,t), choiceind(i,t), welfare(i,t))
                        else if (shock(2,i,t) < movprob) then
                            call adj_func_comp(newhouseholdstate(3:4, i), &
                                transhouseholdstate(:,i), achoiceexpectMov, DchoiceexpectMov,&
                                cchoiceexpectMov, choiceindicatorexpectMov, EVpolMov,&
                                (/.FALSE., .TRUE./), price, consumption(i,t), &
                                rentalind(i,t), choiceind(i,t), welfare(i,t))
                        else
                        call adj_func_comp(newhouseholdstate(3:4, i), &
                            transhouseholdstate(:,i), achoiceexpect, Dchoiceexpect,&
                            cchoiceexpect, choiceindicatorexpect, EVpol,&
                            (/.FALSE., .TRUE./), price, consumption(i,t), &
                            rentalind(i,t), choiceind(i,t), welfare(i,t))
                        end if
                    end if
                ! %>
                ! EXPERIMENTAL: if the policy is anticipated a period in advance,
                ! agent still compares utility of buying a home now or punting
                ! the decision to the policy period. This comparison is already
                ! solved in the DP problem, so just evaluate the function? %<
                else if (period < PolStart) then
                    !TODO: which EV is this?
                    if (shock(2,i,t) < movProbR .AND. currenthouseholdstate(3, i) == 0) then
                        call pol_linworking(transhouseholdstate(:,i),&
                                            achoiceMovR,DchoiceMovR,cchoiceMovR,EVpolMovR,chindMR,&
                        newhouseholdstate(4,i), newhouseholdstate(3,i),&
                        consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))

                    else if (shock(2,i,t) < movprob) then
                        call pol_linworking(transhouseholdstate(:,i),&
                            achoiceexpectMov,DchoiceexpectMov,cchoiceexpectMov,EVpolMov,choiceindicatorexpectMov,&
                            newhouseholdstate(4,i), newhouseholdstate(3,i),&
                            consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                else
                    call pol_linworking(transhouseholdstate(:,i),&
                        achoiceexpect,Dchoiceexpect,cchoiceexpect,EVpol,choiceindicatorexpect,&
                        newhouseholdstate(4,i), newhouseholdstate(3,i),&
                        consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    end if
                ! %>
                else
                    if (shock(2,i,t) < movProbR .AND. currenthouseholdstate(3, i) == 0) then
                        call pol_linworking(transhouseholdstate(:,i),&
                            achoiceMovR,DchoiceMovR,cchoiceMovR,EVMovR,chindMR,&
                        newhouseholdstate(4,i), newhouseholdstate(3,i),&
                        consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    else if (shock(2,i,t) < movprob .AND. currenthouseholdstate(3, i) > 0) then
                        call pol_linworking(transhouseholdstate(:,i),&
                            achoiceMov,DchoiceMov,cchoiceMov,EVMov,choiceindicatorMov,&
                            newhouseholdstate(4,i), newhouseholdstate(3,i),&
                            consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    else
                    call pol_linworking(transhouseholdstate(:,i),&
                        achoice,Dchoice,cchoice,EV,choiceindicator,&
                        newhouseholdstate(4,i), newhouseholdstate(3,i),&
                        consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                    end if
                end if
                ! %>

                ! Transitory income state + policy transfers, encoded in assets
                if ((currenthouseholdstate(7, i) /= 3) .AND. (choiceind(i,t) == 1) .AND.&
                    household_byperiod(2, i, t+1) == 2 .AND. t <= Tretire) then
                    newhouseholdstate(4, i) = newhouseholdstate(4, i) + &
                        (1.0-eta_transfer)*adjTransfer(currenthouseholdstate(1,i), t)
                    subval(i, t) = adjTransfer(currenthouseholdstate(1,i), t)
                else if ((currenthouseholdstate(7, i) /= 3) .AND. (choiceind(i,t) == 1)) then
                    newhouseholdstate(4, i) = newhouseholdstate(4, i) + &
                        (1.0-eta_transfer)*adjTransfer(currenthouseholdstate(1,i), t)
                    subval(i, t) = adjTransfer(currenthouseholdstate(1,i), t)
                end if

                ! To prevent any case of bad interpolation (but wouldn't
                ! pol_linworking catch this?)
                if (newhouseholdstate(3, i)>Dmax) then
                    newhouseholdstate(3, i)=(1.0-1e-7)*Dmax
                end if
                if (newhouseholdstate(4, i)>amax) then
                    newhouseholdstate(4, i)=(1.0-1e-7)*amax
                end if

                incomeholder(i, t)=(currenthouseholdstate(2, i) - 1.0)*&
                                   income(currenthouseholdstate(1, i), t)
                posttaxincome(i, t)=(currenthouseholdstate(2, i) - 1.0)*&
                                    (income(currenthouseholdstate(1, i), t)-&
                                    inctax(currenthouseholdstate(1, i), t))
                if (t < Tdie) then
                    ! Record lags of consumption, durables, assets
                    household_nextperiod(11, i, t+1)=consumption(i,t)
                    household_nextperiod(12, i, t+1)=financialwealth(i,t)
                    ! Welfare recording
                    household_nextperiod(13, i, t+1)=welfare(i,t)
                end if

                ! %< Income shock allocation
                if (t<=Tretire) then
                    if (currenthouseholdstate(2, i) == 2 .AND. &
                        shock(3,i,t) < empprob(t+1)) then
                        ! Remaining employed
                        newhouseholdstate(2, i) = currenthouseholdstate(2, i)
                    else
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
                    newhouseholdstate(2, i)=household_byperiod(2, i, t+1)
                end if

                ! %>

                ! Generate fixed cost shock
                ! Notice we allow arbitrary shifts for assets in the future
                if (choiceind(i,t) == 1.0) newhouseholdstate(4, i) = newhouseholdstate(4,i) - &
                   F*shock(4,i,t)*newhouseholdstate(3, i)*exp(price)

                ! Decision tree begins with unambiguous renting
                if (rentalind(i, t) > .995) then
                    ! %<
                    rentalflow(i, t, period+2)=newhouseholdstate(3, i)
                    newhouseholdstate(3, i)=0
                    fthb_flag(i, 3) = Tdie+1
                    numrent(i, t)=1
                    if (currenthouseholdstate(3, i) > 0) then
                        resaleflow(i, t, period+2) = currenthouseholdstate(3, i)
                        newownerflow(i, t, period+2) = newhouseholdstate(3, i)
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
                    housingstock(i, t, period+2)=currenthouseholdstate(3, i)
                    numowntotal(i, t)=1
                    ! Any house purchase resets timer
                    fthb_flag(i, 2) = Tdie+1
                ! Track agents as eligible repeat buyer if their last purchase old enough
                    if (t >= fthb_flag(i, 3) + EligYrsR) fthb_flag(i, 3) = -1
                    if (choiceind(i, t) == 1.0) then
                ! %< Start the repeat buyer policy countdown whenever a durables purchase is made
                        numtrans(i, t) = 1
                        housingflow(i, t, period+2)=newhouseholdstate(3, i)
                        housingstock(i, t, period+2)=newhouseholdstate(3, i)
                        newownerflow(i, t, period+2) = newhouseholdstate(3, i)
                        if (currenthouseholdstate(3, i) > 0) then
                            resaleflow(i, t, period+2) = currenthouseholdstate(3, i)
                            numresale(i, t) = 1
                        end if
                    ! Insert policy taker characteristics again
                    ! (The final negative indicators are factors for which policy
                    !  was taken since Matlab doesn't handle mixed types well)
                        if ((fthb_flag(i, 1) == Tdie+1)) then
                            if (write_files) then
                            write(73,'(3I6.2,14F19.6,I6.2)') i, t-period, t, newhouseholdstate(3, i), &
                            newhouseholdstate(4, i) - (1-(1-polLevel(3))*theta)*exp(price)*newhouseholdstate(3,i),&
                            rentalflow(i, max(t-1, 1), period+1), fthb_flag(i, 3), &
                            subval(i,t),financialwealth(i,t),consumption(i,t),welfare(i,t),&
                            household_byperiod(1,i,t), householdresultsholder(10, i, t),&
                            household_byperiod(1,i,max(t-1,1)), householdresultsholder(10, i, max(t-1,1)),&
                            household_byperiod(11:12,i,t),-1
                            end if
                        else if ((fthb_flag(i, 3) == -1) .and. (t <= Tdie)) then
                            if (write_files) then
                            write(73,'(3I6.2,14F19.6,I6.2)') i, t-period, t, newhouseholdstate(3, i), &
                            newhouseholdstate(4, i) - (1-(1-polLevel(3))*theta)*exp(price)*newhouseholdstate(3,i),&
                            currenthouseholdstate(3, i), fthb_flag(i, 1), &
                            subval(i,t),financialwealth(i,t),consumption(i,t),welfare(i,t),&
                            household_byperiod(1,i,t), householdresultsholder(10, i, t),&
                            household_byperiod(1,i,max(t-1,1)), householdresultsholder(10, i, max(t-1,1)),&
                            household_byperiod(11:12,i,t),-2
                            end if
                        end if
                        fthb_flag(i, 1) = t
                        fthb_flag(i, 3) = t
                    ! %>
                    end if
                    ! %>
                end if


                ! Debug output file if I need it
                !if ((period < 2) .or. ((period > 10) .and. (period <= 12))) then
                !    if (write_files) then
                !    write(74,'(3I6.2,6F16.6)') i, t-period, t, newhouseholdstate(3, i), rentalind(i, t), newhouseholdstate(4, i), welfare(i, t), currenthouseholdstate(6, i)
                !    write(74,'(3I6.2,6F16.6)') i, t-period, t, currenthouseholdstate(3, i), rentalind(i, t), max(newhouseholdstate(2, i), rentalflow(i, t, period+2)), &
                    !    incomeholder(i, t), currenthouseholdstate(6, i), fthb_flag(i, 1)
                !    end if
                !end if

                ! %< write(*,*) "here 2"
                ! End of life means liquidation of housing
                if (t==TDie) bequestflow(i, t, period+2)=newhouseholdstate(3, i)
                if (t==TDie) numbequest(i, t) = 1
                ! Degree to which durable can be maintained
                newhouseholdstate(3, i) = (1-delta+maint*delta)*newhouseholdstate(3, i)

                ! Carry statistics of living agents into next period
                if (t < TDie) then 
                ! Point here is to shut off income uncertainty adding noise
                ! to the transition results - reuse from simulation code
                household_nextperiod(1, i, t+1)=household_byperiod(1, i, t+1)
                household_nextperiod(2:4, i, t+1)=newhouseholdstate(2:4, i)
                household_nextperiod(5, i, t+1)=household_byperiod(5, i, t+1)

                ! Track housing decisions at end of period (we want to carry these
                ! stats into the next period!)
                household_nextperiod(6:8, i, t+1)=fthb_flag(i, 1:3)
                household_nextperiod(9, i, t+1)=rentalflow(i, t, period+2)
                household_nextperiod(10, i, t)=incomeholder(i,t)
                end if

                ! %>

            end if ! Checking if person is alive
            end do
        !$OMP END PARALLEL DO

            if (write_files) then
                write(71, '(2I6.2, 6F16.6)') period, t, sum(housingstock(:, t, period+2)), sum(housingflow(:, t, period+2)), &
                    sum(numtrans(:, t)), sum(numrent(:, t)), sum(alive(:, t))
            end if

            end do

            do t=period+1,Tdie
                do i=1,numHH
                numadjust(i, t-period, period+2) = numtrans(i, t)
                end do
            end do

            ! this is the housing price elasticity we want
            !write(*,*) t*1.0, (sum(diffc(:, t))/sum(alive(:, t)))/(sum(consumption(:, t))/sum(alive(:, t)))/(hpmin)
            !write(*,*) "elasticity", sum(diffc(:, t)*alive(:, t))/sum(consumption(:, t))/sum(alive(:, t))/0.2


           ! write(*,*) t*1.0, sum(currenthouseholdstate(2, :)), sum(alive(:, t))

            !write(10, 10) (sum((ratioc(:, t)-1))/sum(alive(:, t)))/(exp(hpmin) - 1), (sum(diffc(:, t))/sum(alive(:, t)))/(sum(consumption(:, t))/sum(alive(:, t)))/(hpmin)
            !write(10, 10) (sum((ratioc(:, t)-1))/numHH)/(exp(hpmin) - 1), (sum(diffc(:, t))/numHH)/(sum(consumption(:, t))/numHH)/(exp(hpmin) - 1), ratiocown/numown_t, ratiocrent/numrent_t, numown_t/(numown_t+numrent_t), numsell_t/numHH, numbuy_t/numHH,(sum(diffc(:, t)/hpmin)/numHH), sum(aggregatecontrib(:, t))/numHH

            !write(555, 555) a_t/numHH, aown_t/numown_t, arent_t/numrent_t, leverage_t/numHH, leverageown_t/numown_t, leveragerent_t/numrent_t

        DEALLOCATE(aggregatecontrib, mpc, &
        consumptionmpc, consumptionhpshock, durableconsumption, &
        currenttemp, currentperm, rentalind, choiceind, &
        durableinvestment, actualwealth, financialwealth, &
        housingnet, housinggross, taxrate, welfare, &
        diffc, ratioc, qchg, numtrans)
        DEALLOCATE(numowntotal,numrent,numresale,numbequest)
        DEALLOCATE(currenthouseholdstate,newhouseholdstate)
        if (partial) DEALLOCATE(alive, posttaxincome)

    end subroutine transition ! %>



    ! TODO: The function is a copy of SimSteady right now. The goal:
    ! Change transition to a copy of simulation that accepts a household
    ! state holder argument. Then a transition path is calculated by taking
    ! a loop over the transaction subroutine. But at each step of the loop
    ! we can just call the Brent minimization routine as well.
    REAL(8) FUNCTION SimTransPath(price, numHH, write_files, skipmodel, skipsim)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    This function is called in every iteration of the algorithm used
        !    to find the RE transition path. It solves the DP problem,
        !    runs the simulation and checks market clearing.
        !
        !    Modified: 05/17/2019
        !
        !    PARAMETERS
        !
        !    price: Real: records the price per unit of housing services.
        !
        !    numHH: Marks the number of agents simulated. A total of
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: price
        INTEGER, INTENT(IN) :: numHH
        REAL(8) :: laborshare, construction, taxshare, netrevenue,&
                   lifeinc, lifecons, psi0, avghousing, totrev
        INTEGER :: ind, timestart, timeend, timeres
        LOGICAL :: pol_at_start = .TRUE.
        LOGICAL, INTENT(IN) :: write_files, skipmodel, skipsim
        
        ind = TransTime+2
        hpdelta = price - hpnodes(ind-1)

        write(*,*) "// Price: ", EXP(price)
        write(0,*) NEW_LINE('A') // "// Price: ", EXP(price)
        ! TODO: The generalization for multiple-period policies here is
        ! certainly wrong: revise.
        call bequests((/price/), ind+1)
        call system_clock(timestart, timeres)

        if (.NOT. skipmodel) then
        if (TransTime < PolEnd) then

        if ((PolEnd == 1) .and. (TransTime == 0)) then
            write(0,*) NEW_LINE('A') // "Begin value function iteration with subsidy:"
            call solveworkingproblem((/price, hpnodes(ind+1)/), (/balancer_internal/),&
                achoiceshock(:,:,:,:,:,ind-1:ind-1), Dchoiceshock(:,:,:,:,:,ind-1:ind-1),&
                cchoiceshock(:,:,:,:,:,ind-1:ind-1), choiceindicatorshock(:,:,:,:,:,ind-1:ind-1),& 
                achoiceshockMov(:,:,:,:,:,ind-1:ind-1),dchoiceshockmov(:,:,:,:,:,ind-1:ind-1),&
                cchoiceshockMov(:,:,:,:,:,ind-1:ind-1),choiceindicatorshockMov(:,:,:,:,:,ind-1:ind-1),&
                achoiceshockMovR(:,:,:,:,:,ind-1:ind-1),dchoiceshockmovR(:,:,:,:,:,ind-1:ind-1),&
                cchoiceshockMovR(:,:,:,:,:,ind-1:ind-1),chindshockMR(:,:,:,:,:,ind-1:ind-1),&
                .FALSE., .FALSE., 'None',EV(:,:,:,:,:,ind+1:ind+1),&
                EVMov(:,:,:,:,:,ind+1:ind+1), EVMovR(:,:,:,:,:,ind+1:ind+1), (/2, TDie/))
        else
            ! Any other case involves anticipating the policy in a future period
            ! so refer to the policy array "*expect"
            EVpol(:,:,:,:,TDie+1,ind+1) = EV(:,:,:,:,Tdie+1,ind+1)
            call solveworkingproblem((/price, hpnodes(ind+1)/), (/balancer_internal/),&
                achoiceshock(:,:,:,:,:,ind-1:ind-1),Dchoiceshock(:,:,:,:,:,ind-1:ind-1),&
                cchoiceshock(:,:,:,:,:,ind-1:ind-1),choiceindicatorshock(:,:,:,:,:,ind-1:ind-1),&
                achoiceshockMov(:,:,:,:,:,ind-1:ind-1),DchoiceshockMov(:,:,:,:,:,ind-1:ind-1),&
                cchoiceshockMov(:,:,:,:,:,ind-1:ind-1),choiceindicatorshockMov(:,:,:,:,:,ind-1:ind-1),&
                achoiceshockMovR(:,:,:,:,:,ind-1:ind-1),dchoiceshockmovR(:,:,:,:,:,ind-1:ind-1),&
                cchoiceshockMovR(:,:,:,:,:,ind-1:ind-1),chindshockMR(:,:,:,:,:,ind-1:ind-1),&
                 .FALSE., .TRUE., 'None',EVpol(:,:,:,:,:,ind+1:ind+1),&
                EVpolMov(:,:,:,:,:,ind+1:ind+1),EVpolMovR(:,:,:,:,:,ind+1:ind+1),(/2, TDie/) )

            write(*,*) "Begin value function iteration with anticipated subsidy:"
            write(0,*) NEW_LINE('A') // "Begin value function iteration with anticipated subsidy:"
            EVpol(:,:,3:Dgridsize,:,:,ind+1) = EV(:,:,3:Dgridsize,:,:,ind+1)
            ! Revert as soon as policy DP calculated (variables hardcoded in other
            ! programs)
            noadjTransfer_internal= 0.0
            if (TransTime < PolStart) pol_at_start = .FALSE.

            call solveworkingproblem((/price, hpnodes(ind+1)/), (/balancer_internal/),&
                achoiceexpect(:,:,:,:,:,ind-1:ind-1),Dchoiceexpect(:,:,:,:,:,ind-1:ind-1),&
                cchoiceexpect(:,:,:,:,:,ind-1:ind-1),choiceindicatorexpect(:,:,:,:,:,ind-1:ind-1),&
                achoiceexpectMov(:,:,:,:,:,ind-1:ind-1),DchoiceexpectMov(:,:,:,:,:,ind-1:ind-1),&
                cchoiceexpectMov(:,:,:,:,:,ind-1:ind-1),choiceindicatorexpectMov(:,:,:,:,:,ind-1:ind-1),&
                achoiceexpectMovR(:,:,:,:,:,ind-1:ind-1),DchoiceexpectMovR(:,:,:,:,:,ind-1:ind-1),&
                cchoiceexpectMovR(:,:,:,:,:,ind-1:ind-1),chindexpectMR(:,:,:,:,:,ind-1:ind-1),&
                 pol_at_start, .FALSE., 'first-time', &
                 EVpol(:,:,:,:,:,ind+1:ind+1), EVpolMov(:,:,:,:,:,ind+1:ind+1), EVpolMovR(:,:,:,:,:,ind+1:ind+1))
        end if

            call plotpolicy_constantwealth(1, achoiceshock, Dchoiceshock, cchoiceshock, &
                                           choiceindicatorshock)
        end if

        ! TODO: Rerunning DP for one transition may be superfluous and
        ! counterproductive.
        !if ((TransTime < PolEnd + 5) .OR. (TransTime == hptransLength - 2)) then
            write(0,*) NEW_LINE('A') // "Begin value function iteration:"
            call solveworkingproblem((/price, hpnodes(ind+1)/), (/balancer_internal/),&
                achoice(:,:,:,:,:,ind:ind), Dchoice(:,:,:,:,:,ind:ind),&
                cchoice(:,:,:,:,:,ind:ind), choiceindicator(:,:,:,:,:,ind:ind), &
                achoiceMov(:,:,:,:,:,ind:ind), DchoiceMov(:,:,:,:,:,ind:ind),&
                cchoiceMov(:,:,:,:,:,ind:ind), choiceindicatorMov(:,:,:,:,:,ind:ind), &
                achoiceMovR(:,:,:,:,:,ind:ind), DchoiceMovR(:,:,:,:,:,ind:ind),&
                cchoiceMovR(:,:,:,:,:,ind:ind),chindMR(:,:,:,:,:,ind-1:ind-1),&
                .TRUE., .FALSE., 'None', EV(:,:,:,:,:,ind+1:ind+1), EVmov(:,:,:,:,:,ind+1:ind+1), EVmovR(:,:,:,:,:,ind+1:ind+1))
            write(0,*) "Value function iteration complete"
        !end if
        end if

        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

        write(0,*) "One-period transition starting..."
        housingstock(:,:, ind) = 0
        rentalflow(:,:, ind) = 0
        resaleflow(:,:, ind) = 0
        newownerflow(:,:, ind) = 0
        bequestflow(:,:, ind) = 0
        numadjust(:,:, ind) = 0
        call transition(price,numHH,TransTime,&
            householdtransitionholder(:,:,:,ind-1),&
            householdtransitionholder(:,:,:,ind), shock, .FALSE., write_files)

        if (.NOT. skipsim) then
        ! TODO: Substitute get_mktclear into this?
        call get_mktclear(price, ind, avghousing, construction, totrev)

        ! Balancer calibrated in every iteration, so no need to throw into
        ! objective fn
        ! This is an ad-hoc "smoothing method": to make sure algorithm isn't
        ! all dependent on balanced budget instead of prices, decisions are made
        ! by putting more weight on the last iterated budget transfer value
        balancer_internal = (0.40*balancer(ind)*SUM(alive)+0.60*totrev)/SUM(alive)
        balancer(ind) = (0.25*balancer(ind)*SUM(alive)+0.75*totrev)/SUM(alive)


        !! Calculate new housing construction - does it
        !! 1) make up for depreciation in both the homeowner and rental
        !! industry housing stock;
        !! 2) Balances the amount of housing entering the market from homeowners
        !! entering rental and dead homeowners?
        !avghousing = (SUM(rentalflow(:,:,ind)) - (1-delta)*SUM(rentalflow(:,:,ind-1))) + &
        !             (SUM(newownerflow(:,:,ind)) - SUM(resaleflow(:,:,ind))) - &
        !              SUM(bequestflow(:,:,ind))
        !avghousing = avghousing + delta*SUM(housingstock(:,:,ind-1))
        !if (avghousing < 0) avghousing = 1e-6
        !if (construction < 0) construction = LOG(1e-6)
        !write(0,'(2A30)') "// Initial housing in t", "  |  Final housing in t"
        !write(0,'(2F30.2)') SUM(housingstock(:,:,ind-1)), SUM(housingstock(:,:,ind))
        !write(0,'(2A30)') "// New owners homes in t", "  |  Final rental homes in t"
        !write(0,'(2F30.2)') SUM(newownerflow(:,:,ind)), SUM(rentalflow(:,:,ind))
        !write(0,'(2A30)') "// Resold housing in t", "  |  Bequested housing in t"
        !write(0,'(2F30.2)') SUM(resaleflow(:,:,ind)), SUM(bequestflow(:,:,ind))
        !write(0,'(A40,A20)') "// Demand for housing less resale supply", "  |  New construction"
        !write(0,'(F40.10,F20.10)') avghousing, EXP(construction)
        SimTransPath = (LOG(avghousing) - construction)**2.0
        else
        SimTransPath = 1.0e6  ! Placeholder
        end if
        DEALLOCATE(posttaxincome)
        

    END FUNCTION SimTransPath ! %>

    SUBROUTINE SimTransLeads
        IMPLICIT NONE
        INTEGER :: t, i
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Subroutine outputs one of the output datasets: A dataset
        !    indexed by agent-variable and contains in each column
        !    the variable's observation from 2 periods before to 9 periods
        !    after the policy. Data are then used to make an "Event study graph."
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        OPEN (UNIT=75, FILE="transition_lagsleads.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

!$OMP PARALLEL DO
        do t=1,Tretire+10 ! Recording future paths of working age agents
            do i=1,numhouseholds,10  ! Taking only a 10% sample of all households
        write(75,'(2I6.2, A9, 12F16.6)') i, t, "C", householdresultsholder(6, i, max(t-4,1)), &
            householdresultsholder(6, i, max(t-3,1)),householdresultsholder(6, i, max(t-2,1)),&
            householdresultsholder(6, i, max(t-1,1)),householdresultsholder(6, i, t),&
            householdtransitionholder(11, i, min(t+1,Tdie), 2), householdtransitionholder(11, i, min(t+2,Tdie), 3),& 
            householdtransitionholder(11, i, min(t+3,Tdie), 4), householdtransitionholder(11, i, min(t+4,Tdie), 5),&
            householdtransitionholder(11, i, min(t+5,Tdie), 6), householdtransitionholder(11, i, min(t+6,Tdie), 7),&
            householdtransitionholder(11, i, min(t+10,Tdie), 11)
        write(75,'(2I6.2, A9, 12F16.6)') i, t, "Q", householdresultsholder(7, i, max(t-4,1)), &
            householdresultsholder(7, i, max(t-3,1)),householdresultsholder(7, i, max(t-2,1)),&
            householdresultsholder(7, i, max(t-1,1)),householdresultsholder(7, i, t),&
            householdtransitionholder(4, i, min(t+1,Tdie), 2), householdtransitionholder(4, i, min(t+2,Tdie), 3),& 
            householdtransitionholder(4, i, min(t+3,Tdie), 4), householdtransitionholder(4, i, min(t+4,Tdie), 5),&
            householdtransitionholder(4, i, min(t+5,Tdie), 6), householdtransitionholder(4, i, min(t+6,Tdie), 7),&
            householdtransitionholder(4, i, min(t+10,Tdie), 11)
        write(75,'(2I6.2, A9, 12F16.6)') i, t, "H_rent", householdresultsholder(8,i,max(t-4,1))*householdresultsholder(9,i,max(t-4,1)),&
            householdresultsholder(8,i,max(t-3,1))*householdresultsholder(9,i,max(t-3,1)), householdresultsholder(8,i,max(t-2,1))*householdresultsholder(9,i,max(t-2,1)),&
            householdresultsholder(8,i,max(t-1,1))*householdresultsholder(9,i,max(t-1,1)), householdresultsholder(8,i,t)*householdresultsholder(9,i,t),&
            householdtransitionholder(9, i, min(t+1,Tdie), 2), householdtransitionholder(9, i, min(t+2,Tdie), 3),& 
            householdtransitionholder(9, i, min(t+3,Tdie), 4), householdtransitionholder(9, i, min(t+4,Tdie), 5),&
            householdtransitionholder(9, i, min(t+5,Tdie), 6), householdtransitionholder(9, i, min(t+6,Tdie), 7),&
            householdtransitionholder(9, i, min(t+10,Tdie), 11)
        write(75,'(2I6.2, A9, 12F16.6)') i, t, "H_own", householdresultsholder(8,i,max(t-4,1))*(1.0-householdresultsholder(9,i,max(t-4,1))),&
            householdresultsholder(8,i,max(t-3,1))*(1.0-householdresultsholder(9,i,max(t-3,1))), householdresultsholder(8,i,max(t-2,1))*(1.0-householdresultsholder(9,i,max(t-2,1))),&
            householdresultsholder(8,i,max(t-1,1))*(1.0-householdresultsholder(9,i,max(t-1,1))), householdresultsholder(8,i,t)*(1.0-householdresultsholder(9,i,t)),&
            householdtransitionholder(3, i, min(t+1,Tdie), 2), householdtransitionholder(3, i, min(t+2,Tdie), 3),& 
            householdtransitionholder(3, i, min(t+3,Tdie), 4), householdtransitionholder(3, i, min(t+4,Tdie), 5),&
            householdtransitionholder(3, i, min(t+5,Tdie), 6), householdtransitionholder(3, i, min(t+6,Tdie), 7),&
            householdtransitionholder(3, i, min(t+10,Tdie), 11)
        write(75,'(2I6.2, A9, 12F16.6)') i, t, "Y", householdresultsholder(1,i,max(t-4,1)),&
            householdresultsholder(1,i,max(t-3,1)), householdresultsholder(1,i,max(t-2,1)),&
            householdresultsholder(1,i,max(t-1,1)), householdresultsholder(1,i,t),&
            householdtransitionholder(1, i, t, 1), householdtransitionholder(1, i, min(t+1,Tdie), 2),& 
            householdtransitionholder(1, i, min(t+2,Tdie), 3), householdtransitionholder(1, i, min(t+3,Tdie), 4),&
            householdtransitionholder(1, i, min(t+5,Tdie), 5),householdtransitionholder(1, i, min(t+6,Tdie), 6),&
            householdtransitionholder(1, i, min(t+10,Tdie), 10)
        write(75,'(2I6.2, A9, 12F16.6)') i, t, "V", householdresultsholder(11,i,max(t-4,1)),&
            householdresultsholder(11,i,max(t-3,1)), householdresultsholder(11,i,max(t-2,1)),&
            householdresultsholder(11,i,max(t-1,1)), householdresultsholder(11,i,t),&
            householdtransitionholder(13, i, min(t+1,Tdie), 2), householdtransitionholder(13, i, min(t+2,Tdie), 3),& 
            householdtransitionholder(13, i, min(t+3,Tdie), 4), householdtransitionholder(13, i, min(t+4,Tdie), 5),&
            householdtransitionholder(13, i, min(t+5,Tdie), 6), householdtransitionholder(13, i, min(t+6,Tdie), 7),&
            householdtransitionholder(13, i, min(t+10,Tdie), 11)
            end do
        end do
!$OMP END PARALLEL DO
        close(75)

    END SUBROUTINE SimTransLeads ! %>

end module lifecycle_transition
