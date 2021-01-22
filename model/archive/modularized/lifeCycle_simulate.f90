module lifecycle_simulate
    USE lifecycle_algs
    USE share
    USE OMP_LIB
    IMPLICIT NONE

    REAL(8), dimension(numhouseholds, 5) :: currenthouseholdstate, newhouseholdstate
    REAL(8), dimension(:,:), ALLOCATABLE :: aggregatecontrib, mpc, consumption, consumptionmpc, consumptionhpshock, durableconsumption, currenttemp, currentperm, incomeholder, rentalind, choiceind, durableinvestment, actualwealth, financialwealth, housingnet, housinggross, actualwealthdividedbyincome, totalnetworthdividedbyincome, taxrate, welfare, diffc, ratioc, qchg
    REAL(8), dimension(:,:), ALLOCATABLE :: liquidassetsminusmortgage, changec, changetemp, changeperm, changey, changeyconditional, changed, changedconditional, changetempconditional, changepermconditional, demeanedindicator, adjustmean
    REAL(8), dimension(:,:), ALLOCATABLE :: alive
    REAL(8), dimension(:,:), ALLOCATABLE :: housingstock, housingflow, rentalflow, numrent, numowntotal
    REAL(8), dimension(numhouseholds, 2) :: nearestnode
    REAL(8) :: shock
    integer, dimension(12) :: seedvalue
    REAL(8), dimension(numhouseholds) :: insurancetemp, insuranceperm, cov1, cov2, cov3, cov4, insurancedtemp, insurancedperm, insurancedtempconditional, insurancedpermconditional, numtrans, adjustt
    REAL(8) :: adjust
    REAL(8) :: conditionalindicator
    REAL(8) :: ax, bx, cx, tol
    REAL(8) :: cov1est, cov2est, cov3est, cov4est, minest, constrained
    REAL(8) :: m_eps2, m_eta2, m_cFy, m_yFy, m_epsc, m_etac, m_cIVy, m_yIVy, m_eps, m_eta, m_c, m_y, m_IVy, m_Fy, m_d, m_dFy, m_dIVy, m_epsd, m_etad
    REAL(8) :: pretaxcalcholder, pretaxincome
    REAL(8) :: rep
    REAL(8) :: actualwealthtotal
    REAL(8) :: exanteEVborn, exanteEVoverall, numobswelfare
    REAL(8), dimension(hpgridsize, 1) :: exanteEVbornstate, numobsstate, overallwelfarebystate, overallobsbystate, consumptionbystate, overallwelfarebystateyoung, overallwelfarebystatemiddle, overallwelfarebystateold, overallobsbystateyoung, overallobsbystatemiddle, overallobsbystateold, consumptionbystateyoung, consumptionbystatemiddle, consumptionbystateold
    REAL(8), dimension(hpgridsize, 350) :: overallwelfarebystateCE, overallwelfarebystateyoungCE, overallwelfarebystatemiddleCE, overallwelfarebystateoldCE
    REAL(8), dimension(hpgridsize) :: cutoff
    REAL(8), dimension(350, 1) :: CE
    REAL(8) :: medianincome, ratiocrent, ratiocown, diffcrent, diffcown, numrent_t, numown_t, numbuy_t, numsell_t, fracown
    REAL(8), dimension(8):: numbins ! Initial income quantiles X (Owner, renter)
    REAL(8) :: leverage_t, a_t, leverageown_t, aown_t, leveragerent_t, arent_t
    real(8), dimension(numhouseholds, 3) :: fthb_flag
    ! 1: age when purchase fthb, 2: last period when house sold,
    ! 3: last period when house bought

    CONTAINS

    subroutine simulate(numhouseholds)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: numhouseholds

        do i=1, 350
        CE(i, 1)=.97+(1.0*i)/5000.0
        !write(*,*) CE(i, 1)
        end do

        tol=1.0e-9

        OPEN (UNIT=22, FILE="fthb.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        22 format (I6.2, I6.2, F16.6, F16.6, F16.6, F16.6)

        OPEN (UNIT=1989, FILE="dist_fthb.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        1989 format (I6.2, I6.2, F16.6, F16.6, F16.6, F16.6, I6.2)

        ! Orphaned as of now?
        !OPEN (UNIT=3, FILE="agecoefficients.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        !3 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)
        !
        !4 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

        OPEN (UNIT=7, FILE="lifecycle_adjust.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        7 format (I6.2, I6.2, F16.6)

        OPEN (UNIT=5, FILE="lifecycleprofiles.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        5 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

        OPEN (UNIT=9, FILE="householdresults.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        9 format (I16.6, I6.2, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

        OPEN (UNIT=99, FILE="housingstock.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        99 format (I6.2, F16.6, F16.6, F16.6)

        OPEN (UNIT=555, FILE="leverageprofiles.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        555 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

        OPEN (UNIT=10, FILE="elasticity.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        10  format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)


        OPEN (UNIT=11, FILE="statesandmpc.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        11  format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

        seedvalue(:)=1
        CALL random_seed(put=seedvalue(:))

        ! STATE DEFINITION
        ! state 1: liquid assets a
        ! state 2: durable assets d
        ! state 3: idiosyncratic income shock: z
        ! state 4: house price
        ! state 5: age: t

        do i=1, numhouseholds
            if (real(i)/numhouseholds<Probinitcum(1)) then
                currenthouseholdstate(i, 3)=1
            else
                do j=1, zgridsize-1
                    if (real(i)/numhouseholds > Probinitcum(j) .AND. real(i)/numhouseholds<=Probinitcum(j+1)) then
                        currenthouseholdstate(i, 3)=j+1
                    end if
                end do
            end if
        end do

        fthb_flag = Tdie + 1 ! flag for when purchasing a house (Tdie+1 = not buying a house in lifetime)

        medianincome=income((zgridsize+1)/2, 1)
        numbins=0

        ! These initial values are from the create_initial_conditions Matlab files
        ! Each if block is a quantile of the income distribution, so quartiles for now.
        ! The shock variable determines if agent is house owner. Else agent rents.
        ! Deprecated while playing with common initial asset levels
        do i=1, numhouseholds
            call random_number(shock)
            if (income(currenthouseholdstate(i, 3), 1)/medianincome<.5017) then
                if (shock < .0524) then
                    fthb_flag(i, 1) = 0
                    currenthouseholdstate(i, 2)=1.568
                    currenthouseholdstate(i, 1)=0.841
                    numbins(1)=numbins(1)+1
                else
                    currenthouseholdstate(i, 2)=0
                    currenthouseholdstate(i, 1)=0.0017
                    numbins(2)=numbins(2)+1
                end if
            elseif (income(currenthouseholdstate(i, 3), 1)/medianincome>=.5017 .and. income(currenthouseholdstate(i, 3), 1)/medianincome< 1) then
                if (shock<.0834) then
                    fthb_flag(i, 1) = 0
                    currenthouseholdstate(i, 2)=0.8650
                    currenthouseholdstate(i, 1)=0.3267
                    numbins(3)=numbins(3)+1
                else
                    currenthouseholdstate(i, 2)=0
                    currenthouseholdstate(i, 1)=0.0015
                    numbins(4)=numbins(4)+1
                end if
            elseif (income(currenthouseholdstate(i, 3), 1)/medianincome>= 1 .and. income(currenthouseholdstate(i, 3), 1)/medianincome<1.6931) then
                if (shock<.1209) then
                    fthb_flag(i, 1) = 0
                    currenthouseholdstate(i, 2)=1.8136
                    currenthouseholdstate(i, 1)=0.0041
                    numbins(5)=numbins(5)+1
                else
                    currenthouseholdstate(i, 2)=0
                    currenthouseholdstate(i, 1)=0.007
                    numbins(6)=numbins(6)+1
                end if
            else
                if (shock<.2760) then
                    fthb_flag(i, 1) = 0
                    currenthouseholdstate(i, 2)=1.764
                    currenthouseholdstate(i, 1)=0.012
                    numbins(7)=numbins(7)+1
                else
                    currenthouseholdstate(i, 2)=0
                    currenthouseholdstate(i, 1)=0.0215
                    numbins(8)=numbins(8)+1
                end if
            end if
        end do


        write(*,*) "Distribution of owners/renters by ", size(numbins)/2, "quantiles:"
        write(*,*) numbins/numhouseholds
        write(*,*) sum(numbins)

        currenthouseholdstate(:, 4)=hpgridsize
        newhouseholdstate(:, 4)=hpgridsize

        ! If you wanted to play with changing the initial distribution of wealth

        !currenthouseholdstate(:, 1)=0  ! start with zero assets (could do something here with bequests)
        !currenthouseholdstate(:, 2)=0


        call random_number(shock)


        ! Allocation of automatic arrays for OpenMP (skip this)
        ALLOCATE(aggregatecontrib(numhouseholds, Tdie), mpc(numhouseholds, Tdie), consumption(numhouseholds, Tdie), &
        consumptionmpc(numhouseholds, Tdie), consumptionhpshock(numhouseholds, Tdie), durableconsumption(numhouseholds, Tdie), &
        currenttemp(numhouseholds, Tdie), currentperm(numhouseholds, Tdie), incomeholder(numhouseholds, Tdie), rentalind(numhouseholds, Tdie), &
        choiceind(numhouseholds, Tdie), durableinvestment(numhouseholds, Tdie), actualwealth(numhouseholds, Tdie), financialwealth(numhouseholds, Tdie), &
        housingnet(numhouseholds, Tdie), housinggross(numhouseholds, Tdie), actualwealthdividedbyincome(numhouseholds, Tdie), &
        totalnetworthdividedbyincome(numhouseholds, Tdie), taxrate(numhouseholds, Tdie), welfare(numhouseholds, Tdie), &
        diffc(numhouseholds, Tdie), ratioc(numhouseholds, Tdie), qchg(numhouseholds, Tdie))

        ALLOCATE(liquidassetsminusmortgage(numhouseholds, Tdie), changec(numhouseholds, Tdie), changetemp(numhouseholds, Tdie), &
        changeperm(numhouseholds, Tdie), changey(numhouseholds, Tdie), changeyconditional(numhouseholds, Tdie), &
        changed(numhouseholds, Tdie), changedconditional(numhouseholds, Tdie), changetempconditional(numhouseholds, Tdie), &
        changepermconditional(numhouseholds, Tdie), demeanedindicator(numhouseholds, Tdie), adjustmean(numhouseholds, Tdie))

        ALLOCATE(alive(numhouseholds, Tdie), housingstock(numhouseholds, Tdie), housingflow(numhouseholds, Tdie), &
        rentalflow(numhouseholds, Tdie), numowntotal(numhouseholds, Tdie), numrent(numhouseholds, Tdie))

        numowntotal=0
        numrent=0
        actualwealthtotal=0
        numtrans=0
        adjustt=0

        housingstock=0
        housingflow=0
        rentalflow=0
        alive=1


        write(*,*) "Simulation start"
        do t=1, Tdie
        write(*,*) "Cohort of age", t

        ratiocrent=0
        ratiocown=0
        diffcrent=0
        diffcown=0
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


        !write(*,*) "t", maxval((1+r)*(currenthouseholdstate(:, 1)+theta*Currenthouseholdstate(:, 2))-delta*currenthouseholdstate(:, 2))
        !write(*,*) "maxq", maxval(currenthouseholdstate(:, 1)), maxval(currenthouseholdstate(:, 2))


            ! 5th state is age (or equiv) time
            currenthouseholdstate(:, 5)=t
        !$OMP PARALLEL DO
            do i=1, numhouseholds
                if (alive(i, t)==1) then

                ! first thing is that if you are old, do you die or not
                if (t>Tretire) then
                    if (alive(i, 1)==1) then
                        call random_number(shock)
                        if (shock<deathprob(t-Tretire)) then
                            alive(i, t:Tdie)=0
                        end if
                    end if
                end if




                actualwealth(i, t)=currenthouseholdstate(i, 1)+theta*currenthouseholdstate(i, 2)*exp(hpnodes(currenthouseholdstate(i, 4)))  !total wealth is voluntary equity plus equity in durable
                financialwealth(i, t)=currenthouseholdstate(i, 1)-(1-theta)*currenthouseholdstate(i, 2)*exp(hpnodes(currenthouseholdstate(i, 4)))
                housingnet(i, t)=theta*currenthouseholdstate(i, 2)*exp(hpnodes(currenthouseholdstate(i, 4)))
                housinggross(i, t)=currenthouseholdstate(i, 2)*exp(hpnodes(currenthouseholdstate(i, 4)))
                !qchg(i, t) = financialwealth(i, t) + (1-theta)*currenthouseholdstate(i, 2)*exp(hpnodes(currenthouseholdstate(i, 4)-1))
                qchg(i, t) = currenthouseholdstate(i, 1) + (1-theta)*currenthouseholdstate(i, 2)*(exp(hpnodes(currenthouseholdstate(i, 4)))-exp(hpnodes(currenthouseholdstate(i, 4)+1)))

                !write(*,*) "here"



                !call pol_linworking(qchg(i, t), currenthouseholdstate(i, 2), currenthouseholdstate(i, 3), currenthouseholdstate(i, 4), currenthouseholdstate(i, 5), newhouseholdstate(i, 1), newhouseholdstate(i, 2), consumption(i, t), rentalind(i, t), welfare(i, t))
                !call pol_linworking(currenthouseholdstate(i, 1)-.01, currenthouseholdstate(i, 2), currenthouseholdstate(i, 3), currenthouseholdstate(i, 4), currenthouseholdstate(i, 5), newhouseholdstate(i, 1), newhouseholdstate(i, 2), consumptionmpc(i, t), rentalind(i, t), welfare(i, t))
                call pol_linworking(currenthouseholdstate(i, 1), currenthouseholdstate(i, 2), currenthouseholdstate(i, 3), currenthouseholdstate(i, 4), currenthouseholdstate(i, 5), newhouseholdstate(i, 1), newhouseholdstate(i, 2), consumption(i, t), rentalind(i, t), choiceind(i, t), welfare(i, t))
                mpc(i, t)=(consumptionmpc(i, t)-consumption(i, t))/-.01
                aggregatecontrib(i, t)=mpc(i, t)*currenthouseholdstate(i, 2)*(1-delta)*(1-F)


                ! Track housing decisions going into period
                householdresultsholder(i, t, 10)=fthb_flag(i, 1)
                householdresultsholder(i, t, 11)=fthb_flag(i, 2)
                householdresultsholder(i, t, 12)=fthb_flag(i, 3)
                householdresultsholder(i, t, 13)=rentalflow(i, max(t-1,1))

                if (i<=numhouseholds) then
                    write(22, 22) i, t, newhouseholdstate(i, 2), currenthouseholdstate(i, 2), rentalind(i, t), currenthouseholdstate(i, 1)
                end if

                if (newhouseholdstate(i, 1)<0.0001 .and. t==1) then
                    constrained=constrained+1
                end if

                ! Decision tree begins with unambiguous renting
                if ((rentalind(i, t) > .995) .or.  (rentalind(i, t) /= rentalind(i, t))) then
                    rentalflow(i, t)=newhouseholdstate(i, 2)
                    newhouseholdstate(i, 2)=0
                    fthb_flag(i, 3) = Tdie+1
                ! Starts the FTHB policy timer.
                    if ((fthb_flag(i, 1) /= Tdie+1) .and. (fthb_flag(i, 2) == Tdie+1)) then
                        fthb_flag(i, 2) = t
                    end if
                ! Track agents again as FTHBs after three consecutive years of renting
                    if ((t >= fthb_flag(i, 2) + (EligYrsF - 1)) .and. (t <= Tretire)) then
                        fthb_flag(i, 1) = Tdie+1
                        fthb_flag(i, 2) = Tdie+1
                    end if
                else
                    housingstock(i, t)=currenthouseholdstate(i, 2)
                    fthb_flag(i, 2) = Tdie+1
                ! Track agents as eligible repeat buyer if their last purchase old enough
                        if ((t >= fthb_flag(i, 3) + EligYrsR) .and. (t <= Tretire)) then
                            fthb_flag(i, 3) = 0
                        end if
                ! Start the repeat buyer policy countdown whenever a durables purchase is made
                        if (choiceind(i, t) == 1.0) then
                            numtrans(i) = 1
                            adjustt(i) = adjustt(i) + 1
                            housingflow(i, t)=newhouseholdstate(i, 2)
                            if (currenthouseholdstate(i, 2) > 0) then
                                housingstock(i, t)=newhouseholdstate(i, 2)
                            end if
                ! Record stats for agent if they are eligible for policy
                ! (The final negative indicators are factors for which policy
                !  was taken since Matlab doesn't handle mixed types well)
                            if (fthb_flag(i, 1) == Tdie+1) then
                                fthb_flag(i, 1) = t
                                write(1989, 1989) i, t, newhouseholdstate(i, 2), &
                                rentalflow(i, max(t-1,1)), 0.0, &
                                currenthouseholdstate(i, 2) + income&
                                &(currenthouseholdstate(i, 3), t), -1
                            else if (fthb_flag(i, 3) == 0) then
                                fthb_flag(i, 1) = t
                                housingstock(i, t)=newhouseholdstate(i, 2)
                                write(1989, 1989) i, t, newhouseholdstate(i, 2), &
                                currenthouseholdstate(i, 2), &
                                (currenthouseholdstate(i, 1) - currenthouseholdstate&
                                &(i, 2)*(1-theta))/currenthouseholdstate(i, 2), &
                                currenthouseholdstate(i, 2) + income&
                                &(currenthouseholdstate(i, 3), t), -2
                                ! Stop tracking repeat buyers after first repeat buy
                                fthb_flag(i, 3) = 99
                            end if
                            fthb_flag(i, 3) = t
                        end if
                end if

                if (newhouseholdstate(i, 1)>amax) then
                    newhouseholdstate(i, 1)=.9999999*amax
                end if
                diffc(i, t)=consumptionhpshock(i, t)-consumption(i, t)
                ratioc(i, t) = consumptionhpshock(i, t)/consumption(i, t)

                if (t<=Tretire) then
                incomeholder(i, t)=income(currenthouseholdstate(i, 3), t)
                else
                incomeholder(i, t)=retirementincome(currenthouseholdstate(i, 3))
                end if

                !write(*,*) "here 2"
                householdresultsholder(i, t, 1)=currenthouseholdstate(i, 1)
                householdresultsholder(i, t, 2)=currenthouseholdstate(i, 2)
                householdresultsholder(i, t, 3)=(currenthouseholdstate(i, 3))
                householdresultsholder(i, t, 4)=alive(i, t)
                householdresultsholder(i, t, 5)=consumption(i, t)
                householdresultsholder(i, t, 6)=newhouseholdstate(i, 1)
                householdresultsholder(i, t, 7)=newhouseholdstate(i, 2)
                householdresultsholder(i, t, 8)=rentalind(i, t)
                householdresultsholder(i, t, 9)=incomeholder(i, t)


                durableconsumption(i, t)=newhouseholdstate(i, 2)+.0000001

                !if (newhouseholdstate(i, 2)-currenthouseholdstate(i, 2)>.001 .and. rentalind(i, t)<=.995) then
                !numbuy_t=numbuy_t+1
                !elseif (newhouseholdstate(i, 2)-currenthouseholdstate(i, 2)<-.001) then
                !numsell_t=numsell_t+1
                !end if

                !a_t=a_t+currenthouseholdstate(i, 1)-(1-theta)*currenthouseholdstate(i, 2)
                !if (currenthouseholdstate(i, 2)/(currenthouseholdstate(i, 1)+theta*currenthouseholdstate(i, 2)) .ne. 0/0) then
                !leverage_t=leverage_t+currenthouseholdstate(i, 2)/(currenthouseholdstate(i, 1)+theta*currenthouseholdstate(i, 2))
                !end if


                if (currenthouseholdstate(i, 2)==0) then
                    numrent(i, t)=1
                !    numrent_t=numrent_t+1
                !    ratiocrent=ratiocrent+(ratioc(i, t)-1)/(exp(hpmin)-1)

                !    arent_t=arent_t+currenthouseholdstate(i, 1)-(1-theta)*currenthouseholdstate(i, 2)
                    !if (currenthouseholdstate(i, 2)/(currenthouseholdstate(i, 1)+theta*currenthouseholdstate(i, 2)) .ne. 0/0) then
                !        leveragerent_t=leveragerent_t+currenthouseholdstate(i, 2)/(currenthouseholdstate(i, 1)+theta*currenthouseholdstate(i, 2))
                    !end if


                else
                    numowntotal(i, t)=1
                !    ratiocown=ratiocown+(ratioc(i, t)-1)/(exp(hpmin)-1)

                    !aown_t=aown_t+currenthouseholdstate(i, 1)-(1-theta)*currenthouseholdstate(i, 2)
                    !if (currenthouseholdstate(i, 2)/(currenthouseholdstate(i, 1)+theta*currenthouseholdstate(i, 2)) .ne. 0/0) then
                    !    leverageown_t=leverageown_t+currenthouseholdstate(i, 2)/(currenthouseholdstate(i, 1)+theta*currenthouseholdstate(i, 2))
                    !end if
                end if


                durableinvestment(i, t)=newhouseholdstate(i, 2)-(1-delta)*currenthouseholdstate(i, 2)

                if (t<=Tretire) then
                    call random_number(shock)

                    if (shock < Probzcum(currenthouseholdstate(i, 3), 1)) then
                            newhouseholdstate(i, 3)=1
                    else
                        do j=1, zgridsize-1
                            if (shock > Probzcum(currenthouseholdstate(i, 3), j) .AND. shock<=Probzcum(currenthouseholdstate(i, 3), j+1)) then
                                newhouseholdstate(i, 3)=j+1
                            end if
                        end do
                    end if
                else
                    newhouseholdstate(i, 3)=currenthouseholdstate(i, 3)
                end if

                if (t==3) then
                    write(11, 11) mpc(i, t), currenthouseholdstate(i, 1), currenthouseholdstate(i, 2), currenthouseholdstate(i, 3), currenthouseholdstate(i, 4), aggregatecontrib(i, t), rentalind(i, t)
                end if

                newhouseholdstate(i, 5)= t+1
                end if

            end do
        !$OMP END PARALLEL DO

            write(5, 5) t*1.0, sum(alive(:, t)), sum(alive(:, t)*rentalind(:, t)), sum(alive(:, t)*(1-rentalind(:, t))), &
            sum(consumption(:, t))/sum(alive(:, t)), sum(consumption(:, t)*(1-rentalind(:, t)))/sum(alive(:, t)*(1-rentalind(:, t))), &
            sum(currenthouseholdstate(:, 2)*(1-rentalind(:, t)))/sum(alive(:, t)*(1-rentalind(:, t))), &
            sum(currenthouseholdstate(:, 1))/sum(alive(:, t)), sum(currenthouseholdstate(:, 1)*(1-rentalind(:, t)))/sum(alive(:, t)*(1-rentalind(:, t))), &
            sum(currenthouseholdstate(:, 1)-(1-theta)*currenthouseholdstate(:, 2))/sum(alive(:, t)), &
            sum((currenthouseholdstate(:, 1)-(1-theta)*currenthouseholdstate(:, 2))*(1-rentalind(:, t)))/sum(alive(:, t)*(1-rentalind(:, t))), &
            sum((currenthouseholdstate(:, 2)/(currenthouseholdstate(:, 2)+currenthouseholdstate(:, 1)-(1-theta)*currenthouseholdstate(:, 2)))*(1-rentalind(:, t)))/sum(alive(:, t)*(1-rentalind(:, t))), &
            1-(sum(rentalind(:, t))/sum(alive(:, t)))


            !write(5, 5) sum(financialwealth(:, t)*alive(:, t))/numhouseholds, sum(currenthouseholdstate(:, 2)*alive(:, t))/numhouseholds, sum(consumption(:, t)*alive(:, t))/numhouseholds, sum(actualwealth(:, t)*alive(:, t))/numhouseholds, sum(newhouseholdstate(:, 2)*alive(:, t))/numhouseholds, sum(newhouseholdstate(:, 2)*rentalind(:, t)*alive(:, t))/numhouseholds, sum(rentalind(:, t)*alive(:, t))/numhouseholds, sum(alive(:, t))/numhouseholds, sum(consumption(:, t))/numhouseholds
             currenthouseholdstate=newhouseholdstate


            ! this is the housing price elasticity we want
            write(*,*) t*1.0, (sum(diffc(:, t))/sum(alive(:, t)))/(sum(consumption(:, t))/sum(alive(:, t)))/(hpmin)
            !write(*,*) "elasticity", sum(diffc(:, t)*alive(:, t))/sum(consumption(:, t))/sum(alive(:, t))/0.2


           ! write(*,*) t*1.0, sum(currenthouseholdstate(:, 2)), sum(alive(:, t))

            !write(10, 10) (sum((ratioc(:, t)-1))/sum(alive(:, t)))/(exp(hpmin) - 1), (sum(diffc(:, t))/sum(alive(:, t)))/(sum(consumption(:, t))/sum(alive(:, t)))/(hpmin)
            !write(10, 10) (sum((ratioc(:, t)-1))/numhouseholds)/(exp(hpmin) - 1), (sum(diffc(:, t))/numhouseholds)/(sum(consumption(:, t))/numhouseholds)/(exp(hpmin) - 1), ratiocown/numown_t, ratiocrent/numrent_t, numown_t/(numown_t+numrent_t), numsell_t/numhouseholds, numbuy_t/numhouseholds,(sum(diffc(:, t)/hpmin)/numhouseholds), sum(aggregatecontrib(:, t))/numhouseholds
            write(99, 99) t, sum(housingstock(:, t)), sum(housingflow(:, t)), sum(numtrans(:))
            numtrans = 0

            !write(555, 555) a_t/numhouseholds, aown_t/numown_t, arent_t/numrent_t, leverage_t/numhouseholds, leverageown_t/numown_t, leveragerent_t/numrent_t
        !
        end do

        !$OMP PARALLEL DO
        do i=1, numhouseholds
            do t=1, 60
              write(9, 9) i, t, householdresultsholder(i, t, 1),  householdresultsholder(i, t, 2),  householdresultsholder(i, t, 3),  householdresultsholder(i, t, 4),  householdresultsholder(i, t, 5),  householdresultsholder(i, t, 6),  householdresultsholder(i, t, 7),  householdresultsholder(i, t, 8),   householdresultsholder(i, t, 9)
            end do

            adjustmean(i,:) = sum(incomeholder(i,:))/Tdie ! Does this mean anything?
            incomeholder(i,:) = (incomeholder(i,:) - adjustmean(i,:))**2
            write(7, 7) i, int(adjustt(i)), sum(incomeholder(i,:))/(Tdie-1)
        end do
        !$OMP END PARALLEL DO
        write(*,*) "end"

        DEALLOCATE(aggregatecontrib, mpc, consumption, &
        consumptionmpc, consumptionhpshock, durableconsumption, &
        currenttemp, currentperm, incomeholder, rentalind, &
        choiceind, durableinvestment, actualwealth, financialwealth, &
        housingnet, housinggross, actualwealthdividedbyincome, &
        totalnetworthdividedbyincome, taxrate, welfare, &
        diffc, ratioc, qchg)

        DEALLOCATE(liquidassetsminusmortgage, changec, changetemp, &
        changeperm, changey, changeyconditional, &
        changed, changedconditional, changetempconditional, &
        changepermconditional, demeanedindicator, adjustmean)

        DEALLOCATE(alive, housingflow, &
        rentalflow, numowntotal, numrent)


        close(22)
        close(5)
        close(7)
        close(9)
        close(99)
        close(1989)
        close(10)

    end subroutine simulate

    REAL(8) FUNCTION SimSteady(price, numHH)
        USE lifecycle_solveDP
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: price
        INTEGER, INTENT(IN) :: numHH
        REAL(8) :: TotHousing, psi2
        
        if (price <= 0) then
            SimSteady = 1e10
        else
        psi2 = ConsElas/(1+ConsElas)
        
        call solveworkingproblem(.TRUE., price)
        call output_vfuncs(.TRUE., price)
        call simulate(numHH)
        TotHousing = SUM(housingstock)/(Tdie*numHH)
        DEALLOCATE(housingstock)

        ! Calculate new housing construction - does it make up for depreciation?
        SimSteady = (ABS(delta*TotHousing - psi0*(psi2*psi0*price)&
                    &**((-psi2)/(psi2-1.0))))
        end if
        write(*,*) SimSteady, price
        
    END FUNCTION SimSteady

end module lifecycle_simulate
