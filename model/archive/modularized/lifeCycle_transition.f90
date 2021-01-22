module lifecycle_transition
    USE lifecycle_algs
    USE share
    USE OMP_LIB
    IMPLICIT NONE

    REAL(8), dimension(numhouseholds, 6) :: currenthouseholdstate, newhouseholdstate
    ! the 6th is the flag for fthb
    REAL(8), dimension(:,:), ALLOCATABLE :: aggregatecontrib, mpc, consumption, consumptionmpc, consumptionhpshock, durableconsumption, currenttemp, currentperm, incomeholder, rentalind, choiceind, durableinvestment, actualwealth, financialwealth, housingnet, housinggross, actualwealthdividedbyincome, totalnetworthdividedbyincome, taxrate, welfare, diffc, ratioc, qchg
    REAL(8), dimension(:,:), ALLOCATABLE :: liquidassetsminusmortgage, changec, changetemp, changeperm, changey, changeyconditional, changed, changedconditional, changetempconditional, changepermconditional, demeanedindicator, adjustt, adjustmean
    REAL(8), dimension(:,:), ALLOCATABLE :: alive
    REAL(8), dimension(:,:,:), ALLOCATABLE :: housingstock, housingflow, rentalflow
    REAL(8), dimension(numhouseholds, 2) :: nearestnode
    REAL(8) :: shock
    integer, dimension(12) :: seedvalue
    integer :: tj
    REAL(8), dimension(numhouseholds) :: insurancetemp, insuranceperm, cov1, cov2, cov3, cov4, &
    insurancedtemp, insurancedperm, insurancedtempconditional, insurancedpermconditional, numtrans
    REAL(8), dimension(numhouseholds) :: thetastate, hpdelta_i
    REAL(8) :: adjust
    REAL(8) :: conditionalindicator
    REAL(8) :: ax, bx, cx, tol
    REAL(8) :: cov1est, cov2est, cov3est, cov4est, minest, constrained
    REAL(8) :: m_eps2, m_eta2, m_cFy, m_yFy, m_epsc, m_etac, m_cIVy, m_yIVy, m_eps, &
    m_eta, m_c, m_y, m_IVy, m_Fy, m_d, m_dFy, m_dIVy, m_epsd, m_etad
    REAL(8) :: pretaxcalcholder, pretaxincome
    REAL(8) :: rep
    REAL(8) :: actualwealthtotal
    REAL(8) :: exanteEVborn, exanteEVoverall, numobswelfare
    REAL(8), dimension(hpgridsize, 1) :: exanteEVbornstate, numobsstate, &
    overallwelfarebystate, overallobsbystate, consumptionbystate, &
    overallwelfarebystateyoung, overallwelfarebystatemiddle, overallwelfarebystateold, overallobsbystateyoung, overallobsbystatemiddle, overallobsbystateold, &
    consumptionbystateyoung, consumptionbystatemiddle, consumptionbystateold
    REAL(8), dimension(hpgridsize, 350) :: overallwelfarebystateCE, &
    overallwelfarebystateyoungCE, overallwelfarebystatemiddleCE, overallwelfarebystateoldCE
    REAL(8), dimension(hpgridsize) :: cutoff
    REAL(8), dimension(350, 1) :: CE
    REAL(8) :: medianincome, ratiocrent, ratiocown, diffcrent, diffcown, numrent_t, numbuy_t, numsell_t
    REAL(8), dimension(10):: numbins
    REAL(8) :: leverage_t, a_t, leverageown_t, aown_t, leveragerent_t, arent_t
    real(8), dimension(numhouseholds, 3) :: fthb_flag

    CONTAINS

    subroutine transition
        IMPLICIT NONE


        tol=1.0e-9
        seedvalue(:)=1
        CALL random_seed(put=seedvalue(1:12))

        OPEN (UNIT=999, FILE="housing_transit.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        999 format (I6.2, I6.2, F16.6, F16.6, F16.6)

        OPEN (UNIT=799, FILE="transition_adjust.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        799 format (I6.2, I6.2, I6.2, F16.6)

        OPEN (UNIT=299, FILE="transition_fthb.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        299 format (I6.2, I6.2, I6.2, F16.6, F16.6, F16.6, F16.6, I6.2)

        OPEN (UNIT=499, FILE="transition_debug.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        499 format (I6.2, I6.2, I6.2, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

        write(*, *) "Solving transition dynamics"


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
        changepermconditional(numhouseholds, Tdie), demeanedindicator(numhouseholds, Tdie), adjustt(numhouseholds, Tdie), adjustmean(numhouseholds, Tdie))

        ALLOCATE(alive(numhouseholds, Tdie), housingstock(numhouseholds, Tdie, Tdie), &
        housingflow(numhouseholds, Tdie, Tdie), rentalflow(numhouseholds, Tdie, Tdie))



        do t =1, Tdie

        write(*,*) "Solving for HH at age ", t, " during price drop"
        ! STATE DEFINITION
        ! state 1: liquid assets a
        ! state 2: durable assets d
        ! state 3: idiosyncratic income shock: z
        ! state 4: house price
        ! state 5: age: t
        ! Notice we load in stats from earlier simulation

        currenthouseholdstate(:, 1) = householdresultsholder(:, t, 1)
        currenthouseholdstate(:, 2) = householdresultsholder(:, t, 2)
        currenthouseholdstate(:, 3) = householdresultsholder(:, t, 3)
        currenthouseholdstate(:, 5) = t
        do i = 1, numhouseholds
            ! First-time eligibility: never bought a durable or has
            ! been renting durable for EligYrsF years. The EligYrsF
            ! <= condition ensures policy turned off if that parameter
            ! is large enough.
            if (((householdresultsholder(i, t, 10)>=t) .or. &
                (t >= householdresultsholder(i, t, 11) + EligYrsF)) &
                .and. (EligYrsF <= Tdie) .and. (PolYrDiff >= 0)) then
                currenthouseholdstate(i, 6) = 1
                newhouseholdstate(i, 6) = 1
                currenthouseholdstate(i, 4) = 1
                rentalflow(i, max(t-1,1), t) = householdresultsholder(i, t, 13)
           ! Repeat eligibility: has owned the same durable level
           ! for EligYrsR years.
            else if ((t >= householdresultsholder(i, t, 12) + EligYrsR) &
                     .and. (PolYrDiff <= 0)) then
                currenthouseholdstate(i, 6) = 2
                newhouseholdstate(i, 6) = 2
                currenthouseholdstate(i, 4) = 2
            else
                currenthouseholdstate(i, 6) = 3
                newhouseholdstate(i, 6) = 3
                currenthouseholdstate(i, 4) = 3
            end if
        end do

        fthb_flag(:, 1) = householdresultsholder(:, t, 10)
        fthb_flag(:, 2) = householdresultsholder(:, t, 11)
        fthb_flag(:, 3) = householdresultsholder(:, t, 12)


        call random_number(shock)
        numtrans=0
        alive=1

        do tj=1, t-1
             write(999, 999) Tretire - tj, t, 0.0, 0.0, 0.0
        end do

        do tj=t, Tdie

            ! write(*, *) tj
            ! 5th state is age (or equiv) time
            currenthouseholdstate(:, 5)=tj
        !$OMP PARALLEL DO
            do i=1, numhouseholds

                if (tj + 1 <= t + PolYrs) then
                    thetastate(i) = thetanodes(currenthouseholdstate(i, 6))
                    hpdelta_i(i) = hpdelta(currenthouseholdstate(i, 6))
                end if

                if (alive(i, tj) == 1) then
                    ! first thing is that if you are old, do you die or not
                    if (tj>Tretire) then
                        if (alive(i, 1)==1) then
                            call random_number(shock)
                            if (shock<deathprob(tj-Tretire)) then
                                alive(i, tj:Tdie)=0
                            end if
                        end if
                    end if



                actualwealth(i, tj)=currenthouseholdstate(i, 1)+theta*currenthouseholdstate(i, 2)*exp(hpnodes(currenthouseholdstate(i, 4)))  !total wealth is voluntary equity plus equity in durable
                financialwealth(i, tj)=currenthouseholdstate(i, 1)-(1-theta)*currenthouseholdstate(i, 2)*exp(hpnodes(currenthouseholdstate(i, 4)))
                housingnet(i, tj)=theta*currenthouseholdstate(i, 2)*exp(hpnodes(currenthouseholdstate(i, 4)))
                housinggross(i, tj)=currenthouseholdstate(i, 2)*exp(hpnodes(currenthouseholdstate(i, 4)))


                !qchg(i, t) = financialwealth(i, t) + (1-thetastate(i))*currenthouseholdstate(i, 2)*exp(hpnodes(currenthouseholdstate(i, 4)-1))
                qchg(i, tj) = currenthouseholdstate(i, 1) + (1-theta)*currenthouseholdstate(i, 2)*hpdelta_i(i)
                call pol_linworking(qchg(i, tj), currenthouseholdstate(i, 2), currenthouseholdstate(i, 3), currenthouseholdstate(i, 4), currenthouseholdstate(i, 5), newhouseholdstate(i, 1), newhouseholdstate(i, 2), consumption(i, tj), rentalind(i, tj), choiceind(i, tj), welfare(i, tj))
                !call pol_linworking(currenthouseholdstate(i, 1)-.01, currenthouseholdstate(i, 2), currenthouseholdstate(i, 3), currenthouseholdstate(i, 4), currenthouseholdstate(i, 5), newhouseholdstate(i, 1), newhouseholdstate(i, 2), consumptionmpc(i, tj), rentalind(i, tj), welfare(i, tj))
                !call pol_linworking(currenthouseholdstate(i, 1), currenthouseholdstate(i, 2), currenthouseholdstate(i, 3), currenthouseholdstate(i, 4), currenthouseholdstate(i, 5), newhouseholdstate(i, 1), newhouseholdstate(i, 2), consumption(i, tj), rentalind(i, tj), choiceind(i, tj), welfare(i, tj))

                if (newhouseholdstate(i, 1)>amax) then
                    newhouseholdstate(i, 1)=.9999999*amax
                end if

                if (tj<=Tretire) then
                incomeholder(i, tj)=income(currenthouseholdstate(i, 3), tj)
                else
                incomeholder(i, tj)=retirementincome(currenthouseholdstate(i, 3))
                end if


                ! Decision tree begins with unambiguous renting
                if ((rentalind(i, tj) > .995) .or.  (rentalind(i, tj) /= rentalind(i, tj))) then
                    rentalflow(i, tj, t)=newhouseholdstate(i, 2)
                    newhouseholdstate(i, 2)=0
                    fthb_flag(i, 3) = Tdie+1
                ! Starts the FTHB policy timer.
                    if ((fthb_flag(i, 1) /= Tdie+1) .and. (fthb_flag(i, 2) == Tdie+1)) then
                        fthb_flag(i, 2) = tj
                    end if
                ! Track agents again as FTHBs after three consecutive years of renting
                    if ((tj >= fthb_flag(i, 2) + (EligYrsF - 1)) .and. (tj <= Tretire)) then
                        newhouseholdstate(i, 6)=1
                        fthb_flag(i, 1) = Tdie+1
                        fthb_flag(i, 2) = Tdie+1
                    end if
                else
                    housingstock(i, tj, t)=currenthouseholdstate(i, 2)
                    newhouseholdstate(i, 6)=3
                ! Track agents as eligible repeat buyer if their last purchase old enough
                    if ((tj >= fthb_flag(i, 3) + EligYrsR) .and. (t <= Tretire)) then
                        newhouseholdstate(i, 6)=2
                        fthb_flag(i, 3) = 0
                    end if
                ! Start the repeat buyer policy countdown whenever a durables purchase is made
                    if (choiceind(i, tj) == 1.0) then
                        numtrans(i) = 1
                        adjustt(i, t) = adjustt(i, t) + 1
                        housingflow(i, tj, t)=newhouseholdstate(i, 2)
                        if (currenthouseholdstate(i, 2) > 0) then
                            housingstock(i, tj, t)=newhouseholdstate(i, 2)
                        end if
                        thetastate(i) = thetanodes(currenthouseholdstate(i, 4))
                        hpdelta_i(i) = hpdelta(currenthouseholdstate(i, 4))
                    ! Insert policy taker characteristics again
                    ! (The final negative indicators are factors for which policy
                    !  was taken since Matlab doesn't handle mixed types well)
                        if ((fthb_flag(i, 1) == Tdie+1) .and. (tj <= Tretire)) then
                            fthb_flag(i, 1) = tj
                            write(299, 299) i, t, tj, newhouseholdstate(i, 2), &
                            rentalflow(i, max(tj-1, 1), t), 0.0, &
                            currenthouseholdstate(i, 2) + income&
                            &(currenthouseholdstate(i, 3), t), -1
                        else if ((fthb_flag(i, 3) == 0) .and. (tj <= Tretire)) then
                            fthb_flag(i, 1) = tj
                            write(299, 299) i, t, tj, newhouseholdstate(i, 2), &
                            currenthouseholdstate(i, 2), &
                            (currenthouseholdstate(i, 1) - currenthouseholdstate&
                            &(i, 2)*(1-theta))/currenthouseholdstate(i, 2), &
                            currenthouseholdstate(i, 2) + income&
                            &(currenthouseholdstate(i, 3), t), -2
                            ! Stop tracking repeat buyers after first repeat buy
                            fthb_flag(i, 3) = 99
                        end if
                        fthb_flag(i, 3) = tj
                    end if
                end if


                if ((tj - t < 3) .or. ((tj - t > 8) .and. (tj - t <= 12))) then
                    !write(499, 499) i, t, tj, newhouseholdstate(i, 2), rentalind(i, tj), newhouseholdstate(i, 1), consumption(i, tj), currenthouseholdstate(i, 6)
                    !write(499, 499) i, t, tj, newhouseholdstate(i, 1), choiceind(i, tj), newhouseholdstate(i, 2), currenthouseholdstate(i, 6), newhouseholdstate(i, 6), fthb_flag(i, 1), fthb_flag(i, 3)
                end if

                if (tj<=Tretire) then
                    call random_number(shock)

                    if (shock < Probzcum(currenthouseholdstate(i, 3), 1)) then
                            newhouseholdstate(i, 3)=1
                    else
                        do j=1, zgridsize-1
                            if (shock > Probzcum(currenthouseholdstate(i, 3), j) .AND. &
                                shock<=Probzcum(currenthouseholdstate(i, 3), j+1)) then
                                newhouseholdstate(i, 3)=j+1
                            end if
                        end do
                    end if
                else
                    newhouseholdstate(i, 3)=currenthouseholdstate(i, 3)
                end if

                if ((currenthouseholdstate(i, 4) /= hpgridsize) .and. (tj + 1 >= t + PolYrs)) then
                    newhouseholdstate(i, 4)=hpgridsize ! House prices from next period onwards are high
                else if (tj + 1 >= t + abs(PolYrDiff)) then
                    if (((fthb_flag(i, 1)>=t) .or. (t >= fthb_flag(i, 2) + EligYrsF)) &
                        .and. (EligYrsF <= Tdie) .and. (PolYrDiff >= 0)) then
                        currenthouseholdstate(i, 6) = 1 ! First-time eligibility
                        newhouseholdstate(i, 6) = 1
                        currenthouseholdstate(i, 4) = 1
                    else if ((t >= fthb_flag(i, 3) + EligYrsR) .and. (PolYrDiff <= 0)) then
                        currenthouseholdstate(i, 6) = 2 ! Repeat eligibility
                        newhouseholdstate(i, 6) = 2
                        currenthouseholdstate(i, 4) = 2
                    end if
                end if
                newhouseholdstate(i, 5)= tj+1
                currenthouseholdstate(i, :) = newhouseholdstate(i, :)
                end if

            end do
        !$OMP END PARALLEL DO

            if (tj <= Tretire) then
                write(999, 999) tj - t, t, sum(housingstock(:, tj, t)), sum(housingflow(:, tj, t)), sum(numtrans(:))
            end if

            numtrans = 0
            end do
            !write(5, 5) sum(financialwealth(:, t)*alive(:, t))/numhouseholds, sum(currenthouseholdstate(:, 2)*alive(:, t))/numhouseholds, sum(consumption(:, t)*alive(:, t))/numhouseholds, sum(actualwealth(:, t)*alive(:, t))/numhouseholds, sum(newhouseholdstate(:, 2)*alive(:, t))/numhouseholds, sum(newhouseholdstate(:, 2)*rentalind(:, t)*alive(:, t))/numhouseholds, sum(rentalind(:, t)*alive(:, t))/numhouseholds, sum(alive(:, t))/numhouseholds, sum(consumption(:, t))/numhouseholds
            !currenthouseholdstate=newhouseholdstate


            ! this is the housing price elasticity we want
            !write(*,*) t*1.0, (sum(diffc(:, t))/sum(alive(:, t)))/(sum(consumption(:, t))/sum(alive(:, t)))/(hpmin)
            !write(*,*) "elasticity", sum(diffc(:, t)*alive(:, t))/sum(consumption(:, t))/sum(alive(:, t))/0.2


           ! write(*,*) t*1.0, sum(currenthouseholdstate(:, 2)), sum(alive(:, t))

            !write(10, 10) (sum((ratioc(:, t)-1))/sum(alive(:, t)))/(exp(hpmin) - 1), (sum(diffc(:, t))/sum(alive(:, t)))/(sum(consumption(:, t))/sum(alive(:, t)))/(hpmin)
            !write(10, 10) (sum((ratioc(:, t)-1))/numhouseholds)/(exp(hpmin) - 1), (sum(diffc(:, t))/numhouseholds)/(sum(consumption(:, t))/numhouseholds)/(exp(hpmin) - 1), ratiocown/numown_t, ratiocrent/numrent_t, numown_t/(numown_t+numrent_t), numsell_t/numhouseholds, numbuy_t/numhouseholds,(sum(diffc(:, t)/hpmin)/numhouseholds), sum(aggregatecontrib(:, t))/numhouseholds

            !write(555, 555) a_t/numhouseholds, aown_t/numown_t, arent_t/numrent_t, leverage_t/numhouseholds, leverageown_t/numown_t, leveragerent_t/numrent_t

        do i=1, numhouseholds
        !    do t=1, 60
        !      write(9, 9) i, t, householdresultsholder(i, t, 1),  householdresultsholder(i, t, 2),  householdresultsholder(i, t, 3),  householdresultsholder(i, t, 4),  householdresultsholder(i, t, 5),  householdresultsholder(i, t, 6),  householdresultsholder(i, t, 7),  householdresultsholder(i, t, 8)
        !    end do
            adjustmean(i, :) = sum(incomeholder(i,:))/(Tdie-t)
            incomeholder(i,:) = (incomeholder(i,:) - adjustmean(i,:))**2
            write(799, 799) i, t, int(adjustt(i, t)), sum(incomeholder(i,:))/(Tdie-t-1)
        end do

        end do

        write(*,*) "end"

        close(299)
        close(499)
        close(799)
        close(999)

    end subroutine transition




end module lifecycle_transition
