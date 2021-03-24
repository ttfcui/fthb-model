module lifecycle_models

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE MODELS contain subroutines created to manage the workflow
! of running PE versus GE models. Some output files are also created
! here.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use share
    use lifeCycle_solveDP
    use lifeCycle_simulate
    use lifeCycle_transition
    IMPLICIT NONE
    INTEGER :: timestart, timeend
    INTEGER :: i
    INTEGER :: timeres

    CONTAINS

    subroutine PE_vfunc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Subroutine computes the DP problem under partial equilibrium
        !    (constant prices, no market clearing).
        !
        !    Modified: 08/09/2019
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        write(*,*) "Begin value function iteration:"
        write(0,*) NEW_LINE('A') // "Begin value function iteration:"

        call system_clock(timestart, timeres)
        ! Agents face the competitive price for rental housing
        r_rental = rentPrem + exp(hpnodes)*rent
        r_rental_retire = rentPremRetire + exp(hpnodes)*rent
        call bequests(hpnodes, 1)
        EVMov(:,:,:,:,Tdie+1,:) = EV(:,:,:,:,Tdie+1,:)
        EVMovR(:,:,:,:,Tdie+1,:) = EV(:,:,:,:,Tdie+1,:)
        call solveworkingproblem(hpnodes(1:1), balancer(1:1),&
             achoice, Dchoice, rentchoice, cchoice, choiceindicator,&
             achoiceMov, DchoiceMov, cchoiceMov, choiceindicatorMov,&
             achoiceMovR, DchoiceMovR, cchoiceMovR,chindMR,&
             .TRUE., .TRUE., 'None', EV, EVMov, EVMovR)

        ! Broadcasting DP outputs forward a period for consistency in the
        ! transition subroutine
        achoice(:,:,:,:,:,2) = achoice(:,:,:,:,:,1)
        Dchoice(:,:,:,:,:,2) = Dchoice(:,:,:,:,:,1)
        rentchoice(:,:,:,:,:,2) = rentchoice(:,:,:,:,:,1)
        cchoice(:,:,:,:,:,2) = cchoice(:,:,:,:,:,1)
        choiceindicator(:,:,:,:,:,2) = choiceindicator(:,:,:,:,:,1)
        EV(:,:,:,:,:,2) = EV(:,:,:,:,:,1)
        achoiceMov(:,:,:,:,:,2) = achoiceMov(:,:,:,:,:,1)
        DchoiceMov(:,:,:,:,:,2) = DchoiceMov(:,:,:,:,:,1)
        cchoiceMov(:,:,:,:,:,2) = cchoiceMov(:,:,:,:,:,1)
        choiceindicatorMov(:,:,:,:,:,2) = choiceindicatorMov(:,:,:,:,:,1)
        EVMov(:,:,:,:,:,2) = EVMov(:,:,:,:,:,1)
        achoiceMovR(:,:,:,:,:,2) = achoiceMovR(:,:,:,:,:,1)
        DchoiceMovR(:,:,:,:,:,2) = DchoiceMovR(:,:,:,:,:,1)
        cchoiceMovR(:,:,:,:,:,2) = cchoiceMovR(:,:,:,:,:,1)
        chindMR(:,:,:,:,:,2) = chindMR(:,:,:,:,:,1)
        EVMovR(:,:,:,:,:,2) = EVMovR(:,:,:,:,:,1)

        write(*,*) "Value function iteration complete"
        write(0,*) "Value function iteration complete"
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

    end subroutine PE_vfunc ! %>

    subroutine PE_transition_vfunc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Subroutine computes the policy schedule and a separate
        !    array for policy functions when the temporary subsidy is reached.
        !    (stationary eq. possibly reached, but no price change on the
        !     transition path).
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        LOGICAL :: pol_at_start = .TRUE.
        CHARACTER(11) :: pol_perm = 'None'
        INTEGER, DIMENSION(2) :: pol_vec
        REAL(8) :: transphase_slope, trans_yint

        ! Constructs down payment at market price (since we interpret changes in
        ! hp as subsidies). If downtransfer, down subtracted by value of subsidy.
        call system_clock(timestart, timeres)
        ! New price level (proportionate change from baseline)
        noadjTransfer_internal= noadjTransfer
        ! The percentage/log change in prices from baseline
        write(0,'(2A25)') '% monetary transfer  |', 'Transfer applied to down'
        write(0,'(2F25.5)') adjTransfer(1,1), eta_transfer*adjTransfer(1,1)

        ! %< Generating the policy subsidy schedule
        if (trans_phase_int > 0 .AND. trans_phase_out > 0) then
            transphase_slope = MAXVAL(adjTransfer)/(trans_phase_out - trans_phase_int)
            trans_yint = MAXVAL(adjTransfer) + trans_phase_int*transphase_slope
            WHERE (income >= trans_phase_int .AND. income <= trans_phase_out) &
                adjTransfer = trans_yint - income*transphase_slope
        end if
        if (trans_phasein_int > 0) then
            WHERE (income <= trans_phasein_out) adjTransfer = &
                trans_phasein_int + income*(MAXVAL(adjTransfer) - trans_phasein_int)&
                &/trans_phasein_out
        end if ! %>

        write(*,*) "Begin value function iteration with subsidy:"
        write(0,*) NEW_LINE('A') // "Begin value function iteration with subsidy:"

        ! Set subset of state space over which policy subsidy is applicable
        if ((EligYrsR > TDie) .AND. (EligYrsF <= TDie)) then ! FTHB case
            pol_vec = (/ 2, TDie /)
        else ! Repeat case
            pol_vec = (/ dgridsize, TDie /)
        end if

        ! Generating different kinds of policy arrays in presence of policy
        if (PolEnd == 1 .AND. PolStart == 0) then
            do i=1,size(polParam, 1)
                if (polParam(i) == 1.0) polLevel(i) = polRead(i)
            end do
            call solveworkingproblem(hpnodes(1:1), balancer(1:1),&
                achoiceshock, Dchoiceshock, rentchoiceshock, cchoiceshock, choiceindicatorshock,&
                achoiceshockMov, DchoiceshockMov, cchoiceshockMov, choiceindicatorshockMov,&
                achoiceshockMovR, DchoiceshockMovR, cchoiceshockMovR,chindshockMR,&
                .FALSE., .TRUE., 'None', EVpol, EVpolmov, EVpolmovR, pol_vec)

            EVpol(:,:,:,:,:,2) = EVpol(:,:,:,:,:,1)
            EVpolMov(:,:,:,:,:,2) = EVpolMov(:,:,:,:,:,1)
            EVpolMovR(:,:,:,:,:,2) = EVpolMovR(:,:,:,:,:,1)
            ! Revert as soon as policy DP calculated (variables hardcoded in other
            ! programs)
            noadjTransfer_internal= 0.0
            polLevel = 0.0
        else
            ! Any other case involves anticipating the policy in a future period
            ! so refer to the policy array "*expect"
            EVpol(:,:,:,:,TDie+1,:) = EV(:,:,:,:,Tdie+1,:)
            EVpolMov(:,:,:,:,Tdie+1,:) = EV(:,:,:,:,Tdie+1,:)
            EVpolMovR(:,:,:,:,Tdie+1,:) = EV(:,:,:,:,Tdie+1,:)
            if (PolEnd > Tretire) pol_perm = 'first-time'
            call solveworkingproblem(hpnodes(1:1), balancer(1:1),&
                achoiceshock, Dchoiceshock, rentchoiceshock, cchoiceshock, choiceindicatorshock,&
                achoiceshockMov, DchoiceshockMov, cchoiceshockMov, choiceindicatorshockMov,&
                achoiceshockMovR, DchoiceshockMovR, cchoiceshockMovR,chindshockMR,&
                .FALSE., .TRUE., trim(pol_perm), EVpol, EVpolMov, EVpolMovR, pol_vec)

            write(*,*) "Begin value function iteration with anticipated subsidy:"
            write(0,*) NEW_LINE('A') // "Begin value function iteration with anticipated subsidy:"
            EVpol(:,:,3:Dgridsize,:,:,:) = EV(:,:,3:Dgridsize,:,:,:)
            EVpolMov(:,:,3:Dgridsize,:,:,:) = EVMov(:,:,3:Dgridsize,:,:,:)
            EVpolMovR(:,:,3:Dgridsize,:,:,:) = EVMovR(:,:,3:Dgridsize,:,:,:)
            noadjTransfer_internal= 0.0
            if (PolStart == 0) pol_at_start = .FALSE.

            call solveworkingproblem(hpnodes(1:1), balancer(1:1),&
                achoiceexpect, Dchoiceexpect, rentchoiceexpect, cchoiceexpect, choiceindicatorexpect,&
                achoiceexpectMov, DchoiceexpectMov,cchoiceexpectMov, choiceindicatorexpectMov,&
                achoiceexpectMovR, DchoiceexpectMovR,cchoiceexpectMovR,chindexpectMR,&
                pol_at_start, .FALSE., 'first-time', EVpol, EVpolMov, EVpolMovR, pol_vec)
        end if
        call plotpolicy_constantwealth(1, achoiceshock, Dchoiceshock, cchoiceshock, &
                                       choiceindicatorshock)

        ! Call the aggregate shocks again (because they may be brought up in
        ! simulation)
        do i=1,size(polParam, 1)
            if (polParam(i) == 1.0) polLevel(i) = polRead(i)
        end do

        write(*,*) "Value function iteration complete"
        write(0,*) "Value function iteration complete"
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"
    end subroutine ! %>

    subroutine PE_simulate
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Subroutine computes the simulation in steady-state prior to the policy
        !    and then calls a period-by-period simulation of the transition path.
        !    (stationary eq. possibly reached, but no price change on the
        !     transition path).
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8) :: avghousing, construction, totrev

        !simulate model and check match to some target moments:
        write(*,*) "Simulation starting..."
        write(0,*) "Simulation starting..."
        call system_clock(timestart, timeres)
        ! Plus 1 because first period is calibrated state
        ALLOCATE(incshocks(numhouseholds, Tretire+1))
        call gen_life_shocks(numhouseholds, shock(6,:,:), incshocks)
        call simulate(hpnodes(1),numhouseholds,householdtransitionholder(:,:,:,1),&
                      shock, .TRUE.)
        call get_mktclear(hpnodes(1), 1, avghousing, construction, totrev)
        DEALLOCATE(posttaxincome,alive)


        if (ss_only) then
            close(1)
            RETURN
        end if

        balancer(2) = balancer(1)
        call transition_final(Tretire, .TRUE.)

    end subroutine ! %>

    subroutine transition_final(tot_time, partial)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Subroutine simulates the transition path period by period
        !    (in PE or in Rational Expectations) and writes a file comparing
        !    durable adjustments and policy responses over the transition
        !    simulation's duration.
        !
        !
        !    PARAMETERS
        !
        !    tot_time: Integer, duration of the policy transition period.
        !
        !    partial: Boolean. If TRUE, the policy transition is run in
        !    partial equilibrium, i.e. prices remain constant during the
        !    whole transition path.
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: tot_time
        LOGICAL, INTENT(IN) :: partial
        LOGICAL :: write_files = .TRUE.
        LOGICAL :: partBool
        INTEGER :: i, l, t, tpol, pindex
        REAL :: excess
 
        write(*, *) "Solving transition dynamics"
        agg_policies = .TRUE.
        do l=1,tot_time
            TransTime = l - 1
            write(*,*) "Transition period: ", TransTime
            if (partial) then
                pindex = 1
                partBool = .TRUE.
                call transition(hpnodes(pindex),numhouseholds,TransTime,&
                    householdtransitionholder(:,:,:,l),&
                    householdtransitionholder(:,:,:,l+1),shock,partBool, .TRUE.)
            else
                pindex = l+1
                partBool = .FALSE.
                DEALLOCATE(alive)
                excess = SimTransPath(hpnodes(pindex), numhouseholds, .TRUE., .TRUE., .FALSE.)
            end if
        end do
        close(71)
        close(72)
        close(73)
        close(74)
        if (print_micro .AND. tot_time > 10) call SimTransLeads

        OPEN (UNIT=76, FILE="transition_adjust.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        OPEN (UNIT=77, FILE="transition_difffull.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        OPEN (UNIT=79, FILE="transition_lives.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        if (write_files .AND. print_micro) then
        do i=1,numhouseholds,10
            do t=1, TDie
                write(76, '(3I6.2)') i, t, sum(numadjust(i, t, 2:tot_time+1))
                write(77, *) i, t, &
                            numadjust(i, t, 2:2+(tot_time-1)-max(t+tot_time-1-Tdie, 0))&
                            - numadjust(i,t:min(t+tot_time-1,TDie), 1)
            end do
        end do
        end if

        ! Sample lifecycle housing, consumption paths
        do i=999, numhouseholds/2, 10
            do t=3, Tretire, 5
            if (numadjust(i, t, 2) == 1) then
                do tpol=t,t+tot_time-2
                    if (tpol < Tdie - 1) &
                    write(79, '(3I6.2, 9F16.10, 2I6.2)') i, t, tpol-t, &
                    householdresultsholder(10, i, tpol),&
                    householdresultsholder(8, i, tpol),&
                    householdtransitionholder(3, i, tpol+1, 2+tpol-t),&
                    householdresultsholder(6, i, tpol),&
                    householdtransitionholder(11, i, tpol+1, 2+tpol-t),&
                    householdtransitionholder(12, i, tpol+2, 1),&
                    householdtransitionholder(12, i, tpol+2, 3+tpol-t),&
                    householdresultsholder(9, i, tpol),&
                    householdtransitionholder(9, i, tpol+1, 2+tpol-t),&
                    numadjust(i,tpol, 1), numadjust(i, t, 2+tpol-t)
                end do
            end if
            end do
        end do
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"
        close(76)
        close(77)
       
    end subroutine ! %>

    subroutine GE_steady_state(pstart, yend)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Subroutine is a wrapper for the iterative algorithm finding
        !    the stationary equilibrium in a GE setting.
        !
        !    Modified: 06/22/2017
        !
        !    PARAMETERS
        !
        !    pstart: Vector of length 2: dictates the range in which
        !    Brent's method attempts to find the stationary equilibrium.
        !
        !    yend: Output. Real number indicating value of the objective function
        !    used to find the equilibrium.
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), DIMENSION(2,1), INTENT(INOUT) :: pstart ! One price in GE for the moment
        REAL(8), INTENT(OUT) :: yend

        ! Plus 1 because first period is calibrated state
        ALLOCATE(incshocks(numhouseholds, Tretire+1))
        call gen_life_shocks(numhouseholds, shock(6,:,:), incshocks)

        ! NOTE: the function SimSteady includes the backwards induction calculation
        yend = brentmindist(pstart(1,1), pstart(1,1),pstart(2,1),SimSteady,numhouseholds,&
                            steady_conv_thres,hpnodes(1))
        DEALLOCATE(incshocks)
    end subroutine ! %>

    subroutine GE_transition_vfunc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Subroutine is a wrapper to call solving the DP problem
        !    in GE over the length of the transition period.
        !
        !    Modified: 05/12/2019
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE

        call system_clock(timestart, timeres)
        r_rental = rentPrem + exp(hpnodes)*rent
        r_rental_retire = rentPremRetire + exp(hpnodes)*rent
        ! Recall l=1 is the steady-state price.
        call solveworkingproblem(hpnodes(2:hptransLength+1), balancer(2:hptransLength+1),&
            achoice(:,:,:,:,:,2:hptransLength), Dchoice(:,:,:,:,:,2:hptransLength),&
            rentchoice(:,:,:,:,:,2:hptransLength), cchoice(:,:,:,:,:,2:hptransLength), &
            choiceindicator(:,:,:,:,:,2:hptransLength),&
            achoiceMov(:,:,:,:,:,2:hptransLength), DchoiceMov(:,:,:,:,:,2:hptransLength),&
            cchoiceMov(:,:,:,:,:,2:hptransLength), choiceindicatorMov(:,:,:,:,:,2:hptransLength),&
            achoiceMovR(:,:,:,:,:,2:hptransLength), DchoiceMovR(:,:,:,:,:,2:hptransLength),&
            cchoiceMovR(:,:,:,:,:,2:hptransLength),chindMR(:,:,:,:,:,2:hptransLength),& 
            .TRUE., .TRUE., 'None', EV(:,:,:,:,:,3:hptransLength+1), EVmov(:,:,:,:,:,3:hptransLength+1),&
            EVmovR(:,:,:,:,:,3:hptransLength+1))

        write(*,*) "Value function iteration complete:"
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

     end subroutine ! %>

    REAL(8) FUNCTION SimTransIter(price, numHH)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: price
        INTEGER, INTENT(IN) :: numHH

        SimTransIter = SimTransPath(price, numHH, .FALSE., .FALSE., .FALSE.)
        DEALLOCATE(alive)

    end function


    subroutine GE_transition_clearing(pstart, yend, priceforward)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Subroutine is a wrapper for the iterative algorithm finding
        !    the RE transition path in a GE setting.
        !
        !    Modified: 06/22/2017
        !
        !    PARAMETERS
        !
        !    pstart: Array of 2 x period duration: dictates the range in which
        !    Brent's method attempts to find the stationary equilibrium.
        !
        !    yend: Output. Real vector indicating value of the objective function
        !    used to find the equilibrium.
        !
        !    priceforward: Output. Real vector recording the equilibrium
        !    price vectors found through forward iteration through each period.
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), DIMENSION(2,hptransLength-1), INTENT(INOUT) :: pstart
        REAL(8), DIMENSION(hptranslength-1), INTENT(OUT) :: yend
        REAL(8), DIMENSION(hptransLength-1), INTENT(OUT) :: priceforward
        REAL(8) :: transphase_slope, trans_yint
        INTEGER :: l

        ! New price level (proportionate change from baseline)
        noadjTransfer_internal= noadjTransfer
        ! The percentage/log change in prices from baseline
        write(0,'(2A25)') '% monetary transfer  |', 'Transfer applied to down'
        write(0,'(2F25.5)') adjTransfer, eta_transfer*adjTransfer
        if (trans_phase_int > 0 .AND. trans_phase_out > 0) then
            transphase_slope = MAXVAL(adjTransfer)/(trans_phase_out - trans_phase_int)
            trans_yint = MAXVAL(adjTransfer) + trans_phase_int*transphase_slope
            WHERE (income >= trans_phase_int .AND. income <= trans_phase_out) &
                adjTransfer = trans_yint - income*transphase_slope
        end if
        if (trans_phasein_int > 0) then
            WHERE (income <= trans_phasein_out) adjTransfer = &
                trans_phasein_int + income*(MAXVAL(adjTransfer) - trans_phasein_int)&
                &/trans_phasein_out
        end if
        write(*,*) trans_phase_int, trans_phase_out

        ! TODO: Do I want to run solveworkingproblem for one transition period
        ! in SimTransPath? Could very well be what Stroebel does
        do l=1,hptransLength-1
            TransTime = l - 1
            balancer_internal = balancer(l+1)
            write(*,*) "Transition period:", TransTime
            write(*,*) pstart(:,l)
            yend(l) = brentmindist(hpnodes(l+1), pstart(1,l),pstart(2,l),SimTransIter,&
                numhouseholds, transition_conv_thres,priceforward(l))
        ! Compare the initial price guesses and the recalibrated price
            if ((yend(l) < 100) .AND. (abs(priceforward(l) - hpnodes(l+1)) >= pthres)) then
                hpnodes(l+1) = 0.90*priceforward(l) + 0.10*hpnodes(l+1)
            end if
            write(0,*) "Market clearing in period:", TransTime
            write(0,*) NEW_LINE("A"), "Updated price for period:", priceforward(l), hpnodes(l+1)
        end do

     end subroutine ! %>

end module lifecycle_models
