program durables  !Program shell
    use share
    use lifeCycle_calibrate
    use lifeCycle_models
    USE OMP_LIB
    implicit none
    REAL(8), DIMENSION(hptransLength) :: diffH

    call open_log

    write(0,*) NEW_LINE('A') // "**  Predetermined parameters  **"
    write(0,'(5A20)' ) "F (cost)  |", "  theta (DP)  |", "  alpha  |", "  costlyequity  |", "  rentelasticity"
    write(0,'(3F20.4,L20,F20.4)') F, theta, elasticity2, costlyequity, rentelasticity
    write(0,'(5A20)') "Households  |", "  A grid  |", "  Dur grid  |", "  Z grid  |", "  TransPeriod grid"
    write(0,'(5I20)') numhouseholds, agridsize, Dgridsize, zgridsize, hptransLength
    write(0,'(4A20)') "Minimum A  |", "  Maximum A  |", "  Minimum D  |", "  Maximum D"
    write(0,'(4F20.3)') amin, amax, dmin, dmax
    write(0,*) "**  Max threads:", omp_get_max_threads(), "**"

    call system_clock(timestart, timeres)

    call get_grids
    call income_array

    call system_clock(timeend, timeres)
    write(*,*) "Initial grid setup complete"
    write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"
 
    EV=0
    numiter=1 ! TODO: Change later?
    write(0,*) NEW_LINE('A') // "**  Calibration iteration:", numiter, "**"
    write(0,*) "**  Parameters to be tested  **"
    write(0,'(6A15)') 'Psi  |', '  Rentprem  |', '  RentRetire  |','  Dmin  |','  Ret_wealth  |','  Beta'
    write(0, '(6F15.4)') psi, rentprem, rentPremRetire, dmin, ret_wealth, beta2

    if (pe_start) call PE_subsidies

    if (ge_start) call GE_transition_path(startprice, diffH)


    CONTAINS

    subroutine get_grids ! %<
        IMPLICIT NONE
        INTEGER :: i, Dapprox
        ! FILL IN ALL THE GRID POINTS

        write(0,*) NEW_LINE('A') // "// Asset search grid:"
        
        ! Transform the grid so more mass with small numbers.
        ! The log is a normalizing factor.
        DO i=1, 6
        ! Negative assets have a fixed number of nodes
            anodes(i)=(1.0/5.0)*(1.0*i-1.0)
            anodes(i) = exp(log(0.0-amin+1.0)*anodes(i))+amin-1.0
        end do
        DO i=7, agridsize
            ! Double transformation of uniform grid to concentrate
            ! grid points in the 0-1 range.
            anodes(i)=(1.0/(agridsize-6))*(1.0*(i-5)-1.0)
            anodes(i)=(exp(log(amax+1)*anodes(i))-1.0)/amax
            anodes(i)=exp(log(amax+1)*anodes(i))-1.0
        end do
        write(0,'(10F16.10)') anodes

        write(0,*) NEW_LINE('A') // "// Durable search grid:"
         DO i=1, Dgridsize
            Dnodes(i)=(1.0/(Dgridsize-1))*(1.0*i-1.0)
        end do
        !
        if (Dmin > 0) then
        Dapprox = MINLOC(ABS(Dnodes*Dmax - Dmin), 1)
        write(*,*) Dapprox
        write(*,*) Dnodes(max(Dapprox-3,2):Dgridsize)
        write(*,*) Dmax*Dnodes(1:Dapprox)
        DO i=max(Dapprox-3,2), Dgridsize-6
            Dnodes(i)=exp(log(dmax-Dnodes(max(Dapprox-4, 1))+1)*Dnodes(i))+dmax*Dnodes(max(Dapprox-4, 1))-1.0
        end do
        DO i=Dgridsize-5, Dgridsize
            Dnodes(i) = Dnodes(Dgridsize-6) + (Dmax-Dnodes(Dgridsize-6))*((i-(Dgridsize-6))/6.0)
        END DO
            if (Dapprox-4 > 1) then
                DO i=2, Dapprox-4
                write(*,*) Dnodes(i)
                Dnodes(i) = Dnodes(Dapprox-3)*(i-1.0)/(real(Dapprox-3)-1.0)
                write(*,*) Dnodes(i)
        end do
            end if
        Dapprox = MINLOC(ABS(Dnodes - Dmin), 1)
        Dmin = Dnodes(Dapprox)
        do i=1, Dgridsize
        !    if (i > Dapprox) Dnodes(i) = Dmin*(1.0/(1.0-delta))**i
        end do
        else
            do i=1, Dgridsize
                Dnodes(i)=exp(log(dmax-0.0+1)*Dnodes(i))+0.0-1.0
            end do
        end if
        write(0,'(10F16.10)') Dnodes
        write(0,'(A20, F16.10)') 'Imputed Dmin:', Dmin

        call share_arrays

    end subroutine !%>
   
    subroutine PE_subsidies ! %<
! Main PE model code: Solves the backward induction
        IMPLICIT NONE

        ! Check hptransLength=1 here
        if (hptransLength /= 1) then
            write(*,*) "Partial equilibrium model requires only one price is&
                set for the entire problem. Quitting."
            RETURN
        end if

        ! Agent state holders in each (transition) period
        ALLOCATE(householdtransitionholder(13,numhouseholds,TDie+1,Tretire+1))

        ALLOCATE(housingstock(numhouseholds,Tdie,Tretire+1),&
        housingflow(numhouseholds,Tdie,Tretire+1),rentalflow(numhouseholds,Tdie,Tretire+1),&
        bequestflow(numhouseholds,Tdie,Tretire+1),newownerflow(numhouseholds,Tdie,Tretire+1),&
        resaleflow(numhouseholds,Tdie,Tretire+1),numadjust(numhouseholds,Tdie,Tretire+1))
        ! Fix a flat subsidy/tax value
        balancer_internal = balancer(1)
        
        call PE_vfunc

        if (.NOT. ss_only) call PE_transition_vfunc

        call PE_simulate
        write(0,*) NEW_LINE('A') // "// C sum:", sum(cchoice(:,:,:,:,:,1:hptransLength))

        close(1)
    end subroutine ! %>

    subroutine GE_transition_path(pstart, yend) ! %<
        IMPLICIT NONE
        INTEGER :: i
        REAL(8), DIMENSION(2,hptransLength), INTENT(INOUT) :: pstart
        REAL(8), DIMENSION(hptranslength), INTENT(OUT) :: yend
        REAL(8), DIMENSION(hptransLength-1) :: priceforward
        REAL(8) :: avghousing, construction, totrev, excess

        ! Check hptransLength > 1 here
        if ((hptransLength == 1) .AND. .NOT. ss_only) then
            write(*,*) "General equilibrium transition requires at least one&
                & transition period with endogenous prices away from steady&
                & state. Quitting."
            RETURN
        end if

        ! Agent state holders in each (transition) period
        ALLOCATE(householdtransitionholder(13,numhouseholds,TDie+1,hptransLength))

        ALLOCATE(housingstock(numhouseholds,Tdie,hptransLength),&
        housingflow(numhouseholds,Tdie,hptransLength),rentalflow(numhouseholds,Tdie,hptransLength),&
        bequestflow(numhouseholds,Tdie,hptransLength),newownerflow(numhouseholds,Tdie,hptransLength),&
        resaleflow(numhouseholds,Tdie,hptransLength),numadjust(numhouseholds,Tdie,Tretire+1))

        if (.NOT. steady_state_block) then
        !block: define variable in share
        call system_clock(timestart, timeres)
        write(*,*) "Begin steady-state computation..."
        write(0,*) NEW_LINE('A') // "Begin steady-state computation..."
        call GE_steady_state(pstart(:,1), yend(1))
        write(*,*) "Steady-state agent distribution computed"
        write(0,*) "Steady-state agent distribution computed"
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"
        if (ABS(yend(1)) >= 5e-2) then
            write(*,*) "Market equilibrium failed to converge given the interval of possible&
                       & price values. Quitting."
            RETURN
        end if
        else 
            write(*,*) "Assume partial equilibrium with fixed price level"
            call PE_vfunc
            ! Fix a flat subsidy/tax value
            balancer_internal = balancer(1)
        end if
        
        write(0,*) NEW_LINE('A') // "// Steady-state market clearing price:", hpnodes(1)
        ! end block
        ! TODO: Check this
        hpnodes(hptransLength+1)=hpnodes(1)
        EV(:, :, :, :, :, hptransLength+1) = EV(:, :, :, :, :, 1)
        if (.NOT. loadpricepath) then
            hpnodes(2) = hpnodes(1) + 1.0e-2
            do i=3,hptransLength
                hpnodes(i)= hpnodes(1) - 1.5e-2*(hptransLength+1-i)/(hptransLength-2)
            end do
            pstart(1,2:hptransLength)=pstart(1,2:hptransLength) + hpnodes(2:hptransLength)
            pstart(2,2:hptransLength)=pstart(2,2:hptransLength) + hpnodes(2:hptransLength)
        else
            !col of numbers encoding % change in price
            OPEN(97, FILE='input_data/currentPricePath.txt', STATUS='OLD')
            do i=2, hptransLength
                READ(97,*) hpnodes(i)
            end do
            write(*,*) hpnodes
        end if

        write(*,*) "Calculating full panel in steady-state"
        write(0,*) NEW_LINE('A') // "Calculating full panel in steady-state"
        ALLOCATE(incshocks(numhouseholds, Tretire+1))
        call gen_life_shocks(numhouseholds, shock(5,:,:), incshocks)
        call simulate(hpnodes(1), numhouseholds, householdtransitionholder(:,:,:,1),&
                      shock, .TRUE.)
        call get_mktclear(hpnodes(1), 1, avghousing, construction, totrev)
        ! Blank out the steady state values
        consumption = 0
        mortint = 0
        incomeholder = 0
        DEALLOCATE(mpc, posttaxincome, alive)
        write(0,'(3A30)') "// Housing in SS", "  |  New homeowners in SS", "  |  Resold housing in SS"
        write(0,'(3F30.2)') SUM(housingstock(:,:,1)), SUM(newownerflow(:,:,1)), SUM(resaleflow(:,:,1))
        write(0,'(2A30)') "// Final rental homes in SS", "  |  Bequested housing in SS"
        write(0,'(2F30.2)') SUM(rentalflow(:,:,1)), SUM(bequestflow(:,:,1))
        if (ss_only) then
            close(1)
            RETURN
        end if

        ! If the law of motion for prices pre loaded we don't iterate back and
        ! forth until equilibrium found
        if (loadpricepath) then
            call bequests(hpnodes(2:hptransLength+1), 2)

            call GE_transition_vfunc

            ! Calculation of value function with policy
            TransTime = 0  ! Represents policy period or 0 periods post policy
            excess = SimTransPath(hpnodes(2), numhouseholds, .FALSE., .FALSE., .TRUE.)
        end if

        do while (.NOT. loadpricepath .AND.  &
            maxval(abs(priceforward - hpnodes(2:hptransLength))) >= pthres)

            call bequests(hpnodes(2:hptransLength+1), 2)

            ! Backwards induction: EV grids on the (:, :, :, t, l) dimensions
            ! are updated. Agents face the competitive price for rental housing
            call GE_transition_vfunc

            ! Forward induction: simulate transition path one period at a time
            ! and recalibrate the market clearing price
            call GE_transition_clearing(pstart(:,2:hptranslength),&
                                        yend(2:hptranslength), priceforward)
            write(0,*) NEW_LINE('A') // "Criterion:", maxval(abs(priceforward - hpnodes(2:hptransLength)))

        end do
        
        ! Final GE iteration
        !call bequests(hpnodes(2:hptransLength+1), 2)
        write(0,*) NEW_LINE('A') // "Transition path for prices:"
        do i=2,hptranslength+1
            write(0,'(F16.6)') hpnodes(i)
        end do
        !call GE_transition_vfunc
        call transition_final(hptransLength-1, .FALSE.)
        DEALLOCATE(alive, consumption, incomeholder, mortint)
         
        OPEN (UNIT=78, FILE="transition_aggs.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        write(78,'(I6.2,5F16.6)') 1, hpnodes(1), SUM(housingstock(:,:,1)),&
                                  SUM(rentalflow(:,:,1)),SUM(resaleflow(:,:,1)),&
                                  SUM(newownerflow(:,:,1))
        do i=2,hptranslength
            write(78,'(I6.2,5F16.6)') i, hpnodes(i)/hpnodes(1),&
                                      SUM(housingstock(:,:,i))/SUM(housingstock(:,:,1)),&
                                      SUM(rentalflow(:,:,i))/SUM(rentalflow(:,:,1)),&
                                      SUM(resaleflow(:,:,i))/SUM(resaleflow(:,:,1)),&
                                      SUM(newownerflow(:,:,i))/SUM(newownerflow(:,:,1))
        end do
        close(78) 
        close(1)

    end subroutine ! %>

     
end program durables
