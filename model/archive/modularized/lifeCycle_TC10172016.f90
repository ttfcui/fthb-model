program durables  !Program shell
use share
use lifeCycle_solveDP
use lifeCycle_simulate
use lifeCycle_transition
use lifeCycle_calibrate
USE OMP_LIB
implicit none
integer :: timestart, timeend
REAL(8) :: timeres
REAL(8) :: dd, m2, m2young, m2middle, m2old

REAL(8), DIMENSION(2,hptransLength) :: startprice
REAL(8), DIMENSION(hptransLength) :: diffH

elasticity2=elasticity2O
beta2=beta2O
beta2retire=beta2retireO


write(*,*) "amin", amin

write(*,*) "Max threads:", omp_get_max_threads()

numiter=0


!Policy functions when solving problem
ALLOCATE(achoice(agridsize, Dgridsize, zgridsize, hptransLength+1, Tdie),&
         Dchoice(agridsize, Dgridsize, zgridsize, hptransLength+1, Tdie),&
         cchoice(agridsize, Dgridsize, zgridsize, hptransLength+1, Tdie))
ALLOCATE(rentalindicator(agridsize, Dgridsize, zgridsize, hptransLength+1, Tdie),&
         choiceindicator(agridsize, Dgridsize, zgridsize, hptransLength+1, Tdie))
ALLOCATE(EV(agridsize, Dgridsize, zgridsize, hptransLength+1, Tdie+1))

!Policy functions for storing subsidy policy response
ALLOCATE(achoiceshock(agridsize, Dgridsize, zgridsize, 1, Tdie),&
         Dchoiceshock(agridsize, Dgridsize, zgridsize, 1, Tdie),&
         cchoiceshock(agridsize, Dgridsize, zgridsize, 1, Tdie))
ALLOCATE(rentalindicatorshock(agridsize, Dgridsize, zgridsize, 1, Tdie),&
         choiceindicatorshock(agridsize, Dgridsize, zgridsize, 1, Tdie))

EV=0

diffmoments=1
do while (diffmoments>difftol)
    numiter=numiter+1
    write(*,*) "numiter", numiter
    write(*,*) elasticity2, rentprem, beta2

    !OPEN (UNIT=1, FILE="vfunc.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    !1 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

    call system_clock(timestart, timeres)

    call get_grids
    call tauchen
    call income_array

    call system_clock(timeend, timeres)
    write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

    !call PE_subsidies
    startprice(1,1) = 0.50
    startprice(2,1) = 0.00
    startprice(1,2) = 2.0
    startprice(2,2) = 0.0
    startprice(1,3:hptransLength) = 0.3
    startprice(2,3:hptransLength) = -0.5
    call GE_transition_path(startprice, diffH)
    diffmoments=0.0

end do



write(*,*) "model version"
write(*,*) "F  |  ", "theta  |  ", "rentprem  |  ", "costlyequity  |  ", "rentelasticity  |  "
write(*,*) F, theta, rentprem, costlyequity, rentelasticity
write(*,*) "A grid  |  ", "Dur grid  |  ", "Z grid  |  ", "HPrice grid  |  ", "households"
write(*,*) agridsize, Dgridsize, zgridsize, hptransLength, numhouseholds


    CONTAINS

    subroutine get_grids ! %<
        IMPLICIT NONE
        ! FILL IN ALL THE GRID POINTS

        ! Shocks grid - equal intervals a la Tauchen discretization
        write(*,*) "Shocks grid:"
        do i=1, zgridsize
            znodes(i)=zmin+((zmax-zmin)/(zgridsize-1))*(1.0*i-1.0)
            !write(*,*) znodes(i)
        end do

        write(*,*) "Asset search grid:"
        DO i=1, agridsize
            anodes(i)=(1.0/(agridsize-1))*(1.0*i-1.0)
        end do
        ! Transform the grid so more mass with small numbers.
        ! The log is a normalizing factor.
        do i=1, agridsize
            anodes(i)=exp(log(amax-amin+1)*anodes(i))+amin-1.0
            !write(*,*) anodes(i)
        end do

        write(*,*) "Durable search grid:"
         DO i=1, Dgridsize
            Dnodes(i)=(1.0/(dgridsize-1))*(1.0*i-1.0)
        end do
        do i=1, dgridsize
            Dnodes(i)=exp(log(dmax-dmin+1)*Dnodes(i))+dmin-1.0
            !write(*,*) Dnodes(i)
        end do

        call share_arrays

    end subroutine !%>
   
    subroutine tauchen !%<
        IMPLICIT NONE
        REAL(8) :: w, foundmin, foundmax
   
        ! create transition matrix for log idiosyncratic labor shock using Tauchen 86
        w=znodes(2)-znodes(1)
        do j=1, zgridsize
            Probz(j, 1)=cdfnormal((znodes(1)-rho_z*znodes(j)+w/2)/(sigma_z))
            Probz(j, zgridsize)=1-cdfnormal((znodes(zgridsize)-rho_z*znodes(j)-w/2)/(sigma_z))
            do k=2, zgridsize-1
                Probz(j, k)=cdfnormal((znodes(k)-rho_z*znodes(j)+w/2)/(sigma_z))-cdfnormal((znodes(k)-rho_z*znodes(j)-w/2)/(sigma_z))
            end do
        end do

    !Truncating the transistion probabilities of transitory shock. Want to
    !minimize number of gridpoints have to compute conditional expectation
        minexpectationz=1
        maxexpectationz=zgridsize
        do j=1, zgridsize
            foundmin=0.0
            foundmax=0.0
            do k=1, zgridsize
                if (Probz(j, k)>.01) then
                    foundmin=1.0
                elseif (foundmin==0.0) then
                    minexpectationz(j)=minexpectationz(j)+1
                end if
            end do
            do k=0, zgridsize-1
                if (Probz(j, zgridsize-k)>.01) then
                    foundmax=1.0
                elseif (foundmax==0.0) then
                    maxexpectationz(j)=maxexpectationz(j)-1
                end if
            end do
        end do

        do j=1, zgridsize
            do k=1, zgridsize
                if (k<minexpectationz(j)) then
                    !write(*,*) Probz(j, k)
                    Probz(j, k)=0.0
                end if
                if (k>maxexpectationz(j)) then
                    !write(*,*) Probz(j, k)
                    Probz(j, k)=0.0
                end if
            end do
        end do

        do j=1, zgridsize
            Probz(j,:)=Probz(j,:)/sum(Probz(j,:))
        end do
        do j=1, zgridsize
            do k=1, zgridsize
                Probzcum(j, k)=sum(Probz(j, 1:k))
            end do
        end do

    ! Following line pops out stationary distribution.
        Probinit = ergodic(Probz, zgridsize)
        do j=1, zgridsize
            Probinitcum(j)=sum(Probinit(1:j))
        end do
        write(*,*) Probinitcum
    end subroutine ! %>
   
    ! TODO: Why not have this subroutine contain any calibration loops?
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
        ALLOCATE(householdtransitionholder(numhouseholds,TDie+1,9,Tretire))

        ALLOCATE(housingstock(numhouseholds,Tdie,Tretire),&
        housingflow(numhouseholds,Tdie,Tretire),rentalflow(numhouseholds,Tdie,Tretire),&
        bequestflow(numhouseholds,Tdie,Tretire),resaleflow(numhouseholds,Tdie,Tretire))

        call PE_vfunc

        ! Broadcasting DP outputs forward a period for consistency in the
        ! transition subroutine
        achoice(:,:,:,2,:) = achoice(:,:,:,1,:)
        Dchoice(:,:,:,2,:) = Dchoice(:,:,:,1,:)
        cchoice(:,:,:,2,:) = cchoice(:,:,:,1,:)
        choiceindicator(:,:,:,2,:) = choiceindicator(:,:,:,1,:)
        rentalindicator(:,:,:,2,:) = rentalindicator(:,:,:,1,:)

        call PE_simulate

        write(*,*) "C sum", sum(cchoice(:,:,:,1:hptransLength,:))
        DEALLOCATE(achoice, Dchoice, cchoice, choiceindicator,rentalindicator, EV)

        !parameters are updated, then repeat whole big loop with solution/simulation to target moments.
        close(1)
    end subroutine ! %>

    subroutine PE_vfunc ! %<
        IMPLICIT NONE
        write(*,*) "Begin value function iteration:"

        call system_clock(timestart, timeres)
        thetanodes = theta
        ! Agents face the competitive price for rental housing
        r_rental = rentPrem + exp(hpnodes)*(1.0 - (1-delta)/(1+r))
        call bequests(hpnodes)
        call solveworkingproblem(hpnodes, achoice, Dchoice, cchoice, choiceindicator, &
            rentalindicator, EV, .TRUE.)
        call output_vfuncs(1)
        write(*,*) "Value function iteration complete:"
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"


        ! Constructs down payment at market price (since we interpret changes in
        ! hp as subsidies). If downtransfer, down subtracted by value of subsidy.
        write(*,*) "Begin value function iteration with subsidy:"
        call system_clock(timestart, timeres)
        ! New price level (proportionate change from baseline)
        adjsubsidy_internal   = adjsubsidy 
        rentsubsidy_internal  = rentsubsidy
        adjTransfer_internal  = adjTransfer
        rentTransfer_internal = rentTransfer 
        noadjTransfer_internal= noadjTransfer
        hpshock = exp(adjTransfer_internal - hpnodes(hptransLength)) ! New price level
        thetanodes = theta / hpshock - eta * (1 - hpshock) / hpshock
        ! The percentage/log change in prices from baseline
        write(*,*) adjTransfer, thetanodes

        call solveworkingproblem(hpnodes, achoiceshock, Dchoiceshock, cchoiceshock, &
            choiceindicatorshock, rentalindicatorshock, EV, .FALSE.)
        call plotpolicy_constantwealth(1, achoiceshock, Dchoiceshock, cchoiceshock, &
                                       choiceindicatorshock)
        write(*,*) "Value function iteration complete:"
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"
    end subroutine ! %>

    subroutine PE_simulate ! %<
        IMPLICIT NONE

        !simulate model and check match to some target moments:
        write(*,*) "Simulation starting..."
        call system_clock(timestart, timeres)
        call simulate(hpnodes(1),numhouseholds,householdtransitionholder(:,:,:,1))
        write(*,*) SUM(housingstock(:,:,1)), SUM(rentalflow(:,:,1))

        write(*, *) "Solving transition dynamics"
        do l=1,Tretire
            TransTime = l - 1
            write(*,*) "Transition period: ", TransTime
            call transition(hpnodes(1),numhouseholds,TransTime,&
                householdtransitionholder(:,:,:,l), householdtransitionholder(:,:,:,l+1))
        end do
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

        close(299)
        close(499)
        close(799)
        close(999)

    end subroutine ! %>

    subroutine GE_transition_path(pstart, yend) ! %<
        IMPLICIT NONE
        REAL(8), DIMENSION(2,hptransLength), INTENT(INOUT) :: pstart
        REAL(8), DIMENSION(hptranslength), INTENT(OUT) :: yend
        REAL(8), PARAMETER :: pthres=0.1
        REAL(8), DIMENSION(hptransLength-1) :: priceforward

        ! Agent state holders in each (transition) period
        ALLOCATE(householdtransitionholder(numhouseholds,TDie+1,9,hptransLength))

        ALLOCATE(housingstock(numhouseholds,Tdie,hptransLength),&
        housingflow(numhouseholds,Tdie,hptransLength),rentalflow(numhouseholds,Tdie,hptransLength),&
        bequestflow(numhouseholds,Tdie,hptransLength),resaleflow(numhouseholds,Tdie,hptransLength))

        call system_clock(timestart, timeres)
        call GE_steady_state(startprice(:,1), diffH(1))
        write(*,*) "Steady-state agent distribution computed"
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"
        ! TODO: Check this
        hpnodes(2:hptransLength+1)=hpnodes(1)
        call bequests(hpnodes)
        EV(:, :, :, hptransLength+1, :) = EV(:, :, :, 1, :)

        call simulate(hpnodes(1), numhouseholds, householdtransitionholder(:,:,1:8,1))

        do while (minval(priceforward - hpnodes(2:hptransLength)) < pthres)

            ! Backwards induction: EV grids on the (j, t) dimension are updated.
            ! Agents face the competitive price for rental housing
            call GE_transition_vfunc

            ! Forward induction: simulate transition path one period at a time
            ! and recalibrate the market clearing price
            call GE_transition_clearing(startprice(:,2:hptranslength),&
                                        diffH(2:hptranslength), pthres,&
                                        priceforward)

        end do

    end subroutine ! %>

    subroutine GE_steady_state(pstart, yend) ! %<
        IMPLICIT NONE
        REAL(8), DIMENSION(2,1), INTENT(INOUT) :: pstart ! One price in GE for the moment
        REAL(8), INTENT(OUT) :: yend

        ! NOTE: the function SimSteady includes the backwards induction calculation
        yend = brentmindist(pstart(1,1), pstart(1,1),pstart(2,1),SimSteady,2500,&
                            1e-3,hpnodes(1))
        write(*,*) hpnodes(1)

    end subroutine ! %>

    subroutine GE_transition_vfunc ! %<
        IMPLICIT NONE

        call system_clock(timestart, timeres)
        thetanodes = theta
        r_rental = rentPrem + exp(hpnodes)*(1.0 - (1-delta)/(1+r))
        ! Recall l=1 is the steady-state price.
        call solveworkingproblem(hpnodes(2:hptransLength+1), achoice(:,:,:,2:hptransLength,:), &
            Dchoice(:,:,:,2:hptransLength,:), cchoice(:,:,:,2:hptransLength,:), &
            choiceindicator(:,:,:,2:hptransLength,:), rentalindicator(:,:,:,2:hptransLength,:), &
            EV(:,:,:,2:hptransLength,:), .TRUE.)
        call output_vfuncs(1)
        call output_vfuncs(hptransLength)
        write(*,*) "Value function iteration complete:"
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

        ! Constructs down payment at market price (since we interpret changes in
        ! hp as subsidies). If downtransfer, down subtracted by value of subsidy.
        write(*,*) "Begin value function iteration with subsidy:"
        call system_clock(timestart, timeres)
        ! New price level (proportionate change from baseline)
        adjsubsidy_internal   = adjsubsidy 
        rentsubsidy_internal  = rentsubsidy
        adjTransfer_internal  = adjTransfer
        rentTransfer_internal = rentTransfer 
        noadjTransfer_internal= noadjTransfer
        hpshock = exp(adjTransfer_internal - hpnodes(hptransLength)) ! New price level
        thetanodes = theta / hpshock - eta * (1 - hpshock) / hpshock
        ! The percentage/log change in prices from baseline
        write(*,*) adjTransfer, thetanodes

        call solveworkingproblem(hpnodes(2:2), achoiceshock, Dchoiceshock, cchoiceshock, &
            choiceindicatorshock, rentalindicatorshock, EV, .FALSE.)
        call plotpolicy_constantwealth(1, achoiceshock, Dchoiceshock, cchoiceshock, &
                                       choiceindicatorshock)
        write(*,*) "Value function iteration complete:"
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"
     end subroutine ! %>


    subroutine GE_transition_clearing(pstart, yend, pthres, priceforward) ! %<
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: pthres
        REAL(8), DIMENSION(2,hptransLength-1), INTENT(INOUT) :: pstart
        REAL(8), DIMENSION(hptranslength-1), INTENT(OUT) :: yend
        REAL(8), DIMENSION(hptransLength-1), INTENT(OUT) :: priceforward
        INTEGER :: l

        ! TODO: Do I want to run solveworkingproblem for one transition period
        ! in SimTransPath? Could very well be what Stroebel does
        do l=1,hptransLength-1
            TransTime = l - 1
            yend(l) = brentmindist(hpnodes(l+1), pstart(1,l),pstart(2,l),SimTransPath,&
                numhouseholds,5e-2,priceforward(l))
        ! Compare the initial price guesses and the recalibrated price
            if (priceforward(l) - hpnodes(l+1) >= pthres) then
                hpnodes(l+1) = 0.8*priceforward(l) + 0.2*hpnodes(l+1)
            end if
        end do

     end subroutine ! %>


     
end program durables
