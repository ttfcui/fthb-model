program durables  !Program shell
use share
use lifeCycle_solveDP
use lifeCycle_simulate
use lifeCycle_transition
USE OMP_LIB
implicit none
integer :: timestart, timeend
REAL(8) :: timeres
real(8) :: age
REAL(8) :: w, ll, rr, ww, niter, wholder
REAL(8) :: foundmin, foundmax
REAL(8) :: dd, m2, m2young, m2middle, m2old
REAL(8) :: finwealth

REAL(8), DIMENSION(2,1) :: startprice
REAL(8) :: eqprice, diffH

! EXTERNAL cdfnormal  !Fnc to compute normal cdf
! INTERFACE ! Stationary dist. function for Markov
!     FUNCTION ergodic(mat, states)
!     IMPLICIT NONE
!     INTEGER, intent(in) :: states
!     REAL(8), DIMENSION(states, states), intent(in) :: mat
!     REAL(8), DIMENSION(states) :: ergodic
!     END FUNCTION
! END INTERFACE

r_rental = r_rentalO
r_rental_retire=r_rental_retireO
elasticity2=elasticity2O
beta2=beta2O
beta2retire=beta2retireO


write(*,*) "amin", amin




ALLOCATE(Vnoadjust(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie+1),  Vadjust(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie+1), Vrent(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie+1))  !Conjectured value of not adjusting and adjusting, and new values of not adjusting and adjusting
ALLOCATE(achoiceadjust(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie),  achoicenoadjust(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie), achoicerent(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie),  Dchoiceadjust(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie), Dchoicenoadjust(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie), Dchoicerent(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie))
ALLOCATE(achoice(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie), Dchoice(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie), cchoice(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie))  !Policy functions when solving problem
ALLOCATE(cchoiceadjust(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie), cchoicenoadjust(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie), cchoicerent(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie))
ALLOCATE(rentalindicator(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie), choiceindicator(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie))
ALLOCATE(EV(agridsize, Dgridsize, zgridsize, hpgridsize, Tdie+1))

write(*,*) "Max threads:", omp_get_max_threads()

numiter=0




diffmoments=1
do while (diffmoments>difftol)
    numiter=numiter+1
    write(*,*) "numiter", numiter

    write(*,*) elasticity2, r_rental, beta2


    !OPEN (UNIT=1, FILE="vfunc.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    !1 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

    ! FILL IN ALL THE GRID POINTS

    ! Shocks grid - equal intervals a la Tauchen discretization
    write(*,*) "Shocks grid:"
    do i=1, zgridsize
        znodes(i)=zmin+((zmax-zmin)/(zgridsize-1))*(1.0*i-1.0)
        write(*,*) znodes(i)
    end do


    write(*,*) "Asset search grid:"
    DO i=1, agridsize
        anodes(i)=(1.0/(agridsize-1))*(1.0*i-1.0)
    end do

    ! Transform the grid so more mass with small numbers.
    ! The log is a normalizing factor.
    do i=1, agridsize
        anodes(i)=exp(log(amax-amin+1)*anodes(i))+amin-1.0
        write(*,*) anodes(i)
    end do


    write(*,*) "Durable search grid:"
     DO i=1, Dgridsize
        Dnodes(i)=(1.0/(dgridsize-1))*(1.0*i-1.0)
    end do


    do i=1, dgridsize
        Dnodes(i)=exp(log(dmax-dmin+1)*Dnodes(i))+dmin-1.0
        write(*,*) Dnodes(i)
    end do

    call share_arrays

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
! Not sure what old passage meant.
    Probinit = ergodic(Probz, zgridsize)
!    w=znodes(2)-znodes(1)
!    Probinit(1)=cdfnormal((znodes(1)+w/2)/sigma_eta_init)
!    Probinit(zgridsize)=1-cdfnormal((znodes(zgridsize)-w/2)/sigma_eta_init)
!    do j=2, zgridsize-1
!        Probinit(j)=cdfnormal((znodes(j)+w/2)/sigma_eta_init)-cdfnormal((znodes(j)-w/2)/sigma_eta_init)
!    end do

    do j=1, zgridsize
        Probinitcum(j)=sum(Probinit(1:j))
    end do
write(*,*) Probinitcum

! Chi/permanent portion of earnings proccess
! This is a cubic of log income over age from the
! Blundell, Pistaferri, Preston data, normalized;
! See BPP_calibration.do for more details

ageearnings(1)  = -0.3522492051
ageearnings(2)  = -0.2265817374
ageearnings(3)  = -0.1947079152
ageearnings(4)  = -0.1652790904
ageearnings(5)  = -0.1382162571
ageearnings(6)  = -0.1134403795
ageearnings(7)  = -0.0908724517
ageearnings(8)  = -0.0704334378
ageearnings(9)  = -0.0520443171
ageearnings(10) = -0.0356260762
ageearnings(11) = -0.0210996903
ageearnings(12) = -0.0083861314
ageearnings(13) =  0.0025936160
ageearnings(14) =  0.0119185783
ageearnings(15) =  0.0196677744
ageearnings(16) =  0.0259202309
ageearnings(17) =  0.0307549611
ageearnings(18) =  0.0342509970
ageearnings(19) =  0.0364873558
ageearnings(20) =  0.0375430509
ageearnings(21) =  0.0374971256
ageearnings(22) =  0.0364285782
ageearnings(23) =  0.0344164483
ageearnings(24) =  0.0315397494
ageearnings(25) =  0.0278775059
ageearnings(26) =  0.0235087387
ageearnings(27) =  0.0185124725
ageearnings(28) =  0.0129677244
ageearnings(29) =  0.0069535188
ageearnings(30) =  0.0005488787
ageearnings(31) = -0.0061671734
ageearnings(32) = -0.0131156184
ageearnings(33) = -0.0202174317
ageearnings(34) = -0.0273935925
ageearnings(35) = -0.0345650800
ageearnings(36) = -0.0416528694
ageearnings(37) = -0.0485779382
ageearnings(38) = -0.0552612692
!ageearnings(39) =
!ageearnings(40) =

!    ageearnings=0

! death probabilities (starting at age 60 or 65)
    deathprob(1)=.012
    deathprob(2)=.013
    deathprob(3)=.014
    deathprob(4)=.015
    deathprob(5)=.016
    deathprob(6)=.018
    deathprob(7)=.02
    deathprob(8)=.021
    deathprob(9)=.022
    deathprob(10)=.024
    deathprob(11)=.027
    deathprob(12)=.03
    deathprob(13)=.032
    deathprob(14)=.035
    deathprob(15)=.038
    deathprob(16)=.042
    deathprob(17)=.046
    deathprob(18)=.051
    deathprob(19)=.056
    deathprob(20)=.061
    deathprob(21)=.068
    deathprob(22)=.075
    deathprob(23)=.083
    deathprob(24)=.092
    deathprob(25)=.11

!    deathprob=0



    averagelifetimeearnings=0.0
    !retirement regression needs to first run matlab program simulateearningsprocess
    do j=1, zgridsize
        predictedlifetimeearnings(j)=0.27*(znodes(j))
        predictedlifetimeearnings(j)=exp(predictedlifetimeearnings(j))/exp(averagelifetimeearnings)
    ! This function is from Guvenen and Smith 2014
        if (predictedlifetimeearnings(j)<=0.3) then
            retirementincome(j)=0.9*predictedlifetimeearnings(j)*exp(averagelifetimeearnings)
        elseif (predictedlifetimeearnings(j)>0.3 .and. predictedlifetimeearnings(j)<=2) then
            retirementincome(j)=(0.27+0.32*(predictedlifetimeearnings(j)-0.3))*exp(averagelifetimeearnings)
        elseif (predictedlifetimeearnings(j)>2 .and. predictedlifetimeearnings(j)<=4.1) then
            retirementincome(j)=(0.81+0.15*(predictedlifetimeearnings(j)-2))*exp(averagelifetimeearnings)
        else
            retirementincome(j)=1.13*exp(averagelifetimeearnings)
        end if

      !  retirementincome(j)=exp(znodes(j))
      !  retirementincome(j) = 0.0
        write(*,*) "r i", retirementincome(j)
    end do


    !retirementincome=exp(ageearnings(35))
    !retirementincome=0

    ! Populate income process array
    do j=1, zgridsize
        do t=1, Tretire
            income(j, t)=exp(znodes(j)+ageearnings(t))
        end do
        do t=Tretire+1, Tdie
           income(j, t)=retirementincome(j)
        end do
      !  do t=Tretire+1, Tdie
      !      income(j, t)=exp(znodes(j))
      !  end do
    end do
    do t=1, Tdie
        write(*,*) t, income(1, t), income(13, t)
    end do



   ! write(*,*) "income"
  !  write(*,*) sum(income(:, 1,:)), sum(income(:, 2,:)), sum(income(:, 3,:)), sum(income(:, 4,:)), sum(income(:, 5,:))

    ! Constructs down payment at market price (since we interpret changes in
    ! hp as subsidies). If downtransfer, down subtracted by value of subsidy.
    do i=1, hpgridsize-1
        hpdelta(i) = exp(hpnodes(i) - hpnodes(hpgridsize)) ! New price level
        thetanodes(i) = theta / hpdelta(i) - eta(i) * (1 - hpdelta(i)) / hpdelta(i)
        hpdelta(i) = exp(hpnodes(i)) - exp(hpnodes(hpgridsize)) ! The actual change in prices from baseline
    end do
    hpdelta(hpgridsize) = 0
    thetanodes(hpgridsize) = theta
    write(*,*) hpdelta
    write(*,*) thetanodes




! Put in bequest motive in last period

    call system_clock(timestart, timeres)
    !call solveretirementproblem
    !Vnoadjust(:,:,:,:, 61)=0
    !Vadjust(:,:,:,:, 61)=0
   ! Vrent(:,:,:,:, 61)=0

    do i=1, agridsize
        do j=1, Dgridsize
            do k=1, hpgridsize
                finwealth=anodes(i)-(1-theta)*Dnodes(j)*exp(hpnodes(k))
                wholder=(1+r)*finwealth+(1-delta)*Dnodes(j)*exp(hpnodes(k))*(1-f)

                !wholder=(1+r)*(anodes(i)-(1-theta)*Dnodes(j)*(1-f))+(1-delta)*Dnodes(j)*(1-f)

                !wholder=(1-delta)*Dnodes(j)*exp(hpnodes(k))
                !wholder=(1+r)*(anodes(i)+theta*exp(hpnodes(k))*Dnodes(j)*(1-f))-delta*exp(hpnodes(k))*Dnodes(j)*(1-f)


                if (elasticity .ne. 1.0) then
                !EV(i, j,:, k, Tdie+1)=1.0/(1.0-elasticity)*rent**(elasticity-1)*psi*(income(:, Tdie)/r+wholder)**(1-elasticity)
                EV(i, j,:, k, Tdie+1)=1.0/(1.0-elasticity)*psi*(income(:, Tdie)/r+wholder)**(1-elasticity) ! (still don't know why income is divide by r. Gideon)
                else
                EV(i, j,:, k, Tdie+1)=psi*log(income(:, Tdie)/r+wholder)
                end if
                if (j==1 .and. k==1) then
                write(*,*) EV(i, j, 1, k, Tdie+1)
                end if

                !if (finwealth<0) then
                !    EV(i, j,:, k, Tdie+1)=1/((1-elasticity))*((finwealth*(1+rborrow)+(1-delta)*Dnodes(j)*exp(hpnodes(k))*(1-F)))**(1-elasticity)
                !else
                !    EV(i, j,:, k, Tdie+1)=1/((1-elasticity))*((finwealth*(1+r)+(1-delta)*Dnodes(j)*exp(hpnodes(k))*(1-F)))**(1-elasticity)
                !end if
                !if (anodes(i)-(1-theta)*Dnodes(j)*exp(hpnodes(k))<0) then
                !    EV(i, j,:, k, Tdie+1)=1/(1-beta2)*((1+rborrow)*anodes(i)+Dnodes(j)*(1-F)*(1-delta)*exp(hpnodes(k))-(1+rborrow)*(1-theta)*Dnodes(j)*(1-F)*exp(hpnodes(k)))**(1-elasticity)/(1-elasticity)
                !else
                !    EV(i, j,:, k, Tdie+1)=1/(1-beta2)*((1+r)*anodes(i)+Dnodes(j)*(1-F)*(1-delta)*exp(hpnodes(k))-(1+r)*(1-theta)*Dnodes(j)*(1-F)*exp(hpnodes(k)))**(1-elasticity)/(1-elasticity)
                !end if
            end do
        end do
    end do
    !EV=0
    !EV=max(-10000000000.0, EV)

    write(*,*) "totEV ", sum(EV)

    !EV=0



    call system_clock(timeend, timeres)
    write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

    !call PE_price_subsidy
    startprice(1,1) = 0.2
    startprice(2,1) = 0.19
    call GE_steady_state(startprice, diffH, eqprice)
          
end do


write(*,*) "C sum", sum(cchoice)


write(*,*) "model version"
write(*,*) "F  |  ", "theta  |  ", "r_rental  |  ", "costlyequity  |  ", "rentelasticity  |  ", "hpmin"
write(*,*) F, theta, r_rental, costlyequity, rentelasticity, hpmin
write(*,*) "A grid  |  ", "Dur grid  |  ", "Z grid  |  ", "HPrice grid  |  ", "households"
write(*,*) agridsize, Dgridsize, zgridsize, hpgridsize, numhouseholds


   CONTAINS


! Main model code: Solves the backward induction
    subroutine PE_price_subsidy
        IMPLICIT NONE
        write(*,*) "Begin value function iteration:"
        call system_clock(timestart, timeres)
        call solveworkingproblem(.FALSE., 0.0)
        write(*,*) "Value function iteration complete:"
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

        call plotpolicy_constantwealth(2.0)
        call output_vfuncs(.FALSE., 0.0)

        !simulate model and check match to some target moments:
        write(*,*) "Simulation starting..."
        call system_clock(timestart, timeres)
        call simulate(numhouseholds)

        call transition
        call system_clock(timeend, timeres)
        write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

        diffmoments=0.0
        !parameters are updated, then repeat whole big loop with solution/simulation to target moments.
        close(1)
    end subroutine

    subroutine GE_steady_state(pstart, yend, xend)
        IMPLICIT NONE
        REAL(8), DIMENSION(2,1), INTENT(INOUT) :: pstart ! One price in GE for the moment
        REAL(8), INTENT(OUT) :: yend, xend

        yend = brentmindist(pstart(1,1), pstart(1,1),pstart(2,1),SimSteady,1000,1e-5,xend)
        write(*,*) xend
        diffmoments=0.0

     end subroutine
 

end program durables
