!***********************************************************
!    This file is different from 04032016 because
!    It simulates a policy like the CARS program, not FTHB.
!    Program eligibility is changed, and the code tracks
!    transactions in a period, not "FTHB status."
!    Other unnecessary code from earlier projects have
!    also been removed.
!***********************************************************


module share
    ! Basic model parameters
    REAL(8), parameter :: scalefactor=1.0
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
    REAL(8) :: r_rental, r_rental_retire, elasticity2, beta2, beta2retire

    ! The asset and durable grids
    REAL(8), parameter :: Dmin=0
    REAL(8), parameter :: Dmax=6

    !REAL(8), parameter :: amin=(exp(hpmin)-exp(hpmax))*(1-theta)*Dmax    !minimum asset value (actually voluntary equity if theta!=1).  Note that if theta<1, we can have assets <0 even with amin=0 since a is vol. equity.  But with theta=1 & amin=0, no borrowing.
    !REAL(8), parameter :: amin=-.89
    REAL(8), parameter :: amin=-0.40/(1+r)  !.40 is roughly min y, so this is basically a no default condition
    !REAL(8), parameter :: amax2=8   !max asset
    REAL(8), parameter :: amax=40   !max asset

    integer, parameter :: agridsize = 60
    integer, parameter :: Dgridsize = 45

    !idiosyncratic earnings shocks
    REAL(8), parameter :: sigma_z= .21
    REAL(8), parameter :: rho_z=.91
    REAL(8), parameter :: sigma_eta_init=(sigma_z**2.0/(1-rho_z**2.0))**.5
    REAL(8), parameter :: zmin=-2.5*sigma_eta_init
    REAL(8), parameter :: zmax=2.5*sigma_eta_init

    integer, parameter :: zgridsize=13  !size of ar(1) income shock grid
    REAL(8), parameter :: borrowconstraint=0.0

    !for now, just two point process
    REAL(8), parameter :: hpmin=-0.15
    REAL(8), parameter :: hpmax=0.0
    integer, parameter :: hpgridsize=2 !size of house price grid
    REAL(8), parameter :: downtransfer=1 ! 0-1: Smoothing how much subsidy applies to down
    integer, parameter :: PolYrs=1 ! Number of years subsidy is active
    integer, parameter :: EligYrs=3 ! Years under which agent is in a certain state for policy eligibility
    REAL(8), DIMENSION(hpgridsize) :: hpnodes, thetanodes
    REAL(8), DIMENSION(hpgridsize,hpgridsize) :: Probhp

    !integer, parameter :: Tretire=35
    !integer, parameter :: Tdie=60

    !integer, parameter :: Tretire=100
    integer, parameter :: Tretire=40 !35
    integer, parameter :: Tdie=Tretire+25


    ! Below are the state and policy holders


    REAL(8), DIMENSION(agridsize) :: anodes  !Nodes for idiosyncratic assets (voluntary equity)
    REAL(8), DIMENSION(Dgridsize) :: Dnodes  !Nodes for idiosyncratic durable

    REAL(8), DIMENSION(zgridsize) :: Probinit,Probinitcum



    REAL(8), DIMENSION(zgridsize) :: znodes  !Nodes for idiosyncratic productivity
    REAL(8), DIMENSION(zgridsize,zgridsize) :: Probz,Probzcum  !while probability of being at any z will vary with time, transition probability will not, so max we need to loop over wont
    INTEGER, DIMENSION(zgridsize) :: maxexpectationz,minexpectationz


    REAL(8), DIMENSION(Tretire) :: ageearnings
    REAL(8), DIMENSION(Tdie-Tretire) :: deathprob
    REAL(8), DIMENSION(zgridsize) :: predictedlifetimeearnings,retirementincome
    REAL(8), DIMENSION(zgridsize,Tdie) :: income
    REAL(8) :: averagelifetimeearnings

    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Vnoadjust,  Vadjust,Vrent  !Conjectured value of not adjusting and adjusting, and new values of not adjusting and adjusting
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: achoiceadjust,  achoicenoadjust,achoicerent,  Dchoiceadjust, Dchoicenoadjust,Dchoicerent, achoice, Dchoice,cchoice, cchoiceadjust, cchoicenoadjust, cchoicerent, cchoiceshock,Dchoiceshock  !Policy functions when solving problem
    REAL(8) :: priceshock

    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: rentalindicator, choiceindicator !choice indicator = 1 if adjust, 2 if noadjust, 3 if rent
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: EV


    INTEGER, parameter :: numhouseholds=8000   !Size of household panel
    REAL(8), dimension(numhouseholds,Tdie,11) :: householdresultsholder
    !integer, parameter :: burnin=1000    !Initial periods dropped in simulation


    REAL(8) :: diffmoments,diffmoments1,diffmoments2,diffmoments3
    REAL(8), parameter :: difftol=.01
    REAL(8) :: numiter

    REAL(8), parameter :: costlyequity=1

    REAL(8), parameter :: rentelasticity=1

end module share


program durables  !Program shell
use share
USE OMP_LIB
implicit none
integer :: i,j,k,l,m,n,t, iter
integer :: timestart, timeend
REAL(8) :: timeres
real(8) :: age
REAL(8) :: cdfnormal
REAL(8) :: w,ll,rr,ww,niter,wholder
REAL(8) :: foundmin, foundmax
REAL(8) :: dd,m2,m2young,m2middle,m2old
REAL(8) :: finwealth, hpdelta
REAL(8) :: i2,j2,k2,l2,t2

EXTERNAL cdfnormal  !Fnc to compute normal cdf
INTERFACE ! Stationary dist. function for Markov
    FUNCTION ergodic(mat, states)
    IMPLICIT NONE
    INTEGER, intent(in) :: states
    REAL(8), DIMENSION(states,states), intent(in) :: mat
    REAL(8), DIMENSION(states) :: ergodic
    END FUNCTION
END INTERFACE

r_rental = r_rentalO
r_rental_retire=r_rental_retireO
elasticity2=elasticity2O
beta2=beta2O
beta2retire=beta2retireO


write(*,*) "amin", amin

!Housing shocks
Probhp(1,1)=0.0
Probhp(1,2)=1.0
Probhp(2,1)=0.0
Probhp(2,2)=1.0



ALLOCATE(Vnoadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie+1),  Vadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie+1),Vrent(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie+1))  !Conjectured value of not adjusting and adjusting, and new values of not adjusting and adjusting
ALLOCATE(achoiceadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),  achoicenoadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),achoicerent(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),  Dchoiceadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie), Dchoicenoadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),Dchoicerent(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie))
ALLOCATE(achoice(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie), Dchoice(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),cchoice(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie))  !Policy functions when solving problem
ALLOCATE(cchoiceadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie), cchoicenoadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie), cchoicerent(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie))
ALLOCATE(rentalindicator(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie), choiceindicator(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie))
ALLOCATE(EV(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie+1))

write(*,*) "Max threads:", omp_get_max_threads()

numiter=0




diffmoments=1
do while (diffmoments>difftol)
    numiter=numiter+1
    write(*,*) "numiter", numiter

    write(*,*) elasticity2, r_rental, beta2


    !OPEN (UNIT=1, FILE="vfunc.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    !1 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

    ! FILL IN ALL THE GRID POINTS

    ! Shocks grid - equal intervals a la Tauchen discretization
    write(*,*) "Shocks grid:"
    do i=1,zgridsize
        znodes(i)=zmin+((zmax-zmin)/(zgridsize-1))*(1.0*i-1.0)
        write(*,*) znodes(i)
    end do

    ! House price grid - only two right now, so excessive
    write(*,*) "House price grid:"
    DO i=1,hpgridsize
        hpnodes(i)=hpmin+ ((hpmax -hpmin)/(hpgridsize-1))*(1.0*i-1.0)
    END DO


    write(*,*) "Asset search grid:"
    DO i=1,agridsize
        anodes(i)=(1.0/(agridsize-1))*(1.0*i-1.0)
    end do

    ! Transform the grid so more mass with small numbers.
    ! The log is a normalizing factor.
    do i=1,agridsize
        anodes(i)=exp(log(amax-amin+1)*anodes(i))+amin-1.0
        write(*,*) anodes(i)
    end do



    write(*,*) "Durable search grid:"
     DO i=1,Dgridsize
        Dnodes(i)=(1.0/(dgridsize-1))*(1.0*i-1.0)
    end do


    do i=1,dgridsize
        Dnodes(i)=exp(log(dmax-dmin+1)*Dnodes(i))+dmin-1.0
        write(*,*) Dnodes(i)
    end do


    ! create transition matrix for log idiosyncratic labor shock using Tauchen 86
    w=znodes(2)-znodes(1)
    do j=1,zgridsize
        Probz(j,1)=cdfnormal((znodes(1)-rho_z*znodes(j)+w/2)/(sigma_z))
        Probz(j,zgridsize)=1-cdfnormal((znodes(zgridsize)-rho_z*znodes(j)-w/2)/(sigma_z))
        do k=2,zgridsize-1
            Probz(j,k)=cdfnormal((znodes(k)-rho_z*znodes(j)+w/2)/(sigma_z))-cdfnormal((znodes(k)-rho_z*znodes(j)-w/2)/(sigma_z))
        end do
    end do

!Truncating the transistion probabilities of transitory shock. Want to
!minimize number of gridpoints have to compute conditional expectation
    minexpectationz=1
    maxexpectationz=zgridsize
    do j=1,zgridsize
        foundmin=0.0
        foundmax=0.0
        do k=1,zgridsize
            if (Probz(j,k)>.01) then
                foundmin=1.0
            elseif (foundmin==0.0) then
                minexpectationz(j)=minexpectationz(j)+1
            end if
        end do
        do k=0,zgridsize-1
            if (Probz(j,zgridsize-k)>.01) then
                foundmax=1.0
            elseif (foundmax==0.0) then
                maxexpectationz(j)=maxexpectationz(j)-1
            end if
        end do
    end do



    do j=1,zgridsize
        do k=1,zgridsize
            if (k<minexpectationz(j)) then
                !write(*,*) Probz(j,k)
                Probz(j,k)=0.0
            end if
            if (k>maxexpectationz(j)) then
                !write(*,*) Probz(j,k)
                Probz(j,k)=0.0
            end if
        end do
    end do

    do j=1,zgridsize
        Probz(j,:)=Probz(j,:)/sum(Probz(j,:))
    end do


    do j=1,zgridsize
        do k=1,zgridsize
            Probzcum(j,k)=sum(Probz(j,1:k))
        end do
    end do

! Following line pops out stationary distribution.
! Not sure what old passage meant.
    Probinit = ergodic(Probz, zgridsize)
!    w=znodes(2)-znodes(1)
!    Probinit(1)=cdfnormal((znodes(1)+w/2)/sigma_eta_init)
!    Probinit(zgridsize)=1-cdfnormal((znodes(zgridsize)-w/2)/sigma_eta_init)
!    do j=2,zgridsize-1
!        Probinit(j)=cdfnormal((znodes(j)+w/2)/sigma_eta_init)-cdfnormal((znodes(j)-w/2)/sigma_eta_init)
!    end do

    do j=1,zgridsize
        Probinitcum(j)=sum(Probinit(1:j))
    end do
write(*,*) Probinitcum

! Chi/permanent portion of earnings proccess
! These are binned log income means by age from the Blundell, Pistaferri, Preston
! data, then divided by the mean income from 20 to 60
! (so analogous to the normalization factor for everything else).

ageearnings(1) = -1.1077580452
ageearnings(2) = -0.9300475717
ageearnings(3) = -0.8037734628
ageearnings(4) = -0.6993461251
ageearnings(5) = -0.5463338494
ageearnings(6) = -0.4464326501
ageearnings(7) = -0.3789954782
ageearnings(8) = -0.2823281884
ageearnings(9) = -0.2092328668
ageearnings(10) =-0.1674123406
ageearnings(11) =-0.1048646569
ageearnings(12) =-0.0711455941
ageearnings(13) =-0.0178294759
ageearnings(14) = 0.0032276530
ageearnings(15) = 0.0625872016
ageearnings(16) = 0.0729317069
ageearnings(17) = 0.1144995093
ageearnings(18) = 0.1381124854
ageearnings(19) = 0.1558660865
ageearnings(20) = 0.1640648246
ageearnings(21) = 0.1893905997
ageearnings(22) = 0.1976723075
ageearnings(23) = 0.2042774558
ageearnings(24) = 0.1958660483
ageearnings(25) = 0.2035974860
ageearnings(26) = 0.2287039161
ageearnings(27) = 0.1494678855
ageearnings(28) = 0.1569246650
ageearnings(29) = 0.1292919517
ageearnings(30) = 0.1157287955
ageearnings(31) = 0.1053223014
ageearnings(32) = 0.1010994315
ageearnings(33) = 0.0841526389
ageearnings(34) = 0.0717768073
ageearnings(35) = 0.0287908930
ageearnings(36) = 0.0317749418
ageearnings(37) =-0.0166640859
ageearnings(38) =-0.0598741136
ageearnings(39) =-0.0866485238
ageearnings(40) =-0.1181903481

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
    do j=1,zgridsize
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
    do j=1,zgridsize
        do t=1,Tretire
            income(j,t)=exp(znodes(j)+ageearnings(t))
        end do
        do t=Tretire+1,Tdie
           income(j,t)=retirementincome(j)
        end do
      !  do t=Tretire+1,Tdie
      !      income(j,t)=exp(znodes(j))
      !  end do
    end do
    do t=1,Tdie
        write(*,*) t, income(1,t), income(13,t)
    end do



   ! write(*,*) "income"
  !  write(*,*) sum(income(:,1,:)),sum(income(:,2,:)),sum(income(:,3,:)),sum(income(:,4,:)),sum(income(:,5,:))

    ! Constructs down payment at market price (since we interpret changes in
    ! hp as subsidies). If downtransfer, down subtracted by value of subsidy.
    hpdelta = exp(hpnodes(1) - hpnodes(2))
    if ((downtransfer >= 0.0) .and. (downtransfer <= 1.0)) then
        thetanodes(1) = theta / hpdelta - downtransfer * (1 - hpdelta) / hpdelta
    else
        thetanodes(1) = theta
    end if
    thetanodes(2) = theta
    write(*,*) hpdelta
    write(*,*) thetanodes




! Put in bequest motive in last period

    call system_clock(timestart, timeres)
    !call solveretirementproblem
    !Vnoadjust(:,:,:,:,61)=0
    !Vadjust(:,:,:,:,61)=0
   ! Vrent(:,:,:,:,61)=0

    do i=1,agridsize
        do j=1,Dgridsize
            do k=1,hpgridsize
                finwealth=anodes(i)-(1-thetanodes(k))*Dnodes(j)*exp(hpnodes(k))
                wholder=(1+r)*finwealth+(1-delta)*Dnodes(j)*exp(hpnodes(k))*(1-f)

                !wholder=(1+r)*(anodes(i)-(1-theta)*Dnodes(j)*(1-f))+(1-delta)*Dnodes(j)*(1-f)

                !wholder=(1-delta)*Dnodes(j)*exp(hpnodes(k))
                !wholder=(1+r)*(anodes(i)+theta*exp(hpnodes(k))*Dnodes(j)*(1-f))-delta*exp(hpnodes(k))*Dnodes(j)*(1-f)


                if (elasticity .ne. 1.0) then
                !EV(i,j,:,k,Tdie+1)=1.0/(1.0-elasticity)*rent**(elasticity-1)*psi*(income(:,Tdie)/r+wholder)**(1-elasticity)
                EV(i,j,:,k,Tdie+1)=1.0/(1.0-elasticity)*psi*(income(:,Tdie)/r+wholder)**(1-elasticity) ! (still don't know why income is divide by r. Gideon)
                else
                EV(i,j,:,k,Tdie+1)=psi*log(income(:,Tdie)/r+wholder)
                end if
                if (j==1 .and. k==1) then
                write(*,*) EV(i,j,1,k,Tdie+1)
                end if

                !if (finwealth<0) then
                !    EV(i,j,:,k,Tdie+1)=1/((1-elasticity))*((finwealth*(1+rborrow)+(1-delta)*Dnodes(j)*exp(hpnodes(k))*(1-F)))**(1-elasticity)
                !else
                !    EV(i,j,:,k,Tdie+1)=1/((1-elasticity))*((finwealth*(1+r)+(1-delta)*Dnodes(j)*exp(hpnodes(k))*(1-F)))**(1-elasticity)
                !end if
                !if (anodes(i)-(1-theta)*Dnodes(j)*exp(hpnodes(k))<0) then
                !    EV(i,j,:,k,Tdie+1)=1/(1-beta2)*((1+rborrow)*anodes(i)+Dnodes(j)*(1-F)*(1-delta)*exp(hpnodes(k))-(1+rborrow)*(1-theta)*Dnodes(j)*(1-F)*exp(hpnodes(k)))**(1-elasticity)/(1-elasticity)
                !else
                !    EV(i,j,:,k,Tdie+1)=1/(1-beta2)*((1+r)*anodes(i)+Dnodes(j)*(1-F)*(1-delta)*exp(hpnodes(k))-(1+r)*(1-theta)*Dnodes(j)*(1-F)*exp(hpnodes(k)))**(1-elasticity)/(1-elasticity)
                !end if
            end do
        end do
    end do
    !EV=0
    !EV=max(-10000000000.0,EV)

    write(*,*) "totEV ", sum(EV)

    !EV=0


    OPEN (UNIT=1, FILE="vfunc1.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    1 format (F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12, F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12,F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12)

    OPEN (UNIT=17, FILE="vfunc2.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    17 format (F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12, F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12,F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12)

    OPEN (UNIT=18, FILE="vfunc3.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    18 format (F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12, F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12,F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12)

    OPEN (UNIT=19, FILE="vfunc4.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    19 format (F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12, F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12,F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12)


   ! Main model code: Solves the backward induction
    call system_clock(timeend, timeres)
    write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"
    write(*,*) "Begin value function iteration:"
    call system_clock(timestart, timeres)
    call solveworkingproblem
    write(*,*) "Value function iteration complete:"
    call system_clock(timeend, timeres)
    write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

    call plotpolicy_constantwealth

! outputing the policy functions
    write(*,*) agridsize,Dgridsize,zgridsize,hpgridsize
    t=0
    do i=1,agridsize
        do j=1,Dgridsize
            do k=1,zgridsize
                do l=1,hpgridsize
                t=t+1
                t2=33
                k2=k
                l2=l
                    !call computempc(anodes(i),Dnodes(j),k2,l2,t2,mpc_pol(i,j,k,l,33))
                    write(1,1) EV(i,j,k,l,33), achoice(i,j,k,l,33), dchoice(i,j,k,l,33), cchoice(i,j,k,l,33),rentalindicator(i,j,k,l,33), anodes(i), Dnodes(j), znodes(k),hpnodes(l), income(k,1) !, mpc_pol(i,j,k,l,33) , Vadjust(i,j,k,l,1), Vnoadjust(i,j,k,l,1)
                end do
            end do
        end do
  end do

  t=0
    do i=1,agridsize
        do j=1,Dgridsize
            do k=1,zgridsize
                do l=1,hpgridsize
                t=t+1
                    write(17,17) EV(i,j,k,l,34), achoice(i,j,k,l,34), dchoice(i,j,k,l,34), cchoice(i,j,k,l,34),rentalindicator(i,j,k,l,34), anodes(i), Dnodes(j), znodes(k),hpnodes(l), income(k,1) !, Vadjust(i,j,k,l,1), Vnoadjust(i,j,k,l,1)
                end do
            end do
        end do
  end do

  t=0
    do i=1,agridsize
        do j=1,Dgridsize
            do k=1,zgridsize
                do l=1,hpgridsize
                t=t+1
                    write(18,18) EV(i,j,k,l,35), achoice(i,j,k,l,35), dchoice(i,j,k,l,35), cchoice(i,j,k,l,35),rentalindicator(i,j,k,l,35), anodes(i), Dnodes(j), znodes(k),hpnodes(l), income(k,1) !, Vadjust(i,j,k,l,1), Vnoadjust(i,j,k,l,1)
                end do
            end do
        end do
  end do

  t=0
    do i=1,agridsize
        do j=1,Dgridsize
            do k=1,zgridsize
                do l=1,hpgridsize
                t=t+1
                    write(19,19) EV(i,j,k,l,36), achoice(i,j,k,l,36), dchoice(i,j,k,l,36), cchoice(i,j,k,l,36),rentalindicator(i,j,k,l,36), anodes(i), Dnodes(j), znodes(k),hpnodes(l), income(k,1) !, Vadjust(i,j,k,l,1), Vnoadjust(i,j,k,l,1)
                end do
            end do
        end do
  end do



  write(*,*) t

    !simulate model and check match to some target moments:
    write(*,*) "Simulation starting..."
    call system_clock(timestart, timeres)
    call simulate

    call transition
    call system_clock(timeend, timeres)
    write(*,*) "Time elapsed:", (timeend - timestart) / timeres, "seconds"

    diffmoments=0
    !parameters are updated, then repeat whole big loop with solution/simulation to target moments.
    close(1)
end do


write(*,*) "C sum", sum(cchoice)


write(*,*) "model version"
write(*,*) "F  |  ", "theta  |  ", "r_rental  |  ", "costlyequity  |  ", "rentelasticity  |  ", "hpmin"
write(*,*) F, theta, r_rental, costlyequity, rentelasticity, hpmin
write(*,*) "A grid  |  ", "Dur grid  |  ", "Z grid  |  ", "HPrice grid  |  ", "households"
write(*,*) agridsize, Dgridsize, zgridsize, hpgridsize, numhouseholds

end program durables




subroutine simulate
USE nrtype; USE nrutil
USE nr
USE share
USE OMP_LIB
IMPLICIT NONE

REAL(8), dimension(numhouseholds,5) :: currenthouseholdstate, newhouseholdstate
REAL(8), dimension(:,:), ALLOCATABLE :: aggregatecontrib,mpc,consumption,consumptionmpc,consumptionhpshock,durableconsumption, currenttemp,currentperm,incomeholder,rentalind,choiceind,durableinvestment,actualwealth,financialwealth,housingnet,housinggross,actualwealthdividedbyincome,totalnetworthdividedbyincome,taxrate,welfare,diffc, ratioc, qchg
REAL(8), dimension(:,:), ALLOCATABLE :: liquidassetsminusmortgage,changec,changetemp,changeperm,changey,changeyconditional,changed,changedconditional,changetempconditional,changepermconditional,demeanedindicator
REAL(8), dimension(:,:), ALLOCATABLE :: alive
REAL(8), dimension(:,:), ALLOCATABLE :: housingstock, housingflow, rentalstock, numrent, numowntotal
REAL(8), dimension(numhouseholds,2) :: nearestnode
REAL(8) :: shock
integer, dimension(12) :: seedvalue
integer :: j,t,i,k
REAL(8), dimension(numhouseholds) :: insurancetemp, insuranceperm, cov1, cov2, cov3, cov4, insurancedtemp,insurancedperm,numtrans, insurancedtempconditional,insurancedpermconditional
REAL(8) :: adjust
REAL(8) :: conditionalindicator
EXTERNAL brentmindist
REAL(8) :: brentmindist
EXTERNAL calcpretax
REAL(8) :: calcpretax
EXTERNAL calcpretaxss
REAL(8) :: calcpretaxss
REAL(8) :: ax,bx,cx,tol
REAL(8) :: cov1est, cov2est, cov3est, cov4est,minest, constrained
REAL(8) :: m_eps2, m_eta2, m_cFy, m_yFy, m_epsc, m_etac, m_cIVy, m_yIVy, m_eps, m_eta, m_c, m_y, m_IVy, m_Fy, m_d, m_dFy, m_dIVy, m_epsd, m_etad
REAL(8) :: pretaxcalcholder,pretaxincome
REAL(8) :: tot1,tot2,tot3,adjustt,changedconditionalholder,changeyconditionalholder
REAL(8) :: rep
REAL(8) :: actualwealthtotal
REAL(8) :: exanteEVborn,exanteEVoverall,numobswelfare
REAL(8), dimension(hpgridsize,1) :: exanteEVbornstate, numobsstate, overallwelfarebystate,overallobsbystate, consumptionbystate, overallwelfarebystateyoung,overallwelfarebystatemiddle,overallwelfarebystateold,overallobsbystateyoung,overallobsbystatemiddle,overallobsbystateold, consumptionbystateyoung,consumptionbystatemiddle,consumptionbystateold
REAL(8), dimension(hpgridsize,350) :: overallwelfarebystateCE, overallwelfarebystateyoungCE,overallwelfarebystatemiddleCE,overallwelfarebystateoldCE
REAL(8), dimension(hpgridsize) :: cutoff
REAL(8), dimension(350,1) :: CE
REAL(8) :: medianincome, ratiocrent, ratiocown,diffcrent,diffcown,numrent_t,numown_t,numbuy_t,numsell_t, fracown
REAL(8), dimension(8):: numbins ! Initial income quantiles X (Owner, renter)
REAL(8) :: leverage_t, a_t, leverageown_t,aown_t,leveragerent_t,arent_t
real(8), dimension(numhouseholds, 5) :: fthb_flag ! 1: age when purchase fthb, 2: size of fth, 3: wealth while purchasing fth, 4: income while purchasing fth, 5: last period when house sold

REAL(8), dimension(:,:), ALLOCATABLE :: stateindicator1,stateindicator2,stateindicator5


do i=1,350
CE(i,1)=.97+(1.0*i)/5000.0
!write(*,*) CE(i,1)
end do

tol=1.0e-9
OPEN (UNIT=2, FILE="singlehouseholdsim.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
2 format (I6.2, I6.2, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=22, FILE="fthb.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
22 format (I6.2, I6.2, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=1989, FILE="dist_fthb.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
1989 format (I6.2, F6.2, F16.6, F16.6, F16.6)

! Orphaned as of now?
!OPEN (UNIT=3, FILE="agecoefficients.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
!3 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)
!
!4 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=5, FILE="lifecycleprofiles.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
5 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=9, FILE="householdresults.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
9 format (I16.6, I6.2, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6)

OPEN (UNIT=99, FILE="housingstock.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
99 format (I6.2, F16.6, F16.6, F16.6)

OPEN (UNIT=555, FILE="leverageprofiles.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
555 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=10, FILE="elasticity.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
10  format (F16.6, F16.6, F16.6,F16.6, F16.6, F16.6,F16.6, F16.6, F16.6)


OPEN (UNIT=11, FILE="statesandmpc.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
11  format (F16.6, F16.6, F16.6,F16.6, F16.6, F16.6,F16.6, F16.6, F16.6,F16.6, F16.6, F16.6)

seedvalue(:)=1
CALL random_seed(put=seedvalue(:))

! STATE DEFINITION
! state 1: liquid assets a
! state 2: durable assets d
! state 3: idiosyncratic income shock: z
! state 4: house price
! state 5: age: t

do i=1,numhouseholds
    if (real(i)/numhouseholds<Probinitcum(1)) then
        currenthouseholdstate(i,3)=1
    else
        do j=1,zgridsize-1
            if (real(i)/numhouseholds > Probinitcum(j) .AND. real(i)/numhouseholds<=Probinitcum(j+1)) then
                currenthouseholdstate(i,3)=j+1
            end if
        end do
    end if
end do

fthb_flag = 0 !flag for when purchasing a house (Tdie+1 = not buying a house in lifetime)

medianincome=income((zgridsize+1)/2,1)
numbins=0

! These initial values are from the create_initial_conditions Matlab files
! Each if block is a quantile of the income distribution, so quartiles for now.
! The shock variable determines if agent is house owner. Else agent rents.
! Deprecated while playing with common initial asset levels
do i=1,numhouseholds
    call random_number(shock)
    if (income(currenthouseholdstate(i,3),1)/medianincome<.4799) then
        if (shock < .0601) then
            fthb_flag(i,1) = 0
            currenthouseholdstate(i,2)=1.591
            currenthouseholdstate(i,1)=1.3336
            numbins(1)=numbins(1)+1
        else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.0017
            numbins(2)=numbins(2)+1
        end if
    elseif (income(currenthouseholdstate(i,3),1)/medianincome>=.4799 .and. income(currenthouseholdstate(i,3),1)/medianincome< 1) then
        if (shock<.1103) then
            fthb_flag(i,1) = 0
            currenthouseholdstate(i,2)=1.0458
            currenthouseholdstate(i,1)=0.3959
            numbins(3)=numbins(3)+1
        else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.00105
            numbins(4)=numbins(4)+1
        end if
    elseif (income(currenthouseholdstate(i,3),1)/medianincome>= 1 .and. income(currenthouseholdstate(i,3),1)/medianincome<1.8929) then
        if (shock<.0464) then
            fthb_flag(i,1) = 0
            currenthouseholdstate(i,2)=0.3978
            currenthouseholdstate(i,1)=0.0042
            numbins(5)=numbins(5)+1
        else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.0118
            numbins(6)=numbins(6)+1
        end if
    else
        if (shock<.2881) then
            fthb_flag(i,1) = 0
            currenthouseholdstate(i,2)=1.6706
            currenthouseholdstate(i,1)=0.0502
            numbins(7)=numbins(7)+1
        else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.0134
            numbins(8)=numbins(8)+1
        end if
    end if
end do


write(*,*) "Distribution of owners/renters by ", size(numbins)/2, "quantiles:"
write(*,*) numbins/numhouseholds
write(*,*) sum(numbins)

currenthouseholdstate(:,4)=2  !high house price state
newhouseholdstate(:,4)=2

! If you wanted to play with changing the initial distribution of wealth

currenthouseholdstate(:,1)=0  ! start with zero assets (could do something here with bequests)
currenthouseholdstate(:,2)=0


call random_number(shock)


! Allocation of automatic arrays for OpenMP (skip this)
ALLOCATE(aggregatecontrib(numhouseholds,Tdie),mpc(numhouseholds,Tdie),consumption(numhouseholds,Tdie), &
consumptionmpc(numhouseholds,Tdie),consumptionhpshock(numhouseholds,Tdie),durableconsumption(numhouseholds,Tdie), &
currenttemp(numhouseholds,Tdie),currentperm(numhouseholds,Tdie),incomeholder(numhouseholds,Tdie),rentalind(numhouseholds,Tdie), &
choiceind(numhouseholds,Tdie),durableinvestment(numhouseholds,Tdie),actualwealth(numhouseholds,Tdie),financialwealth(numhouseholds,Tdie), &
housingnet(numhouseholds,Tdie),housinggross(numhouseholds,Tdie),actualwealthdividedbyincome(numhouseholds,Tdie), &
totalnetworthdividedbyincome(numhouseholds,Tdie),taxrate(numhouseholds,Tdie),welfare(numhouseholds,Tdie), &
diffc(numhouseholds,Tdie), ratioc(numhouseholds,Tdie), qchg(numhouseholds,Tdie))

ALLOCATE(liquidassetsminusmortgage(numhouseholds,Tdie),changec(numhouseholds,Tdie),changetemp(numhouseholds,Tdie), &
changeperm(numhouseholds,Tdie),changey(numhouseholds,Tdie),changeyconditional(numhouseholds,Tdie), &
changed(numhouseholds,Tdie),changedconditional(numhouseholds,Tdie),changetempconditional(numhouseholds,Tdie), &
changepermconditional(numhouseholds,Tdie),demeanedindicator(numhouseholds,Tdie))

ALLOCATE(alive(numhouseholds,Tdie),housingstock(numhouseholds,Tdie),housingflow(numhouseholds,Tdie), &
rentalstock(numhouseholds,Tdie), numowntotal(numhouseholds,Tdie), numrent(numhouseholds,Tdie))
ALLOCATE(stateindicator1(numhouseholds,Tdie),stateindicator2(numhouseholds,Tdie), &
         stateindicator5(numhouseholds,Tdie))

numowntotal=0
numrent=0
actualwealthtotal=0
numtrans=0

housingstock=0
housingflow=0
rentalstock=0
alive=1


write(*,*) "Simulation start"
do t=1,Tdie
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


!write(*,*) "t", maxval((1+r)*(currenthouseholdstate(:,1)+theta*Currenthouseholdstate(:,2))-delta*currenthouseholdstate(:,2))
!write(*,*) "maxq", maxval(currenthouseholdstate(:,1)), maxval(currenthouseholdstate(:,2))
    stateindicator1=0
    stateindicator2=0
    stateindicator5=0


    ! 5th state is age (or equiv) time
    currenthouseholdstate(:,5)=t
!$OMP PARALLEL DO
    do i=1,numhouseholds
        if (alive(i,t)==1) then

        ! first thing is that if you are old, do you die or not
        if (t>Tretire) then
            if (alive(i,1)==1) then
                call random_number(shock)
                if (shock<deathprob(t-Tretire)) then
                    alive(i,t:Tdie)=0
                end if
            end if
        end if




        actualwealth(i,t)=currenthouseholdstate(i,1)+theta*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))  !total wealth is voluntary equity plus equity in durable
        financialwealth(i,t)=currenthouseholdstate(i,1)-(1-theta)*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))
        housingnet(i,t)=theta*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))
        housinggross(i,t)=currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))
        !qchg(i,t) = financialwealth(i,t) + (1-theta)*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)-1))
        qchg(i,t) = currenthouseholdstate(i,1) + (1-theta)*currenthouseholdstate(i,2)*(exp(hpnodes(currenthouseholdstate(i,4)))-exp(hpnodes(currenthouseholdstate(i,4)+1)))

        !write(*,*) "here"



        !call pol_linworking(qchg(i,t),currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumption(i,t),rentalind(i,t),welfare(i,t))
        !call pol_linworking(currenthouseholdstate(i,1)-.01,currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumptionmpc(i,t),rentalind(i,t),welfare(i,t))
        call pol_linworking(currenthouseholdstate(i,1),currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumption(i,t),rentalind(i,t),choiceind(i,t),welfare(i,t))
        mpc(i,t)=(consumptionmpc(i,t)-consumption(i,t))/-.01
        aggregatecontrib(i,t)=mpc(i,t)*currenthouseholdstate(i,2)*(1-delta)*(1-F)


        ! Track housing decisions going into period
        householdresultsholder(i,t,10)=fthb_flag(i, 1)
        householdresultsholder(i,t,11)=fthb_flag(i, 5)


        ! Decision tree begins with unambiguous renting
        if ((rentalind(i,t) > .995) .or.  (rentalind(i,t) /= rentalind(i,t))) then
            rentalstock(i,t)=newhouseholdstate(i,2)
            newhouseholdstate(i,2)=0
            fthb_flag(i,5) = Tdie+1
        else
            housingstock(i,t)=newhouseholdstate(i,2)
        ! Track agents as CARS eligible if their car is old enough
                if ((t >= fthb_flag(i, 5) + EligYrs) .and. (t <= Tretire)) then
                    fthb_flag(i, 1) = Tdie+1
                    fthb_flag(i, 5) = Tdie+1
                end if
        ! Start the CARS policy countdown whenever a durables purchase is made
                if (choiceind(i,t) == 1.0) then
                    numtrans(i) = 1
                    housingflow(i,t)=newhouseholdstate(i,2)
                    fthb_flag(i,5) = t
        ! Record stats for agent if they are eligible for policy
                    if (fthb_flag(i, 1) == Tdie+1) then
                        fthb_flag(i, 1) = t
                        fthb_flag(i, 2) = newhouseholdstate(i,2)
                        fthb_flag(i, 3) = currenthouseholdstate(i, 1)
                        fthb_flag(i, 4) = income(currenthouseholdstate(i, 3), t)
                        write(1989, 1989) i, fthb_flag(i, 1), fthb_flag(i, 2), & 
                        fthb_flag(i, 3), fthb_flag(i, 4)
                    end if
                end if
        end if

        if (newhouseholdstate(i,1)>amax) then
            newhouseholdstate(i,1)=.9999999*amax
        end if
        diffc(i,t)=consumptionhpshock(i,t)-consumption(i,t)
        ratioc(i,t) = consumptionhpshock(i,t)/consumption(i,t)

        if (t<=Tretire) then
        incomeholder(i,t)=income(currenthouseholdstate(i,3),t)
        else
        incomeholder(i,t)=retirementincome(currenthouseholdstate(i,3))
        end if

        !write(*,*) "here 2"
        householdresultsholder(i,t,1)=currenthouseholdstate(i,1)
        householdresultsholder(i,t,2)=currenthouseholdstate(i,2)
        householdresultsholder(i,t,3)=(currenthouseholdstate(i,3))
        householdresultsholder(i,t,4)=alive(i,t)
        householdresultsholder(i,t,5)=consumption(i,t)
        householdresultsholder(i,t,6)=newhouseholdstate(i,1)
        householdresultsholder(i,t,7)=newhouseholdstate(i,2)
        householdresultsholder(i,t,8)=rentalind(i,t)
        householdresultsholder(i,t,9)=incomeholder(i,t)

        if (i<=2500) then
            write(2,2) i, t, alive(i,t), currenthouseholdstate(i,1), currenthouseholdstate(i,2),newhouseholdstate(i,2),consumption(i,t),rentalind(i,t), currenthouseholdstate(i,3), currenthouseholdstate(i,4), currenthouseholdstate(i,1)-(1-theta)*currenthouseholdstate(i,2)
        end if

         if (i<=numhouseholds) then
            write(22,22) i, t, newhouseholdstate(i,2), currenthouseholdstate(i,2),rentalind(i,t), currenthouseholdstate(i, 1)
        end if


        durableinvestment(i,t)=newhouseholdstate(i,2)-(1-delta)*currenthouseholdstate(i,2)

        if (t<=Tretire) then
            call random_number(shock)

            if (shock < Probzcum(currenthouseholdstate(i,3),1)) then
                    newhouseholdstate(i,3)=1
            else
                do j=1,zgridsize-1
                    if (shock > Probzcum(currenthouseholdstate(i,3),j) .AND. shock<=Probzcum(currenthouseholdstate(i,3),j+1)) then
                        newhouseholdstate(i,3)=j+1
                    end if
                end do
            end if
        else
            newhouseholdstate(i,3)=currenthouseholdstate(i,3)
        end if

        if (t==3) then
            write(11,11) mpc(i,t), currenthouseholdstate(i,1), currenthouseholdstate(i,2), currenthouseholdstate(i,3), currenthouseholdstate(i,4), aggregatecontrib(i,t), rentalind(i,t)
        end if

        newhouseholdstate(i,5)= t+1
        end if

    end do
!$OMP END PARALLEL DO

    write(5,5) t*1.0, sum(alive(:,t)), sum(alive(:,t)*rentalind(:,t)), sum(alive(:,t)*(1-rentalind(:,t))), &
    sum(consumption(:,t))/sum(alive(:,t)), sum(consumption(:,t)*(1-rentalind(:,t)))/sum(alive(:,t)*(1-rentalind(:,t))), &
    sum(currenthouseholdstate(:,2)*(1-rentalind(:,t)))/sum(alive(:,t)*(1-rentalind(:,t))), &
    sum(currenthouseholdstate(:,1))/sum(alive(:,t)), sum(currenthouseholdstate(:,1)*(1-rentalind(:,t)))/sum(alive(:,t)*(1-rentalind(:,t))), &
    sum(currenthouseholdstate(:,1)-(1-theta)*currenthouseholdstate(:,2))/sum(alive(:,t)), &
    sum((currenthouseholdstate(:,1)-(1-theta)*currenthouseholdstate(:,2))*(1-rentalind(:,t)))/sum(alive(:,t)*(1-rentalind(:,t))), &
    sum((currenthouseholdstate(:,2)/(currenthouseholdstate(:,2)+currenthouseholdstate(:,1)-(1-theta)*currenthouseholdstate(:,2)))*(1-rentalind(:,t)))/sum(alive(:,t)*(1-rentalind(:,t))), &
    1-(sum(rentalind(:,t))/sum(alive(:,t)))


    !write(5,5) sum(financialwealth(:,t)*alive(:,t))/numhouseholds,sum(currenthouseholdstate(:,2)*alive(:,t))/numhouseholds, sum(consumption(:,t)*alive(:,t))/numhouseholds, sum(actualwealth(:,t)*alive(:,t))/numhouseholds, sum(newhouseholdstate(:,2)*alive(:,t))/numhouseholds, sum(newhouseholdstate(:,2)*rentalind(:,t)*alive(:,t))/numhouseholds, sum(rentalind(:,t)*alive(:,t))/numhouseholds, sum(alive(:,t))/numhouseholds, sum(consumption(:,t))/numhouseholds
     currenthouseholdstate=newhouseholdstate


    ! this is the housing price elasticity we want
    write(*,*) t*1.0, (sum(diffc(:,t))/sum(alive(:,t)))/(sum(consumption(:,t))/sum(alive(:,t)))/(hpmin)
    !write(*,*) "elasticity", sum(diffc(:,t)*alive(:,t))/sum(consumption(:,t))/sum(alive(:,t))/0.2


   ! write(*,*) t*1.0,sum(currenthouseholdstate(:,2)), sum(alive(:,t))

    !write(10,10) (sum((ratioc(:,t)-1))/sum(alive(:,t)))/(exp(hpmin) - 1), (sum(diffc(:,t))/sum(alive(:,t)))/(sum(consumption(:,t))/sum(alive(:,t)))/(hpmin)
    !write(10,10) (sum((ratioc(:,t)-1))/numhouseholds)/(exp(hpmin) - 1), (sum(diffc(:,t))/numhouseholds)/(sum(consumption(:,t))/numhouseholds)/(exp(hpmin) - 1), ratiocown/numown_t,ratiocrent/numrent_t,numown_t/(numown_t+numrent_t),numsell_t/numhouseholds,numbuy_t/numhouseholds,(sum(diffc(:,t)/hpmin)/numhouseholds),sum(aggregatecontrib(:,t))/numhouseholds
    write(99,99) t, sum(housingstock(:,t)), sum(housingflow(:,t)), sum(numtrans(:))
    numtrans = 0

    !write(555,555) a_t/numhouseholds, aown_t/numown_t, arent_t/numrent_t, leverage_t/numhouseholds, leverageown_t/numown_t, leveragerent_t/numrent_t
!
end do

!$OMP PARALLEL DO
do i=1,numhouseholds
    do t=1,60
      write(9,9) i,t, householdresultsholder(i,t,1),  householdresultsholder(i,t,2),  householdresultsholder(i,t,3),  householdresultsholder(i,t,4),  householdresultsholder(i,t,5),  householdresultsholder(i,t,6),  householdresultsholder(i,t,7),  householdresultsholder(i,t,8),   householdresultsholder(i,t,9)
    end do
end do
!$OMP END PARALLEL DO
write(*,*) "end"

close(2)
close(22)
close(5)
close(9)
close(99)
close(1989)
close(10)

end subroutine simulate



subroutine transition
USE nrtype; USE nrutil
USE nr
USE share
USE OMP_LIB
IMPLICIT NONE

REAL(8), dimension(numhouseholds,6) :: currenthouseholdstate, newhouseholdstate
! the 6th is the flag for fthb
REAL(8), dimension(:,:), ALLOCATABLE :: aggregatecontrib,mpc,consumption,consumptionmpc,consumptionhpshock,durableconsumption, currenttemp,currentperm,incomeholder,rentalind,choiceind,durableinvestment,actualwealth,financialwealth,housingnet,housinggross,actualwealthdividedbyincome,totalnetworthdividedbyincome,taxrate,welfare,diffc, ratioc, qchg
REAL(8), dimension(:,:), ALLOCATABLE :: liquidassetsminusmortgage,changec,changetemp,changeperm,changey,changeyconditional,changed,changedconditional,changetempconditional,changepermconditional,demeanedindicator
REAL(8), dimension(:,:), ALLOCATABLE :: alive
REAL(8), dimension(:,:,:), ALLOCATABLE :: housingstock, housingflow, rentalstock
REAL(8), dimension(numhouseholds,2) :: nearestnode
REAL(8) :: shock
REAL(8) :: hpdelta
integer, dimension(12) :: seedvalue
integer :: t,j,tj,i,k
REAL(8), dimension(numhouseholds) :: insurancetemp, insuranceperm, cov1, cov2, cov3, cov4, &
insurancedtemp,insurancedperm,numtrans,insurancedtempconditional,insurancedpermconditional
REAL(8), dimension(numhouseholds) :: thetastate
REAL(8) :: adjust
REAL(8) :: conditionalindicator
EXTERNAL brentmindist
REAL(8) :: brentmindist
EXTERNAL calcpretax
REAL(8) :: calcpretax
EXTERNAL calcpretaxss
REAL(8) :: calcpretaxss
REAL(8) :: ax,bx,cx,tol
REAL(8) :: cov1est, cov2est, cov3est, cov4est,minest, constrained
REAL(8) :: m_eps2, m_eta2, m_cFy, m_yFy, m_epsc, m_etac, m_cIVy, m_yIVy, m_eps, &
m_eta, m_c, m_y, m_IVy, m_Fy, m_d, m_dFy, m_dIVy, m_epsd, m_etad
REAL(8) :: pretaxcalcholder,pretaxincome
REAL(8) :: tot1,tot2,tot3,adjustt,changedconditionalholder,changeyconditionalholder
REAL(8) :: rep
REAL(8) :: actualwealthtotal
REAL(8) :: exanteEVborn,exanteEVoverall,numobswelfare
REAL(8), dimension(hpgridsize,1) :: exanteEVbornstate, numobsstate, &
overallwelfarebystate,overallobsbystate, consumptionbystate, &
overallwelfarebystateyoung,overallwelfarebystatemiddle,overallwelfarebystateold,overallobsbystateyoung,overallobsbystatemiddle,overallobsbystateold, &
consumptionbystateyoung,consumptionbystatemiddle,consumptionbystateold
REAL(8), dimension(hpgridsize,350) :: overallwelfarebystateCE, &
overallwelfarebystateyoungCE,overallwelfarebystatemiddleCE,overallwelfarebystateoldCE
REAL(8), dimension(hpgridsize) :: cutoff
REAL(8), dimension(350,1) :: CE
REAL(8) :: medianincome, ratiocrent,ratiocown,diffcrent,diffcown,numrent_t,numbuy_t,numsell_t
REAL(8), dimension(10):: numbins
REAL(8) :: leverage_t, a_t, leverageown_t,aown_t,leveragerent_t,arent_t
real(8), dimension(numhouseholds, 2) :: fthb_flag




tol=1.0e-9


seedvalue(:)=1
CALL random_seed(put=seedvalue(1:12))

OPEN (UNIT=999, FILE="housing_transit.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
999 format (I6.2, I6.2, F16.6, F16.6, F16.6)

OPEN (UNIT=299, FILE="transition_fthb.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
299 format (I6.2, I6.2, I6.2, F16.6, F16.6, F16.6)

OPEN (UNIT=499, FILE="transition_debug.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
499 format (I6.2, I6.2, I6.2, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

write(*, *) "Solving transition dynamics"


! Allocation of automatic arrays for OpenMP (skip this)
ALLOCATE(aggregatecontrib(numhouseholds,Tdie),mpc(numhouseholds,Tdie),consumption(numhouseholds,Tdie), &
consumptionmpc(numhouseholds,Tdie),consumptionhpshock(numhouseholds,Tdie),durableconsumption(numhouseholds,Tdie), &
currenttemp(numhouseholds,Tdie),currentperm(numhouseholds,Tdie),incomeholder(numhouseholds,Tdie),rentalind(numhouseholds,Tdie), &
choiceind(numhouseholds,Tdie),durableinvestment(numhouseholds,Tdie),actualwealth(numhouseholds,Tdie),financialwealth(numhouseholds,Tdie), &
housingnet(numhouseholds,Tdie),housinggross(numhouseholds,Tdie),actualwealthdividedbyincome(numhouseholds,Tdie), &
totalnetworthdividedbyincome(numhouseholds,Tdie),taxrate(numhouseholds,Tdie),welfare(numhouseholds,Tdie), &
diffc(numhouseholds,Tdie), ratioc(numhouseholds,Tdie), qchg(numhouseholds,Tdie))

ALLOCATE(liquidassetsminusmortgage(numhouseholds,Tdie),changec(numhouseholds,Tdie),changetemp(numhouseholds,Tdie), &
changeperm(numhouseholds,Tdie),changey(numhouseholds,Tdie),changeyconditional(numhouseholds,Tdie), &
changed(numhouseholds,Tdie),changedconditional(numhouseholds,Tdie),changetempconditional(numhouseholds,Tdie), &
changepermconditional(numhouseholds,Tdie),demeanedindicator(numhouseholds,Tdie))

ALLOCATE(alive(numhouseholds,Tdie),housingstock(numhouseholds,Tdie,Tdie), housingflow(numhouseholds,Tdie,Tdie), &
rentalstock(numhouseholds,Tdie,Tdie))



do t =1,Tdie

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
    if (t < householdresultsholder(i, t, 11) + EligYrs) then
        currenthouseholdstate(i, 6) = 0 ! Not CARS status
        newhouseholdstate(i, 6) = 0
        currenthouseholdstate(i, 4) = 2
    else
        currenthouseholdstate(i, 6) = 1 ! CARS status (bought durable more than 3 years ago)
        newhouseholdstate(i, 6) = 1
        currenthouseholdstate(i, 4) = 1
    end if
end do

fthb_flag(:, 1) = householdresultsholder(:, t, 10)
fthb_flag(:, 2) = householdresultsholder(:, t, 11)


call random_number(shock)
numtrans=0
alive=1

do tj=1,t-1
     write(999,999) Tretire - tj, t, 0.0, 0.0, 0.0
end do

do tj=t,Tdie

    ! write(*, *) tj
    ! 5th state is age (or equiv) time
    currenthouseholdstate(:,5)=tj
!$OMP PARALLEL DO
    do i=1,numhouseholds

        if (tj + 1 <= t + PolYrs) then
            thetastate(i) = thetanodes(currenthouseholdstate(i, 4))
        end if

        if (alive(i,tj) == 1) then
            ! first thing is that if you are old, do you die or not
            if (tj>Tretire) then
                if (alive(i,1)==1) then
                    call random_number(shock)
                    if (shock<deathprob(tj-Tretire)) then
                        alive(i,tj:Tdie)=0
                    end if
                end if
            end if



        actualwealth(i,tj)=currenthouseholdstate(i,1)+thetastate(i)*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))  !total wealth is voluntary equity plus equity in durable
        financialwealth(i,tj)=currenthouseholdstate(i,1)-(1-thetastate(i))*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))
        housingnet(i,tj)=thetastate(i)*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))
        housinggross(i,tj)=currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))


        if (fthb_flag(i, 1) < t + PolYrs) then
        !qchg(i,t) = financialwealth(i,t) + (1-thetastate(i))*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)-1))
        qchg(i,tj) = currenthouseholdstate(i,1) + (1-thetastate(i))*currenthouseholdstate(i,2)*hpdelta
        call pol_linworking(qchg(i,tj),currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumption(i,tj),rentalind(i,tj),choiceind(i,tj),welfare(i,tj))
        else
        !call pol_linworking(currenthouseholdstate(i,1)-.01,currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumptionmpc(i,tj),rentalind(i,tj),welfare(i,tj))
        call pol_linworking(currenthouseholdstate(i,1),currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumption(i,tj),rentalind(i,tj),choiceind(i,tj),welfare(i,tj))
        end if


        if (newhouseholdstate(i,1)>amax) then
            newhouseholdstate(i,1)=.9999999*amax
        end if

        if (tj<=Tretire) then
        incomeholder(i,tj)=income(currenthouseholdstate(i,3),tj)
        else
        incomeholder(i,tj)=retirementincome(currenthouseholdstate(i,3))
        end if


        ! Decision tree begins with unambiguous renting
        if ((rentalind(i,tj) > .995) .or.  (rentalind(i,tj) /= rentalind(i,tj))) then
        !   rentalstock(i,t)=newhouseholdstate(i,2)
            newhouseholdstate(i,2)=0
            fthb_flag(i,2) = Tdie+1
        else
            housingstock(i,tj,t)=newhouseholdstate(i,2)
            newhouseholdstate(i,6)=0
        ! Track agents as CARS eligible if their car is old enough
            if ((tj >= fthb_flag(i,2) + EligYrs) .and. (t <= Tretire)) then
                newhouseholdstate(i,6)=1
                fthb_flag(i, 1) = Tdie+1
                fthb_flag(i, 2) = Tdie+1
            end if
        ! Indicator of new durable purchase. Start the CARS policy countdown.
            if (choiceind(i,tj) == 1.0) then
                numtrans(i) = 1
                housingflow(i,tj,t)=newhouseholdstate(i,2)
                thetastate(i) = thetanodes(currenthouseholdstate(i, 4))
                fthb_flag(i,2) = tj
            ! Insert fthb characteristics again
                if ((fthb_flag(i, 1) == Tdie+1) .and. (tj <= Tretire)) then
                    fthb_flag(i, 1) = tj
                    write(299, 299) i, t, tj, newhouseholdstate(i, 2), &
                    currenthouseholdstate(i, 1), incomeholder(i, tj)
                end if
            end if
        end if


        if ((tj - t < 3) .or. ((tj - t > 8) .and. (tj - t <= 12))) then
            !write(499,499) i, t, tj, newhouseholdstate(i,2), rentalind(i,tj), newhouseholdstate(i, 1), consumption(i,tj), currenthouseholdstate(i,6)
            write(499,499) i, t, tj, newhouseholdstate(i,1), choiceind(i,tj), newhouseholdstate(i, 2), currenthouseholdstate(i,6), newhouseholdstate(i,6), fthb_flag(i,1), fthb_flag(i,2)
        end if

        if (tj<=Tretire) then
            call random_number(shock)

            if (shock < Probzcum(currenthouseholdstate(i,3),1)) then
                    newhouseholdstate(i,3)=1
            else
                do j=1,zgridsize-1
                    if (shock > Probzcum(currenthouseholdstate(i,3),j) .AND. &
                        shock<=Probzcum(currenthouseholdstate(i,3),j+1)) then
                        newhouseholdstate(i,3)=j+1
                    end if
                end do
            end if
        else
            newhouseholdstate(i,3)=currenthouseholdstate(i,3)
        end if

        if (tj + 1 >= t + PolYrs) then
            newhouseholdstate(i,4)=2 ! House prices from next period onwards are high
        else
            newhouseholdstate(i,4)=1
        end if
        newhouseholdstate(i,5)= tj+1
        currenthouseholdstate(i, :) = newhouseholdstate(i, :)
        end if

    end do
!$OMP END PARALLEL DO

    if (tj <= Tretire) then
        write(999,999) tj - t, t, sum(housingstock(:,tj,t)), sum(housingflow(:,tj,t)), sum(numtrans(:))
    end if

    numtrans = 0
    end do
    !write(5,5) sum(financialwealth(:,t)*alive(:,t))/numhouseholds,sum(currenthouseholdstate(:,2)*alive(:,t))/numhouseholds, sum(consumption(:,t)*alive(:,t))/numhouseholds, sum(actualwealth(:,t)*alive(:,t))/numhouseholds, sum(newhouseholdstate(:,2)*alive(:,t))/numhouseholds, sum(newhouseholdstate(:,2)*rentalind(:,t)*alive(:,t))/numhouseholds, sum(rentalind(:,t)*alive(:,t))/numhouseholds, sum(alive(:,t))/numhouseholds, sum(consumption(:,t))/numhouseholds
    !currenthouseholdstate=newhouseholdstate


    ! this is the housing price elasticity we want
    !write(*,*) t*1.0, (sum(diffc(:,t))/sum(alive(:,t)))/(sum(consumption(:,t))/sum(alive(:,t)))/(hpmin)
    !write(*,*) "elasticity", sum(diffc(:,t)*alive(:,t))/sum(consumption(:,t))/sum(alive(:,t))/0.2


   ! write(*,*) t*1.0,sum(currenthouseholdstate(:,2)), sum(alive(:,t))

    !write(10,10) (sum((ratioc(:,t)-1))/sum(alive(:,t)))/(exp(hpmin) - 1), (sum(diffc(:,t))/sum(alive(:,t)))/(sum(consumption(:,t))/sum(alive(:,t)))/(hpmin)
    !write(10,10) (sum((ratioc(:,t)-1))/numhouseholds)/(exp(hpmin) - 1), (sum(diffc(:,t))/numhouseholds)/(sum(consumption(:,t))/numhouseholds)/(exp(hpmin) - 1), ratiocown/numown_t,ratiocrent/numrent_t,numown_t/(numown_t+numrent_t),numsell_t/numhouseholds,numbuy_t/numhouseholds,(sum(diffc(:,t)/hpmin)/numhouseholds),sum(aggregatecontrib(:,t))/numhouseholds

    !write(555,555) a_t/numhouseholds, aown_t/numown_t, arent_t/numrent_t, leverage_t/numhouseholds, leverageown_t/numown_t, leveragerent_t/numrent_t

end do

!do i=1,numhouseholds
!    do t=1,60
!      write(9,9) i,t, householdresultsholder(i,t,1),  householdresultsholder(i,t,2),  householdresultsholder(i,t,3),  householdresultsholder(i,t,4),  householdresultsholder(i,t,5),  householdresultsholder(i,t,6),  householdresultsholder(i,t,7),  householdresultsholder(i,t,8)
!    end do
!end do
write(*,*) "end"

close(299)
close(499)
close(999)

end subroutine transition



subroutine pol_linworking(astate,Dstate,zstate,hpstate,t,achoicelin,Dchoicelin,cchoicelin,rentallin,choicelin,welfare)
    USE nrtype; USE nrutil
    USE nr
    USE share
    USE OMP_LIB
    IMPLICIT NONE
    REAL(8) :: weightal, weightDl
    REAL(8) :: nearestanode, nearestDnode
    REAL(8) :: al,ah, Dl,Dh
    REAL(8) :: astate,Dstate,zstate,hpstate,t,achoicelin,Dchoicelin,cchoicelin,rentallin,choicelin,welfare

    if (astate>amax) then
    astate=amax-.00000000001
    elseif (astate<amin) then
    astate=amin+.000000000001
    end if

    nearestanode=minloc(abs(astate-anodes),1)
    if (astate==anodes(nearestanode)) then
        if (nearestanode<agridsize) then
            al=nearestanode
            ah=nearestanode+1
            weightal=1.0
        else
            al=nearestanode-1
            ah=nearestanode
            weightal=0.0
        end if
    else
        if (astate-anodes(nearestanode)>0) then
            al=nearestanode
            ah=nearestanode+1
            if (ah>agridsize) then
            write(*,*) "error", ah, al, astate,nearestanode,anodes(nearestanode)
            end if
            weightal=1-(astate-anodes(al))/(anodes(ah)-anodes(al))
        else
            al=nearestanode-1
            if (al<1) then
            write(*,*) "error", ah, al, astate,nearestanode,anodes(nearestanode)
            end if
            ah=nearestanode
            weightal=1-(astate-anodes(al))/(anodes(ah)-anodes(al))
        end if
    end if


    nearestDnode=minloc(abs(Dstate-Dnodes),1)
    if (Dstate==Dnodes(nearestDnode)) then
        if (nearestDnode<Dgridsize) then
            Dl=nearestDnode
            Dh=nearestDnode+1
            weightDl=1.0
        else
            Dl=nearestDnode-1
            Dh=nearestDnode
            weightDl=0.0
        end if
    else
        if (Dstate-Dnodes(nearestDnode)>0) then
            Dl=nearestDnode
            Dh=nearestDnode+1
            weightDl=1-(Dstate-Dnodes(Dl))/(Dnodes(Dh)-Dnodes(Dl))
        else
            Dl=nearestDnode-1
            Dh=nearestDnode
            weightDl=1-(Dstate-Dnodes(Dl))/(Dnodes(Dh)-Dnodes(Dl))
        end if
    end if

    !write(*,*) weightal, weightDl,al,Dl,estate,Dchoice(al,Dl,estate)

    !write(*,*) "r  ", "r  ", al,ah,Dl,Dh,zstate,epsstate,t

   ! Explicitly coding rental indicator as one or zero instead of interpolating,
   ! experimental.
    if ((weightal>0.5) .AND. (weightDl>0.5)) then
        choicelin=choiceindicator(ah, Dh, zstate, hpstate, t)

        if (choiceindicator(ah, Dh, zstate, hpstate, t)==1) then
            !achoicelin=weightal*weightDl*achoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoiceadjust(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoiceadjust(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoiceadjust(ah,Dh,zstate,hpstate,t)
            rentallin=0

        elseif  (choiceindicator(ah, Dh, zstate, hpstate, t)==2) then
            !achoicelin=weightal*weightDl*achoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoicenoadjust(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoicenoadjust(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoicenoadjust(ah,Dh,zstate,hpstate,t)
            rentallin=0

        else
            !achoicelin=weightal*weightDl*achoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoicerent(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoicerent(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoicerent(ah,Dh,zstate,hpstate,t)
            rentallin=1
        endif

    elseif ((weightal<=0.5) .AND. (weightDl>0.5)) then
        choicelin=choiceindicator(al, Dh, zstate, hpstate, t)

        if (choiceindicator(al, Dh, zstate, hpstate, t)==1) then
            !achoicelin=weightal*weightDl*achoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoiceadjust(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoiceadjust(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoiceadjust(ah,Dh,zstate,hpstate,t)
            rentallin=0

        elseif  (choiceindicator(al, Dh, zstate, hpstate, t)==2) then
            !achoicelin=weightal*weightDl*achoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoicenoadjust(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoicenoadjust(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoicenoadjust(ah,Dh,zstate,hpstate,t)
            rentallin=0

        else
            !achoicelin=weightal*weightDl*achoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoicerent(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoicerent(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoicerent(ah,Dh,zstate,hpstate,t)
            rentallin=1
        endif

    elseif ((weightal>0.5) .and. (weightDl<=0.5)) then
        choicelin=choiceindicator(ah, Dl, zstate, hpstate, t)

        if (choiceindicator(ah, Dl, zstate, hpstate, t)==1) then
            !achoicelin=weightal*weightDl*achoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoiceadjust(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoiceadjust(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoiceadjust(ah,Dh,zstate,hpstate,t)
            rentallin=0

        elseif  (choiceindicator(ah, Dl, zstate, hpstate, t)==2) then
            !achoicelin=weightal*weightDl*achoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoicenoadjust(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoicenoadjust(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoicenoadjust(ah,Dh,zstate,hpstate,t)
            rentallin=0

        else
            !achoicelin=weightal*weightDl*achoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoicerent(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoicerent(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoicerent(ah,Dh,zstate,hpstate,t)
            rentallin=1
        endif

    else
        choicelin=choiceindicator(al, Dl, zstate, hpstate, t)

        if (choiceindicator(al, Dl, zstate, hpstate, t)==1) then
            !achoicelin=weightal*weightDl*achoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoiceadjust(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoiceadjust(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoiceadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoiceadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoiceadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoiceadjust(ah,Dh,zstate,hpstate,t)
            rentallin=0

        elseif  (choiceindicator(al, Dl, zstate, hpstate, t)==2) then
            !achoicelin=weightal*weightDl*achoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoicenoadjust(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoicenoadjust(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoicenoadjust(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoicenoadjust(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoicenoadjust(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoicenoadjust(ah,Dh,zstate,hpstate,t)
            rentallin=0

        else
            !achoicelin=weightal*weightDl*achoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoicerent(ah,Dh,zstate,hpstate,t)
            !Dchoicelin=weightal*weightDl*Dchoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoicerent(ah,Dh,zstate,hpstate,t)
            !cchoicelin=weightal*weightDl*cchoicerent(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoicerent(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoicerent(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoicerent(ah,Dh,zstate,hpstate,t)
            rentallin=1
        endif

    endif


    achoicelin=weightal*weightDl*achoice(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoice(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoice(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoice(ah,Dh,zstate,hpstate,t)
    Dchoicelin=weightal*weightDl*Dchoice(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoice(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoice(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoice(ah,Dh,zstate,hpstate,t)
    cchoicelin=weightal*weightDl*cchoice(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoice(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoice(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoice(ah,Dh,zstate,hpstate,t)
   ! rentallin=weightal*weightDl*rentalindicator(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*rentalindicator(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*rentalindicator(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*rentalindicator(ah,Dh,zstate,hpstate,t)
    welfare=weightal*weightDl*EV(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*EV(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*EV(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*EV(ah,Dh,zstate,hpstate,t)

end subroutine pol_linworking




subroutine plotpolicy_constantwealth
    USE nrtype; USE nrutil
    USE nr
    USE share
    USE OMP_LIB
    IMPLICIT NONE
    REAL(8) :: i,j,k,l,m,n,t, iter,p
    REAL(8) :: wealth, h, q, achoicelin, Dchoicelin, cchoicelin, rentallin, choicelin, welfare

    OPEN (UNIT=266, FILE="policy_constantwealth.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    266 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

    t=3.0
    p=2.0
    do k=1,zgridsize
    do i=1,50
        wealth=0+.1*i
        do j=1,50
            h=0+(j-1)*.1
            q=wealth-theta*h
            call pol_linworking(q,h,k,p,t,achoicelin,Dchoicelin,cchoicelin,rentallin,choicelin,welfare)
            write(266,266) k, wealth, h, Dchoicelin, cchoicelin, achoicelin
        end do
    end do
    end do



end subroutine plotpolicy_constantwealth

subroutine computempc(astate,Dstate,zstate,hpstate,t,mpcholder)
    USE nrtype; USE nrutil
    USE nr
    USE share
    USE OMP_LIB
    IMPLICIT NONE
    REAL(8) :: weightal, weightDl
    REAL(8) :: nearestanode, nearestDnode
    REAL(8) :: al,ah, Dl,Dh
    REAL(8), intent(in) :: astate,Dstate,zstate,hpstate,t
    REAL(8) :: achoicelin,Dchoicelin,cchoicelin1,rentallin,welfare, mpcholder,cchoicelin2, astate2,Dstate2

    astate2=astate
    if (astate2>amax) then
    astate2=amax
    elseif (astate2<amin) then
    astate2=amin
    end if

    nearestanode=minloc(abs(astate2-anodes),1)
    if (astate2==anodes(nearestanode)) then
        if (nearestanode<agridsize) then
            al=nearestanode
            ah=nearestanode+1
            weightal=1.0
        else
            al=nearestanode-1
            ah=nearestanode
            weightal=0.0
        end if
    else
        if (astate2-anodes(nearestanode)>0) then
            al=nearestanode
            ah=nearestanode+1
            if (ah>agridsize) then
            write(*,*) "error", ah, al, astate,nearestanode,anodes(nearestanode)
            end if
            weightal=1-(astate2-anodes(al))/(anodes(ah)-anodes(al))
        else
            al=nearestanode-1
            ah=nearestanode
            weightal=1-(astate2-anodes(al))/(anodes(ah)-anodes(al))
        end if
    end if


    nearestDnode=minloc(abs(Dstate-Dnodes),1)
    if (Dstate==Dnodes(nearestDnode)) then
        if (nearestDnode<Dgridsize) then
            Dl=nearestDnode
            Dh=nearestDnode+1
            weightDl=1.0
        else
            Dl=nearestDnode-1
            Dh=nearestDnode
            weightDl=0.0
        end if
    else
        if (Dstate-Dnodes(nearestDnode)>0) then
            Dl=nearestDnode
            Dh=nearestDnode+1
            weightDl=1-(Dstate-Dnodes(Dl))/(Dnodes(Dh)-Dnodes(Dl))
        else
            Dl=nearestDnode-1
            Dh=nearestDnode
            weightDl=1-(Dstate-Dnodes(Dl))/(Dnodes(Dh)-Dnodes(Dl))
        end if
    end if

    cchoicelin1=weightal*weightDl*cchoice(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoice(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoice(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoice(ah,Dh,zstate,hpstate,t)

    astate2=astate+.01
    !astate2=astate2+.01
    if (astate2>amax) then
    astate2=amax-.000000001
    elseif (astate2<amin) then
    astate2=amin
    end if

   nearestanode=minloc(abs(astate2-anodes),1)
    if (astate2==anodes(nearestanode)) then
        if (nearestanode<agridsize) then
            al=nearestanode
            ah=nearestanode+1
            weightal=1.0
        else
            al=nearestanode-1
            ah=nearestanode
            weightal=0.0
        end if
    else
        if (astate2-anodes(nearestanode)>0) then
            al=nearestanode
            ah=nearestanode+1
            if (ah>agridsize) then
            write(*,*) "error", ah, al, astate2,nearestanode,anodes(nearestanode)
            end if
            weightal=1-(astate2-anodes(al))/(anodes(ah)-anodes(al))
        else
            al=nearestanode-1
            ah=nearestanode
            weightal=1-(astate2-anodes(al))/(anodes(ah)-anodes(al))
        end if
    end if


    nearestDnode=minloc(abs(Dstate-Dnodes),1)
    if (Dstate==Dnodes(nearestDnode)) then
        if (nearestDnode<Dgridsize) then
            Dl=nearestDnode
            Dh=nearestDnode+1
            weightDl=1.0
        else
            Dl=nearestDnode-1
            Dh=nearestDnode
            weightDl=0.0
        end if
    else
        if (Dstate-Dnodes(nearestDnode)>0) then
            Dl=nearestDnode
            Dh=nearestDnode+1
            weightDl=1-(Dstate-Dnodes(Dl))/(Dnodes(Dh)-Dnodes(Dl))
        else
            Dl=nearestDnode-1
            Dh=nearestDnode
            weightDl=1-(Dstate-Dnodes(Dl))/(Dnodes(Dh)-Dnodes(Dl))
        end if
    end if


    cchoicelin2=weightal*weightDl*cchoice(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoice(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoice(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoice(ah,Dh,zstate,hpstate,t)
    mpcholder=(cchoicelin2-cchoicelin1)/.01

end subroutine computempc













subroutine solveworkingproblem
    USE nrtype; USE nrutil
    USE nr
    USE share
    USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: i,j,k,l,t
    INTEGER :: thread_count
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,hpgridsize,2) :: optpolicynoadjust,optpolicyadjust
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,hpgridsize,3,2) :: pstartadjust
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,hpgridsize,3) :: ystartadjust
    REAL(8), dimension(:,:,:,:), ALLOCATABLE :: ax,bx,cx,adjust
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,hpgridsize,5) :: state

    REAL(8), DIMENSION(2) :: p2

    EXTERNAL valfuncnoadjust
    REAL(8) :: valfuncnoadjust

    REAL(8) valfuncadjust
    EXTERNAL valfuncadjust

    EXTERNAL valfuncadjust2
    REAL(8) valfuncadjust2

    REAL(8) valfuncrent
    EXTERNAL valfuncrent

    EXTERNAL valfuncrent2
    REAL(8) valfuncrent2

    EXTERNAL valfuncrent3
    REAL(8) valfuncrent3

    EXTERNAL brentnew
    REAL(8) brentnew
    REAL(8) ftol



    ftol=5e-9


! Allocation of automatic arrays for OpenMP (skip this)
ALLOCATE(ax(agridsize,Dgridsize,zgridsize,hpgridsize),bx(agridsize,Dgridsize,zgridsize,hpgridsize),cx(agridsize,Dgridsize,zgridsize,hpgridsize),adjust(agridsize,Dgridsize,zgridsize,hpgridsize))


    adjust=0

    do t=Tdie,1,-1
        write(*,*) "Processing year state:", t
!$OMP PARALLEL
!$OMP DO
        do i=1,agridsize
            do j=1,Dgridsize
                do k=1,zgridsize
                    do l=1,hpgridsize
                        pstartadjust(i,j,k,l,1,1)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,1,2)=.8*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,1)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,2)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,1)=.5*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,2)=.2*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))


                        state(i,j,k,l,1)=anodes(i)
                        state(i,j,k,l,2)=Dnodes(j)
                        state(i,j,k,l,3)=k
                        state(i,j,k,l,4)=l
                        state(i,j,k,l,5)=t

                        ystartadjust(i,j,k,l,1)=valfuncadjust2(pstartadjust(i,j,k,l,1,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,2)=valfuncadjust2(pstartadjust(i,j,k,l,2,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,3)=valfuncadjust2(pstartadjust(i,j,k,l,3,:),state(i,j,k,l,:))


                        call amoeba(state(i,j,k,l,:),pstartadjust(i,j,k,l,:,:),ystartadjust(i,j,k,l,:),ftol,valfuncadjust)
                        Vadjust(i,j,k,l,t)=ystartadjust(i,j,k,l,1)
                        achoiceadjust(i,j,k,l,t)=pstartadjust(i,j,k,l,1,1)
                        Dchoiceadjust(i,j,k,l,t)=pstartadjust(i,j,k,l,1,2)


                        pstartadjust(i,j,k,l,1,1)=0.0
                        pstartadjust(i,j,k,l,1,2)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,1)=.05*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,2)=.1*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,1)=.4*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,2)=.21*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))

                        ystartadjust(i,j,k,l,1)=valfuncadjust2(pstartadjust(i,j,k,l,1,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,2)=valfuncadjust2(pstartadjust(i,j,k,l,2,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,3)=valfuncadjust2(pstartadjust(i,j,k,l,3,:),state(i,j,k,l,:))

                        call amoeba(state(i,j,k,l,:),pstartadjust(i,j,k,l,:,:),ystartadjust(i,j,k,l,:),ftol,valfuncadjust)
                        if (ystartadjust(i,j,k,l,1)<Vadjust(i,j,k,l,t)) then  !(again we're minimizing)
                            Vadjust(i,j,k,l,t)=ystartadjust(i,j,k,l,1)
                            achoiceadjust(i,j,k,l,t)=pstartadjust(i,j,k,l,1,1)
                            Dchoiceadjust(i,j,k,l,t)=pstartadjust(i,j,k,l,1,2)
                        end if










                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        pstartadjust(i,j,k,l,1,1)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,1,2)=.8*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,1)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,2)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,1)=.5*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,2)=.2*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))

                        state(i,j,k,l,1)=anodes(i)
                        state(i,j,k,l,2)=Dnodes(j)
                        state(i,j,k,l,3)=k
                        state(i,j,k,l,4)=l
                        state(i,j,k,l,5)=t

                        ystartadjust(i,j,k,l,1)=valfuncrent2(pstartadjust(i,j,k,l,1,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,2)=valfuncrent2(pstartadjust(i,j,k,l,2,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,3)=valfuncrent2(pstartadjust(i,j,k,l,3,:),state(i,j,k,l,:))

                        call amoeba(state(i,j,k,l,:),pstartadjust(i,j,k,l,:,:),ystartadjust(i,j,k,l,:),ftol,valfuncrent)


                        Vrent(i,j,k,l,t)=ystartadjust(i,j,k,l,1)
                        achoicerent(i,j,k,l,t)=pstartadjust(i,j,k,l,1,1)
                        Dchoicerent(i,j,k,l,t)=pstartadjust(i,j,k,l,1,2)


                        pstartadjust(i,j,k,l,1,1)=0.0
                        pstartadjust(i,j,k,l,1,2)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,1)=.05*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,2)=.1*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,1)=.4*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,2)=.21*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))


                        ystartadjust(i,j,k,l,1)=valfuncrent2(pstartadjust(i,j,k,l,1,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,2)=valfuncrent2(pstartadjust(i,j,k,l,2,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,3)=valfuncrent2(pstartadjust(i,j,k,l,3,:),state(i,j,k,l,:))

                        call amoeba(state(i,j,k,l,:),pstartadjust(i,j,k,l,:,:),ystartadjust(i,j,k,l,:),ftol,valfuncrent)


                        if (ystartadjust(i,j,k,l,1)<Vrent(i,j,k,l,t)) then  !(again we're minimizing)
                            Vrent(i,j,k,l,t)=ystartadjust(i,j,k,l,1)
                            achoicerent(i,j,k,l,t)=pstartadjust(i,j,k,l,1,1)
                            Dchoicerent(i,j,k,l,t)=pstartadjust(i,j,k,l,1,2)
                        end if


                        if (costlyequity==1) then

                            if (anodes(i)>=(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))) then
                                ax(i,j,k,l)=(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                            else
                                ax(i,j,k,l)=anodes(i)
                            end if
                        bx(i,j,k,l)=ax(i,j,k,l)
                        else
                        ax(i,j,k,l)=0.0
                        bx(i,j,k,l)=borrowconstraint
                        end if
                        cx(i,j,k,l)=amax

                        state(i,j,k,l,2)=j
                        Vnoadjust(i,j,k,l,t)=brentnew(ax(i,j,k,l),bx(i,j,k,l),cx(i,j,k,l),valfuncnoadjust,ftol,state(i,j,k,l,1),state(i,j,k,l,2),state(i,j,k,l,3),state(i,j,k,l,4),state(i,j,k,l,5),achoicenoadjust(i,j,k,l,t))
                        Dchoicenoadjust(i,j,k,l,t)=Dnodes(j)

                        ! Gideon: computing cchoiceadjust, cchoicenoadjust and cchoicerent
                        if (anodes(i)-(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))<0) then
                                cchoiceadjust(i,j,k,l,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoiceadjust(i,j,k,l,t)*exp(hpnodes(l))-achoiceadjust(i,j,k,l,t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                cchoicenoadjust(i,j,k,l,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoicenoadjust(i,j,k,l,t)*exp(hpnodes(l))-achoicenoadjust(i,j,k,l,t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                cchoicerent(i,j,k,l,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoicerent(i,j,k,l,t)*exp(rentelasticity*hpnodes(l))-achoicerent(i,j,k,l,t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                        else
                                cchoiceadjust(i,j,k,l,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoicenoadjust(i,j,k,l,t)*exp(hpnodes(l))-achoicenoadjust(i,j,k,l,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                cchoicenoadjust(i,j,k,l,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoicenoadjust(i,j,k,l,t)*exp(hpnodes(l))-achoicenoadjust(i,j,k,l,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                cchoicerent(i,j,k,l,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoicerent(i,j,k,l,t)*exp(rentelasticity*hpnodes(l))-achoicerent(i,j,k,l,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                        end if

                        if (Vadjust(i,j,k,l,t)<Vnoadjust(i,j,k,l,t) .and. Vadjust(i,j,k,l,t)<Vrent(i,j,k,l,t)) then  ! since V = - V from minimization
                            achoice(i,j,k,l,t)=achoiceadjust(i,j,k,l,t)
                            Dchoice(i,j,k,l,t)=Dchoiceadjust(i,j,k,l,t)
                            if (anodes(i)-(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))<0) then
                                cchoice(i,j,k,l,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                !write(*,*) "diff1", income(k,t)+(1+rborrow*exp(hpnodes(l))**)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**)*(1-thetanodes(l))*(1-F)*Dnodes(j)*exp(hpnodes(l))-(income(k,t)+(1+rborrow*exp(hpnodes(l))**(-0.0))*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**(-0.0))*(1-thetanodes(l))*(1-F)*Dnodes(j)*exp(hpnodes(l)))
                            else
                                cchoice(i,j,k,l,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                            end if
                            rentalindicator(i,j,k,l,t)=0
                            choiceindicator(i,j,k,l,t)=1
                        elseif (Vnoadjust(i,j,k,l,t)<Vadjust(i,j,k,l,t) .and. Vnoadjust(i,j,k,l,t)<Vrent(i,j,k,l,t)) then
                            achoice(i,j,k,l,t)=achoicenoadjust(i,j,k,l,t)
                            Dchoice(i,j,k,l,t)=Dchoicenoadjust(i,j,k,l,t)


                            if (anodes(i)-(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))<0) then
                                cchoice(i,j,k,l,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                !write(*,*) "diff2", income(k,t)+(1+rborrow*exp(hpnodes(l))**)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))-(income(k,t)+(1+rborrow*exp(hpnodes(l))**(-0.0))*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**(-0.0))*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l)))
                            else
                                cchoice(i,j,k,l,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                            end if
                            rentalindicator(i,j,k,l,t)=0
                            choiceindicator(i,j,k,l,t)=2
                        else
                            achoice(i,j,k,l,t)=achoicerent(i,j,k,l,t)
                            Dchoice(i,j,k,l,t)=Dchoicerent(i,j,k,l,t)
                            if (anodes(i)-(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))<0) then
                                cchoice(i,j,k,l,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoice(i,j,k,l,t)*exp(rentelasticity*hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                !write(*,*) "diff3", income(k,t)+(1+rborrow*exp(hpnodes(l))**)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**)*(1-thetanodes(l))*(1-F)*Dnodes(j)*exp(hpnodes(l)) - (income(k,t)+(1+rborrow*exp(hpnodes(l))**(-0.0))*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**(-0.0))*(1-thetanodes(l))*(1-F)*Dnodes(j)*exp(hpnodes(l)))
                            else
                                cchoice(i,j,k,l,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoice(i,j,k,l,t)*exp(rentelasticity*hpnodes(l))-achoice(i,j,k,l,t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                            end if
                            !write(*,*) l,(1+rborrow*exp(hpnodes(l))**)
                            rentalindicator(i,j,k,l,t)=1
                            choiceindicator(i,j,k,l,t)=3
                        end if
                       ! write(*,*) exp(hpnodes(l)), (1+rborrow*exp(hpnodes(l))**), income(6,l,1)
                    end do
                end do
            end do
        end do
!$OMP END DO
!$OMP END PARALLEL
        Vnoadjust(:,:,:,:,t)=-Vnoadjust(:,:,:,:,t)
        Vadjust(:,:,:,:,t)=-Vadjust(:,:,:,:,t)
        Vrent(:,:,:,:,t)=-Vrent(:,:,:,:,t)

        EV(:,:,:,:,t)=max(Vrent(:,:,:,:,t),max(Vnoadjust(:,:,:,:,t),Vadjust(:,:,:,:,t)))

       ! write(*,*) EV(1,1,2,1,2), EV(2,1,2,1,2)

        if (t==1) then
        write(*,*) "policy functions:"
        write(*,*) achoicerent(1,1,2,1,1),dchoicerent(1,1,2,1,1)
        write(*,*) achoiceadjust(1,1,2,1,1),dchoiceadjust(1,1,2,1,1)
        write(*,*) achoicenoadjust(1,1,2,1,1),dchoicenoadjust(1,1,2,1,1)
        write(*,*) achoice(1,1,2,1,1),dchoice(1,1,2,1,1),cchoice(1,1,2,1,1)

                         state(1,1,1,1,1)=anodes(1)
                        state(1,1,1,1,2)=Dnodes(1)
                        state(1,1,1,1,3)=2
                        state(1,1,1,1,4)=1
                        state(1,1,1,1,5)=1
                        pstartadjust(1,1,1,1,1,1)=achoicerent(1,1,2,1,1)
                        pstartadjust(1,1,1,1,1,2)=dchoicerent(1,1,2,1,1)

                        p2(1)=pstartadjust(1,1,2,1,1,1)
                        p2(2)=pstartadjust(1,1,2,1,1,2)
        ! write(*,*) valfuncrent3(p2,state(1,1,1,1,:))
        ! p2(1)=achoicerent(1,1,2,1,1)+.000001
        ! write(*,*) valfuncrent3(p2,state(1,1,1,1,:))
        end if

        if (priceshock==1) then
        cchoiceshock=cchoice
        Dchoiceshock=Dchoice
        end if
end do

DEALLOCATE(ax,bx,cx,adjust)

end subroutine solveworkingproblem










REAL(8) FUNCTION valfuncadjust(p,state)  !p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
USE SHARE
IMPLICIT NONE
REAL(8), DIMENSION(2) :: p
REAL(8), DIMENSION(5) :: state
REAL(8) :: a, D, z, hp, t, zindex,hpindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder,hpholder,thetaholder

!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)

hpholder=hpnodes(hpindex)
thetaholder=thetanodes(hpindex)

currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dnodes(Dgridsize) .AND. Dcurrent>=0) then
    nearestDnode=minloc(abs(Dprime-Dnodes),1)
    if (Dprime==Dnodes(nearestDnode)) then
        if (nearestDnode<Dgridsize) then
            Dprimel=nearestDnode
            Dprimeh=nearestDnode+1
            weightDprimel=1.0
        else
            Dprimel=nearestDnode-1
            Dprimeh=nearestDnode
            weightDprimel=0.0
        end if
    else
        if (Dprime-Dnodes(nearestDnode)>0) then
            Dprimel=nearestDnode
            Dprimeh=nearestDnode+1
            weightDprimel=1-(Dprime-Dnodes(Dprimel))/(Dnodes(Dprimeh)-Dnodes(Dprimel))
        else
            Dprimel=nearestDnode-1
            Dprimeh=nearestDnode
            weightDprimel=1-(Dprime-Dnodes(Dprimel))/(Dnodes(Dprimeh)-Dnodes(Dprimel))
        end if
    end if
    if ((a-D*(1-thetaholder)*exp(hpholder))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)*exp(hpholder)-thetaholder*Dcurrent*exp(hpholder)-aprimeholder-(1+rborrow)*(1-thetaholder)*D*exp(hpholder))
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)*exp(hpholder)-thetaholder*Dcurrent*exp(hpholder)-aprimeholder-(1+r)*(1-thetaholder)*D*exp(hpholder))
    end if

    !need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.  Find nearest grid points and linearly interpolate.
    if (consumption>0) then
        if (elasticity .ne. 1) then
            valfuncadjust=((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))**(1-elasticity)/(1-elasticity)
        else
            valfuncadjust=log((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))
        end if
        EVholder=0
        EVholder2=0
        if (t<=Tretire) then
            do j=1,hpgridsize
                nearestanode=minloc(abs(aprime-anodes),1)
                if (aprime==anodes(nearestanode)) then
                    if (nearestanode<agridsize) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1.0
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=0.0
                    end if
                else
                    if (aprime-anodes(nearestanode)>0) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    end if
                end if
                do i=minexpectationz(zindex),maxexpectationz(zindex)

                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,i,j,t+1)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,i,j,t+1)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,i,j,t+1)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,i,j,t+1)
                end do
            end do
            valfuncadjust=valfuncadjust+beta2*EVholder
        else
            do j=1,hpgridsize
                nearestanode=minloc(abs(aprime-anodes),1)
                if (aprime==anodes(nearestanode)) then
                    if (nearestanode<agridsize) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1.0
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=0.0
                    end if
                else
                    if (aprime-anodes(nearestanode)>0) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    end if
                end if

                    EVholder=EVholder+Probhp(hpindex,j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,j,t+1)
                    EVholder=EVholder+Probhp(hpindex,j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,j,t+1)
                    EVholder=EVholder+Probhp(hpindex,j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,j,t+1)
                    EVholder=EVholder+Probhp(hpindex,j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,j,t+1)

                    EVholder2=EVholder2+Probhp(hpindex,j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,j,Tdie+1)
                    EVholder2=EVholder2+Probhp(hpindex,j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,j,Tdie+1)
                    EVholder2=EVholder2+Probhp(hpindex,j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,j,Tdie+1)
                    EVholder2=EVholder2+Probhp(hpindex,j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,j,Tdie+1)
            end do
            valfuncadjust=valfuncadjust+beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if


    else
        valfuncadjust=-1.0e13
    end if
else
     valfuncadjust=-1.0e13
end if



valfuncadjust=-valfuncadjust  !powell minimizes



END FUNCTION valfuncadjust


REAL(8) FUNCTION valfuncrent(p,state)  !p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
USE SHARE
IMPLICIT NONE
REAL(8), DIMENSION(2) :: p
REAL(8), DIMENSION(5) :: state
REAL(8) :: a, D, z, hp, t, zindex,hpindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder,hpholder,thetaholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)

hpholder=hpnodes(hpindex)
thetaholder=thetanodes(hpindex)

currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dnodes(Dgridsize) .AND. Dcurrent>=0) then
    if ((a-D*(1-thetaholder)*exp(hpholder))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)*exp(hpholder)-r_rental*Dcurrent*exp(rentelasticity*hpholder)-aprimeholder-(1+rborrow)*(1-thetaholder)*D*exp(hpholder))
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)*exp(hpholder)-r_rental*Dcurrent*exp(rentelasticity*hpholder)-aprimeholder-(1+r)*(1-thetaholder)*D*exp(hpholder))
    end if
    !need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.  Find nearest grid points and linearly interpolate.
    if (consumption>0) then
        if (elasticity .ne. 1) then
            valfuncrent=((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))**(1-elasticity)/(1-elasticity)
        else
            valfuncrent=log((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))
        end if
        EVholder=0
        EVholder2=0
        if (t<=Tretire) then
            do j=1,hpgridsize
                nearestanode=minloc(abs(aprime-anodes),1)
                if (aprime==anodes(nearestanode)) then
                    if (nearestanode<agridsize) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1.0
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=0.0
                    end if
                else
                    if (aprime-anodes(nearestanode)>0) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    end if
                end if

                do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*weightaprimel*EV(aprimel,1,i,j,t+1)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,1,i,j,t+1)
                end do
            end do
            valfuncrent=valfuncrent+beta2*EVholder
        else
            do j=1,hpgridsize
                nearestanode=minloc(abs(aprime-anodes),1)
                if (aprime==anodes(nearestanode)) then
                    if (nearestanode<agridsize) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1.0
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=0.0
                    end if
                else
                    if (aprime-anodes(nearestanode)>0) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    end if
                end if

                EVholder=EVholder+Probhp(hpindex,j)*weightaprimel*EV(aprimel,1,zindex,j,t+1)
                EVholder=EVholder+Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,1,zindex,j,t+1)

                EVholder2=EVholder2+Probhp(hpindex,j)*weightaprimel*EV(aprimel,1,zindex,j,Tdie+1)
                EVholder2=EVholder2+Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,1,zindex,j,Tdie+1)
            end do
            valfuncrent=valfuncrent+beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if


    else
        valfuncrent=-1.0e13
    end if
else
     valfuncrent=-1.0e13
end if



valfuncrent=-valfuncrent  !powell minimizes



END FUNCTION valfuncrent



REAL(8) FUNCTION valfuncadjust2(p,state)  !p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
USE SHARE
IMPLICIT NONE
REAL(8), DIMENSION(2) :: p
REAL(8), DIMENSION(5) :: state
REAL(8) :: a, D, z, hp, t, zindex,hpindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder,hpholder,thetaholder

!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)

hpholder=hpnodes(hpindex)
thetaholder=thetanodes(hpindex)

currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dnodes(Dgridsize) .AND. Dcurrent>=0) then
    nearestDnode=minloc(abs(Dprime-Dnodes),1)
    if (Dprime==Dnodes(nearestDnode)) then
        if (nearestDnode<Dgridsize) then
            Dprimel=nearestDnode
            Dprimeh=nearestDnode+1
            weightDprimel=1.0
        else
            Dprimel=nearestDnode-1
            Dprimeh=nearestDnode
            weightDprimel=0.0
        end if
    else
        if (Dprime-Dnodes(nearestDnode)>0) then
            Dprimel=nearestDnode
            Dprimeh=nearestDnode+1
            weightDprimel=1-(Dprime-Dnodes(Dprimel))/(Dnodes(Dprimeh)-Dnodes(Dprimel))
        else
            Dprimel=nearestDnode-1
            Dprimeh=nearestDnode
            weightDprimel=1-(Dprime-Dnodes(Dprimel))/(Dnodes(Dprimeh)-Dnodes(Dprimel))
        end if
    end if
    if ((a-D*(1-thetaholder)*exp(hpholder))<0) then
    consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)*exp(hpholder)-thetaholder*Dcurrent*exp(hpholder)-aprimeholder-(1+rborrow)*(1-thetaholder)*D*exp(hpholder))
    else
    consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)*exp(hpholder)-thetaholder*Dcurrent*exp(hpholder)-aprimeholder-(1+r)*(1-thetaholder)*D*exp(hpholder))
    end if

    !need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.  Find nearest grid points and linearly interpolate.
    if (consumption>0) then
        if (elasticity .ne. 1) then
            valfuncadjust2=((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))**(1-elasticity)/(1-elasticity)
        else
            valfuncadjust2=log((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))
        end if
        EVholder=0
        EVholder2=0
        if (t<=Tretire) then
            do j=1,hpgridsize
                nearestanode=minloc(abs(aprime-anodes),1)
                if (aprime==anodes(nearestanode)) then
                    if (nearestanode<agridsize) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1.0
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=0.0
                    end if
                else
                    if (aprime-anodes(nearestanode)>0) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    end if
                end if

                do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,i,j,t+1)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,i,j,t+1)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,i,j,t+1)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,i,j,t+1)
                end do
            end do
            valfuncadjust2=valfuncadjust2+beta2*EVholder
        else
            do j=1,hpgridsize
                nearestanode=minloc(abs(aprime-anodes),1)
                if (aprime==anodes(nearestanode)) then
                    if (nearestanode<agridsize) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1.0
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=0.0
                    end if
                else
                    if (aprime-anodes(nearestanode)>0) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    end if
                end if

                EVholder=EVholder+Probhp(hpindex,j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,j,t+1)
                EVholder=EVholder+Probhp(hpindex,j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,j,t+1)
                EVholder=EVholder+Probhp(hpindex,j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,j,t+1)
                EVholder=EVholder+Probhp(hpindex,j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,j,t+1)

                EVholder2=EVholder2+Probhp(hpindex,j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,j,Tdie+1)
                EVholder2=EVholder2+Probhp(hpindex,j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,j,Tdie+1)
                EVholder2=EVholder2+Probhp(hpindex,j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,j,Tdie+1)
                EVholder2=EVholder2+Probhp(hpindex,j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,j,Tdie+1)
            end do
        valfuncadjust2=valfuncadjust2+beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
    end if


else
    valfuncadjust2=-1.0e13
end if
else
 valfuncadjust2=-1.0e13
end if



valfuncadjust2=-valfuncadjust2  !powell minimizes



END FUNCTION valfuncadjust2


REAL(8) FUNCTION valfuncrent2(p,state)  !p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
USE SHARE
IMPLICIT NONE
REAL(8), DIMENSION(2) :: p
REAL(8), DIMENSION(5) :: state
REAL(8) :: a, D, z, hp, t, zindex,hpindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder,hpholder,thetaholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)

hpholder=hpnodes(hpindex)
thetaholder=thetanodes(hpindex)

currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dnodes(Dgridsize) .AND. Dcurrent>=0) then
    if ((a-D*(1-thetaholder)*exp(hpholder))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)*exp(hpholder)-r_rental*Dcurrent*exp(rentelasticity*hpholder)-aprimeholder-(1+rborrow)*(1-thetaholder)*D*exp(hpholder))
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)*exp(hpholder)-r_rental*Dcurrent*exp(rentelasticity*hpholder)-aprimeholder-(1+r)*(1-thetaholder)*D*exp(hpholder))
    end if

    !need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.  Find nearest grid points and linearly interpolate.
    if (consumption>0) then

        if (elasticity .ne. 1) then
            valfuncrent2=((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))**(1-elasticity)/(1-elasticity)
        else
            valfuncrent2=log((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))
        end if
        EVholder=0
        EVholder2=0
        if (t<=Tretire) then
            do j=1,hpgridsize
                nearestanode=minloc(abs(aprime-anodes),1)
                if (aprime==anodes(nearestanode)) then
                    if (nearestanode<agridsize) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1.0
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=0.0
                    end if
                else
                    if (aprime-anodes(nearestanode)>0) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    end if
                end if

                do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*weightaprimel*EV(aprimel,1,i,j,t+1)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,1,i,j,t+1)
                end do
            end do
            valfuncrent2=valfuncrent2+beta2*EVholder
        else
            do j=1,hpgridsize
                nearestanode=minloc(abs(aprime-anodes),1)
                if (aprime==anodes(nearestanode)) then
                    if (nearestanode<agridsize) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1.0
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=0.0
                    end if
                else
                    if (aprime-anodes(nearestanode)>0) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    end if
                end if

                EVholder=EVholder+Probhp(hpindex,j)*weightaprimel*EV(aprimel,1,zindex,j,t+1)
                EVholder=EVholder+Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,1,zindex,j,t+1)

                EVholder2=EVholder2+Probhp(hpindex,j)*weightaprimel*EV(aprimel,1,zindex,j,Tdie+1)
                EVholder2=EVholder2+Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,1,zindex,j,Tdie+1)
            end do
        valfuncrent2=valfuncrent2+beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if


    else
        valfuncrent2=-1.0e13
    end if
else
     valfuncrent2=-1.0e13
end if



valfuncrent2=-valfuncrent2  !powell minimizes



END FUNCTION valfuncrent2






REAL(8) FUNCTION valfuncrent3(p,state)  !p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
USE SHARE
IMPLICIT NONE
REAL(8), DIMENSION(2) :: p
REAL(8), DIMENSION(5) :: state
REAL(8) :: a, D, z, hp, t, zindex,hpindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder,hpholder,thetaholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)

hpholder=hpnodes(hpindex)
thetaholder=thetanodes(hpindex)

currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dnodes(Dgridsize) .AND. Dcurrent>=0) then
    if ((a-D*(1-thetaholder)*exp(hpholder))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)*exp(hpholder)-r_rental*Dcurrent*exp(rentelasticity*hpholder)-aprimeholder-(1+rborrow)*(1-thetaholder)*D*exp(hpholder))
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)*exp(hpholder)-r_rental*Dcurrent*exp(rentelasticity*hpholder)-aprimeholder-(1+r)*(1-thetaholder)*D*exp(hpholder))
    end if

    !need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.  Find nearest grid points and linearly interpolate.
    if (consumption>0) then
        if (elasticity .ne. 1) then
            valfuncrent3=((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))**(1-elasticity)/(1-elasticity)
        else
            valfuncrent3=log((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))
        end if
        EVholder=0
        EVholder2=0
        if (t<=Tretire) then
            do j=1,hpgridsize
                nearestanode=minloc(abs(aprime-anodes),1)
                if (aprime==anodes(nearestanode)) then
                    if (nearestanode<agridsize) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1.0
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=0.0
                    end if
                else
                    if (aprime-anodes(nearestanode)>0) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    end if
                end if

                do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*weightaprimel*EV(aprimel,1,i,j,t+1)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,1,i,j,t+1)
                end do
            end do
            valfuncrent3=valfuncrent3+beta2*EVholder
        else
                do j=1,hpgridsize
                    nearestanode=minloc(abs(aprime-anodes),1)
                    if (aprime==anodes(nearestanode)) then
                        if (nearestanode<agridsize) then
                            aprimel=nearestanode
                            aprimeh=nearestanode+1
                            weightaprimel=1.0
                        else
                            aprimel=nearestanode-1
                            aprimeh=nearestanode
                            weightaprimel=0.0
                        end if
                    else
                        if (aprime-anodes(nearestanode)>0) then
                            aprimel=nearestanode
                            aprimeh=nearestanode+1
                            weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                        else
                            aprimel=nearestanode-1
                            aprimeh=nearestanode
                            weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                        end if
                    end if

                    EVholder=EVholder+Probhp(hpindex,j)*weightaprimel*EV(aprimel,1,zindex,j,t+1)
                    EVholder=EVholder+Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,1,zindex,j,t+1)

                    EVholder2=EVholder2+Probhp(hpindex,j)*weightaprimel*EV(aprimel,1,zindex,j,Tdie+1)
                    EVholder2=EVholder2+Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,1,zindex,j,Tdie+1)
                end do
        valfuncrent3=valfuncrent3+beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if


    else
        valfuncrent3=-1.0e13
    end if
else
     valfuncrent3=-1.0e13
end if

write(*,*) "A:s"
write(*,*) aprimel,aprimeh,weightaprimel



valfuncrent3=-valfuncrent3  !powell minimizes



END FUNCTION valfuncrent3





REAL(8) FUNCTION valfuncnoadjust(aprime,a,Dindex,zindex,hpindex,t)  !p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
USE SHARE
IMPLICIT NONE
REAL(8), DIMENSION(2) :: p
REAL(8) :: a, z, hp, t, zindex,hpindex,Dindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel, EVholder, EVholder2
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: Vadjustholder,Vnoadjustholder
REAL(8) :: Vadjustretireholder,Vnoadjustretireholder
REAL(8) :: hpholder,thetaholder
!write(*,*) state

!write(*,*) state

currentincome=income(zindex,t)

hpholder=hpnodes(hpindex)
thetaholder=thetanodes(hpindex)

Dcurrent=Dnodes(Dindex)
Dprime=Dnodes(Dindex)
aprimeholder=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize)) then

    !using that Dprime=D
    if ((a-Dcurrent*(1-thetaholder)*exp(hpholder))<0) then
        consumption=(currentincome+(1+rborrow)*a+Dprime*(1-delta)*exp(hpholder)-thetaholder*Dprime*exp(hpholder)-aprimeholder-(1+rborrow)*(1-thetaholder)*Dprime*exp(hpholder))
    else
        consumption=(currentincome+(1+r)*a+Dprime*(1-delta)*exp(hpholder)-thetaholder*Dprime*exp(hpholder)-aprimeholder-(1+r)*(1-thetaholder)*Dprime*exp(hpholder))
    end if

    !need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.  Find nearest grid points and linearly interpolate.
    if (consumption>0) then
        if (elasticity .ne. 1) then
            valfuncnoadjust=((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))**(1-elasticity)/(1-elasticity)
        else
            valfuncnoadjust=log((consumption+0)**elasticity2*(Dcurrent+0)**(1-elasticity2))
        end if

        EVholder=0
        EVholder2=0
        if (t<=Tretire) then
            do j=1,hpgridsize
                nearestanode=minloc(abs(aprime-anodes),1)
                if (aprime==anodes(nearestanode)) then
                    if (nearestanode<agridsize) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1.0
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=0.0
                    end if
                else
                    if (aprime-anodes(nearestanode)>0) then
                        aprimel=nearestanode
                        aprimeh=nearestanode+1
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    else
                        aprimel=nearestanode-1
                        aprimeh=nearestanode
                        weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                    end if
                end if

                do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*weightaprimel*EV(aprimel,Dindex,i,j,t+1)
                    EVholder=EVholder+Probz(zindex,i)*Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,Dindex,i,j,t+1)
                end do
            end do
        valfuncnoadjust=valfuncnoadjust+beta2*EVholder
        else
                do j=1,hpgridsize
                    nearestanode=minloc(abs(aprime-anodes),1)
                    if (aprime==anodes(nearestanode)) then
                        if (nearestanode<agridsize) then
                            aprimel=nearestanode
                            aprimeh=nearestanode+1
                            weightaprimel=1.0
                        else
                            aprimel=nearestanode-1
                            aprimeh=nearestanode
                            weightaprimel=0.0
                        end if
                    else
                        if (aprime-anodes(nearestanode)>0) then
                            aprimel=nearestanode
                            aprimeh=nearestanode+1
                            weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                        else
                            aprimel=nearestanode-1
                            aprimeh=nearestanode
                            weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
                        end if
                    end if

                    EVholder=EVholder+Probhp(hpindex,j)*weightaprimel*EV(aprimel,Dindex,zindex,j,t+1)
                    EVholder=EVholder+Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,Dindex,zindex,j,t+1)

                    EVholder2=EVholder2+Probhp(hpindex,j)*weightaprimel*EV(aprimel,Dindex,zindex,j,Tdie+1)
                    EVholder2=EVholder2+Probhp(hpindex,j)*(1-weightaprimel)*EV(aprimeh,Dindex,zindex,j,Tdie+1)
                end do
        valfuncnoadjust=valfuncnoadjust+beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if



    else
        valfuncnoadjust=-1.0e13
    end if
else
     valfuncnoadjust=-1.0e13
end if



valfuncnoadjust=-valfuncnoadjust  !powell minimizes



END FUNCTION valfuncnoadjust











SUBROUTINE amoeba(stateholder,p,y,ftol,valfuncadjust)
    USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
    IMPLICIT NONE

    INTERFACE
        FUNCTION valfuncadjust(x,statehold)
        USE nrtype
        IMPLICIT NONE
        REAL(SP), DIMENSION(2), INTENT(IN) :: x
        REAL(8), DIMENSION(5), INTENT(IN) :: statehold
        REAL(SP) :: valfuncadjust
        END FUNCTION valfuncadjust
    END INTERFACE
    INTEGER(I4B) :: iter
    REAL(SP), INTENT(IN) :: ftol
    REAL(SP), DIMENSION(3), INTENT(INOUT) :: y
    REAL(SP), DIMENSION(3,2), INTENT(INOUT) :: p
    REAL(SP), DIMENSION(3) :: yholder
    REAL(SP), DIMENSION(3,2) :: pholder
    REAL(SP), DIMENSION(5), INTENT(IN) :: stateholder
    INTEGER(I4B), PARAMETER :: ITMAX=5000000
    REAL(SP), PARAMETER :: TINY=1.0e-10
    INTEGER(I4B) :: ihi,ndim
    REAL(SP), DIMENSION(size(p,2)) :: psum
    call amoeba_private
    CONTAINS
!BL
    SUBROUTINE amoeba_private
    IMPLICIT NONE
    INTEGER(I4B) :: i,ilo,inhi
    REAL(SP) :: rtol,ysave,ytry,ytmp
    ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
    iter=0
    psum(:)=sum(p(:,:),dim=1)

    pholder=p
    yholder=y

    !write(*,*) "11",stateholder
    do
        ilo=iminloc(y(:))
        ihi=imaxloc(y(:))
        ytmp=y(ihi)
        y(ihi)=y(ilo)
        inhi=imaxloc(y(:))
        y(ihi)=ytmp

        !write(*,*) y(ihi), y(ilo)

        rtol=2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
        if (rtol < ftol) then
            call swap(y(1),y(ilo))
            call swap(p(1,:),p(ilo,:))
            RETURN
        end if
        if (iter >= ITMAX) then
        write(*,*) "stateholder", stateholder
        write(*,*) "p", p
        write(*,*) "pholder", pholder
        write(*,*) "y", y
        write(*,*) "yholder", yholder
        call nrerror('ITMAX exceeded in amoeba')
        end if
        !write(*,*) "2",stateholder

        ytry=amotry(-1.0_dp,stateholder)
        iter=iter+1
        if (ytry <= y(ilo)) then

        !write(*,*) stateholder

            ytry=amotry(2.0_dp,stateholder)
            iter=iter+1
        else if (ytry >= y(inhi)) then
            ysave=y(ihi)

            !write(*,*) "3", stateholder

            ytry=amotry(0.5_dp,stateholder)
            iter=iter+1
            if (ytry >= ysave) then
                p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
                do i=1,ndim+1

                !write(*,*) "4", stateholder

                    if (i /= ilo) y(i)=valfuncadjust(p(i,:),stateholder)
                    !write(*,*) y(i)
                end do
                iter=iter+ndim
                psum(:)=sum(p(:,:),dim=1)
            end if
        end if
    end do
    END SUBROUTINE amoeba_private
!BL
    FUNCTION amotry(fac,stateholder)
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: fac
    REAL(SP) :: amotry
    REAL(SP) :: fac1,fac2,ytry
    REAL(8), DIMENSION(5) :: stateholder
    REAL(SP), DIMENSION(size(p,2)) :: ptry
    fac1=(1.0_dp-fac)/ndim
    fac2=fac1-fac
    ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
    !write(*,*) "private", stateholder
    ytry=valfuncadjust(ptry,stateholder)
    !write(*,*) "done?"
    if (ytry < y(ihi)) then
        y(ihi)=ytry
        psum(:)=psum(:)-p(ihi,:)+ptry(:)
        p(ihi,:)=ptry(:)
    end if
    amotry=ytry
    END FUNCTION amotry
END SUBROUTINE amoeba





FUNCTION brentnew(ax,bx,cx,value,tol,aholder,Dholder,zholder,hpholder,tholder,xmin)
USE nrtype; USE nrutil, ONLY : nrerror
USE share
IMPLICIT NONE
REAL(SP), INTENT(IN) :: ax,bx,cx,tol,aholder,Dholder,zholder,hpholder,tholder
REAL(SP), INTENT(OUT) :: xmin
REAL(SP) :: brentnew
EXTERNAL :: value
REAL(SP) :: value



INTEGER(I4B), PARAMETER :: ITMAX=1000
REAL(SP), PARAMETER :: CGOLD=0.3819660_dp,ZEPS=1.0e-3_dp*epsilon(ax)
INTEGER(I4B) :: iter
REAL(SP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r2,tol1,tol2,u,v,w,x,xm
a=min(ax,cx)
b=max(ax,cx)
v=bx
w=v
x=v
e=0.0
fx=value(x,aholder,Dholder,zholder,hpholder,tholder)
fv=fx
fw=fx
do iter=1,ITMAX
    xm=0.5_dp*(a+b)
    tol1=tol*abs(x)+ZEPS
    tol2=2.0_dp*tol1
    if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
        xmin=x
        brentnew=fx
        RETURN
    end if
    if (abs(e) > tol1) then
        r2=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r2
        q=2.0_dp*(q-r2)
        if (q > 0.0) p=-p
        q=abs(q)
        etemp=e
        e=d
        if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
            p <= q*(a-x) .or. p >= q*(b-x)) then
            e=merge(a-x,b-x, x >= xm )
            d=CGOLD*e
        else
            d=p/q
            u=x+d
            if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
        end if
    else
        e=merge(a-x,b-x, x >= xm )
        d=CGOLD*e
    end if
    u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
    fu=value(u,aholder,Dholder,zholder,hpholder,tholder)
    if (fu <= fx) then
        if (u >= x) then
            a=x
        else
            b=x
        end if
        call shft(v,w,x,u)
        call shft(fv,fw,fx,fu)
    else
        if (u < x) then
            a=u
        else
            b=u
        end if
        if (fu <= fw .or. w == x) then
            v=w
            fv=fw
            w=u
            fw=fu
        else if (fu <= fv .or. v == x .or. v == w) then
            v=u
            fv=fu
        end if
    end if
end do
call nrerror('brentnew: exceed maximum iterations')
CONTAINS
!BL
SUBROUTINE shft(a,b,c,d)
REAL(SP), INTENT(OUT) :: a
REAL(SP), INTENT(INOUT) :: b,c
REAL(SP), INTENT(IN) :: d
a=b
b=c
c=d
END SUBROUTINE shft
END FUNCTION brentnew




FUNCTION brentmindist(ax,bx,cx,value,aftertax,tol,xmin)
USE nrtype; USE nrutil, ONLY : nrerror
USE share
IMPLICIT NONE
REAL(SP), INTENT(IN) :: ax,bx,cx,tol
REAL(8), DIMENSION(numhouseholds) :: aftertax
REAL(SP), INTENT(OUT) :: xmin
REAL(SP) :: brentmindist
EXTERNAL :: value
REAL(SP) :: value



INTEGER(I4B), PARAMETER :: ITMAX=200
REAL(SP), PARAMETER :: CGOLD=0.3819660_dp,ZEPS=1.0e-3_dp*epsilon(ax)
INTEGER(I4B) :: iter
REAL(SP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r2,tol1,tol2,u,v,w,x,xm
a=min(ax,cx)
b=max(ax,cx)
v=bx
w=v
x=v
e=0.0
fx=value(x,aftertax)
fv=fx
fw=fx
do iter=1,ITMAX
    xm=0.5_dp*(a+b)
    tol1=tol*abs(x)+ZEPS
    tol2=2.0_dp*tol1
    if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
        xmin=x
        brentmindist=fx
        RETURN
    end if
    if (abs(e) > tol1) then
        r2=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r2
        q=2.0_dp*(q-r2)
        if (q > 0.0) p=-p
        q=abs(q)
        etemp=e
        e=d
        if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
            p <= q*(a-x) .or. p >= q*(b-x)) then
            e=merge(a-x,b-x, x >= xm )
            d=CGOLD*e
        else
            d=p/q
            u=x+d
            if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
        end if
    else
        e=merge(a-x,b-x, x >= xm )
        d=CGOLD*e
    end if
    u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
    fu=value(u,aftertax)
    if (fu <= fx) then
        if (u >= x) then
            a=x
        else
            b=x
        end if
        call shft(v,w,x,u)
        call shft(fv,fw,fx,fu)
    else
        if (u < x) then
            a=u
        else
            b=u
        end if
        if (fu <= fw .or. w == x) then
            v=w
            fv=fw
            w=u
            fw=fu
        else if (fu <= fv .or. v == x .or. v == w) then
            v=u
            fv=fu
        end if
    end if
end do
call nrerror('brentnew: exceed maximum iterations')
CONTAINS
!BL
SUBROUTINE shft(a,b,c,d)
REAL(SP), INTENT(OUT) :: a
REAL(SP), INTENT(INOUT) :: b,c
REAL(SP), INTENT(IN) :: d
a=b
b=c
c=d
END SUBROUTINE shft
END FUNCTION brentmindist


! Computes the stationary distribution of an ergodic matrix
! using state space reduction. We write the transition
! matrix as P = [(A C) / (B D)], A being n-1xn-1, then
! iterate downwards. Got this from the R. Feres notes
FUNCTION ergodic(mat, states)
    IMPLICIT NONE
    INTEGER, intent(in) :: states
    INTEGER :: iter, j, j2, dim1, dim2
    REAL(8), DIMENSION(states,states), intent(in) :: mat
    REAL(8), DIMENSION(states,states) :: P
    REAL(8) :: inverse, total
    REAL(8), DIMENSION(states) :: ergodic

    iter  = SIZE(mat, 1)
    P = mat

    do while (iter > 1)
        dim1 = iter - 1
        dim2 = iter - 1
        inverse = SUM(P(iter, 1:dim1))
        ! Storing resulting vectors in P
        P(1:dim1, iter) = P(1:dim1, iter)/inverse
        ! Using the censored Markov identity P_n = A + B(1-D)C
        do while (dim2 > 0)
            P(1:dim1,dim2) = P(1:dim1,dim2) + &
            P(iter,dim2)*P(1:dim1,iter)
            dim2 = dim2 - 1
        end do
        iter = iter - 1
    end do

    ergodic(1) = 1
    j=2
    !Back out the stationary probabilities
    do while (j <= size(P, 1))
        j2 = j - 1
        !write(*,*) DOT_PRODUCT(ergodic(1:j2),P(1:j2,j))
        ergodic(j) = DOT_PRODUCT(ergodic(1:j2),P(1:j2,j))
        j = j + 1
    end do
    total = sum(ergodic)
    ergodic = ergodic/total
    RETURN

END FUNCTION



   REAL(8) FUNCTION cdfnormal (x)
!
!*******************************************************************************
!
!! cdfnormal evaluates the Normal 01 CDF.
!
!
!  Reference:
!
!    A G Adams,
!    Areas Under the Normal Curve,
!    Algorithm 39,
!    Computer j.,
!    Volume 12, pages 197-198, 1969.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Output, real CDF, the value of the CDF.
!
  implicit none
!
  real, parameter :: a1 = 0.398942280444E+00
  real, parameter :: a2 = 0.399903438504E+00
  real, parameter :: a3 = 5.75885480458E+00
  real, parameter :: a4 = 29.8213557808E+00
  real, parameter :: a5 = 2.6813121679E+00
  real, parameter :: a6 = 48.6959930692E+00
  real, parameter :: a7 = 5.92885724438E+00
  real, parameter :: b0 = 0.398942280385E+00
  real, parameter :: b1 = 3.8052E-08
  real, parameter :: b2 = 1.00000615302E+00
  real, parameter :: b3 = 3.98064794E-04
  real, parameter :: b4 = 1.98615381364E+00
  real, parameter :: b5 = 0.151679116635E+00
  real, parameter :: b6 = 5.29330324926E+00
  real, parameter :: b7 = 4.8385912808E+00
  real, parameter :: b8 = 15.1508972451E+00
  real, parameter :: b9 = 0.742380924027E+00
  real, parameter :: b10 = 30.789933034E+00
  real, parameter :: b11 = 3.99019417011E+00
  !real cdfnormal
  real q
  real(8), intent(in) :: x
  real y
!
!  |X| <= 1.28.
!

  if ( abs ( x ) <= 1.28 ) then

    y = 0.5E+00 * x**2

    q = 0.5E+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )

!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7E+00 ) then

    y = 0.5E+00 * x**2

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0E+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0E+00 ) then
    cdfnormal = q
  else
    cdfnormal = 1.0E+00 - q
  end if

  return
end
