
module share
    REAL(8), PARAMETER :: delta=.022 !depreciation rate of durable  was .03
    REAL(8), parameter :: elasticity= 2  !elasticity of substitution

    REAL(8), parameter :: r=.024
    REAL(8), parameter :: rborrow=.024
    REAL(8), parameter :: thetamatlab=.25
    REAL(8), parameter :: theta=1-(1-thetamatlab)*(1-delta)/(1+r) !down payment
    REAL(8) :: relutility
    REAL(8), parameter :: F = 0.05 ! fixed cost
    

    REAL(8), parameter :: rent=1-(1-delta)/(1+r)
    REAL(8) :: beta2   !Quarterly discount factor
    REAL(8):: beta2retire
    REAL(8):: elasticity2 !non-durable share
    REAL(8) :: psi
    REAL(8) :: r_rental
    REAL(8) :: r_rental_retire
    REAL(8) :: ret_wealth
    
    REAL(8) :: const

    !for now, just two point process
    REAL(8), parameter :: hpmin=0.0
    REAL(8), parameter :: hpmax=0.15 
    integer, parameter :: hpgridsize=2 !size of house price grid
    REAL(8), DIMENSION(hpgridsize) :: hpnodes
    REAL(8), DIMENSION(hpgridsize,hpgridsize) :: Probhp

    !idiosyncratic earnings shocks
    REAL(8), parameter :: sigma_z=.15
    REAL(8), parameter :: rho_z=.91
    REAL(8), parameter :: sigma_eta_init=(sigma_z**2.0/(1-rho_z**2.0))**.5
    REAL(8), parameter :: zmin=-2.5*(sigma_z**2.0/(1-rho_z**2.0))**.5
    REAL(8), parameter :: zmax=2.5*(sigma_z**2.0/(1-rho_z**2.0))**.5

 
    
    REAL(8), parameter :: Dmin=0    !minimum durable choice
    REAL(8), parameter :: Dmax=8  !maximum durable choice
    integer, parameter :: agridsize = 80 ! size of asset grid
    integer, parameter :: Dgridsize = 60 ! size of durable grid

    !REAL(8), parameter :: amin=(exp(hpmin)-exp(hpmax))*(1-theta)*Dmax    !minimum asset value (actually voluntary equity if theta!=1).  Note that if theta<1, we can have assets <0 even with amin=0 since a is vol. equity.  But with theta=1 & amin=0, no borrowing.
    !REAL(8), parameter :: amin=-.89
    REAL(8), parameter :: amin=-0.40/(1+r)  !.40 is roughly min y, so this is basically a no default condition
    !REAL(8), parameter :: amax2=8   !max asset
    REAL(8), parameter :: amax=40   !max asset

    integer, parameter :: zgridsize=5  !size of ar(1) income shock grid
    REAL(8), parameter :: borrowconstraint=amin
    

    !integer, parameter :: Tretire=35
    !integer, parameter :: Tdie=60

    !integer, parameter :: Tretire=100
    integer, parameter :: Tretire=35
    integer, parameter :: Tdie=Tretire+25

    REAL(8), dimension(Tdie) :: r_rental_t
    
    
    
    REAL(8), DIMENSION(agridsize) :: anodes  !Nodes for idiosyncratic assets (voluntary equity)
    REAL(8), DIMENSION(Dgridsize) :: Dnodes  !Nodes for idiosyncratic durable

    REAL(8), DIMENSION(zgridsize) :: Probinit,Probinitcum

    
    
    REAL(8), DIMENSION(zgridsize) :: znodes  !Nodes for idiosyncratic productivity
    REAL(8), DIMENSION(zgridsize,zgridsize) :: Probz,Probzcum  !while probability of being at any z will vary with time, transition probability will not, so max we need to loop over wont
    REAL(8), DIMENSION(zgridsize) :: maxexpectationz,minexpectationz  

    
    REAL(8), DIMENSION(Tretire) :: ageearnings
    REAL(8), DIMENSION(Tdie-Tretire) :: deathprob
    REAL(8), DIMENSION(zgridsize) :: predictedlifetimeearnings,retirementincome
    REAL(8), DIMENSION(zgridsize,Tdie) :: income
    REAL(8) :: averagelifetimeearnings
    
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Vnoadjust,  Vadjust,Vrent  !Conjectured value of not adjusting and adjusting, and new values of not adjusting and adjusting
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: achoiceadjust,  achoicenoadjust,achoicerent,  Dchoiceadjust, Dchoicenoadjust,Dchoicerent, achoice, Dchoice,cchoice,mpc_pol  !Policy functions when solving problem
    REAL(8) :: priceshock
   
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: rentalindicator
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: EV
   
    
    INTEGER, parameter :: numhouseholds=15000    !Size of household panel
    !integer, parameter :: burnin=1000    !Initial periods dropped in simulation
    
    
    REAL(8) :: diffmoments,diffmoments1,diffmoments2,diffmoments3
    REAL(8), parameter :: difftol=.01
    REAL(8) :: numiter

    REAL(8), parameter :: costlyequity=0

    REAL(8), parameter :: rentelasticity=1

    REAL(8), DIMENSION(Tdie+1) :: Ps, gs,P_exp
    REAL(8) :: g

    integer :: totalbubblelength=2
    REAL(8), parameter :: rho=0.5

    REAL(8) :: P1
    REAL(8) :: P2

    INTEGER :: bubbletime

    REAL(8) :: age

    integer :: ageatbubblestart

    !REAL(8), DIMENSION(Tdie,Tdie+totalbubblelength-1):: hlevelnobubble,hlevelwithbubble,clevelnobubble,clevelwithbubble,ribubble,rinobubble,Bmatrix
    REAL(8), DIMENSION(Tdie,Tdie)::Bmatrix
    REAL(8), DIMENSION(Tdie,2*Tdie) :: hlevelwithbubble,clevelwithbubble,sufficientstatwithbubble,hlevelnobubble,clevelnobubble,rinobubble
    REAL(8), DIMENSION(Tdie,1):: irfh, irfc, irfri, totalhlevelwithbubble,totalhlevelnobubble,totalclevelwithbubble,totalclevelnobubble,totalsufficientstatwithbubble,totalrilevelnobubble
    
    REAL(8) :: switch

end module share


program durables  !Program shell
use share
USE OMP_LIB
implicit none
integer :: i,j,k,l,m,n,t, iter
EXTERNAL cdfnormal  !Fnc to compute normal cdf
REAL(8) :: cdfnormal
REAL(8) :: w,ll,rr,ww,niter,wholder
REAL(8) :: foundmin, foundmax
REAL(8) :: timestart, timeend
REAL(8) :: dd,m2,m2young,m2middle,m2old
REAL(8) :: finwealth
REAL(8) :: i2,j2,k2,l2,t2
REAL(8) :: param1, param2,param3,param4,param5,param6
REAL(8) giter


ALLOCATE(Vnoadjust(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie+1),  Vadjust(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie+1),Vrent(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie+1))  !Conjectured value of not adjusting and adjusting, and new values of not adjusting and adjusting
ALLOCATE(mpc_pol(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie),achoiceadjust(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie),  achoicenoadjust(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie),achoicerent(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie),  Dchoiceadjust(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie), Dchoicenoadjust(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie),Dchoicerent(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie))
ALLOCATE(achoice(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie), Dchoice(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie),cchoice(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie))  !Policy functions when solving problem
ALLOCATE(rentalindicator(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie))
ALLOCATE(EV(agridsize,Dgridsize,zgridsize,totalbubblelength,Tdie+1))


do giter=1,3



param1=1
param2=3
param3=1
param4=1
param5=3
param6=3

beta2=.92+.002*param1
beta2retire=beta2

elasticity2=.83+.01*param2
psi=2100+200*param3
r_rental=.059+.0015*param4
r_rental_retire=.0442+.001*param5
ret_wealth=.75+.5*param6


const = (elasticity2**elasticity2)*((1.0-elasticity2)**(1.0-elasticity2))

r_rental_t(1:Tretire)=r_rental
r_rental_t(Tretire+1:Tdie)=r_rental_retire


Ps=1.0
gs(1)=0
gs(2)=0-(.1)**(4*giter)


hlevelnobubble=0
hlevelwithbubble=0
clevelnobubble=0
clevelwithbubble=0
sufficientstatwithbubble=0
rinobubble=0


Probhp(1,1)=1.0
Probhp(1,2)=0.0
Probhp(2,1)=0.0
Probhp(2,2)=1.0




!CALL OMP_SET_NUM_THREADS(12)

numiter=0




diffmoments=1
do while (diffmoments>difftol)
    numiter=numiter+1
    !write(*,*) "numiter", numiter

    !write(*,*) elasticity2, r_rental, beta2



    !OPEN (UNIT=1, FILE="vfunc.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    !1 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

    ! FILL IN ALL THE GRID POINTS

    do i=1,zgridsize
        znodes(i)=zmin+((zmax-zmin)/(zgridsize-1))*(1.0*i-1.0)
        !write(*,*) znodes(i)
    end do
 
    DO i=1,hpgridsize
        hpnodes(i)=hpmin+ ((hpmax -hpmin)/(hpgridsize-1))*(1.0*i-1.0)
    END DO


    DO i=1,agridsize
        anodes(i)=(1.0/(agridsize-1))*(1.0*i-1.0)
    end do
    

    do i=1,agridsize
        anodes(i)=exp(log(amax-amin+1)*anodes(i))+amin-1.0
    end do


   ! write(*,*) "anodes2",anodes



     DO i=1,Dgridsize
        dnodes(i)=(1.0/(dgridsize-1))*(1.0*i-1.0)
    end do
    

    do i=1,dgridsize
        dnodes(i)=exp(log(dmax-dmin+1)*dnodes(i))+dmin-1.0
    end do

  !  write(*,*) "dnodes2",dnodes
  
    ! create transition matrix for log idiosyncratic labor shock using Tauchen 86
    w=znodes(2)-znodes(1)
    do j=1,zgridsize
        Probz(j,1)=cdfnormal((znodes(1)-rho_z*znodes(j)+w/2)/(sigma_z))
        Probz(j,zgridsize)=1-cdfnormal((znodes(zgridsize)-rho_z*znodes(j)-w/2)/(sigma_z))
        do k=2,zgridsize-1
            Probz(j,k)=cdfnormal((znodes(k)-rho_z*znodes(j)+w/2)/(sigma_z))-cdfnormal((znodes(k)-rho_z*znodes(j)-w/2)/(sigma_z))
        end do
    end do

!Truncating the transistion probabilities of transitory shock. Want to minimize number of gridpoints have to compute conditional expectation
    minexpectationz=1.0
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

    w=znodes(2)-znodes(1)
    Probinit(1)=cdfnormal((znodes(1)+w/2)/sigma_eta_init)
    Probinit(zgridsize)=1-cdfnormal((znodes(zgridsize)-w/2)/sigma_eta_init)
    do j=2,zgridsize-1
        Probinit(j)=cdfnormal((znodes(j)+w/2)/sigma_eta_init)-cdfnormal((znodes(j)-w/2)/sigma_eta_init)
    end do
    
    do j=1,zgridsize
        Probinitcum(j)=sum(Probinit(1:j))
    end do

   
! Kai portion of earnings proccess
 ageearnings(1)=-0.520437830949535
 ageearnings(2)=-0.445723729949535
 ageearnings(3)=-0.378341829949534
 ageearnings(4)=-0.317627929949534
 ageearnings(5)=-0.262959589949535
 ageearnings(6)=-0.213756129949534
 ageearnings(7)=-0.169478629949534
 ageearnings(8)=-0.129629929949534
 ageearnings(9)=-0.0937546299495352
 ageearnings(10)=-0.0614390899495344
 ageearnings(11)=-0.0323114299495341
 ageearnings(12)=-0.00604152994953496
 ageearnings(13)=0.0176589700504652
 ageearnings(14)=0.0390366700504658
 ageearnings(15)=0.0582964100504647
 ageearnings(16)=0.0756012700504658
 ageearnings(17)=0.0910725700504656
 ageearnings(18)=0.104789870050465
 ageearnings(19)=0.116790970050464
 ageearnings(20)=0.127071910050466
 ageearnings(21)=0.135586970050465
 ageearnings(22)=0.142248670050465
 ageearnings(23)=0.146927770050466
 ageearnings(24)=0.149453270050464
 ageearnings(25)=0.149612410050466
 ageearnings(26)=0.147150670050465
 ageearnings(27)=0.141771770050466
 ageearnings(28)=0.133137670050466
 ageearnings(29)=0.120868570050465
 ageearnings(30)=0.104542910050465
 ageearnings(31)=0.0836973700504653
 ageearnings(32)=0.0578268700504655
 ageearnings(33)=0.0263845700504652
 ageearnings(34)=-0.0112181299495353
 ageearnings(35)=-0.0556115899495339


   ! ageearnings=0



! death probabilityes (starting at age 60 or 65)
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

    deathprob=0



    averagelifetimeearnings=0.0
    !retirement regression needs to first run matlab program simulateearningsprocess
    do j=1,zgridsize
        predictedlifetimeearnings(j)=0.3083*(znodes(j))
        predictedlifetimeearnings(j)=exp(predictedlifetimeearnings(j))/exp(averagelifetimeearnings)
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
      !  write(*,*) "r i", retirementincome(j)
    end do


    !retirementincome=exp(ageearnings(35))
    !retirementincome=0

    do j=1,zgridsize
        do t=1,Tretire
            income(j,t)=exp(znodes(j)+ageearnings(t))
            income(j,t)=income(j,t)/sum(Probinit*exp(znodes))
        end do
        do t=Tretire+1,Tdie
           income(j,t)=retirementincome(j)
        end do
      !  do t=Tretire+1,Tdie
      !      income(j,t)=exp(znodes(j))
      !  end do
    end do
    income(:,Tretire+1)=income(:,Tretire+1)+income(:,Tretire)*ret_wealth

    

  !  write(*,*) "income"
  !  write(*,*) sum(income(:,1,:)),sum(income(:,2,:)),sum(income(:,3,:)),sum(income(:,4,:)),sum(income(:,5,:))




! Put in bequest motive in last period

    call cpu_time(timestart)
    !call solveretirementproblem
    !Vnoadjust(:,:,:,:,61)=0
    !Vadjust(:,:,:,:,61)=0
   ! Vrent(:,:,:,:,61)=0

    do i=1,agridsize
        do j=1,Dgridsize
            do k=1,totalbubblelength
             
                finwealth=anodes(i)-(1-theta)*Dnodes(j)*Ps(k)
                wholder=(1+r)*finwealth+(1-delta)*Dnodes(j)*Ps(k)*(1-f)


                
                
                if (elasticity .ne. 1.0) then
                EV(i,j,:,k,Tdie+1)=1.0/(1.0-elasticity)*const**(1-elasticity)*rent**((1-elasticity2)*(elasticity-1))*psi*(income(:,Tdie)/r+wholder)**(1-elasticity)
                else
                EV(i,j,:,k,Tdie+1)=psi*log(income(:,Tdie)/r+wholder)
                end if
                if (j==1 .and. k==1) then
                !write(*,*) EV(i,j,1,k,Tdie+1)
                end if
                
            end do
        end do
    end do
    !EV=0
    !EV=max(-10000000000.0,EV)

   ! write(*,*) "totEV ", sum(EV)

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
    call cpu_time(timeend)
   ! write(*,*) (timeend-timestart)/8.0
    call cpu_time(timestart)
    call solveworkingproblem
    call cpu_time(timeend)
  !  write(*,*) (timeend-timestart)/8.0

   ! write(*,*) "MINIMUM", minval(Dchoice)


    diffmoments=0
   
    !write(*,*) agridsize,Dgridsize,zgridsize,hpgridsize
    t=0
    do i=1,agridsize
        do j=1,Dgridsize
            do k=1,zgridsize
                do l=1,totalbubblelength
                t=t+1
                t2=13
                k2=k
                l2=l
                   write(1,1) EV(i,j,k,l,t2), achoice(i,j,k,l,t2), dchoice(i,j,k,l,t2), cchoice(i,j,k,l,t2),rentalindicator(i,j,k,l,t2), anodes(i), dnodes(j), znodes(k),real(l), income(k,1), achoiceadjust(i,j,k,l,t2), dchoiceadjust(i,j,k,l,t2), achoicerent(i,j,k,l,t2), dchoicerent(i,j,k,l,t2), achoicenoadjust(i,j,k,l,t2),dchoicenoadjust(i,j,k,l,t2), Vadjust(i,j,k,l,t2), Vnoadjust(i,j,k,l,t2), Vrent(i,j,k,l,t2)

                t2=33
                   write(17,17) EV(i,j,k,l,t2), achoice(i,j,k,l,t2), dchoice(i,j,k,l,t2), cchoice(i,j,k,l,t2),rentalindicator(i,j,k,l,t2), anodes(i), dnodes(j), znodes(k),real(l), income(k,1), achoiceadjust(i,j,k,l,t2), dchoiceadjust(i,j,k,l,t2), achoicerent(i,j,k,l,t2), dchoicerent(i,j,k,l,t2), achoicenoadjust(i,j,k,l,t2),dchoicenoadjust(i,j,k,l,t2), Vadjust(i,j,k,l,t2), Vnoadjust(i,j,k,l,t2), Vrent(i,j,k,l,t2)
                end do
            end do
        end do
  end do


    !parameters are updated, then repeat whole big loop with solution/simulation to target moments.
    close(1)
end do

do ageatbubblestart=1,Tdie+totalbubblelength-1
call simulate
end do

do i=Tdie+totalbubblelength,2*Tdie
hlevelwithbubble(:,i)=hlevelwithbubble(:,Tdie+totalbubblelength-1)
clevelwithbubble(:,i)=clevelwithbubble(:,Tdie+totalbubblelength-1)
sufficientstatwithbubble(:,i)=sufficientstatwithbubble(:,Tdie+totalbubblelength-1)
end do






do t=1,60
Bmatrix=0
do i=t,60
    do j=1,61-t
        if (i+1-t==j) then
            Bmatrix(i,j)=1
        else
            Bmatrix(i,j)=0
        end if
    end do
end do
totalhlevelwithbubble(t,1)=sum((hlevelwithbubble(1:Tdie,1:Tdie)*Bmatrix))
totalhlevelnobubble(t,1)=sum(hlevelnobubble(:,1))

totalclevelwithbubble(t,1)=sum((clevelwithbubble(1:Tdie,1:Tdie)*Bmatrix))
totalclevelnobubble(t,1)=sum(clevelnobubble(:,1))

totalsufficientstatwithbubble(t,1)=sum((sufficientstatwithbubble(1:Tdie,1:Tdie)*Bmatrix))


end do


do t=2,60
    do i=t-1,1,-1
    totalhlevelwithbubble(t,1)=totalhlevelwithbubble(t,1)+hlevelwithbubble(i,Tdie+t-i)
    totalclevelwithbubble(t,1)=totalclevelwithbubble(t,1)+clevelwithbubble(i,Tdie+t-i)
    totalsufficientstatwithbubble(t,1)=totalsufficientstatwithbubble(t,1)+sufficientstatwithbubble(i,Tdie+t-i)
    end do
end do

write(*,*) "Housing IRF", gs(2)
do t=1,Tdie,5
write(*,*) (totalhlevelwithbubble(t,1)-totalhlevelnobubble(t,1))/totalhlevelnobubble(t,1),totalhlevelwithbubble(t,1),totalhlevelnobubble(t,1)
end do

write(*,*) "Consumption IRF", gs(2)
do t=1,Tdie,5
write(*,*) (totalclevelwithbubble(t,1)-totalclevelnobubble(t,1))/totalclevelnobubble(t,1),totalclevelwithbubble(t,1),totalclevelnobubble(t,1)
end do

write(*,*) "Residential Investment IRF"
t=1
write(*,*) (totalhlevelwithbubble(t,1)-(1-delta)*totalhlevelnobubble(t,1))/(totalhlevelnobubble(t,1)-(1-delta)*totalhlevelnobubble(t,1)),(totalhlevelwithbubble(t,1)-(1-delta)*totalhlevelnobubble(t,1)),(totalhlevelnobubble(t,1)-(1-delta)*totalhlevelnobubble(t,1))
do t=2,Tdie,5
write(*,*) (totalhlevelwithbubble(t,1)-(1-delta)*totalhlevelwithbubble(t-1,1))/(totalhlevelnobubble(t,1)-(1-delta)*totalhlevelnobubble(t-1,1)), (totalhlevelwithbubble(t,1)-(1-delta)*totalhlevelwithbubble(t-1,1)), (totalhlevelnobubble(t,1)-(1-delta)*totalhlevelnobubble(t-1,1))
!write(*,*) (totalrilevelwithbubble(t,1)-totalrilevelnobubble(t,1))/totalrilevelnobubble(t,1),totalrilevelwithbubble(t,1),totalrilevelnobubble(t,1)
end do

write(*,*) "Residential Investment Rate"
t=1
write(*,*) (totalhlevelwithbubble(t,1)-(1-delta)*totalhlevelnobubble(t,1))/totalhlevelnobubble(t,1)
do t=2,Tdie,5
write(*,*) (totalhlevelwithbubble(t,1)-(1-delta)*totalhlevelwithbubble(t-1,1))/totalhlevelnobubble(t,1)
!write(*,*) (totalrilevelwithbubble(t,1)-totalrilevelnobubble(t,1))/totalrilevelnobubble(t,1),totalrilevelwithbubble(t,1),totalrilevelnobubble(t,1)
end do





end do

!write(*,*) "Check last term", hlevelwithbubble(60,1), hlevelnobubble(60,1)



!write(*,*) "C sum", sum(cchoice)

end program durables






subroutine simulate
USE nrtype; USE nrutil
USE nr
USE share
USE OMP_LIB
IMPLICIT NONE

REAL(8), dimension(numhouseholds,5) :: currenthouseholdstate, newhouseholdstate, currenthouseholdstatewithbubble,newhouseholdstatewithbubble

REAL(8), dimension(numhouseholds,60) :: currenthouseholdstate_agebubblestate, newhouseholdstate_agebubblestate

REAL(8), dimension(numhouseholds,1) :: hstockwithrental,hstockwithrentalwithbubble

REAL(8), dimension(numhouseholds,Tdie) :: aggregatecontrib,mpc,consumption,consumptionwithbubble,consumptionwithbubblempc,consumptionhpshock,durableconsumption, currenttemp,currentperm,incomeholder,rentalind,rentalindwithbubble,durableinvestment,durableinvestmentwithbubble,actualwealth,financialwealth,housingnet,housinggross,actualwealthdividedbyincome,totalnetworthdividedbyincome,taxrate,welfare,welfarewithbubble,diffc, ratioc, qchg
REAL(8), dimension(numhouseholds,Tdie) :: changec,changetemp,changeperm,changey,changeyconditional,changed,changedconditional,changetempconditional,changepermconditional,demeanedindicator
REAL(8), dimension(numhouseholds,Tdie) :: alive
REAL(8), dimension(numhouseholds,Tdie) :: housingstock, rentalstock
REAL(8), dimension(numhouseholds,2) :: nearestnode
REAL(8), dimension(numhouseholds,Tdie,8) :: householdresultsholder
REAL(8) :: shock
integer, dimension(12) :: seedvalue
integer :: j,i,k,t
REAL(8) :: adjust
REAL(8) :: conditionalindicator
REAL(8) :: numrent
REAL(8) :: pretaxcalcholder,pretaxincome
REAL(8) :: tot1,tot2,tot3,adjustt,changedconditionalholder,changeyconditionalholder
REAL(8) :: numowntotal
REAL(8) :: rep
REAL(8) :: actualwealthtotal
REAL(8) :: exanteEVborn,exanteEVoverall,numobswelfare
REAL(8) :: medianincome, ratiocrent, ratiocown,diffcrent,diffcown,numrent_t,numown_t,numbuy_t,numsell_t
REAL(8), dimension(10):: numbins
REAL(8) :: leverage_t, a_t, leverageown_t,aown_t,leveragerent_t,arent_t

REAL(8) :: totalhousing, totalincome,totalfinancialassets,totalincomeowners


REAL(8), dimension(numhouseholds,Tdie) :: stateindicator1,stateindicator2,stateindicator5

REAL(8) :: bubbletimeholder
integer :: actualageatbubblestart

totalhousing=0
totalincome=0
totalfinancialassets=0
totalincomeowners=0



seedvalue=1
CALL random_seed(put=seedvalue(1:12))

! STATE DEFINITION
! state 1: liquid assets a
! state 2: durable assets d
! state 3: idiosyncratic income shock: z
! state 4: house price
! state 5: age: t

! If you wanted to play with changing the initial distribution of wealth
!currenthouseholdstate(:,1)=0  ! start with zero assets (could do something here with bequests)
!currenthouseholdstate(:,2)=0


do i=1,numhouseholds
    if (i/numhouseholds<Probinitcum(1)) then
        currenthouseholdstate(i,3)=1
    else
        do j=1,zgridsize-1
            if (i/numhouseholds > Probinitcum(j) .AND. i/numhouseholds<=Probinitcum(j+1)) then
                currenthouseholdstate(i,3)=j+1
            end if
        end do
    end if
end do

medianincome=income((zgridsize+1)/2,1)
numbins=0
do i=1,numhouseholds
    call random_number(shock)
    if (income(currenthouseholdstate(i,3),1)/medianincome<.4632) then 
            if (shock<.0957) then
            currenthouseholdstate(i,2)=2.9545
            currenthouseholdstate(i,1)=0.4003
            else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.0
            end if
            
            

            numbins(1)=numbins(1)+1
    elseif (income(currenthouseholdstate(i,3),1)/medianincome>=.4632 .and. income(currenthouseholdstate(i,3),1)/medianincome<1) then
        if (shock<.1111) then
            currenthouseholdstate(i,2)=1.2727
            currenthouseholdstate(i,1)=1.3628
            numbins(2)=numbins(2)+1
        else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.0023
            numbins(3)=numbins(3)+1
        end if
    elseif (income(currenthouseholdstate(i,3),1)/medianincome>=1 .and. income(currenthouseholdstate(i,3),1)/medianincome<1.5455) then
        if (shock<.2250) then
            currenthouseholdstate(i,2)=4.0909
            currenthouseholdstate(i,1)=2.734
            numbins(4)=numbins(4)+1
        else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.0139
            numbins(5)=numbins(5)+1
        end if
    else
        if (shock<.4286) then
            currenthouseholdstate(i,2)=3.0455
            currenthouseholdstate(i,1)=0.0909
            numbins(6)=numbins(6)+1
        else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.0998
            numbins(7)=numbins(7)+1
        end if
    end if
    !currenthouseholdstate(:,1)=0
end do
currenthouseholdstate(:,1)=currenthouseholdstate(:,1)*exp(ageearnings(1))/sum(exp(ageearnings)/Tretire)
currenthouseholdstate(:,2)=currenthouseholdstate(:,2)*exp(ageearnings(1))/sum(exp(ageearnings)/Tretire)
currenthouseholdstate(:,1)=0
currenthouseholdstate(:,2)=0


!write(*,*) numbins/numhouseholds




currenthouseholdstate(:,4)=1  !no bubble in baseline
newhouseholdstate(:,4)=1
currenthouseholdstatewithbubble=currenthouseholdstate
newhouseholdstatewithbubble=newhouseholdstate



call random_number(shock)
alive=1






numowntotal=0
numrent=0
actualwealthtotal=0

housingstock=0
rentalstock=0


!write(*,*) "ageatbubble start", ageatbubblestart


do t=1,Tdie

!write(*,*) "min 2", minval(currenthouseholdstate(:,2))
!write(*,*) "min 3", minval(currenthouseholdstatewithbubble(:,2))

if (ageatbubblestart<=Tdie) then
bubbletimeholder=min(max(t-ageatbubblestart+1,1),totalbubblelength)
currenthouseholdstatewithbubble(:,4)=bubbletimeholder
actualageatbubblestart=ageatbubblestart
else
actualageatbubblestart=Tdie+1-ageatbubblestart  !age is 60 - age for age > 60, need to account for households being born after bubble start
bubbletimeholder=min(2-actualageatbubblestart,totalbubblelength)
currenthouseholdstatewithbubble(:,4)=bubbletimeholder
end if


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


!write(*,*) "t", maxval((1+r)*(currenthouseholdstate(:,1)+theta*currenthouseholdstate(:,2))-delta*currenthouseholdstate(:,2))
!write(*,*) "maxq", maxval(currenthouseholdstate(:,1)), maxval(currenthouseholdstate(:,2))
    stateindicator1=0
    stateindicator2=0
    stateindicator5=0
    

    ! 5th state is age (or equiv) time
    currenthouseholdstate(:,5)=t
    currenthouseholdstatewithbubble(:,5)=t
    do i=1,numhouseholds
        if (alive(i,t)==1) then
        
        end if

        ! first thing is that if you are old, do you die or not
        if (t>Tretire) then
            if (alive(i,1)==1) then
                call random_number(shock)
                if (shock<deathprob(t-Tretire)) then
                    alive(i,t:Tdie)=0
                end if
            end if
        end if


   

      !  actualwealth(i,t)=currenthouseholdstate(i,1)+theta*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))  !total wealth is voluntary equity plus equity in durable
      !  financialwealth(i,t)=currenthouseholdstate(i,1)-(1-theta)*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))
      !  housingnet(i,t)=theta*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))
      !  housinggross(i,t)=currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)))
        !qchg(i,t) = financialwealth(i,t) + (1-theta)*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)-1))
        

        

        if (t>actualageatbubblestart .and. t<=actualageatbubblestart+totalbubblelength) then
        !write(*,*) bubbletimeholder, t, ageatbubblestart
        currenthouseholdstatewithbubble(i,1) = currenthouseholdstatewithbubble(i,1) + (1-theta)*currenthouseholdstatewithbubble(i,2)*(Ps(bubbletimeholder)-Ps(bubbletimeholder-1))
        end if
   

        call pol_linworking(currenthouseholdstatewithbubble(i,1)-.01,currenthouseholdstatewithbubble(i,2),currenthouseholdstatewithbubble(i,3),currenthouseholdstatewithbubble(i,4),currenthouseholdstatewithbubble(i,5),newhouseholdstatewithbubble(i,1),newhouseholdstatewithbubble(i,2),consumptionwithbubblempc(i,t),rentalindwithbubble(i,t),welfarewithbubble(i,t))
        call pol_linworking(currenthouseholdstatewithbubble(i,1),currenthouseholdstatewithbubble(i,2),currenthouseholdstatewithbubble(i,3),currenthouseholdstatewithbubble(i,4),currenthouseholdstatewithbubble(i,5),newhouseholdstatewithbubble(i,1),newhouseholdstatewithbubble(i,2),consumptionwithbubble(i,t),rentalindwithbubble(i,t),welfarewithbubble(i,t))
        mpc(i,t)=(consumptionwithbubblempc(i,t)-consumptionwithbubble(i,t))/-.01
        aggregatecontrib(i,t)=mpc(i,t)*currenthouseholdstatewithbubble(i,2)*(1-delta)*(1-F)
        call pol_linworking(currenthouseholdstate(i,1),currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumption(i,t),rentalind(i,t),welfare(i,t))

        if (newhouseholdstate(i,1)>amax) then
            newhouseholdstate(i,1)=.9999999*amax
        end if
         if (newhouseholdstatewithbubble(i,1)>amax) then
            newhouseholdstatewithbubble(i,1)=.9999999*amax
        end if

        hstockwithrental(i,1)=newhouseholdstate(i,2)
        hstockwithrentalwithbubble(i,1)=newhouseholdstatewithbubble(i,2)

        if (rentalind(i,t)>.995) then
            newhouseholdstate(i,2)=0
        end if
        if (rentalindwithbubble(i,t)>.995) then
            newhouseholdstatewithbubble(i,2)=0
        end if
            
        durableinvestment(i,t)=newhouseholdstate(i,2)-(1-delta)*currenthouseholdstate(i,2)
        durableinvestmentwithbubble(i,t)=newhouseholdstatewithbubble(i,2)-(1-delta)*currenthouseholdstatewithbubble(i,2)
        
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
        newhouseholdstatewithbubble(i,3)=newhouseholdstate(i,3) 
    end do


    !hlevelnobubble(t,ageatbubblestart)=sum(newhouseholdstate(:,2))
    !hlevelwithbubble(t,ageatbubblestart)=sum(newhouseholdstatewithbubble(:,2))
    hlevelnobubble(t,ageatbubblestart)=sum(hstockwithrental(:,1))
    hlevelwithbubble(t,ageatbubblestart)=sum(hstockwithrentalwithbubble(:,1))
    clevelnobubble(t,ageatbubblestart)=sum(consumption(:,t))/numhouseholds
    clevelwithbubble(t,ageatbubblestart)=sum(consumptionwithbubble(:,t))/numhouseholds
    rinobubble(t,ageatbubblestart)=sum(newhouseholdstate(:,2))-(1-delta)*sum(currenthouseholdstate(:,2))
    sufficientstatwithbubble(t,ageatbubblestart)=sum(aggregatecontrib(:,t))/numhouseholds
    

    !irfh(t,1)=sum(hlevelwithbubble(t,:))/sum(hlevelnobubble(t,:))-1
    !irfc(t,1)=sum(clevelwithbubble(t,:))/sum(clevelnobubble(t,:))-1
    !irfri(t,1)=sum(ribubble(t,:))/sum(rinobubble(t,:))-1

    !irfh(t,ageatbubblestart)=sum(newhouseholdstatewithbubble(:,2))/sum(newhouseholdstate(:,2))-1
    !irfc(t,ageatbubblestart)=sum(consumptionwithbubble(:,t))/sum(consumption(:,t))-1

    !write(*,*) t
    !write(*,*) irfh(t,1),irfri(t,1),irfc(t,1)
  !  write(*,*) "accuracy"
  !  write(*,*) t,ageatbubblestart,sufficientstatwithbubble(t,ageatbubblestart)

    !write(5,5) t*1.0,sum(currenthouseholdstate(:,1))/numhouseholds, sum(currenthouseholdstate(:,2))/numhouseholds,sum((1+r)*(currenthouseholdstate(:,1)+theta*currenthouseholdstate(:,2))-delta*currenthouseholdstate(:,2))/numhouseholds, sum(newhouseholdstate(:,2))/numhouseholds, sum(consumption(:,t))/numhouseholds, sum(rentalind(:,t))/numhouseholds,sum(currenthouseholdstate(:,1)-(1-theta)*currenthouseholdstate(:,2))/numhouseholds,sum(currenthouseholdstate(:,2)/(currenthouseholdstate(:,2)+currenthouseholdstate(:,1)-(1-theta)*currenthouseholdstate(:,2)))/numhouseholds

    



    !write(5,5) sum(financialwealth(:,t)*alive(:,t))/numhouseholds,sum(currenthouseholdstate(:,2)*alive(:,t))/numhouseholds, sum(consumption(:,t)*alive(:,t))/numhouseholds, sum(actualwealth(:,t)*alive(:,t))/numhouseholds, sum(newhouseholdstate(:,2)*alive(:,t))/numhouseholds, sum(newhouseholdstate(:,2)*rentalind(:,t)*alive(:,t))/numhouseholds, sum(rentalind(:,t)*alive(:,t))/numhouseholds, sum(alive(:,t))/numhouseholds, sum(consumption(:,t))/numhouseholds
   !  currenthouseholdstate=newhouseholdstate

    
    ! this is the housing price elasticity we want
   ! write(*,*) t*1.0, (sum(diffc(:,t))/sum(alive(:,t)))/(sum(consumption(:,t))/sum(alive(:,t)))/(hpmin)
    !write(*,*) "elasticity", sum(diffc(:,t)*alive(:,t))/sum(consumption(:,t))/sum(alive(:,t))/0.2


   ! write(*,*) t*1.0,sum(currenthouseholdstate(:,2)), sum(alive(:,t))
   
    !write(10,10) (sum((ratioc(:,t)-1))/sum(alive(:,t)))/(exp(hpmin) - 1), (sum(diffc(:,t))/sum(alive(:,t)))/(sum(consumption(:,t))/sum(alive(:,t)))/(hpmin)
   ! write(10,10) (sum((ratioc(:,t)-1))/numhouseholds)/(exp(hpmin) - 1), (sum(diffc(:,t))/numhouseholds)/(sum(consumption(:,t))/numhouseholds)/(exp(hpmin) - 1), ratiocown/numown_t,ratiocrent/numrent_t,numown_t/(numown_t+numrent_t),numsell_t/numhouseholds,numbuy_t/numhouseholds,(sum(diffc(:,t)/hpmin)/numhouseholds),sum(aggregatecontrib(:,t))/numhouseholds,sum(mpc(:,t))/numhouseholds,sum(currenthouseholdstate(:,2)*(1-delta)*(1-F))/numhouseholds

   ! write(555,555) a_t/numhouseholds, aown_t/numown_t, arent_t/numrent_t, leverage_t/numhouseholds, leverageown_t/numown_t, leveragerent_t/numrent_t
    

    currenthouseholdstate=newhouseholdstate
    currenthouseholdstatewithbubble=newhouseholdstatewithbubble
end do




end subroutine simulate


















  

subroutine pol_linworking(astate,Dstate,zstate,bubbletimestate,t,achoicelin,Dchoicelin,cchoicelin,rentallin,welfare)
    USE nrtype; USE nrutil
    USE nr
    USE share
    USE OMP_LIB
    IMPLICIT NONE
    REAL(8) :: weightal, weightDl
    REAL(8) :: nearestanode, nearestDnode
    REAL(8) :: al,ah, Dl,Dh
    REAL(8) :: astate,Dstate,zstate,bubbletimestate,t,achoicelin,Dchoicelin,cchoicelin,rentallin,welfare
    
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
            if (Dh>Dgridsize) then
            write(*,*) "error", Dh, Dl, Dstate,nearestDnode,Dnodes(nearestDnode)
            end if
            weightDl=1-(Dstate-Dnodes(Dl))/(Dnodes(Dh)-Dnodes(Dl))
        else
            Dl=nearestDnode-1
            Dh=nearestDnode
            if (Dl<1) then
            write(*,*) "error", Dh, Dl, Dstate,nearestDnode,Dnodes(nearestDnode), astate, bubbletimestate,t
            end if
            weightDl=1-(Dstate-Dnodes(Dl))/(Dnodes(Dh)-Dnodes(Dl))
        end if
    end if
    
    !write(*,*) weightal, weightDl,al,Dl,estate,Dchoice(al,Dl,estate)
    
    !write(*,*) "r  ", "r  ", al,ah,Dl,Dh,zstate,epsstate,t

    achoicelin=weightal*weightDl*achoice(al,Dl,zstate,bubbletimestate,t)+(1-weightal)*weightDl*achoice(ah,Dl,zstate,bubbletimestate,t)+weightal*(1-weightDl)*achoice(al,Dh,zstate,bubbletimestate,t)+(1-weightal)*(1-weightDl)*achoice(ah,Dh,zstate,bubbletimestate,t)
    Dchoicelin=weightal*weightDl*Dchoice(al,Dl,zstate,bubbletimestate,t)+(1-weightal)*weightDl*Dchoice(ah,Dl,zstate,bubbletimestate,t)+weightal*(1-weightDl)*Dchoice(al,Dh,zstate,bubbletimestate,t)+(1-weightal)*(1-weightDl)*Dchoice(ah,Dh,zstate,bubbletimestate,t)
    cchoicelin=weightal*weightDl*cchoice(al,Dl,zstate,bubbletimestate,t)+(1-weightal)*weightDl*cchoice(ah,Dl,zstate,bubbletimestate,t)+weightal*(1-weightDl)*cchoice(al,Dh,zstate,bubbletimestate,t)+(1-weightal)*(1-weightDl)*cchoice(ah,Dh,zstate,bubbletimestate,t)
    rentallin=weightal*weightDl*rentalindicator(al,Dl,zstate,bubbletimestate,t)+(1-weightal)*weightDl*rentalindicator(ah,Dl,zstate,bubbletimestate,t)+weightal*(1-weightDl)*rentalindicator(al,Dh,zstate,bubbletimestate,t)+(1-weightal)*(1-weightDl)*rentalindicator(ah,Dh,zstate,bubbletimestate,t)
    welfare=weightal*weightDl*EV(al,Dl,zstate,bubbletimestate,t)+(1-weightal)*weightDl*EV(ah,Dl,zstate,bubbletimestate,t)+weightal*(1-weightDl)*EV(al,Dh,zstate,bubbletimestate,t)+(1-weightal)*(1-weightDl)*EV(ah,Dh,zstate,bubbletimestate,t)
    
end subroutine pol_linworking





subroutine computempc(astate,Dstate,zstate,bubbletimestate,t,mpcholder)
    USE nrtype; USE nrutil
    USE nr
    USE share
    USE OMP_LIB
    IMPLICIT NONE
    REAL(8) :: weightal, weightDl
    REAL(8) :: nearestanode, nearestDnode
    REAL(8) :: al,ah, Dl,Dh
    REAL(8), intent(in) :: astate,Dstate,zstate,bubbletimestate,t
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
 
    cchoicelin1=weightal*weightDl*cchoice(al,Dl,zstate,bubbletimestate,t)+(1-weightal)*weightDl*cchoice(ah,Dl,zstate,bubbletimestate,t)+weightal*(1-weightDl)*cchoice(al,Dh,zstate,bubbletimestate,t)+(1-weightal)*(1-weightDl)*cchoice(ah,Dh,zstate,bubbletimestate,t)
    
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


    cchoicelin2=weightal*weightDl*cchoice(al,Dl,zstate,bubbletimestate,t)+(1-weightal)*weightDl*cchoice(ah,Dl,zstate,bubbletimestate,t)+weightal*(1-weightDl)*cchoice(al,Dh,zstate,bubbletimestate,t)+(1-weightal)*(1-weightDl)*cchoice(ah,Dh,zstate,bubbletimestate,t)
    mpcholder=(cchoicelin2-cchoicelin1)/.01

end subroutine computempc







subroutine solveworkingproblem
    USE nrtype; USE nrutil
    USE nr
    USE share
    USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: i,j,k,m,t
    REAL(8) :: timestart2, timeend2
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,2) :: optpolicynoadjust,optpolicyadjust
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,3,2) :: pstartadjust
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,3) :: ystartadjust
    REAL(8), dimension(agridsize,Dgridsize,zgridsize) :: ax,bx,cx, adjust
    INTEGER, dimension(agridsize,Dgridsize,zgridsize) :: iter
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,4) :: state

  
    
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
    iter=0
    
    adjust=0
    
    do bubbletime=1,totalbubblelength
    write(*,*) "Price trajectory version:", bubbletime
    P1 = Ps(bubbletime);
    g = gs(bubbletime);
    do t=1,Tdie+1
    P_exp(t)=P1*exp(g*(1-rho**(t-1))/(1-rho))
   ! write(*,*) "P_exp", t, P_exp(t)
    end do

    

    do age=1,Tdie
    
    !write(*,*) bubbletime, age
    P_exp(Tdie+2-age)=P_exp(Tdie-age+1)
    
    do t=Tdie,age,-1
    

        
        !$OMP PARALLEL
        !$OMP DO
        do i=1,agridsize
            do j=1,Dgridsize
                do k=1,zgridsize
                        pstartadjust(i,j,k,1,1)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,1,2)=.8*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,2,1)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,2,2)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,3,1)=.5*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,3,2)=.2*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))

                        
                    
                        state(i,j,k,1)=anodes(i)
                        state(i,j,k,2)=Dnodes(j)
                        state(i,j,k,3)=k
                        state(i,j,k,4)=t
                    
                        ystartadjust(i,j,k,1)=valfuncadjust2(pstartadjust(i,j,k,1,:),state(i,j,k,:))
                        ystartadjust(i,j,k,2)=valfuncadjust2(pstartadjust(i,j,k,2,:),state(i,j,k,:))
                        ystartadjust(i,j,k,3)=valfuncadjust2(pstartadjust(i,j,k,3,:),state(i,j,k,:))

                    
                        call amoeba(state(i,j,k,:),pstartadjust(i,j,k,:,:),ystartadjust(i,j,k,:),ftol,valfuncadjust)
                        Vadjust(i,j,k,bubbletime,t)=ystartadjust(i,j,k,1)
                        achoiceadjust(i,j,k,bubbletime,t)=pstartadjust(i,j,k,1,1)
                        Dchoiceadjust(i,j,k,bubbletime,t)=pstartadjust(i,j,k,1,2)
                    
                    
                        pstartadjust(i,j,k,1,1)=0.0
                        pstartadjust(i,j,k,1,2)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,2,1)=.05*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,2,2)=.1*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,3,1)=.4*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,3,2)=.21*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        
                        ystartadjust(i,j,k,1)=valfuncadjust2(pstartadjust(i,j,k,1,:),state(i,j,k,:))
                        ystartadjust(i,j,k,2)=valfuncadjust2(pstartadjust(i,j,k,2,:),state(i,j,k,:))
                        ystartadjust(i,j,k,3)=valfuncadjust2(pstartadjust(i,j,k,3,:),state(i,j,k,:))

                        call amoeba(state(i,j,k,:),pstartadjust(i,j,k,:,:),ystartadjust(i,j,k,:),ftol,valfuncadjust)
                        if (ystartadjust(i,j,k,1)<Vadjust(i,j,k,bubbletime,t)) then  !(again we're minimizing)
                            Vadjust(i,j,k,bubbletime,t)=ystartadjust(i,j,k,1)
                            achoiceadjust(i,j,k,bubbletime,t)=pstartadjust(i,j,k,1,1)
                            Dchoiceadjust(i,j,k,bubbletime,t)=pstartadjust(i,j,k,1,2)
                        end if


                        







                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        pstartadjust(i,j,k,1,1)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,1,2)=.8*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,2,1)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,2,2)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,3,1)=.5*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,3,2)=.2*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                    
                        state(i,j,k,1)=anodes(i)
                        state(i,j,k,2)=Dnodes(j)
                        state(i,j,k,3)=k
                        state(i,j,k,4)=t
                    
                        ystartadjust(i,j,k,1)=valfuncrent2(pstartadjust(i,j,k,1,:),state(i,j,k,:))
                        ystartadjust(i,j,k,2)=valfuncrent2(pstartadjust(i,j,k,2,:),state(i,j,k,:))
                        ystartadjust(i,j,k,3)=valfuncrent2(pstartadjust(i,j,k,3,:),state(i,j,k,:))
                    
                        call amoeba(state(i,j,k,:),pstartadjust(i,j,k,:,:),ystartadjust(i,j,k,:),ftol,valfuncrent)
                    
                    
                        Vrent(i,j,k,bubbletime,t)=ystartadjust(i,j,k,1)
                        achoicerent(i,j,k,bubbletime,t)=pstartadjust(i,j,k,1,1)
                        Dchoicerent(i,j,k,bubbletime,t)=pstartadjust(i,j,k,1,2)
                    
                    
                        pstartadjust(i,j,k,1,1)=0.0
                        pstartadjust(i,j,k,1,2)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,2,1)=.05*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,2,2)=.1*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,3,1)=.4*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,3,2)=.21*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
 
                    
                        ystartadjust(i,j,k,1)=valfuncrent2(pstartadjust(i,j,k,1,:),state(i,j,k,:))
                        ystartadjust(i,j,k,2)=valfuncrent2(pstartadjust(i,j,k,2,:),state(i,j,k,:))
                        ystartadjust(i,j,k,3)=valfuncrent2(pstartadjust(i,j,k,3,:),state(i,j,k,:))

                        call amoeba(state(i,j,k,:),pstartadjust(i,j,k,:,:),ystartadjust(i,j,k,:),ftol,valfuncrent)
                    
                    
                        if (ystartadjust(i,j,k,1)<Vrent(i,j,k,bubbletime,t)) then  !(again we're minimizing)
                            Vrent(i,j,k,bubbletime,t)=ystartadjust(i,j,k,1)
                            achoicerent(i,j,k,bubbletime,t)=pstartadjust(i,j,k,1,1)
                            Dchoicerent(i,j,k,bubbletime,t)=pstartadjust(i,j,k,1,2)
                        end if
                        

                        if (costlyequity==1) then

                            if (anodes(i)>=(1-theta)*Dnodes(j)*P_exp(t-age+1)) then
                                ax(i,j,k)=(1-theta)*Dnodes(j)*P_exp(t-age+1)
                            else
                                ax(i,j,k)=anodes(i)
                            end if
                        bx(i,j,k)=ax(i,j,k)
                        else
                        ax(i,j,k)=0.0
                        bx(i,j,k)=borrowconstraint
                        end if
                        cx(i,j,k)=amax
           
                        state(i,j,k,2)=j
                        Vnoadjust(i,j,k,bubbletime,t)=brentnew(ax(i,j,k),bx(i,j,k),cx(i,j,k),valfuncnoadjust,ftol,state(i,j,k,1),state(i,j,k,2),state(i,j,k,3),state(i,j,k,4),achoicenoadjust(i,j,k,bubbletime,t))
                        !if (bubbletime==2) then
                        !    write(*,*), ax(i,j,k), cx(i,j,k), Vnoadjust(i,j,k,bubbletime,t), achoicenoadjust(i,j,k,bubbletime,t)
                        !end if
                        Dchoicenoadjust(i,j,k,bubbletime,t)=Dnodes(j)

                        if (Vadjust(i,j,k,bubbletime,t)<Vnoadjust(i,j,k,bubbletime,t) .and. Vadjust(i,j,k,bubbletime,t)<Vrent(i,j,k,bubbletime,t)) then  ! since V = - V from minimization
                            achoice(i,j,k,bubbletime,t)=achoiceadjust(i,j,k,bubbletime,t)
                            Dchoice(i,j,k,bubbletime,t)=Dchoiceadjust(i,j,k,bubbletime,t)
                            if (anodes(i)-(1-theta)*Dnodes(j)*P_exp(t-age+1)<0) then
                                cchoice(i,j,k,bubbletime,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*P_exp(t-age+1)-theta*Dchoice(i,j,k,bubbletime,t)*P_exp(t-age+1)-achoice(i,j,k,bubbletime,t)-(1+rborrow)*(1-theta)*(1-F)*Dnodes(j)*P_exp(t-age+1)
                                !write(*,*) "diff1", income(k,t)+(1+rborrow*exp(hpnodes(l))**)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-theta*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**)*(1-theta)*(1-F)*Dnodes(j)*exp(hpnodes(l))-(income(k,t)+(1+rborrow*exp(hpnodes(l))**(-0.0))*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-theta*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**(-0.0))*(1-theta)*(1-F)*Dnodes(j)*exp(hpnodes(l)))
                            else
                                cchoice(i,j,k,bubbletime,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*P_exp(t-age+1)-theta*Dchoice(i,j,k,bubbletime,t)*P_exp(t-age+1)-achoice(i,j,k,bubbletime,t)-(1+r)*(1-theta)*Dnodes(j)*(1-F)*P_exp(t-age+1)
                            end if
                            rentalindicator(i,j,k,bubbletime,t)=0
                        elseif (Vnoadjust(i,j,k,bubbletime,t)<Vadjust(i,j,k,bubbletime,t) .and. Vnoadjust(i,j,k,bubbletime,t)<Vrent(i,j,k,bubbletime,t)) then
                            achoice(i,j,k,bubbletime,t)=achoicenoadjust(i,j,k,bubbletime,t)
                            Dchoice(i,j,k,bubbletime,t)=Dchoicenoadjust(i,j,k,bubbletime,t)


                            if (anodes(i)-(1-theta)*Dnodes(j)*P_exp(t-age+1)<0) then
                                cchoice(i,j,k,bubbletime,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*P_exp(t-age+1)-theta*Dchoice(i,j,k,bubbletime,t)*P_exp(t-age+1)-achoice(i,j,k,bubbletime,t)-(1+rborrow)*(1-theta)*Dnodes(j)*P_exp(t-age+1)
                                !write(*,*) "diff2", income(k,t)+(1+rborrow*exp(hpnodes(l))**)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-theta*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**)*(1-theta)*Dnodes(j)*exp(hpnodes(l))-(income(k,t)+(1+rborrow*exp(hpnodes(l))**(-0.0))*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-theta*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**(-0.0))*(1-theta)*Dnodes(j)*exp(hpnodes(l)))
                            else
                                cchoice(i,j,k,bubbletime,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*P_exp(t-age+1)-theta*Dchoice(i,j,k,bubbletime,t)*P_exp(t-age+1)-achoice(i,j,k,bubbletime,t)-(1+r)*(1-theta)*Dnodes(j)*P_exp(t-age+1)
                            end if
                            rentalindicator(i,j,k,bubbletime,t)=0
                        else
                            achoice(i,j,k,bubbletime,t)=achoicerent(i,j,k,bubbletime,t)
                            Dchoice(i,j,k,bubbletime,t)=Dchoicerent(i,j,k,bubbletime,t)
                            if (anodes(i)-(1-theta)*Dnodes(j)*P_exp(t-age+1)<0) then
                                cchoice(i,j,k,bubbletime,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*P_exp(t-age+1)-r_rental_t(t)*Dchoice(i,j,k,bubbletime,t)*P_exp(t-age+1)-achoice(i,j,k,bubbletime,t)-(1+rborrow)*(1-theta)*(1-F)*Dnodes(j)*P_exp(t-age+1)
                                !write(*,*) "diff3", income(k,t)+(1+rborrow*exp(hpnodes(l))**)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental_t(t)*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**)*(1-theta)*(1-F)*Dnodes(j)*exp(hpnodes(l)) - (income(k,t)+(1+rborrow*exp(hpnodes(l))**(-0.0))*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental_t(t)*Dchoice(i,j,k,l,t)*exp(hpnodes(l))-achoice(i,j,k,l,t)-(1+rborrow*exp(hpnodes(l))**(-0.0))*(1-theta)*(1-F)*Dnodes(j)*exp(hpnodes(l)))
                            else
                                cchoice(i,j,k,bubbletime,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*P_exp(t-age+1)-r_rental_t(t)*Dchoice(i,j,k,bubbletime,t)*P_exp(t-age+1)-achoice(i,j,k,bubbletime,t)-(1+r)*(1-theta)*(1-F)*Dnodes(j)*P_exp(t-age+1)
                            end if
                            !write(*,*) l,(1+rborrow*exp(hpnodes(l))**)
                            rentalindicator(i,j,k,bubbletime,t)=1
                        end if
                       ! write(*,*) exp(hpnodes(l)), (1+rborrow*exp(hpnodes(l))**), income(6,l,1)
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        Vnoadjust(:,:,:,bubbletime,t)=-Vnoadjust(:,:,:,bubbletime,t)
        Vadjust(:,:,:,bubbletime,t)=-Vadjust(:,:,:,bubbletime,t)
        Vrent(:,:,:,bubbletime,t)=-Vrent(:,:,:,bubbletime,t)

        EV(:,:,:,bubbletime,t)=max(Vrent(:,:,:,bubbletime,t),max(Vnoadjust(:,:,:,bubbletime,t),Vadjust(:,:,:,bubbletime,t)))

       ! write(*,*) EV(1,1,2,1,2), EV(2,1,2,1,2)

        
end do
end do
end do
     
end subroutine solveworkingproblem










REAL(8) FUNCTION valfuncadjust(p,state)  !p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
USE SHARE
IMPLICIT NONE
REAL(8), DIMENSION(2) :: p
REAL(8), DIMENSION(4) :: state
REAL(8) :: a, D, z, hp, t, zindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
t=state(4)


currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
aprime=aprime+(1-theta)*Dcurrent*(P_exp(t-age+1+1)-P_exp(t-age+1))

if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dnodes(Dgridsize) .AND. Dcurrent>=0) then
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
    if ((a-D*(1-theta)*P_exp(t-age+1))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)*P_exp(t-age+1)-theta*Dcurrent*P_exp(t-age+1)-aprimeholder-(1+rborrow)*(1-theta)*D*(1-F)*P_exp(t-age+1))
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)*P_exp(t-age+1)-theta*Dcurrent*P_exp(t-age+1)-aprimeholder-(1+r)*(1-theta)*D*(1-F)*P_exp(t-age+1))
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
            do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,i,bubbletime,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,i,bubbletime,t+1)
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,i,bubbletime,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,i,bubbletime,t+1)
            end do
            valfuncadjust=valfuncadjust+beta2*EVholder
        else
                    EVholder=EVholder+weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,bubbletime,t+1)
                    EVholder=EVholder+(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,bubbletime,t+1)
                    EVholder=EVholder+weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,bubbletime,t+1)
                    EVholder=EVholder+(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,bubbletime,t+1)

                    EVholder2=EVholder2+weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,bubbletime,Tdie+1)
                    EVholder2=EVholder2+(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,bubbletime,Tdie+1)
                    EVholder2=EVholder2+weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,bubbletime,Tdie+1)
                    EVholder2=EVholder2+(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,bubbletime,Tdie+1)
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
REAL(8), DIMENSION(4) :: state
REAL(8) :: a, D, z, hp, t, zindex,hpindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
t=state(4)



currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
aprime=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dnodes(Dgridsize) .AND. Dcurrent>=0) then
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
    if ((a-D*(1-theta)*P_exp(t-age+1))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)*P_exp(t-age+1)-r_rental_t(t)*Dcurrent*P_exp(t-age+1)-aprimeholder-(1+rborrow)*(1-theta)*D*(1-F)*P_exp(t-age+1))
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)*P_exp(t-age+1)-r_rental_t(t)*Dcurrent*P_exp(t-age+1)-aprimeholder-(1+r)*(1-theta)*D*(1-F)*P_exp(t-age+1))
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
            do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*EV(aprimel,1,i,bubbletime,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*EV(aprimeh,1,i,bubbletime,t+1)
            end do
            valfuncrent=valfuncrent+beta2*EVholder
        else
                    EVholder=EVholder+weightaprimel*EV(aprimel,1,zindex,bubbletime,t+1)
                    EVholder=EVholder+(1-weightaprimel)*EV(aprimeh,1,zindex,bubbletime,t+1)

                    EVholder2=EVholder2+weightaprimel*EV(aprimel,1,zindex,bubbletime,Tdie+1)
                    EVholder2=EVholder2+(1-weightaprimel)*EV(aprimeh,1,zindex,bubbletime,Tdie+1)
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
REAL(8), DIMENSION(4) :: state
REAL(8) :: a, D, z, hp, t, zindex,hpindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
t=state(4)



currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
aprime=aprime+(1-theta)*Dcurrent*(P_exp(t-age+1+1)-P_exp(t-age+1))
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dnodes(Dgridsize) .AND. Dcurrent>=0) then
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

            if (aprimel==0) then
            write(*,*) "aprimel==0", nearestanode, aprime, anodes(nearestanode)
            end if
            weightaprimel=1-(aprime-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
        end if
    end if
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
    if ((a-D*(1-theta)*P_exp(t-age+1))<0) then
    consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)*P_exp(t-age+1)-theta*Dcurrent*P_exp(t-age+1)-aprimeholder-(1+rborrow)*(1-theta)*D*(1-F)*P_exp(t-age+1))
    else
    consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)*P_exp(t-age+1)-theta*Dcurrent*P_exp(t-age+1)-aprimeholder-(1+r)*(1-theta)*D*(1-F)*P_exp(t-age+1))
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
            do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,i,bubbletime,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,i,bubbletime,t+1)
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,i,bubbletime,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,i,bubbletime,t+1)
            end do
            valfuncadjust2=valfuncadjust2+beta2*EVholder
        else
                    EVholder=EVholder+weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,bubbletime,t+1)
                    EVholder=EVholder+(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,bubbletime,t+1)
                    EVholder=EVholder+weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,bubbletime,t+1)
                    EVholder=EVholder+(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,bubbletime,t+1)

                    EVholder2=EVholder2+weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,bubbletime,Tdie+1)
                    EVholder2=EVholder2+(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,bubbletime,Tdie+1)
                    EVholder2=EVholder2+weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,bubbletime,Tdie+1)
                    EVholder2=EVholder2+(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,bubbletime,Tdie+1)
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
REAL(8), DIMENSION(4) :: state
REAL(8) :: a, D, z, hp, t, zindex,hpindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder

!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
t=state(4)



currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
aprime=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dnodes(Dgridsize) .AND. Dcurrent>=0) then
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
    if ((a-D*(1-theta)*P_exp(t-age+1))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)*P_exp(t-age+1)-r_rental_t(t)*Dcurrent*P_exp(t-age+1)-aprimeholder-(1+rborrow)*(1-theta)*D*(1-F)*P_exp(t-age+1))
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)*P_exp(t-age+1)-r_rental_t(t)*Dcurrent*P_exp(t-age+1)-aprimeholder-(1+r)*(1-theta)*D*(1-F)*P_exp(t-age+1))
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
            do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*EV(aprimel,1,i,bubbletime,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*EV(aprimeh,1,i,bubbletime,t+1)
            end do
            valfuncrent2=valfuncrent2+beta2*EVholder
        else
                    EVholder=EVholder+weightaprimel*EV(aprimel,1,zindex,bubbletime,t+1)
                    EVholder=EVholder+(1-weightaprimel)*EV(aprimeh,1,zindex,bubbletime,t+1)

                    EVholder2=EVholder2+weightaprimel*EV(aprimel,1,zindex,bubbletime,Tdie+1)
                    EVholder2=EVholder2+(1-weightaprimel)*EV(aprimeh,1,zindex,bubbletime,Tdie+1)
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
REAL(8) :: pholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)



currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
aprime=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if



if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dnodes(Dgridsize) .AND. Dcurrent>=0) then
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
    if ((a-D*(1-theta)*P_exp(t-age+1))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)*P_exp(t-age+1)-r_rental_t(t)*Dcurrent*P_exp(t-age+1)-aprimeholder-(1+rborrow)*(1-theta)*D*(1-F)*P_exp(t-age+1))
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)*P_exp(t-age+1)-r_rental_t(t)*Dcurrent*P_exp(t-age+1)-aprimeholder-(1+r)*(1-theta)*D*(1-F)*P_exp(t-age+1))
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
            do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*EV(aprimel,1,i,bubbletime,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*EV(aprimeh,1,i,bubbletime,t+1)
            end do
            valfuncrent3=valfuncrent3+beta2*EVholder
        else
                    EVholder=EVholder+weightaprimel*EV(aprimel,1,zindex,bubbletime,t+1)
                    EVholder=EVholder+(1-weightaprimel)*EV(aprimeh,1,zindex,bubbletime,t+1)

                    EVholder2=EVholder2+weightaprimel*EV(aprimel,1,zindex,bubbletime,Tdie+1)
                    EVholder2=EVholder2+(1-weightaprimel)*EV(aprimeh,1,zindex,bubbletime,Tdie+1)
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





REAL(8) FUNCTION valfuncnoadjust(aprime,a,Dindex,zindex,t)  !p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
USE SHARE
IMPLICIT NONE
REAL(8), DIMENSION(2) :: p
REAL(8) :: a, z, hp, t, zindex,Dindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel, EVholder, EVholder2
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: Vadjustholder,Vnoadjustholder
REAL(8) :: Vadjustretireholder,Vnoadjustretireholder
REAL(8) :: aprime_adj

!write(*,*) state

!write(*,*) state

currentincome=income(zindex,t)


Dcurrent=Dnodes(Dindex)
Dprime=Dnodes(Dindex)
aprimeholder=aprime

if (switch==2) then
write(*,*) "other stuff", aprime, Dcurrent, (P_exp(t-age+1+1)-P_exp(t-age+1))
end if

aprime_adj=max(aprime+ &
    (1-theta)*Dcurrent*(P_exp(t-age+1+1)-P_exp(t-age+1)), borrowconstraint)
if (switch==2) then
write(*,*) "other stuff 2", aprime
end if
if (aprime_adj>anodes(agridsize)) then
aprime_adj=anodes(agridsize)
end if



if (aprime_adj>=borrowconstraint .AND. aprime_adj<=anodes(agridsize)) then
    nearestanode=minloc(abs(aprime-anodes),1)
    if (aprime_adj==anodes(nearestanode)) then
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
        if (aprime_adj-anodes(nearestanode)>0) then
            aprimel=nearestanode
            aprimeh=nearestanode+1
            weightaprimel=1-(aprime_adj-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
        else
            aprimel=nearestanode-1
            aprimeh=nearestanode
            weightaprimel=1-(aprime_adj-anodes(aprimel))/(anodes(aprimeh)-anodes(aprimel))
        end if
    end if
    
    !using that Dprime=D
    if ((a-Dcurrent*(1-theta)*P_exp(t-age+1))<0) then
        consumption=(currentincome+(1+rborrow)*a+Dprime*(1-delta)*P_exp(t-age+1)-theta*Dprime*P_exp(t-age+1)-aprimeholder-(1+rborrow)*(1-theta)*Dprime*P_exp(t-age+1))
    else
        consumption=(currentincome+(1+r)*a+Dprime*(1-delta)*P_exp(t-age+1)-theta*Dprime*P_exp(t-age+1)-aprimeholder-(1+r)*(1-theta)*Dprime*P_exp(t-age+1))
    end if

    if (switch==2) then
    write(*,*) "C stuff"
    write(*,*) consumption, currentincome,(1+r)*a,Dprime*(1-delta)*P_exp(t-age+1),theta*Dprime*P_exp(t-age+1),aprimeholder,(1+r)*(1-theta)*Dprime*P_exp(t-age+1)
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
            do i=minexpectationz(zindex),maxexpectationz(zindex)
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*EV(aprimel,Dindex,i,bubbletime,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*EV(aprimeh,Dindex,i,bubbletime,t+1)
            end do
        valfuncnoadjust=valfuncnoadjust+beta2*EVholder
        else
                    EVholder=EVholder+weightaprimel*EV(aprimel,Dindex,zindex,bubbletime,t+1)
                    EVholder=EVholder+(1-weightaprimel)*EV(aprimeh,Dindex,zindex,bubbletime,t+1)

                    EVholder2=EVholder2+weightaprimel*EV(aprimel,Dindex,zindex,bubbletime,Tdie+1)
                    EVholder2=EVholder2+(1-weightaprimel)*EV(aprimeh,Dindex,zindex,bubbletime,Tdie+1)
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
	INTEGER(I4B) :: iter
	REAL(SP), INTENT(IN) :: ftol
	REAL(SP), DIMENSION(3), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(3,2), INTENT(INOUT) :: p
    REAL(SP), DIMENSION(3) :: yholder
	REAL(SP), DIMENSION(3,2) :: pholder
	REAL(SP), DIMENSION(4), INTENT(IN) :: stateholder
	INTERFACE
		FUNCTION valfuncadjust(x,stateholder)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(2), INTENT(IN) :: x
		REAL(8), DIMENSION(4), INTENT(IN) :: stateholder
		REAL(SP) :: valfuncadjust
		END FUNCTION valfuncadjust
	END INTERFACE
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
	REAL(8), DIMENSION(4) :: stateholder
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
	




FUNCTION brentnew(ax,bx,cx,value,tol,aholder,Dholder,zholder,tholder,xmin)
USE nrtype; USE nrutil, ONLY : nrerror
USE share
IMPLICIT NONE
REAL(SP), INTENT(IN) :: ax,bx,cx,tol,aholder,Dholder,zholder,tholder
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
fx=value(x,aholder,Dholder,zholder,tholder)
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
    fu=value(u,aholder,Dholder,zholder,tholder)
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

write(*,*) "stuff" , aholder,Dholder,zholder,tholder,xmin,bubbletime,age
switch=2
fu=value(u,aholder,Dholder,zholder,tholder)
fx=value(x,aholder,Dholder,zholder,tholder)
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
