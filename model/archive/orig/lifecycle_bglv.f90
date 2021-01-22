
module share
    REAL(8), PARAMETER :: scalefactor=1.0
    REAL(8), parameter :: delta=.022 !depreciation rate of durable  was .03
    REAL(8), parameter :: deltak=.05 !depreciation rate of capital
    REAL(8), parameter :: r=.024
    REAL(8), parameter :: rborrow=r
    REAL(8), parameter :: elasticity= 2  !elasticity of substitution
    REAL(8), parameter :: thetamatlab=.2
    REAL(8), parameter :: theta=1-(1-thetamatlab)*(1-delta)/(1+r) !down payment
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

    
    REAL(8), parameter :: hpshockgridsize=5
    REAL(8), DIMENSION(hpshockgridsize) :: hpshock
    REAL(8), parameter :: hpshockmin=-.125
    REAL(8), parameter :: hpshockmax=0.125
    REAL(8), DIMENSION(hpshockgridsize) :: Probhpshock,Probhpshockcum
    
    

    !for now, just two point process
    REAL(8), parameter :: hpmin=0.0
    REAL(8), parameter :: hpmax=0.0
    integer, parameter :: hpgridsize=2 !size of house price grid
    REAL(8), DIMENSION(hpgridsize) :: hpnodes
    REAL(8), DIMENSION(hpgridsize,hpgridsize) :: Probhp,Probhpcum
    REAL(8), DIMENSION(hpgridsize) :: maxexpectationhp,minexpectationhp 
    REAL(8), parameter :: rho_hp=1.0
    REAL(8), parameter :: sigma_hp=.0459

    !idiosyncratic earnings shocks
    REAL(8), parameter :: sigma_z=.21
    REAL(8), parameter :: rho_z=.91
    REAL(8), parameter :: sigma_eta_init=(sigma_z**2.0/(1-rho_z**2.0))**.5
    REAL(8), parameter :: zmin=-2.5*(sigma_z**2.0/(1-rho_z**2.0))**.5
    REAL(8), parameter :: zmax=2.5*(sigma_z**2.0/(1-rho_z**2.0))**.5
integer, parameter :: zgridsize=13  !size of ar(1) income shock grid
 
    
    REAL(8), parameter :: Dmin=0    !minimum durable choice
    REAL(8), parameter :: Dmax=12   !maximum durable choice
    integer, parameter :: agridsize = 120 ! size of asset grid
    integer, parameter :: Dgridsize = 40 ! size of durable grid

    !REAL(8), parameter :: amin=(exp(hpmin)-exp(hpmax))*(1-theta)*Dmax    !minimum asset value (actually voluntary equity if theta!=1).  Note that if theta<1, we can have assets <0 even with amin=0 since a is vol. equity.  But with theta=1 & amin=0, no borrowing.
    !REAL(8), parameter :: amin=-.89
    REAL(8), parameter :: amin=-0.9/(1+r)  !.40 is roughly min y, so this is basically a no default condition
    !REAL(8), parameter :: amax2=8   !max asset
    REAL(8), parameter :: amax=30   !max asset

    
    REAL(8), parameter :: borrowconstraint=0.0
    

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
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: achoiceadjust,  achoicenoadjust,achoicerent,  Dchoiceadjust, Dchoicenoadjust,Dchoicerent, achoice, Dchoice,cchoice,cchoiceshock,Dchoiceshock,mpc_pol  !Policy functions when solving problem
    REAL(8) :: priceshock
   
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: rentalindicator
    REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: EV
   
    
    REAL(8), parameter :: numhouseholds=100000    !Size of household panel
    !integer, parameter :: burnin=1000    !Initial periods dropped in simulation
    
    
    REAL(8) :: diffmoments,diffmoments1,diffmoments2,diffmoments3
    REAL(8), parameter :: difftol=.01
    REAL(8) :: numiter

    REAL(8), parameter :: costlyequity=0

    REAL(8), parameter :: rentelasticity=1

    REAL(8), DIMENSION(11,6) :: moments

    REAL(8), parameter :: shocksize=2.0*sigma_hp
	
	REAL(8), parameter :: allowunderwater=0
    
    REAL(8) :: xfactor
    
    REAL(8), parameter :: mu=.012
    !REAL(8), parameter :: mu=0.0
    
    !REAL(8) :: amort_rate=.969
    REAL(8) :: amort_rate=1
	
	integer :: param_iter
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

OPEN (UNIT=901, FILE="momentsgridsearch.txt", ACTION="WRITE", POSITION="REWIND")
901 format (F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12, F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12,F40.12, F40.12, F40.12, F40.12, F40.12, F40.12,F40.12)


param_iter=0



ALLOCATE(Vnoadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie+1),  Vadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie+1),Vrent(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie+1))  !Conjectured value of not adjusting and adjusting, and new values of not adjusting and adjusting
ALLOCATE(mpc_pol(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),achoiceadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),  achoicenoadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),achoicerent(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),  Dchoiceadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie), Dchoicenoadjust(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),Dchoicerent(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie))
ALLOCATE(achoice(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie), Dchoice(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie),cchoice(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie))  !Policy functions when solving problem
ALLOCATE(rentalindicator(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie))
ALLOCATE(EV(agridsize,Dgridsize,zgridsize,hpgridsize,Tdie+1))

CALL OMP_SET_NUM_THREADS(28)


do param1=1,3
do param2=1,3
do param3=1,3
do param4=1,3
do param5=1,3
do param6=1,3

param_iter=param_iter+1

!param1=1
!param2=3
!param3=1
!param4=1
!param5=3
!param6=3

beta2=.9125+.005*param1
beta2retire=beta2

!elasticity2=.83+.01*param2
elasticity2=.871+.0055*param2
psi=2325+75*param3
r_rental=.045+.005*param4
r_rental_retire=.022+.004*param5
ret_wealth=1.0+.2*param6

write(*,*) "parameters", beta2, elasticity2, psi, r_rental, r_rental_retire, ret_wealth


const = (elasticity2**elasticity2)*((1.0-elasticity2)**(1.0-elasticity2))



r_rental_t(1:Tretire)=r_rental
r_rental_t(Tretire+1:Tdie)=r_rental_retire







numiter=0





    numiter=numiter+1
    write(*,*) "numiter", numiter

    write(*,*) elasticity2, r_rental, beta2



    !OPEN (UNIT=1, FILE="vfunc.txt", ACTION="WRITE", POSITION="REWIND")
    !1 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

    ! FILL IN ALL THE GRID POINTS
    
    DO i=1,hpshockgridsize
        hpshock(i)=hpshockmin+ ((hpshockmax -hpshockmin)/(hpshockgridsize-1))*(1.0*i-1.0)
        write(*,*) "hpshock", hpshock(i)
    END DO
    
    
    do i=1,zgridsize
        znodes(i)=zmin+((zmax-zmin)/(zgridsize-1))*(1.0*i-1.0)
        write(*,*) znodes(i)
    end do
 
    DO i=1,hpgridsize
        hpnodes(i)=hpmin+ ((hpmax -hpmin)/(hpgridsize-1))*(1.0*i-1.0)
    END DO

    
    
    
    w=hpshock(2)-hpshock(1)
    Probhpshock(1)=cdfnormal((hpshock(1)+w/2)/sigma_hp)
    Probhpshock(hpshockgridsize)=1-cdfnormal((hpshock(hpshockgridsize)-w/2)/sigma_hp)
    do j=2,hpshockgridsize-1
        Probhpshock(j)=cdfnormal((hpshock(j)+w/2)/sigma_hp)-cdfnormal((hpshock(j)-w/2)/sigma_hp)
    end do
    
    do j=1,hpshockgridsize
        Probhpshockcum(j)=sum(Probhpshock(1:j))
        write(*,*) "Probhpshock", Probhpshockcum(j)
    end do
    
    xfactor=0
    do j=1,hpshockgridsize
        xfactor=xfactor+Probhpshock(j)*((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))
    end do
    write(*,*) "xfactor", xfactor
    
    minexpectationhp=1.0
    maxexpectationhp=hpgridsize
    do j=1,hpgridsize
        foundmin=0.0
        foundmax=0.0
        do k=1,hpgridsize
            if (Probhp(j,k)>.01) then
                foundmin=1.0
            elseif (foundmin==0.0) then
                minexpectationhp(j)=minexpectationhp(j)+1
            end if
        end do
        do k=0,hpgridsize-1
            if (Probhp(j,hpgridsize-k)>.01) then
                foundmax=1.0
            elseif (foundmax==0.0) then
                maxexpectationhp(j)=maxexpectationhp(j)-1
            end if
        end do
    end do



    do j=1,hpgridsize
        do k=1,hpgridsize
            if (k<minexpectationhp(j)) then
                !write(*,*) Probz(j,k)
                Probhp(j,k)=0.0
            end if
            if (k>maxexpectationhp(j)) then
                !write(*,*) Probz(j,k)
                Probhp(j,k)=0.0
            end if
        end do
    end do

    do j=1,hpgridsize
        Probhp(j,:)=Probhp(j,:)/sum(Probhp(j,:))
    end do
    do j=1,hpgridsize
        do k=1,hpgridsize
            Probhpcum(j,k)=sum(Probhp(j,1:k))  
        end do
    end do


    DO i=1,agridsize
        anodes(i)=(1.0/(agridsize-1))*(1.0*i-1.0)
    end do
    

    do i=1,agridsize
        anodes(i)=exp(log(amax-amin+1)*anodes(i))+amin-1.0
    end do


    write(*,*) "anodes2",anodes



     DO i=1,Dgridsize
        dnodes(i)=(1.0/(dgridsize-1))*(1.0*i-1.0)
    end do
    

    do i=1,dgridsize
        dnodes(i)=exp(log(dmax-dmin+1)*dnodes(i))+dmin-1.0
    end do

    write(*,*) "dnodes2",dnodes
  
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
        write(*,*) "r i", retirementincome(j)
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

    do t=1,Tdie
        write(*,*) t, income(1,t), income(13,t)
    end do

    

    write(*,*) "income"
  !  write(*,*) sum(income(:,1,:)),sum(income(:,2,:)),sum(income(:,3,:)),sum(income(:,4,:)),sum(income(:,5,:))




! Put in bequest motive in last period

    call cpu_time(timestart)
    !call solveretirementproblem
    !Vnoadjust(:,:,:,:,61)=0
    !Vadjust(:,:,:,:,61)=0
   ! Vrent(:,:,:,:,61)=0

    do i=1,agridsize
        do j=1,Dgridsize
            do k=1,hpgridsize
                finwealth=anodes(i)-(1-theta)*Dnodes(j)
                wholder=(1+r)*finwealth+(1-delta)*Dnodes(j)*(1-f)
  
                if (elasticity .ne. 1.0) then
                EV(i,j,:,k,Tdie+1)=1.0/(1.0-elasticity)*const**(1-elasticity)*rent**((1-elasticity2)*(elasticity-1))*psi*(income(:,Tdie)/r+wholder)**(1-elasticity)
                else
                EV(i,j,:,k,Tdie+1)=psi*log(income(:,Tdie)/r+wholder)
                end if
                if (j==1 .and. k==1) then
                write(*,*) EV(i,j,1,k,Tdie+1)
                end if
                
                !if (finwealth<0) then
                !    EV(i,j,:,k,Tdie+1)=1/((1-elasticity))*((finwealth*(1+rborrow)+(1-delta)*dnodes(j)*exp(hpnodes(k))*(1-F)))**(1-elasticity)
                !else
                !    EV(i,j,:,k,Tdie+1)=1/((1-elasticity))*((finwealth*(1+r)+(1-delta)*dnodes(j)*exp(hpnodes(k))*(1-F)))**(1-elasticity)
                !end if
                !if (anodes(i)-(1-theta)*Dnodes(j)*exp(hpnodes(k))<0) then
                !    EV(i,j,:,k,Tdie+1)=1/(1-beta2)*((1+rborrow)*anodes(i)+dnodes(j)*(1-F)*(1-delta)*exp(hpnodes(k))-(1+rborrow)*(1-theta)*dnodes(j)*(1-F)*exp(hpnodes(k)))**(1-elasticity)/(1-elasticity)
                !else
                !    EV(i,j,:,k,Tdie+1)=1/(1-beta2)*((1+r)*anodes(i)+dnodes(j)*(1-F)*(1-delta)*exp(hpnodes(k))-(1+r)*(1-theta)*dnodes(j)*(1-F)*exp(hpnodes(k)))**(1-elasticity)/(1-elasticity)
                !end if
            end do
        end do
    end do
    !EV=0
    !EV=max(-10000000000.0,EV)

    write(*,*) "totEV ", sum(EV)

    !EV=0

   

   ! Main model code: Solves the backward induction
    call cpu_time(timeend)
    write(*,*) (timeend-timestart)/8.0
    call cpu_time(timestart)
    call solveworkingproblem
    call cpu_time(timeend)
    write(*,*) (timeend-timestart)/8.0


    call simulate
    write(901,901) param_iter,beta2,elasticity2,psi,r_rental,r_rental_retire,ret_wealth
    do t=1,11
    write(901,901), moments(t,:)
    end do



end do
end do
end do
end do
end do
end do


end program durables




subroutine simulate
USE nrtype; USE nrutil
USE nr
USE share
USE OMP_LIB
IMPLICIT NONE

REAL(8), dimension(numhouseholds,5) :: currenthouseholdstate, newhouseholdstate
REAL(8), dimension(numhouseholds,1) :: h_state_wrent
REAL(8), dimension(numhouseholds,Tdie) :: pholder,aggregatecontrib,mpc,consumption,consumptionmpc,consumptionhpshock,durableconsumption, currenttemp,currentperm,incomeholder,rentalind,durableinvestment,actualwealth,financialwealth,housingnet,housinggross,actualwealthdividedbyincome,totalnetworthdividedbyincome,taxrate,welfare,diffc, ratioc, qchg
REAL(8), dimension(numhouseholds,Tdie) :: changec,changetemp,changeperm,changey,changeyconditional,changed,changedconditional,changetempconditional,changepermconditional,demeanedindicator
REAL(8), dimension(numhouseholds,Tdie) :: alive
REAL(8), dimension(numhouseholds,Tdie) :: housingstock, rentalstock
REAL(8), dimension(numhouseholds,2) :: nearestnode
REAL(8), dimension(numhouseholds,Tdie,8) :: householdresultsholder
REAL(8), dimension(20) :: aggregatecontrib_a,aggregatecontrib_h,mpc_a,mpc_h,num_a,num_h,consumption_a,consumption_h,diffc_a,diffc_h
REAL(8), dimension(zgridsize) :: aggregatecontrib_y,mpc_y,num_y,consumption_y,diffc_y
REAL(8) :: shock
integer, dimension(2) :: seedvalue
integer :: j,t,i,k
REAL(8), dimension(numhouseholds) :: insurancetemp, insuranceperm, cov1, cov2, cov3, cov4, insurancedtemp,insurancedperm,numown, insurancedtempconditional,insurancedpermconditional
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
REAL(8) :: numrent
REAL(8) :: pretaxcalcholder,pretaxincome
REAL(8) :: tot1,tot2,tot3,adjustt,changedconditionalholder,changeyconditionalholder
REAL(8) :: numowntotal
REAL(8) :: rep
REAL(8) :: actualwealthtotal
REAL(8) :: exanteEVborn,exanteEVoverall,numobswelfare
REAL(8), dimension(hpgridsize,1) :: exanteEVbornstate, numobsstate, overallwelfarebystate,overallobsbystate, consumptionbystate, overallwelfarebystateyoung,overallwelfarebystatemiddle,overallwelfarebystateold,overallobsbystateyoung,overallobsbystatemiddle,overallobsbystateold, consumptionbystateyoung,consumptionbystatemiddle,consumptionbystateold
REAL(8), dimension(hpgridsize,350) :: overallwelfarebystateCE, overallwelfarebystateyoungCE,overallwelfarebystatemiddleCE,overallwelfarebystateoldCE
REAL(8), dimension(hpgridsize) :: cutoff
REAL(8), dimension(350,1) :: CE
REAL(8) :: medianincome, ratiocrent, ratiocown,diffcrent,diffcown,numrent_t,numown_t,numbuy_t,numsell_t
REAL(8), dimension(10):: numbins
REAL(8) :: leverage_t, a_t, leverageown_t,aown_t,leveragerent_t,arent_t

REAL(8) :: totalhousing, totalincome,totalfinancialassets,totalincomeowners

REAL(8), dimension(Tdie,6) :: lifecyclemomentsholder


REAL(8), dimension(numhouseholds,Tdie) :: stateindicator1,stateindicator2,stateindicator5

REAL(8) :: housepos, numberown,liquid,liquidown,liquidrent
integer ( kind = 4 ) seed
REAL(8) :: rnorm
real(8) :: r8_normal_ab
real(8), parameter :: mean=0

character(len=6) :: char_i     ! use your maximum expected len
character(len=60) :: filename





seed = 123456789

pholder(:,1)=1

totalhousing=0
totalincome=0
totalfinancialassets=0
totalincomeowners=0

aggregatecontrib_a=0
aggregatecontrib_h=0
aggregatecontrib_y=0
mpc_a=0
mpc_h=0
mpc_y=0
consumption_a=0
consumption_h=0
consumption_y=0
num_a=0
num_h=0
num_y=0
diffc_a=0
diffc_h=0
diffc_y=0



do i=1,350
CE(i,1)=.97+(1.0*i)/5000.0
!write(*,*) CE(i,1)
end do

tol=.000000001

OPEN (UNIT=2, FILE="singlehouseholdsim.txt", ACTION="WRITE", POSITION="REWIND")
2 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=3, FILE="agecoefficients.txt", ACTION="WRITE", POSITION="REWIND")
3 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

4 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=5, FILE="lifecycleprofiles.txt", ACTION="WRITE", POSITION="REWIND")
5 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=9, FILE="householdresults.txt", ACTION="WRITE", POSITION="REWIND")
9 format (I16.6, I16.6, F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=555, FILE="leverageprofiles.txt", ACTION="WRITE", POSITION="REWIND")
555 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=10, FILE="elasticity.txt", ACTION="WRITE", POSITION="REWIND")
10 format (F16.6, F16.6, F16.6,F16.6, F16.6, F16.6,F16.6, F16.6, F16.6,F16.6, F16.6, F16.6)

write(char_i, '(I5)') param_iter        ! convert integer to char
write(filename, '("mpc_results/statesandmpc_param_iter", A, ".txt")') trim(adjustl(char_i))
   
!open(11, file='mpc_results/statesandmpc_param_iter'//trim(str(param_iter))//'.txt')
open(11, file=''//filename//'.txt')
!write (11, *) i
!close (11)

!OPEN (UNIT=11, FILE="statesandmpc.txt", ACTION="WRITE", POSITION="REWIND")
11  format (F16.6, F16.6, F16.6,F16.6, F16.6, F16.6,F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6)


OPEN (UNIT=111, FILE="lifecyclemoments.txt", ACTION="WRITE", POSITION="REWIND")
111 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

OPEN (UNIT=222, FILE="elasticity_bystate.txt", ACTION="WRITE", POSITION="REWIND")
222 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

!OPEN (UNIT=12, FILE="allhouseholdstates.txt", ACTION="WRITE", POSITION="REWIND")
!12  format (F16.6, F16.6, F16.6,F16.6, F16.6, F16.6,F16.6, F16.6, F16.6,F16.6, F16.6, F16.6)

seedvalue(1)=1
seedvalue(2)=1
CALL random_seed(put=seedvalue(1:2))

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
    !if (income(currenthouseholdstate(i,3),1)/medianincome<.4632) then 
     if (Probinitcum(currenthouseholdstate(i,3))<.25) then    
            if (shock<.1156) then
            currenthouseholdstate(i,2)=2.3214
            currenthouseholdstate(i,1)=.3145
            else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.0
            end if
            
            

            numbins(1)=numbins(1)+1
    !elseif (income(currenthouseholdstate(i,3),1)/medianincome>=.4632 .and. income(currenthouseholdstate(i,3),1)/medianincome<1) then
    elseif (Probinitcum(currenthouseholdstate(i,3))>=.25 .and. Probinitcum(currenthouseholdstate(i,3))<0.5) then  
        if (shock<.2433) then
            currenthouseholdstate(i,2)=1.7857
            currenthouseholdstate(i,1)=.2556
            numbins(2)=numbins(2)+1
        else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.0075
            numbins(3)=numbins(3)+1
        end if
    !elseif (income(currenthouseholdstate(i,3),1)/medianincome>=1 .and. income(currenthouseholdstate(i,3),1)/medianincome<1.5455) then
    elseif (Probinitcum(currenthouseholdstate(i,3))>=.5 .and. Probinitcum(currenthouseholdstate(i,3))<0.75) then  
        if (shock<.2658) then
            currenthouseholdstate(i,2)=2.8571
            currenthouseholdstate(i,1)=.1611
            numbins(4)=numbins(4)+1
        else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.0289
            numbins(5)=numbins(5)+1
        end if
    else
        if (shock<.5406) then
            currenthouseholdstate(i,2)=4.0357
            currenthouseholdstate(i,1)=0.4902
            numbins(6)=numbins(6)+1
        else
            currenthouseholdstate(i,2)=0
            currenthouseholdstate(i,1)=0.1689
            numbins(7)=numbins(7)+1
        end if
    end if
    !currenthouseholdstate(:,1)=0
end do
currenthouseholdstate(:,1)=currenthouseholdstate(:,1)*exp(ageearnings(1))/sum(exp(ageearnings)/Tretire)
currenthouseholdstate(:,2)=currenthouseholdstate(:,2)*exp(ageearnings(1))/sum(exp(ageearnings)/Tretire)



write(*,*) numbins/numhouseholds




!currenthouseholdstate(:,4)=(hpgridsize+1)/2  !high house price state
!newhouseholdstate(:,4)=(hpgridsize+1)/2

currenthouseholdstate(:,4)=1
newhouseholdstate(:,4)=1



call random_number(shock)
alive=1






numowntotal=0
numrent=0
actualwealthtotal=0

housingstock=0
rentalstock=0



write(*,*) "start"
do t=1,Tdie

ratiocrent=0
ratiocown=0
diffcrent=0
diffcown=0
numrent_t=0
numown_t=0
numbuy_t=0
numsell_t=0

housepos=0 
numberown=0
liquidown=0
liquidrent=0
liquid=0


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


		actualwealth(i,t)=currenthouseholdstate(i,1)+theta*currenthouseholdstate(i,2)  !total wealth is voluntary equity plus equity in durable
        financialwealth(i,t)=currenthouseholdstate(i,1)-(1-theta)*currenthouseholdstate(i,2)
        housingnet(i,t)=theta*currenthouseholdstate(i,2)
        housinggross(i,t)=currenthouseholdstate(i,2)
		
		
        if (currenthouseholdstate(i,2)>0) then
            numberown=numberown+1
            housepos=housepos+currenthouseholdstate(i,2)
            liquidown=liquidown+financialwealth(i,t)
            liquid=liquid+financialwealth(i,t)
        else
            liquidrent=liquidrent+financialwealth(i,t)
            liquid=liquid+financialwealth(i,t)
        end if

        
        !qchg(i,t) = financialwealth(i,t) + (1-theta)*currenthouseholdstate(i,2)*exp(hpnodes(currenthouseholdstate(i,4)-1))
        
        
        
        

        qchg(i,t) = currenthouseholdstate(i,1) + (1-theta)*currenthouseholdstate(i,2)*shocksize
            
        if (t==1) then
            call pol_linworking(currenthouseholdstate(i,1),currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumption(i,t),rentalind(i,t),welfare(i,t))
            h_state_wrent(i,1)=newhouseholdstate(i,2)
        end if
        
        call pol_linworking(qchg(i,t),currenthouseholdstate(i,2)*(1+shocksize),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumptionhpshock(i,t),rentalind(i,t),welfare(i,t))
        !call pol_linworking(currenthouseholdstate(i,1)+.01,currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumptionmpc(i,t),rentalind(i,t),welfare(i,t))
		currenthouseholdstate(i,4)=2
        call pol_linworking(currenthouseholdstate(i,1),currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumptionmpc(i,t),rentalind(i,t),welfare(i,t))
		currenthouseholdstate(i,4)=1
        !call pol_linworking(max(borrowconstraint,currenthouseholdstate(i,1))+h_state_wrent(i,1)*.1,currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumptionmpc(i,t),rentalind(i,t),welfare(i,t))
		

        call pol_linworking(currenthouseholdstate(i,1),currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5),newhouseholdstate(i,1),newhouseholdstate(i,2),consumption(i,t),rentalind(i,t),welfare(i,t))
        
        if (newhouseholdstate(i,1)<0) then
            !write(*,*) "Some neg should exist"
        end if
        
        !mpc(i,t)=(consumptionmpc(i,t)-consumption(i,t))/(min(0.0,currenthouseholdstate(i,1)-borrowconstraint)+h_state_wrent(i,1)*0.1)
        !mpc(i,t)=(consumptionmpc(i,t)-consumption(i,t))/(.01)
        if (currenthouseholdstate(i,2)>0) then
            mpc(i,t)=(consumptionmpc(i,t)-consumption(i,t))/(currenthouseholdstate(i,2)*shocksize)
        else
            mpc(i,t)=0
        end if
        aggregatecontrib(i,t)=mpc(i,t)*currenthouseholdstate(i,2)*(1-delta)*(1-F)
		if (newhouseholdstate(i,1)>amax) then
            newhouseholdstate(i,1)=.9999999*amax
        end if
        diffc(i,t)=consumptionhpshock(i,t)-consumption(i,t)
        ratioc(i,t) = consumptionhpshock(i,t)/consumption(i,t)
		
		
		do j=1,20
			if (t<Tretire) then
				!if (currenthouseholdstate(i,1)<=amin+j*(amax*.5-amin)/20 .and. currenthouseholdstate(i,1)>amin+(j-1)*(amax*5-amin)/20) then
				if (currenthouseholdstate(i,1)<=amin+j*.5 .and. currenthouseholdstate(i,1)>amin+(j-1)*.5) then
					if (isnan(aggregatecontrib(i,t))) then
						write(*,*) "bad", i,t, mpc(i,t),currenthouseholdstate(i,2)*(1-delta)*(1-F),currenthouseholdstate(i,2),currenthouseholdstate(i,3),currenthouseholdstate(i,4),currenthouseholdstate(i,5)
					end if
					aggregatecontrib_a(j)=aggregatecontrib_a(j)+aggregatecontrib(i,t)
					mpc_a(j)=mpc_a(j)+mpc(i,t)
					consumption_a(j)=consumption_a(j)+consumption(i,t)
					num_a(j)=num_a(j)+1
					diffc_a(j)=diffc_a(j)+diffc(i,t)
				end if
				
				if (currenthouseholdstate(i,2)<=dmin+j*.5 .and. currenthouseholdstate(i,2)>dmin+(j-1)*.5) then
					if (isnan(aggregatecontrib(i,t))) then
						write(*,*) "bad2", i,t, mpc(i,t),currenthouseholdstate(i,2)*(1-delta)*(1-F)
					end if
					aggregatecontrib_h(j)=aggregatecontrib_h(j)+aggregatecontrib(i,t)
					mpc_h(j)=mpc_h(j)+mpc(i,t)
					consumption_h(j)=consumption_h(j)+consumption(i,t)
					num_h(j)=num_h(j)+1
					diffc_h(j)=diffc_h(j)+diffc(i,t)
				end if
			
			end if
		end do
		do j=1,zgridsize
			if (t<Tretire) then
				if (currenthouseholdstate(i,3)==j) then
					aggregatecontrib_y(j)=aggregatecontrib_y(j)+aggregatecontrib(i,t)
					mpc_y(j)=mpc_y(j)+mpc(i,t)
					consumption_y(j)=consumption_y(j)+consumption(i,t)
					num_y(j)=num_y(j)+1
					diffc_y(j)=diffc_y(j)+diffc(i,t)
				end if
			end if
		end do

        

        !if (i==numhouseholds) then
       ! write(*,*) "ind diffs"
        !write(*,*) t,sum(diffc(:,t)),sum(ratioc(:,t))
       ! end if

        !write(*,*) "here 2"
        householdresultsholder(i,t,1)=currenthouseholdstate(i,1)
        householdresultsholder(i,t,2)=currenthouseholdstate(i,2)
        householdresultsholder(i,t,3)=znodes(currenthouseholdstate(i,3))
        householdresultsholder(i,t,4)=alive(i,t)
        householdresultsholder(i,t,5)=consumption(i,t)
        householdresultsholder(i,t,6)=newhouseholdstate(i,1)
        householdresultsholder(i,t,7)=newhouseholdstate(i,2)
        householdresultsholder(i,t,8)=rentalind(i,t)
        if (t<=Tretire) then
        incomeholder(i,t)=income(currenthouseholdstate(i,3),t)
        else
        incomeholder(i,t)=retirementincome(currenthouseholdstate(i,3))
        end if

        if (i==1) then
            write(2,2) currenthouseholdstate(1,1), currenthouseholdstate(1,2),newhouseholdstate(1,2),consumption(1,t),rentalind(1,t), currenthouseholdstate(1,3), currenthouseholdstate(1,4), currenthouseholdstate(1,1)-(1-theta)*currenthouseholdstate(1,2)
        end if
        
        if (newhouseholdstate(i,1)<0.0001 .and. t==1) then
            constrained=constrained+1
        end if
        
        durableconsumption(i,t)=newhouseholdstate(i,2)+.0000001
        
        if (newhouseholdstate(i,2)-currenthouseholdstate(i,2)>.001 .and. rentalind(i,t)<=.995) then
        numbuy_t=numbuy_t+1
        elseif (newhouseholdstate(i,2)-currenthouseholdstate(i,2)<-.001) then
        numsell_t=numsell_t+1
        end if

        a_t=a_t+currenthouseholdstate(i,1)-(1-theta)*currenthouseholdstate(i,2)
        !if (currenthouseholdstate(i,2)/(currenthouseholdstate(i,1)+theta*currenthouseholdstate(i,2)) .ne. 0/0) then
        leverage_t=leverage_t+currenthouseholdstate(i,2)/(currenthouseholdstate(i,1)+theta*currenthouseholdstate(i,2))
        !end if

        
        if (currenthouseholdstate(i,2)==0) then
            numrent=numrent+1
            numrent_t=numrent_t+1
            ratiocrent=ratiocrent+(ratioc(i,t)-1)/(exp(hpmin)-1)

            arent_t=arent_t+currenthouseholdstate(i,1)-(1-theta)*currenthouseholdstate(i,2)
            !if (currenthouseholdstate(i,2)/(currenthouseholdstate(i,1)+theta*currenthouseholdstate(i,2)) .ne. 0/0) then
                leveragerent_t=leveragerent_t+currenthouseholdstate(i,2)/(currenthouseholdstate(i,1)+theta*currenthouseholdstate(i,2))
            !end if


        else
            numowntotal=numowntotal+1
            numown_t=numown_t+1
            ratiocown=ratiocown+(ratioc(i,t)-1)/(exp(hpmin)-1)

            aown_t=aown_t+currenthouseholdstate(i,1)-(1-theta)*currenthouseholdstate(i,2)
            !if (currenthouseholdstate(i,2)/(currenthouseholdstate(i,1)+theta*currenthouseholdstate(i,2)) .ne. 0/0) then
                leverageown_t=leverageown_t+currenthouseholdstate(i,2)/(currenthouseholdstate(i,1)+theta*currenthouseholdstate(i,2))
            !end if
        end if

         totalfinancialassets=totalfinancialassets+financialwealth(i,t)
         totalincome=totalincome+incomeholder(i,t)
         
         
        if (rentalind(i,t)>.995) then
            rentalstock(i,t)=newhouseholdstate(i,2)
			h_state_wrent(i,1)=newhouseholdstate(i,2)
            newhouseholdstate(i,2)=0
        end if

        call random_number(shock)
        if (shock <= Probhpshockcum(1)) then
                newhouseholdstate(i,1)=newhouseholdstate(i,1) + (1-theta)*(newhouseholdstate(i,2)/(1.0-mu+hpshock(1))-newhouseholdstate(i,2))
                newhouseholdstate(i,2)=newhouseholdstate(i,2)/(1.0-mu+hpshock(1))
                if (t>1) then
                pholder(i,t)=pholder(i,t-1)/(1.0-mu+hpshock(1))
                end if
        else
            do j=1,hpshockgridsize-1
                if (shock > Probhpshockcum(j) .AND. shock<=Probhpshockcum(j+1)) then
                   newhouseholdstate(i,1)=newhouseholdstate(i,1) + (1-theta)*(newhouseholdstate(i,2)/(1.0-mu+hpshock(j+1))-newhouseholdstate(i,2))
                    newhouseholdstate(i,2)=newhouseholdstate(i,2)/(1.0-mu+hpshock(j+1))
                    if (t>1) then
                        pholder(i,t)=pholder(i,t-1)/(1.0-mu+hpshock(j+1))
                    end if
                end if
            end do
        end if
        

        

        if (rentalind(i,t)<=.995) then
            housingstock(i,t)=newhouseholdstate(i,2)
			h_state_wrent(i,1)=newhouseholdstate(i,2)
            totalhousing=totalhousing+housingstock(i,t)
            totalincomeowners=totalincomeowners+incomeholder(i,t)
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

		call random_number(shock)
        if (t<Tretire .and. shock<.01) then
            write(11,11) mpc(i,t),aggregatecontrib(i,t),diffc(i,t), currenthouseholdstate(i,1), currenthouseholdstate(i,2), currenthouseholdstate(i,3), hpnodes(currenthouseholdstate(i,4)),currenthouseholdstate(i,5), rentalind(i,t), consumption(i,t),income(currenthouseholdstate(i,3),currenthouseholdstate(i,5))
			if (t<2 .and. shock<.01) then
				!write(*,*) mpc(i,t),aggregatecontrib(i,t),diffc(i,t),.1, currenthouseholdstate(i,1), currenthouseholdstate(i,2), currenthouseholdstate(i,3), hpnodes(currenthouseholdstate(i,4)),currenthouseholdstate(i,5), rentalind(i,t), consumption(i,t),income(currenthouseholdstate(i,3),currenthouseholdstate(i,5))
			end if
		end if

        !write(12,12) financialwealth(i,t), incomeholder(i,t)
        
    end do

    write(*,*) "MIN q vals", minval(currenthouseholdstate(:,1)), minval(newhouseholdstate(:,1))

    write(5,5) t*1.0,sum(currenthouseholdstate(:,1))/numhouseholds, sum(currenthouseholdstate(:,2))/numhouseholds,sum((1+r)*(currenthouseholdstate(:,1)+theta*currenthouseholdstate(:,2))-delta*currenthouseholdstate(:,2))/numhouseholds, sum(newhouseholdstate(:,2))/numhouseholds, sum(consumption(:,t))/numhouseholds, sum(rentalind(:,t))/numhouseholds,sum(currenthouseholdstate(:,1)-(1-theta)*currenthouseholdstate(:,2))/numhouseholds,sum(currenthouseholdstate(:,2)/(currenthouseholdstate(:,2)+currenthouseholdstate(:,1)-(1-theta)*currenthouseholdstate(:,2)))/numhouseholds

    write(111,111) t*1.0, numberown/numhouseholds, liquidown/numberown, liquidrent/(numhouseholds-numberown), liquid/numhouseholds,housepos/numberown,housepos/numhouseholds
    lifecyclemomentsholder(t,1)=numberown/numhouseholds
    lifecyclemomentsholder(t,2)=liquidown/numberown
    lifecyclemomentsholder(t,3)=liquidrent/(numhouseholds-numberown)
    lifecyclemomentsholder(t,4)=liquid/numhouseholds
    lifecyclemomentsholder(t,5)=housepos/numberown
    lifecyclemomentsholder(t,6)=housepos/numhouseholds

    !write(5,5) sum(financialwealth(:,t)*alive(:,t))/numhouseholds,sum(currenthouseholdstate(:,2)*alive(:,t))/numhouseholds, sum(consumption(:,t)*alive(:,t))/numhouseholds, sum(actualwealth(:,t)*alive(:,t))/numhouseholds, sum(newhouseholdstate(:,2)*alive(:,t))/numhouseholds, sum(newhouseholdstate(:,2)*rentalind(:,t)*alive(:,t))/numhouseholds, sum(rentalind(:,t)*alive(:,t))/numhouseholds, sum(alive(:,t))/numhouseholds, sum(consumption(:,t))/numhouseholds
     currenthouseholdstate=newhouseholdstate

    
    ! this is the housing price elasticity we want
    write(*,*) t*1.0, (sum(diffc(:,t))/sum(alive(:,t)))/(sum(consumption(:,t))/sum(alive(:,t)))/shocksize
    !write(*,*) "elasticity", sum(diffc(:,t)*alive(:,t))/sum(consumption(:,t))/sum(alive(:,t))/0.2


   ! write(*,*) t*1.0,sum(currenthouseholdstate(:,2)), sum(alive(:,t))
   
    !write(10,10) (sum((ratioc(:,t)-1))/sum(alive(:,t)))/(exp(hpmin) - 1), (sum(diffc(:,t))/sum(alive(:,t)))/(sum(consumption(:,t))/sum(alive(:,t)))/(hpmin)
    write(10,10) (sum((ratioc(:,t)-1))/numhouseholds)/shocksize, (sum(diffc(:,t))/numhouseholds)/(sum(consumption(:,t))/numhouseholds)/shocksize, ratiocown/numown_t,ratiocrent/numrent_t,numown_t/(numown_t+numrent_t),numsell_t/numhouseholds,numbuy_t/numhouseholds,(sum(diffc(:,t)/shocksize)/numhouseholds),sum(aggregatecontrib(:,t))/numhouseholds,sum(mpc(:,t))/numhouseholds,sum(currenthouseholdstate(:,2)*(1-delta)*(1-F))/numhouseholds

    write(555,555) a_t/numhouseholds, aown_t/numown_t, arent_t/numrent_t, leverage_t/numhouseholds, leverageown_t/numown_t, leveragerent_t/numrent_t
    
end do


do j=1,20
	if (j<=zgridsize) then
		write(222,222) shocksize,aggregatecontrib_a(j), mpc_a(j), consumption_a(j),num_a(j),diffc_a(j),aggregatecontrib_h(j), mpc_h(j), consumption_h(j),num_h(j),diffc_h(j),aggregatecontrib_y(j), mpc_y(j), consumption_y(j),num_y(j),diffc_y(j)
	else
		write(222,222) shocksize,aggregatecontrib_a(j), mpc_a(j), consumption_a(j),num_a(j),diffc_a(j),aggregatecontrib_h(j), mpc_h(j), consumption_h(j),num_h(j),diffc_h(j) 
	end if
end do


do t=1,11
    do i=1,6
    moments(t,i)=sum(lifecyclemomentsholder((t-1)*5+1:5*t,i))/5
    end do
end do



!do i=1,numhouseholds  
!    do t=1,60 
!      write(9,9) i,t, householdresultsholder(i,t,1),  householdresultsholder(i,t,2),  householdresultsholder(i,t,3),  householdresultsholder(i,t,4),  householdresultsholder(i,t,5),  householdresultsholder(i,t,6),  householdresultsholder(i,t,7),  householdresultsholder(i,t,8)
!    end do
!end do
write(*,*) "end"


write(*,*) "housing stock", sum(housingstock*alive)/sum(alive)
write(*,*) "rental stock", sum(rentalstock*alive)/sum(alive)
write(*,*) "total residential investment", sum(housingstock*delta*alive+rentalstock*(r_rental-r)*alive)/sum(alive)
write(*,*) "total non dur consumption", sum(consumption*alive)/sum(alive)

write(*,*) "H/(H+W)", sum(housingnet*alive)/(sum(actualwealth*alive)), sum(housinggross*alive)/(sum(actualwealth*alive)),sum(housinggross*alive)/(sum(actualwealth*alive-housingnet*alive+housinggross*alive))


write(*,*) "wealth over income", sum(actualwealth*alive)/sum(incomeholder*alive)  !Target = 1.52 (median net worth / earnings in 2001 SCF (see Kaplan and Violante table 2)
write(*,*) "financial wealth over income", sum(financialwealth*alive)/sum(incomeholder*alive)  !Target = -.35
write(*,*) "housing (owners) over income", sum(housingstock*alive)/sum(incomeholder*alive*numowntotal/(numowntotal+numrent)) !target 2.23


write(*,*) "C / I_d", sum(consumption*alive)/sum(housingstock*delta*alive+rentalstock*(r_rental-r)*alive)   !Target = 15 from BEA Nondura+Services / Residential investment 1999-2012 chained dollars table 1.1.6
write(*,*) "fracown", numowntotal/(numowntotal+numrent)  !Target 69% from SCF 98 (luengo-prado and diaz)

write(*,*) beta2, r_rental, elasticity2

write(*,*) "housing over GDP", sum(housingstock*alive+rentalstock*alive)/sum(incomeholder*alive)
write(*,*) "assets (debt if neg) over GDP", sum(financialwealth*alive)/sum(incomeholder*alive)

write(*,*) "fracown target .69", numowntotal/(numowntotal+numrent)  !Target 69% from SCF 98 (luengo-prado and diaz)
write(*,*) "total housing to total income (owners) target 2.56", totalhousing/totalincomeowners
write(*,*) "financial assets to total income (everyone target -0.3", totalfinancialassets/totalincome



diffmoments1=log(.69/(numowntotal/(numowntotal+numrent)))**2
diffmoments2=log(1.52/(sum(actualwealth*alive)/sum(incomeholder*alive)))**2
diffmoments3=log((sum(consumption*alive)/sum(housingstock*delta*alive+rentalstock*(r_rental-r)*alive))/15)**2
diffmoments=max(max(diffmoments2,diffmoments3),diffmoments1)
write(*,*) "diff", diffmoments,diffmoments1,diffmoments2,diffmoments3
diffmoments = 0.0




    


close(2)
close(3)
close(4)
close(5)
close(9)
close(10)
close(555)
close(11)
close(111)
close(222)

end subroutine simulate


















  

subroutine pol_linworking(astate,Dstate,zstate,hpstate,t,achoicelin,Dchoicelin,cchoicelin,rentallin,welfare)
    USE nrtype; USE nrutil
    USE nr
    USE share
    USE OMP_LIB
    IMPLICIT NONE
    REAL(8) :: weightal, weightDl
    REAL(8) :: nearestanode, nearestDnode
    REAL(8) :: al,ah, Dl,Dh
    REAL(8) :: astate,Dstate,zstate,hpstate,t,achoicelin,Dchoicelin,cchoicelin,rentallin,welfare
    
    if (astate>amax) then
    astate=amax-.00000000001
    elseif (astate<amin) then
    astate=amin+.000000000001
    end if
    
    if (Dstate>Dmax) then
    Dstate=Dmax-.00000000001
    elseif (Dstate<Dmin) then
    Dstate=Dmin
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

    achoicelin=weightal*weightDl*achoice(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*achoice(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*achoice(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*achoice(ah,Dh,zstate,hpstate,t)
    Dchoicelin=weightal*weightDl*Dchoice(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*Dchoice(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*Dchoice(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*Dchoice(ah,Dh,zstate,hpstate,t)
    cchoicelin=weightal*weightDl*cchoice(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*cchoice(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*cchoice(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*cchoice(ah,Dh,zstate,hpstate,t)
    rentallin=weightal*weightDl*rentalindicator(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*rentalindicator(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*rentalindicator(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*rentalindicator(ah,Dh,zstate,hpstate,t)
    welfare=weightal*weightDl*EV(al,Dl,zstate,hpstate,t)+(1-weightal)*weightDl*EV(ah,Dl,zstate,hpstate,t)+weightal*(1-weightDl)*EV(al,Dh,zstate,hpstate,t)+(1-weightal)*(1-weightDl)*EV(ah,Dh,zstate,hpstate,t)
    
end subroutine pol_linworking




subroutine plotpolicy_constantwealth
    USE nrtype; USE nrutil
    USE nr
    USE share
    USE OMP_LIB
    IMPLICIT NONE
    REAL(8) :: i,j,k,l,m,n,t, iter,p
    REAL(8) :: wealth, h, q, achoicelin, Dchoicelin, cchoicelin, rentallin, welfare

    OPEN (UNIT=266, FILE="policy_constantwealth.txt", ACTION="WRITE", POSITION="REWIND")
    266 format (F16.6, F16.6, F16.6, F16.6, F16.6, F16.6,F16.6, F16.6, F16.6, F16.6, F16.6, F16.6)

    t=3.0
    p=2.0
    do k=1,zgridsize
    do i=1,50
        wealth=0+.1*i
        do j=1,50
            h=0+(j-1)*.1
            q=wealth-theta*h
            call pol_linworking(q,h,k,p,t,achoicelin,Dchoicelin,cchoicelin,rentallin,welfare)
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
    INTEGER :: i,j,k,l,m,t
    REAL(8) :: timestart2, timeend2
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,hpgridsize,2) :: optpolicynoadjust,optpolicyadjust
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,hpgridsize,3,2) :: pstartadjust
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,hpgridsize,3) :: ystartadjust
    REAL(8), dimension(agridsize,Dgridsize,zgridsize,hpgridsize) :: ax,bx,cx, adjust
    INTEGER, dimension(agridsize,Dgridsize,zgridsize,hpgridsize) :: iter, counts_neg_all,counts_neg_noadjust
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
    
    

    ftol=.000000001
    iter=0
    
    adjust=0
    
    do t=Tdie,1,-1
        write(*,*) t
        counts_neg_all=0
        counts_neg_noadjust=0
        !$OMP PARALLEL
        !$OMP DO
        do i=1,agridsize
            do j=1,Dgridsize
                do k=1,zgridsize
                    do l=1,hpgridsize
                        pstartadjust(i,j,k,l,1,1)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,1,2)=.8*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,1)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,2)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,1)=.5*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,2)=.2*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))

                        
                    
                        state(i,j,k,l,1)=anodes(i)
                        state(i,j,k,l,2)=Dnodes(j)
                        state(i,j,k,l,3)=k
                        state(i,j,k,l,4)=l
                        state(i,j,k,l,5)=t
                    
                        ystartadjust(i,j,k,l,1)=valfuncadjust2(pstartadjust(i,j,k,l,1,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,2)=valfuncadjust2(pstartadjust(i,j,k,l,2,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,3)=valfuncadjust2(pstartadjust(i,j,k,l,3,:),state(i,j,k,l,:))

                    
                        call amoeba(pstartadjust(i,j,k,l,:,:),ystartadjust(i,j,k,l,:),ftol,valfuncadjust,iter(i,j,k,l),state(i,j,k,l,:))
                        Vadjust(i,j,k,l,t)=ystartadjust(i,j,k,l,1)
                        achoiceadjust(i,j,k,l,t)=pstartadjust(i,j,k,l,1,1)
                        Dchoiceadjust(i,j,k,l,t)=pstartadjust(i,j,k,l,1,2)
                    
                    
                        pstartadjust(i,j,k,l,1,1)=0.0
                        pstartadjust(i,j,k,l,1,2)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,1)=.05*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,2)=.1*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,1)=.4*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,2)=.21*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        
                        ystartadjust(i,j,k,l,1)=valfuncadjust2(pstartadjust(i,j,k,l,1,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,2)=valfuncadjust2(pstartadjust(i,j,k,l,2,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,3)=valfuncadjust2(pstartadjust(i,j,k,l,3,:),state(i,j,k,l,:))

                        call amoeba(pstartadjust(i,j,k,l,:,:),ystartadjust(i,j,k,l,:),ftol,valfuncadjust,iter(i,j,k,l),state(i,j,k,l,:))
                        if (ystartadjust(i,j,k,l,1)<Vadjust(i,j,k,l,t)) then  !(again we're minimizing)
                            Vadjust(i,j,k,l,t)=ystartadjust(i,j,k,l,1)
                            achoiceadjust(i,j,k,l,t)=pstartadjust(i,j,k,l,1,1)
                            Dchoiceadjust(i,j,k,l,t)=pstartadjust(i,j,k,l,1,2)
                        end if


                        







                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        pstartadjust(i,j,k,l,1,1)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,1,2)=.8*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,1)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,2)=.25*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,1)=.5*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,2)=.2*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                    
                        state(i,j,k,l,1)=anodes(i)
                        state(i,j,k,l,2)=Dnodes(j)
                        state(i,j,k,l,3)=k
                        state(i,j,k,l,4)=l
                        state(i,j,k,l,5)=t
                    
                        ystartadjust(i,j,k,l,1)=valfuncrent2(pstartadjust(i,j,k,l,1,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,2)=valfuncrent2(pstartadjust(i,j,k,l,2,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,3)=valfuncrent2(pstartadjust(i,j,k,l,3,:),state(i,j,k,l,:))
                    
                        call amoeba(pstartadjust(i,j,k,l,:,:),ystartadjust(i,j,k,l,:),ftol,valfuncrent,iter(i,j,k,l),state(i,j,k,l,:))
                    
                    
                        Vrent(i,j,k,l,t)=ystartadjust(i,j,k,l,1)
                        achoicerent(i,j,k,l,t)=pstartadjust(i,j,k,l,1,1)
                        Dchoicerent(i,j,k,l,t)=pstartadjust(i,j,k,l,1,2)
                    
                    
                        pstartadjust(i,j,k,l,1,1)=0.0
                        pstartadjust(i,j,k,l,1,2)=.01*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,1)=.05*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,2,2)=.1*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,1)=.4*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                        pstartadjust(i,j,k,l,3,2)=.21*((1+r)*anodes(i)+income(k,t)-(1+r)*(1-theta)*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
 
                    
                        ystartadjust(i,j,k,l,1)=valfuncrent2(pstartadjust(i,j,k,l,1,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,2)=valfuncrent2(pstartadjust(i,j,k,l,2,:),state(i,j,k,l,:))
                        ystartadjust(i,j,k,l,3)=valfuncrent2(pstartadjust(i,j,k,l,3,:),state(i,j,k,l,:))

                        call amoeba(pstartadjust(i,j,k,l,:,:),ystartadjust(i,j,k,l,:),ftol,valfuncrent,iter(i,j,k,l),state(i,j,k,l,:))
                    
                    
                        if (ystartadjust(i,j,k,l,1)<Vrent(i,j,k,l,t)) then  !(again we're minimizing)
                            Vrent(i,j,k,l,t)=ystartadjust(i,j,k,l,1)
                            achoicerent(i,j,k,l,t)=pstartadjust(i,j,k,l,1,1)
                            Dchoicerent(i,j,k,l,t)=pstartadjust(i,j,k,l,1,2)
                        end if
                        

                        if (costlyequity==1) then

                            if (anodes(i)>=(1-theta)*Dnodes(j)) then
                                ax(i,j,k,l)=(1-theta)*Dnodes(j)
                            else
                                ax(i,j,k,l)=amort_rate*anodes(i)+(1-theta)*(1-amort_rate)*Dnodes(j)  !this comes from setting a'=amort_rate*a' and using q=(1-theta)h+a
                            end if
							bx(i,j,k,l)=ax(i,j,k,l)
                        else
							ax(i,j,k,l)=min(state(i,j,k,l,1),borrowconstraint)
							bx(i,j,k,l)=borrowconstraint
                        end if
                        cx(i,j,k,l)=amax
           
                        state(i,j,k,l,2)=j
                        Vnoadjust(i,j,k,l,t)=brentnew(ax(i,j,k,l),bx(i,j,k,l),cx(i,j,k,l),valfuncnoadjust,ftol,state(i,j,k,l,1),state(i,j,k,l,2),state(i,j,k,l,3),state(i,j,k,l,4),state(i,j,k,l,5),achoicenoadjust(i,j,k,l,t))
                        Dchoicenoadjust(i,j,k,l,t)=Dnodes(j)
                        
                        if (achoicenoadjust(i,j,k,l,t)<-.05) then
                            counts_neg_noadjust(i,j,k,l)=counts_neg_noadjust(i,j,k,l)+1
                        end if
                        
        

                        if (Vadjust(i,j,k,l,t)<Vnoadjust(i,j,k,l,t) .and. Vadjust(i,j,k,l,t)<Vrent(i,j,k,l,t)) then  ! since V = - V from minimization
                            achoice(i,j,k,l,t)=achoiceadjust(i,j,k,l,t)
                            Dchoice(i,j,k,l,t)=Dchoiceadjust(i,j,k,l,t)
                            if (anodes(i)-(1-theta)*Dnodes(j)<0) then
                                cchoice(i,j,k,l,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)-theta*Dchoice(i,j,k,l,t)-achoice(i,j,k,l,t)-(1+rborrow)*(1-theta)*Dnodes(j)+(l-1)*Dnodes(j)*(1-F)*shocksize
                            else
                                cchoice(i,j,k,l,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)-theta*Dchoice(i,j,k,l,t)-achoice(i,j,k,l,t)-(1+r)*(1-theta)*Dnodes(j)+(l-1)*Dnodes(j)*(1-F)*shocksize
                            end if
                            rentalindicator(i,j,k,l,t)=0
                        elseif (Vnoadjust(i,j,k,l,t)<Vadjust(i,j,k,l,t) .and. Vnoadjust(i,j,k,l,t)<Vrent(i,j,k,l,t)) then
                            achoice(i,j,k,l,t)=achoicenoadjust(i,j,k,l,t)
                            Dchoice(i,j,k,l,t)=Dchoicenoadjust(i,j,k,l,t)


                            if (anodes(i)-(1-theta)*Dnodes(j)<0) then
                                cchoice(i,j,k,l,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)-theta*Dchoice(i,j,k,l,t)-achoice(i,j,k,l,t)-(1+rborrow)*(1-theta)*Dnodes(j)+(l-1)*Dnodes(j)*(1-F)*shocksize
                            else
                                cchoice(i,j,k,l,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)-theta*Dchoice(i,j,k,l,t)-achoice(i,j,k,l,t)-(1+r)*(1-theta)*Dnodes(j)+(l-1)*Dnodes(j)*(1-F)*shocksize
                            end if
                            rentalindicator(i,j,k,l,t)=0
                        else
                            achoice(i,j,k,l,t)=achoicerent(i,j,k,l,t)
                            Dchoice(i,j,k,l,t)=Dchoicerent(i,j,k,l,t)
                            if (anodes(i)-(1-theta)*Dnodes(j)<0) then
                                cchoice(i,j,k,l,t)=income(k,t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)-r_rental_t(t)*Dchoice(i,j,k,l,t)-achoice(i,j,k,l,t)-(1+rborrow)*(1-theta)*Dnodes(j)+(l-1)*Dnodes(j)*(1-F)*shocksize
                            else
                                cchoice(i,j,k,l,t)=income(k,t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)-r_rental_t(t)*Dchoice(i,j,k,l,t)-achoice(i,j,k,l,t)-(1+r)*(1-theta)*Dnodes(j)+(l-1)*Dnodes(j)*(1-F)*shocksize
                            end if
                            rentalindicator(i,j,k,l,t)=1
                        end if
                        
                        if (achoice(i,j,k,l,t)<-.05) then
                            counts_neg_all(i,j,k,l)=counts_neg_all(i,j,k,l)+1
                        end if
                        
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

       !write(*,*) "Policy functions"
       write(*,*) t
       !write(*,*) sum(cchoice(:,:,:,1,t)), sum(achoice(:,:,:,1,t)), sum(dchoice(:,:,:,1,t))
       !write(*,*) sum(counts_neg_all),sum(counts_neg_noadjust), minval(achoice)
       write(*,*) sum(cchoice(:,:,:,1,t)),sum(cchoice(:,:,:,2,t))
end do

     
end subroutine solveworkingproblem










REAL(8) FUNCTION valfuncadjust(p,state)  !p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
USE SHARE
IMPLICIT NONE
REAL(8), DIMENSION(2) :: p
REAL(8), DIMENSION(5) :: state
REAL(8) :: a, D, z, hp, t, zindex,hpindex
REAL(8) :: currentincome
REAL(8) :: aprime, Dprime, consumption, Dcurrent,aprimeholder2
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder,hpholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)

hpholder=hpnodes(hpindex)

currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if
aprimeholder2=aprime


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
    if ((a-D*(1-theta))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)-theta*Dcurrent-aprimeholder-(1+rborrow)*(1-theta)*D)+(hpindex-1)*D*(1-F)*shocksize
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)-theta*Dcurrent-aprimeholder-(1+r)*(1-theta)*D)+(hpindex-1)*D*(1-F)*shocksize
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
             do j=1,hpshockgridsize
                Dprime=Dcurrent*1.0/(1.0-mu+hpshock(j))
                if (Dprime>Dnodes(Dgridsize)) then
                    Dprime=Dnodes(Dgridsize)
                elseif (Dprime<Dnodes(1)) then
                    Dprime=Dnodes(1)
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
                aprime=aprimeholder2+(1-theta)*(Dprime-Dcurrent)
                if (aprime>anodes(agridsize)) then
                    aprime=anodes(agridsize)
                elseif (aprime<anodes(1)) then
                    aprime=anodes(1)
                end if
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
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,i,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,i,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,i,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,i,1,t+1)
                end do
            end do
            valfuncadjust=valfuncadjust+beta2*EVholder
        else
            do j=1,hpshockgridsize
                Dprime=Dcurrent*1.0/(1.0-mu+hpshock(j))
                if (Dprime>Dnodes(Dgridsize)) then
                    Dprime=Dnodes(Dgridsize)
                elseif (Dprime<Dnodes(1)) then
                    Dprime=Dnodes(1)
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
                aprime=aprimeholder2+(1-theta)*(Dprime-Dcurrent)
                if (aprime>anodes(agridsize)) then
                    aprime=anodes(agridsize)
                elseif (aprime<anodes(1)) then
                    aprime=anodes(1)
                end if
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



                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,1,t+1)

                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,1,Tdie+1)
                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,1,Tdie+1)
                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,1,Tdie+1)
                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,1,Tdie+1)
            end do
            valfuncadjust=valfuncadjust+beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if

        
    else
        valfuncadjust=-10000000000000
    end if
else
     valfuncadjust=-10000000000000
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
REAL(8) :: pholder,hpholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)

hpholder=hpnodes(hpindex)

currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
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
    if ((a-D*(1-theta))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)-r_rental_t(t)*Dcurrent-aprimeholder-(1+rborrow)*(1-theta)*D)+(hpindex-1)*D*(1-F)*shocksize
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)-r_rental_t(t)*Dcurrent-aprimeholder-(1+r)*(1-theta)*D)+(hpindex-1)*D*(1-F)*shocksize
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
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*EV(aprimel,1,i,1,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*EV(aprimeh,1,i,1,t+1)
            end do
            valfuncrent=valfuncrent+beta2*EVholder
        else
                EVholder=EVholder+weightaprimel*EV(aprimel,1,zindex,1,t+1)
                EVholder=EVholder+(1-weightaprimel)*EV(aprimeh,1,zindex,1,t+1)

                EVholder2=EVholder2+weightaprimel*EV(aprimel,1,zindex,1,Tdie+1)
                EVholder2=EVholder2+(1-weightaprimel)*EV(aprimeh,1,zindex,1,Tdie+1)
            valfuncrent=valfuncrent+xfactor*beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if
        
        
    else
        valfuncrent=-10000000000000
    end if
else
     valfuncrent=-10000000000000
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
REAL(8) :: aprime, Dprime, consumption, Dcurrent,aprimeholder2
REAL(8) :: nearestDnode, nearestanode, Dprimel, Dprimeh, aprimel, aprimeh, weightDprimel, weightaprimel
INTEGER i,j
REAL(8) :: laborsupply
REAL(8) :: aprimeholder
REAL(8) :: EVholder, EVholder2
REAL(8) :: Vadjustretireholder, Vnoadjustretireholder, Vretirerentholder
REAL(8) :: pholder,hpholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)

hpholder=hpnodes(hpindex)

currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if
aprimeholder2=aprime


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
    if ((a-D*(1-theta))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)-theta*Dcurrent-aprimeholder-(1+rborrow)*(1-theta)*D)+(hpindex-1)*D*(1-F)*shocksize
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)-theta*Dcurrent-aprimeholder-(1+r)*(1-theta)*D)+(hpindex-1)*D*(1-F)*shocksize
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
             do j=1,hpshockgridsize
                Dprime=Dcurrent*1.0/(1.0-mu+hpshock(j))
                if (Dprime>Dnodes(Dgridsize)) then
                    Dprime=Dnodes(Dgridsize)
                elseif (Dprime<Dnodes(1)) then
                    Dprime=Dnodes(1)
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
                aprime=aprimeholder2+(1-theta)*(Dprime-Dcurrent)
                if (aprime>anodes(agridsize)) then
                    aprime=anodes(agridsize)
                elseif (aprime<anodes(1)) then
                    aprime=anodes(1)
                end if
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
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,i,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,i,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,i,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,i,1,t+1)
                end do
            end do
            valfuncadjust2=valfuncadjust2+beta2*EVholder
        else
            do j=1,hpshockgridsize
                Dprime=Dcurrent*1.0/(1.0-mu+hpshock(j))
                if (Dprime>Dnodes(Dgridsize)) then
                    Dprime=Dnodes(Dgridsize)
                elseif (Dprime<Dnodes(1)) then
                    Dprime=Dnodes(1)
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
                aprime=aprimeholder2+(1-theta)*(Dprime-Dcurrent)
                if (aprime>anodes(agridsize)) then
                    aprime=anodes(agridsize)
                elseif (aprime<anodes(1)) then
                    aprime=anodes(1)
                end if
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



                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,1,t+1)

                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,1,Tdie+1)
                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,1,Tdie+1)
                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,1,Tdie+1)
                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,1,Tdie+1)
            end do
            valfuncadjust2=valfuncadjust2+beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if

        
    else
        valfuncadjust2=-10000000000000
    end if
else
     valfuncadjust2=-10000000000000
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
REAL(8) :: pholder,hpholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)

hpholder=hpnodes(hpindex)

currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
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
    if ((a-D*(1-theta))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)-r_rental_t(t)*Dcurrent-aprimeholder-(1+rborrow)*(1-theta)*D)+(hpindex-1)*D*(1-F)*shocksize
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)-r_rental_t(t)*Dcurrent-aprimeholder-(1+r)*(1-theta)*D)+(hpindex-1)*D*(1-F)*shocksize
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
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*EV(aprimel,1,i,1,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*EV(aprimeh,1,i,1,t+1)
            end do
            valfuncrent2=valfuncrent2+beta2*EVholder
        else
                EVholder=EVholder+weightaprimel*EV(aprimel,1,zindex,1,t+1)
                EVholder=EVholder+(1-weightaprimel)*EV(aprimeh,1,zindex,1,t+1)

                EVholder2=EVholder2+weightaprimel*EV(aprimel,1,zindex,1,Tdie+1)
                EVholder2=EVholder2+(1-weightaprimel)*EV(aprimeh,1,zindex,1,Tdie+1)
            valfuncrent2=valfuncrent2+xfactor*beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if
        
        
    else
        valfuncrent2=-10000000000000
    end if
else
     valfuncrent2=-10000000000000
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
REAL(8) :: pholder,hpholder
!write(*,*) state

!write(*,*) state

a=state(1)
D=state(2)
zindex=state(3)
hpindex=state(4)
t=state(5)

hpholder=hpnodes(hpindex)

currentincome=income(zindex,t)

aprime=p(1)
Dcurrent=p(2)
Dprime=Dcurrent
aprimeholder=aprime
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
    if ((a-D*(1-theta))<0) then
        consumption=(currentincome+(1+rborrow)*a+D*(1-F)*(1-delta)-r_rental_t(t)*Dcurrent-aprimeholder-(1+rborrow)*(1-theta)*D)+(hpindex-1)*D*(1-F)*shocksize
    else
        consumption=(currentincome+(1+r)*a+D*(1-F)*(1-delta)-r_rental_t(t)*Dcurrent-aprimeholder-(1+r)*(1-theta)*D)+(hpindex-1)*D*(1-F)*shocksize
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
                    EVholder=EVholder+Probz(zindex,i)*weightaprimel*EV(aprimel,1,i,1,t+1)
                    EVholder=EVholder+Probz(zindex,i)*(1-weightaprimel)*EV(aprimeh,1,i,1,t+1)
            end do
            valfuncrent3=valfuncrent3+beta2*EVholder
        else
                EVholder=EVholder+weightaprimel*EV(aprimel,1,zindex,1,t+1)
                EVholder=EVholder+(1-weightaprimel)*EV(aprimeh,1,zindex,1,t+1)

                EVholder2=EVholder2+weightaprimel*EV(aprimel,1,zindex,1,Tdie+1)
                EVholder2=EVholder2+(1-weightaprimel)*EV(aprimeh,1,zindex,1,Tdie+1)
            valfuncrent3=valfuncrent3+xfactor*beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if
        
        
    else
        valfuncrent3=-10000000000000
    end if
else
     valfuncrent3=-10000000000000
end if



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
REAL(8) :: hpholder
REAL(8) :: aprimeholder2
REAL(8) :: constraint
!write(*,*) state

!write(*,*) state

currentincome=income(zindex,t)

hpholder=hpnodes(hpindex)

Dcurrent=Dnodes(Dindex)
Dprime=Dnodes(Dindex)
aprimeholder=aprime
if (aprime>anodes(agridsize)) then
aprime=anodes(agridsize)
end if
aprimeholder2=aprime


if (allowunderwater==1) then
constraint=min(a,borrowconstraint)
else
constraint=borrowconstraint
end if

if (aprime>=constraint .AND. aprime<=anodes(agridsize)) then
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
    
    !using that Dprime=D
    if ((a-Dcurrent*(1-theta))<0) then
        consumption=(currentincome+(1+rborrow)*a+Dprime*(1-delta)-theta*Dprime-aprimeholder-(1+rborrow)*(1-theta)*Dprime)+(hpindex-1)*Dprime*(1-F)*shocksize
    else
        consumption=(currentincome+(1+r)*a+Dprime*(1-delta)-theta*Dprime-aprimeholder-(1+r)*(1-theta)*Dprime)+(hpindex-1)*Dprime*(1-F)*shocksize
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
             do j=1,hpshockgridsize
                Dprime=Dcurrent*1.0/(1.0-mu+hpshock(j))
                if (Dprime>Dnodes(Dgridsize)) then
                    Dprime=Dnodes(Dgridsize)
                elseif (Dprime<Dnodes(1)) then
                    Dprime=Dnodes(1)
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
                aprime=aprimeholder2+(1-theta)*(Dprime-Dcurrent)
                if (aprime>anodes(agridsize)) then
                    aprime=anodes(agridsize)
                elseif (aprime<anodes(1)) then
                    aprime=anodes(1)
                end if
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
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,i,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,i,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,i,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probz(zindex,i)*Probhpshock(j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,i,1,t+1)
                end do
            end do
            valfuncnoadjust=valfuncnoadjust+beta2*EVholder
        else
            do j=1,hpshockgridsize
                Dprime=Dcurrent*1.0/(1.0-mu+hpshock(j))
                if (Dprime>Dnodes(Dgridsize)) then
                    Dprime=Dnodes(Dgridsize)
                elseif (Dprime<Dnodes(1)) then
                    Dprime=Dnodes(1)
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
                aprime=aprimeholder2+(1-theta)*(Dprime-Dcurrent)
                if (aprime>anodes(agridsize)) then
                    aprime=anodes(agridsize)
                elseif (aprime<anodes(1)) then
                    aprime=anodes(1)
                end if
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



                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,1,t+1)
                    EVholder=EVholder+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,1,t+1)

                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*weightDprimel*EV(aprimel,Dprimel,zindex,1,Tdie+1)
                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*weightDprimel*EV(aprimeh,Dprimel,zindex,1,Tdie+1)
                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*weightaprimel*(1-weightDprimel)*EV(aprimel,Dprimeh,zindex,1,Tdie+1)
                    EVholder2=EVholder2+((1.0+mu-hpshock(j)))**(-(1-elasticity)*(1-elasticity2))*Probhpshock(j)*(1-weightaprimel)*(1-weightDprimel)*EV(aprimeh,Dprimeh,zindex,1,Tdie+1)
            end do
            valfuncnoadjust=valfuncnoadjust+beta2retire*(1-deathprob(t-Tretire))*EVholder+beta2retire*deathprob(t-Tretire)*EVholder2
        end if

        
        
    else
        valfuncnoadjust=-10000000000000
    end if
else
     valfuncnoadjust=-10000000000000
end if



valfuncnoadjust=-valfuncnoadjust  !powell minimizes

aprime=aprimeholder

END FUNCTION valfuncnoadjust











SUBROUTINE amoeba(p,y,ftol,valfuncadjust,iter,stateholder)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(SP), INTENT(IN) :: ftol
	REAL(SP), DIMENSION(3), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(3,2), INTENT(INOUT) :: p
    REAL(SP), DIMENSION(3) :: yholder
	REAL(SP), DIMENSION(3,2) :: pholder
	REAL(SP), DIMENSION(5), INTENT(IN) :: stateholder
	INTERFACE
		FUNCTION valfuncadjust(x,stateholder)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(2), INTENT(IN) :: x
		REAL(8), DIMENSION(5), INTENT(IN) :: stateholder
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
		
		rtol=2.0_sp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
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
		
		ytry=amotry(-1.0_sp,stateholder)
		iter=iter+1
		if (ytry <= y(ilo)) then
		
		!write(*,*) stateholder
		
			ytry=amotry(2.0_sp,stateholder)
			iter=iter+1
		else if (ytry >= y(inhi)) then
			ysave=y(ihi)
			
			!write(*,*) "3", stateholder
			
			ytry=amotry(0.5_sp,stateholder)
			iter=iter+1
			if (ytry >= ysave) then
				p(:,:)=0.5_sp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
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
	fac1=(1.0_sp-fac)/ndim
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



INTEGER(I4B), PARAMETER :: ITMAX=200
REAL(SP), PARAMETER :: CGOLD=0.3819660_sp,ZEPS=1.0e-3_sp*epsilon(ax)
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
    xm=0.5_sp*(a+b)
    tol1=tol*abs(x)+ZEPS
    tol2=2.0_sp*tol1
    if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) then
        xmin=x
        brentnew=fx
        RETURN
    end if
    if (abs(e) > tol1) then
        r2=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r2
        q=2.0_sp*(q-r2)
        if (q > 0.0) p=-p
        q=abs(q)
        etemp=e
        e=d
        if (abs(p) >= abs(0.5_sp*q*etemp) .or. &
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
REAL(SP), PARAMETER :: CGOLD=0.3819660_sp,ZEPS=1.0e-3_sp*epsilon(ax)
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
    xm=0.5_sp*(a+b)
    tol1=tol*abs(x)+ZEPS
    tol2=2.0_sp*tol1
    if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) then
        xmin=x
        brentmindist=fx
        RETURN
    end if
    if (abs(e) > tol1) then
        r2=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r2
        q=2.0_sp*(q-r2)
        if (q > 0.0) p=-p
        q=abs(q)
        etemp=e
        e=d
        if (abs(p) >= abs(0.5_sp*q*etemp) .or. &
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
    
    
character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str