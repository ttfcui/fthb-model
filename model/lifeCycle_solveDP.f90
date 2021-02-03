module lifecycle_solveDP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE SOLVEDP contains routines that solve the main dynamic programming
! problem, as well as certain debug/test subroutines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    USE share
    USE lifecycle_algs
    USE lifecycle_vfuncs
    USE OMP_LIB
    IMPLICIT NONE

    CONTAINS

    subroutine solveworkingproblem(price, balancer, achoice, Dchoice, rentchoice, cchoice, &
                                   choiceindicator, achoiceMov, dchoiceMov, cchoiceMov,&
                                   choiceindicatorMov,achoiceMovR, dchoiceMovR, cchoiceMovR,&
                                   choiceindicatorMovR, noShock, repl, expected, &
                                   EV, EVMov, EVMovR, state_end)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    This is the core of the module, a subroutine which solves
        !    for optimal housing choice over the lifecycle through backwards
        !    induction.
        !
        !    Modified: 08/09/2019
        !
        !    PARAMETERS
        !
        !    price: The guess for the housing price that balances supply and
        !    demand in stationary equilibrium, or along the transition path.
        !    Vector-valued. Set a vector of length 1 for PE.
        !
        !    balancer: The level of government transfers needed to balance
        !    the residual of government expenditure after taxes.
        !
        !    *choice: The actual choice and next-period asset policy functions
        !    which maximizes utility over all three possible housing decisions.
        !    Choose arrays defined outside of the module that store these results.
        !    
        !    *choiceMov: Same as *choice, but maximized assuming an adjustment
        !    must be made (when household encounters a moving shock).
        !
        !    *choiceMovR: Same as *choice, but assuming forced home purchase
        !    (Renters may be forced into buying)
        !
        !    choiceindicator/choiceindicatorMov: Indicates which housing
        !    decision was taken. Categorical variable taking three values.
        !    
        !    noShock: if TRUE, does nothing. If FALSE, indicates that a
        !    temporary subsidy is in effect during the current period.
        !
        !    repl: short for "Replace": if TRUE, writes value functions computed
        !    in this period into memory for usage in the DP problem. The one time
        !    when *not* writing it into memory is acceptable is in the case of
        !    policy transition assuming prices are held constant throughout
        !    the path.
        !
        !    expected: String. If one of accepted strings, indicates if DP
        !    problem accounts for expectations of a continued policy in the
        !    period after. Different strings indicate which agents retain
        !    eligibility and needs to expect policy in future.
        !
        !    EV: Optimal values over different states for continuation value    
        !    calculations.
        !    
        !    state_end: Optional vector of integers that specifies if the
        !    DP program only needs solving over a subset of the state space
        !    (i.e. when crunching choice under temporary policies)
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        INTEGER :: i, i2, j, k, t, l
        INTEGER :: pricelen, p1, p2, d_end, t_end, ue_start
        LOGICAL, INTENT(IN) :: noShock, repl
        CHARACTER(LEN=*), INTENT(IN) :: expected
        REAL(8), dimension(:), INTENT(IN) :: price, balancer
        REAL(8), dimension(:,:,:,:,:,:), INTENT(INOUT) :: cchoice, achoice, dchoice, &
            rentchoice, choiceindicator, cchoiceMov, achoiceMov, dchoiceMov, rentchoiceMov, &
            choiceindicatorMov, cchoiceMovR, achoiceMovR, dchoiceMovR, rentchoiceMov,&
            choiceindicatorMovR
        REAL(8), dimension(:,:,:,:,:,:), INTENT(INOUT) :: EV, EVMov, EVMovR
        INTEGER, dimension(2), INTENT(IN), OPTIONAL :: state_end
        REAL(8), dimension(:,:,:,:,:,:), ALLOCATABLE :: achoiceadjust, achoicenoadjust, &
            achoicerent, Dchoiceadjust, Dchoicenoadjust, Dchoicerent, cchoiceadjust, &
            cchoicenoadjust, cchoicerent, Vnoadjust, Vadjust, Vrent
        REAL(8), dimension(:,:,:,:,:,:,:), ALLOCATABLE :: pstartadjust
        REAL(8), dimension(:,:,:,:,:,:), ALLOCATABLE :: ystartadjust
        REAL(8), dimension(:,:,:,:,:), ALLOCATABLE :: ax, bx, cx, adjust
        REAL(8), dimension(:,:,:,:,:,:), ALLOCATABLE :: state, assetHolder, wealthHolder, incomeHolder, costholder, costDenom, myUsercost
        REAL(8), dimension(:,:,:,:,:,:,:), ALLOCATABLE :: transferholder
        REAL(8), dimension(:,:,:,:), ALLOCATABLE :: perturb
        REAL(8), dimension(:,:,:,:,:,:,:), ALLOCATABLE :: mySimplex
        REAL(8), dimension(:,:,:,:,:,:,:), ALLOCATABLE :: mySimplex_rent
        LOGICAL, DIMENSION(3,2) :: valfuncpol

        pricelen=MAX(SIZE(price,1)-1,1)

        ! %< Allocation of automatic arrays for OpenMP (skip this)
        ALLOCATE(achoiceadjust(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 Dchoiceadjust(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 cchoiceadjust(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 achoicenoadjust(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 Dchoicenoadjust(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 cchoicenoadjust(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 achoicerent(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 Dchoicerent(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 cchoicerent(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 Vnoadjust(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 Vadjust(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 Vrent(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen))
        ALLOCATE(pstartadjust(vars+1, vars, zgridsize, 2, Dgridsize, agridsize, pricelen),&
                 ystartadjust(vars+1, zgridsize, 2, Dgridsize, agridsize, pricelen),&
                 state(9, zgridsize, 2, Dgridsize, agridsize, pricelen),&
                 incomeHolder(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 assetHolder(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 wealthHolder(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 costHolder(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 costDenom(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 transferHolder(2, zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen),&
                 myUsercost(zgridsize, 2, Dgridsize, agridsize, Tdie, pricelen))
        ALLOCATE(ax(zgridsize, 2, Dgridsize, agridsize, pricelen),&
                 bx(zgridsize, 2, Dgridsize, agridsize, pricelen),&
                 cx(zgridsize, 2, Dgridsize, agridsize, pricelen),&
                 adjust(zgridsize, 2, Dgridsize, agridsize, pricelen),&
                 perturb(8, zgridsize, Dgridsize, agridsize))
        ALLOCATE(mySimplex(vars, vars+1, zgridsize, 2, Dgridsize, agridsize, pricelen))
        ALLOCATE(mySimplex_rent(vars, vars+1, zgridsize, 2, Dgridsize, agridsize,pricelen))

        adjust=0
        ! %>

        ! %< PARAMETERS FOR SCOPE OF DP EXECUTION
        ! Sets policies to be activated for the specific DP program
        ! Different string arguments turn on different kinds of policies.
        valfuncpol = .FALSE.
        if (noshock .EQV. .FALSE.) valfuncpol(1,1) = .TRUE.
        if (expected == 'continuing') then
            valfuncpol(1,2) = .TRUE.
        else if (expected == 'first-time') then
            valfuncpol(2,2) = .TRUE.
        else if (expected == 'repeat') then
            valfuncpol(3,2) = .TRUE.
        end if
        ! write(*,*) valfuncpol ! Testing anticipated policy markers

        ! Limits the state spaces needed for calculation
        ! (E.g. arrays holding FTHB policy response only need to be run
        !  over cases with zero durables held)
        if (present(state_end)) then
            d_end = state_end(1)
            t_end = state_end(2)
        else
            d_end = Dgridsize
            t_end = TDie
        end if ! %>

        ! PLAN FOR THE STATE SPACES

        ! i=current perm. income shock.
        ! i2=Unemployment or not
        ! a=current assets.
        ! d=current housing services.
        ! l=number of transition periods.
        ! t= time.

        ! Pricelen==1 here either means
        ! 1) calculating equilibrium price for the initial steady state:
        ! 2) calculating the PE dynamic programming scenario.
        do l=pricelen, 1,-1
            if (pricelen > 1) then
                write(*,*) "Transition period:", l
                write(0,*) "Transition period:", l
            end if

            ! Anticipated price changes in durable turned off in PE simulation
            if (SIZE(price,1) == 1) then
                hpdelta = 0
            else
                hpdelta = price(l+1) - price(l)
            end if

        do t=t_end, 1,-1

            CALL random_number(perturb)
            ! No risk of unemployment in retirement
            ! In addition, unemployment can just be turned off
            if (t > Tretire .OR. empprob(Tretire) == 1.0) then
                ue_start = 2
            else
                ue_start = 1
            end if
!$OMP PARALLEL
!$OMP DO
            do j=1, d_end
            do k=agridsize,1,-1
                do i2=ue_start,2
                    do i=1, zgridsize

            ! %< State vector
            state(1, i, i2, j, k, l)=i
            state(2, i, i2, j, k, l)=i2
            state(3, i, i2, j, k, l)=Dnodes(j)
            state(4, i, i2, j, k, l)=anodes(k)
            state(5, i, i2, j, k, l)=price(l)
            if (pricelen == 1) then
                state(6, i, i2, j, k, l)=l-1
            else
                state(6, i, i2, j, k, l)=l+1
            end if
            state(7, i, i2, j, k, l)=t
            state(8, i, i2, j, k, l)=balancer(l)
            state(9, i, i2, j, k, l)=0.0  !placeholder
            ! %>

            ! %< Functions of state space (simplifies code later)
            incomeHolder(i,i2,j,k,t,l) = incomeCall(state(:, i, i2, j, k, l)) + polLevel(1) - polLevel(4) !currentincome in valfuncadjust
            ! This is the part of assetReturns that excludes future interest
            assetHolder(i,i2,j,k,t,l) = anodes(k)-&
                (1-theta)*((1.0-scrapped)*((1.0-delta)/(1.0+rborrow))*Dnodes(j)*exp(price(l)) -&
                scrapped*scrapvalue)
            if (assetHolder(i,i2,j,k,t,l) < 0) then
                assetHolder(i,i2,j,k,t,l) = (1+rborrow)*assetHolder(i,i2,j,k,t,l)&
                + MIN(-himdflag*rborrow*(anodes(k)-(1-theta)*Dnodes(j)*exp(price(l))), inctax(i,t))
                myUsercost(i,i2,j,k,t,l) = (1.0-delta-dtau)/(1.0+rborrow)  
            else
                assetHolder(i,i2,j,k,t,l) = (1+r)*assetHolder(i,i2,j,k,t,l)
                myUsercost(i,i2,j,k,t,l) = (1.0-delta-dtau)/(1.0+r)
            end if
                wealthHolder(i,i2,j,k,t,l)=wealthCall(state(3:5, i, i2, j, k, l),&
                                                      incomeHolder(i,i2,j,k,t,l),&
                                                      assetHolder(i,i2,j,k,t,l))
            if (t > Tretire) then
                ! Because retirees don't pay taxes
                wealthHolder(i,i2,j,k,t,l)=wealthHolder(i,i2,j,k,t,l)+inctax(i,t)
            end if

            ! %< Optimal (a', D') assuming entry into rental sector

            ! Set up initial 2-simplex and find local minimum
            ! The initial simplex is identical to the adjustment initialization, since
            ! selling the house has costs too and someone still in rental has D=0
            
            call get_simplex_rent(mySimplex_rent(:, :, i, i2, j, k, l), state(:, i, i2, j, k, l), incomeHolder(i,i2,j,k,t,l), assetHolder(i,i2,j,k,t,l))
            do p1=1,3
                do p2=1,2
                    pstartadjust(p1, p2, i, i2, j, k, l)=mySimplex_rent(p2, p1, i, i2, j, k, l)
                end do
                ystartadjust(p1, i, i2, j, k, l)=valfuncrent(&
                    pstartadjust(p1,:,i,i2,j,k,l), state(:, i, i2, j, k, l),&
                    valfuncpol(2,:))
            end do

            call amoeba(state(:, i, i2, j, k, l), pstartadjust(:,:,i, i2, j, k, l), &
                        ystartadjust(:,i,i2,j,k,l), ftol, valfuncpol(2,:), &
                        valfuncrent)
            ! p1 dimension doesn't matter because vertices should converge
            Vrent(i, i2, j, k, t, l)=ystartadjust(1, i, i2, j, k, l)
            achoicerent(i, i2, j, k, t, l)=min(max(pstartadjust(1, 1, i, i2, j, k, l),borrowconstraint),amax)
            Dchoicerent(i, i2, j, k, t, l)=min(max(pstartadjust(1, 2, i, i2, j, k, l),0.0),Dmax)

            ! Having found Dchoicerent we assign it to be fed into adjust
            ! function
            state(9, i, i2, j, k, l)=Dchoicerent(i, i2, j, k, t, l)
            ! %>

            ! %< Optimal (a', D') assuming adjustment / paid trans. costs

            ! Set up initial 2-simplex and find local minimum
                !create fn that takes init guess for policy fn. check if downpayment for min house size > income? if F, dont call
                !amoeba (this is outside loop). if T, find optimal a', D'  
            call get_simplex(mySimplex(:, :, i, i2, j, k, l), state(:, i, i2, j, k, l), incomeHolder(i,i2,j,k,t,l), assetHolder(i,i2,j,k,t,l), myUsercost(i,i2,j,k,t,l), valfuncpol(1,:))
            
            
            do p1=1,3
                do p2=1,2
                    pstartadjust(p1, p2, i, i2, j, k, l)=mySimplex(p2, p1, i, i2, j, k, l)
            !        pstartadjust(p1, p2, i, i2, j, k, l)=amoebaGrid(p2, p1, 1) * max(&
            !            wealthHolder(i,i2,j,k,t,l), (1+F/2.0)*Dmin*exp(price(l)) + 1e-2)
                    !should be this
                    !pstartadjust(p1, p2, i, i2, j, k, l)=amoebaGrid(p2, p1, 1) * wealthHolder(i,i2,j,k,t,l)
                end do
                ystartadjust(p1, i, i2, j, k, l)=valfuncadjust(&
                    pstartadjust(p1,:,i,i2,j,k,l), state(:, i, i2, j, k, l),&
                    valfuncpol(1,:))
            end do

            call amoeba(state(:, i, i2, j, k, l), pstartadjust(:,:,i, i2, j, k, l), &
                        ystartadjust(:,i,i2,j,k,l), ftol, valfuncpol(1,:), &
                        valfuncadjust)
            ! p1 dimension doesn't matter because vertices should converge

            Vadjust(i, i2, j, k, t,l)=ystartadjust(1, i, i2, j, k, l)
            achoiceadjust(i, i2, j, k, t,l)=min(max(pstartadjust(1, 1, i, i2, j, k, l),borrowconstraint),amax)
            Dchoiceadjust(i, i2, j, k, t,l)=min(max(pstartadjust(1, 2, i, i2, j, k, l),0.0),Dmax)
            ! %>

            ! %< Optimal (a') assuming no adjustment (D' = D)
            ax(i, i2, j, k, l)=0.0
            if (costlyequity) then
                if (anodes(k)>=(1.0-dtau)*(1-theta)*Dnodes(j)*exp(price(l))) then
                    bx(i, i2, j, k, l)=(1-theta)*Dnodes(j)*exp(price(l))
                else
                    bx(i, i2, j, k, l)=anodes(k)
                end if
            else
            bx(i, i2, j, k, l)=borrowconstraint
            end if
            cx(i, i2, j, k, l)=(1-1e-7)*amax

            state(3, i, i2, j, k, l)=j ! Need this to retrieve EV
            ! %>

            ! %< Optimization in one variable, hence use Brent's method for
            ! optimization without derivatives
            Vnoadjust(i, i2, j, k, t, l)=brentnew(&
                ax(i, i2, j, k, l), bx(i, i2, j, k, l), cx(i, i2, j, k, l), valfuncpol(3,2), valfuncnoadjust,&
                ftol, state(:, i, i2, j, k, l), achoicenoadjust(i, i2, j, k, t, l))
            Dchoicenoadjust(i, i2, j, k, t, l)=Dnodes(j)
            ! Idea is that aprime_adj in function must've eliminated negative
            if (achoicenoadjust(i,i2,j,k,t,l) < 0) achoicenoadjust(i,i2,j,k,t,l) = borrowconstraint
            !>

            ! %>

            ! %< Computing cchoiceadjust, cchoicenoadjust and cchoicerent
            cchoicerent(i,i2,j,k,t,l)=presCons(&
                wealthHolder(i,i2,j,k,t,l) - achoicerent(i,i2,j,k,t,l), 0.0,&
                price(l), 0.0, renttransfer_internal - &
                (r_rental(l)*Dchoicerent(i,i2,j,k,t,l)*exp(rentelasticity*price(l))))

            cchoicenoadjust(i,i2,j,k,t,l)=presCons(&
                wealthHolder(i,i2,j,k,t,l) -  achoicenoadjust(i,i2,j,k,t,l),&
                Dnodes(j), price(l), theta, 0.0)

            ! If adjusting, we need to include transaction costs
            ! TODO: get downpay in this
            costDenom(i,i2,j,k,t,l) = ((Dchoiceadjust(i,i2,j,k,t,l)+Dchoicerent(i,i2,j,k,t,l)))
            costholder(i,i2,j,k,t,l) = &
                F2*(&
                (Dchoiceadjust(i,i2,j,k,t,l)-0.0)/costDenom(i,i2,j,k,t,l))**(2.0)*&
                costDenom(i,i2,j,k,t,l) + F*Dnodes(j)*(1-maint*delta-dtau)*exp(price(l))

            if (valfuncpol(1,1)) then
                ! The first term of transferholder is if the policy applies on BC
                ! The second is if it's applied as a price reduction on down
                ! payment only
                transferholder(:,i,i2,j,k,t,l) = eta_transfer*adjtransfer(i,t)*(/ &
                    (1.0-downflag)*(1.0-discountflag),&
                    (1.0-discountflag)/(exp(price(l))*Dchoiceadjust(i,i2,j,k,t,l)) /)
                if (pctageflag(i,t)) then
                    transferholder(:,i,i2,j,k,t,l) = transferholder(:,i,i2,j,k,t,l)*&
                    &(exp(price(l))*Dchoiceadjust(i,i2,j,k,t,l)) 
                end if
            end if
           cchoiceadjust(i,i2,j,k,t,l) = presCons(&
               wealthHolder(i,i2,j,k,t,l) - achoiceadjust(i,i2,j,k,t,l),&
               Dchoiceadjust(i,i2,j,k,t,l), price(l),&
               (1.0-polLevel(3))*theta * 1.0 * (1.0 - transferholder(2,i,i2,j,k,t,l)/((1-polLevel(3)*theta))), &
               - costholder(i,i2,j,k,t,l) + transferholder(1,i,i2,j,k,t,l)) 
            ! %>

            ! %< Finding the best housing decision and recording it to array 
            rentchoice(i, i2, j, k, t, l)=Dchoicerent(i, i2, j, k, t, l)
            !this chunk is fine, empirically verified. but what's polLevel(3)?
            if (Vadjust(i, i2, j, k, t, l)<Vnoadjust(i, i2, j, k, t, l) .and. &
                Vadjust(i, i2, j, k, t, l)<Vrent(i, i2, j, k, t, l)) then  ! since V = - V from minimization
                    achoice(i, i2, j, k, t, l)=achoiceadjust(i, i2, j, k, t, l)
                    if (polLevel(3) > 0) then
                        achoice(i, i2, j, k, t, l) = (1-delta+maint*delta)*Dchoiceadjust(i, i2, j, k, t, l)*(polLevel(3)*theta)*exp(price(l))
                    end if
                    Dchoice(i, i2, j, k, t, l)=Dchoiceadjust(i, i2, j, k, t, l)
                    cchoice(i, i2, j, k, t, l)=cchoiceadjust(i, i2, j, k, t, l)
                    choiceindicator(i, i2, j, k, t, l)=1
            elseif (Vnoadjust(i, i2, j, k, t, l)<Vadjust(i, i2, j, k, t, l) .and. &
                    Vnoadjust(i, i2, j, k, t, l)<Vrent(i, i2, j, k, t, l)) then
                    achoice(i, i2, j, k, t, l)=achoicenoadjust(i, i2, j, k, t, l)
                    Dchoice(i, i2, j, k, t, l)=Dchoicenoadjust(i, i2, j, k, t, l)
                    cchoice(i, i2, j, k, t, l)=cchoicenoadjust(i, i2, j, k, t, l)
                    choiceindicator(i, i2, j, k, t, l)=2
            else
                    achoice(i, i2, j, k, t, l)=achoicerent(i, i2, j, k, t, l)
                    Dchoice(i, i2, j, k, t, l)=Dchoicerent(i, i2, j, k, t, l)
                    cchoice(i, i2, j, k, t, l)=cchoicerent(i, i2, j, k, t, l)
                    choiceindicator(i, i2, j, k, t, l)=3
            end if

            ! Forced move: either adjust or rent. Here adjust is more valuable
            !where does this come in?
            if (movprob > 0.0 .AND. Vadjust(i, i2, j, k, t, l)<Vrent(i, i2, j, k, t, l)) then
                    achoiceMov(i, i2, j, k, t, l)=achoiceadjust(i, i2, j, k, t, l)
                    DchoiceMov(i, i2, j, k, t, l)=Dchoiceadjust(i, i2, j, k, t, l)
                    cchoiceMov(i, i2, j, k, t, l)=cchoiceadjust(i, i2, j, k, t, l)
                    choiceindicatorMov(i, i2, j, k, t, l)=1
            else if (movprob > 0.0) then
                    achoiceMov(i, i2, j, k, t, l)=achoicerent(i, i2, j, k, t, l)
                    DchoiceMov(i, i2, j, k, t, l)=Dchoicerent(i, i2, j, k, t, l)
                    cchoiceMov(i, i2, j, k, t, l)=cchoicerent(i, i2, j, k, t, l)
                    choiceindicatorMov(i, i2, j, k, t, l)=3
            end if
             
            ! Forced move for renters *if* their BC not violated by owning with
            ! a min house size present. If BC violated, you still rent.
            if (movprobR > 0.0 .AND. Vadjust(i, i2, j, k, t, l) < 1e4 ) then !.and. anodes(k) > 5e-3) then
                    achoiceMovR(i, i2, j, k, t, l)=achoiceadjust(i, i2, j, k, t, l)
                    DchoiceMovR(i, i2, j, k, t, l)=Dchoiceadjust(i, i2, j, k, t, l)
                    cchoiceMovR(i, i2, j, k, t, l)=cchoiceadjust(i, i2, j, k, t, l)
                    choiceindicatorMovR(i, i2, j, k, t, l)=1
                    ! The coef in front of Vadjust represents how "expected" the
                    ! renter moving shock is (0 = MIT)
                    EVMovR(i,i2,j,k, t,l)= -(1.00*Vadjust(i, i2, j, k, t, l) + (1.0-1.00)*Vrent(i, i2, j, k, t, l))
            else if (movprobR > 0.0) then
                    achoiceMovR(i, i2, j, k, t, l)=achoicerent(i, i2, j, k, t, l)
                    DchoiceMovR(i, i2, j, k, t, l)=Dchoicerent(i, i2, j, k, t, l)
                    cchoiceMovR(i, i2, j, k, t, l)=cchoicerent(i, i2, j, k, t, l)
                    choiceindicatorMovR(i, i2, j, k, t, l)=3
                    EVMovR(i,i2,j,k, t,l)= -Vrent(i, i2, j, k, t, l)
            end if


            ! %< Policy choice printout for certain states
            if ((i==zgridsize/2) .AND. (((scrapped == 0) .AND. (j==20) .AND. (k==6)) .OR. &
                ((scrapped == 1) .AND. (j==6) .AND. (k==6)) .OR. ((j==1) .AND. (k==30)) &
                   ) .AND. show_pol) THEN
                !((t == t_end) .OR. (mod(t, 15)==1)) .AND. show_pol) THEN
                write(0,*) "Check optimal policy in period", t
                write(0,*) "Sample for D-q state:", dnodes(j), anodes(k)
                write(0,'(4A20)') "  |", "  With adjustment  |", "  With rental  |", "  Without adjustment"
                write(0,'(A20,3F20.10)') "Vol. equity  |", &
                    achoiceadjust(i,i2,j,k,t,l), achoicerent(i,i2,j,k,t,l), achoicenoadjust(i,i2,j,k,t,l)
                write(0,'(A20,3F20.10)') "Durables  |", &
                    exp(price(l))*Dchoiceadjust(i,i2,j,k,t,l), exp(price(l))*Dchoicerent(i,i2,j,k,t,l), exp(price(l))*Dchoicenoadjust(i,i2,j,k,t,l)
                write(0,'(A20,3F20.10)') "Consumption  |", &
                    cchoiceadjust(i,i2,j,k,t,l), cchoicerent(i,i2,j,k,t,l), cchoicenoadjust(i,i2,j,k,t,l)
                write(0,'(A20,3F20.10)') "Expected Value  |", &
                    -Vadjust(i,i2,j,k,t,l), -Vrent(i,i2,j,k,t,l), -Vnoadjust(i,i2,j,k,t,l)
                write(0,*) "// Choice  :", choiceindicator(i,i2,j,k,t,l)
                write(0,*) "// EV  :", -min(Vrent(i,i2,j,k, t,l),min(Vnoadjust(i,i2,j,k, t,l), Vadjust(i,i2,j,k, t,l)))
            end if ! %>

                    end do
                end do
            end do
        end do
!$OMP END DO
!$OMP END PARALLEL
            ! Broadcast value functions for next age cohort iteration 
            Vnoadjust(:,:,:,:, t,l)=-Vnoadjust(:,:,:,:, t,l)
            Vadjust(:,:,:,:, t,l)=-Vadjust(:,:,:,:, t,l)
            Vrent(:,:,:,:, t,l)=-Vrent(:,:,:,:, t,l)
            if (repl) THEN
                EV(:,:,:,:, t,l)=max(Vrent(:,:,:,:, t,l),max(&
                                   Vnoadjust(:,:,:,:, t,l), Vadjust(:,:,:,:, t,l)))
                if (movprob > 0.0) EVMov(:,:,:,:, t,l)=max(Vrent(:,:,:,:, t,l), Vadjust(:,:,:,:, t,l))
                if (movprobR > 0.0) then
                    !EVMovR(:,:,:,:, t,l)= Vadjust(:,:,:,:, t,l)
                    !where(choiceindicatormovR(:,:,:,:, t,l) == 3) EVMovR(:,:,:,:, t,l) = Vrent(:,:,:,:, t,l)
                end if
            end if
    end do
    end do

    if (noshock .EQV. .TRUE.) call output_vfuncs( (/ 3, 8, 37, 38, 39, 40, 47 /), &
        l+1, EV, Vadjust, Vrent, Vnoadjust, achoice, dchoiceadjust, dchoice, cchoiceadjust, cchoice)

    DEALLOCATE(ax, bx, cx, adjust)

    end subroutine solveworkingproblem ! %>

    subroutine plotpolicy_constantwealth(p, aArray, DArray, cArray,choiceArray)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Legacy module from BGLV code that plots the policy response over
        !    (a, D) values holding one's net wealth constant. Not used much nowadays.
        !    Accepts arbitrary policy arrays to test out different policy responses.
        !
        !    Modified: 07/06/2019
        !
        !    PARAMETERS
        !
        !    p: Integer. Indicates the policy time period to look at, with
        !    the corresponding aggregate house price level.
        !
        !    *Array: Policy arrays and adjustment indicator considered with this
        !    function.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE OMP_LIB
        IMPLICIT NONE
        INTEGER :: i, j, k, t
        INTEGER, INTENT(IN) :: p
        REAL(8), dimension(:,:,:,:,:,:), INTENT(IN) :: aArray, DArray, cArray, choiceArray
        REAL(8) :: wealth, h, q, achoicelin, Dchoicelin, cchoicelin, rentallin, choicelin, welfare

        OPEN (UNIT=25, FILE="policy_constantwealth.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        t=3
        do k=1, zgridsize
        do i=1, 50
            wealth=0+.1*i
            do j=1, 50
                h=0+(j-1)*.1
                q=wealth-theta*h
                ! State space vector follows same order as in solveworkingproblem
                ! NOTE: p here is not the value of the price state, but its coordinate!
                call pol_linworking((/ REAL(k), 2.0,  h, q, REAL(t), REAL(p) /),&
                    aArray,Darray,cArray,EV,choiceArray,&
                    achoicelin, Dchoicelin, cchoicelin, rentallin, choicelin, welfare)
                write(25, '(I6.2,99F16.6)') k, wealth, h, Dchoicelin, cchoicelin, achoicelin
            end do
        end do
        end do

    end subroutine plotpolicy_constantwealth ! %>

    ! TODO: Wrap a loop that accepts arbitrary values of t.
    subroutine output_vfuncs(agelist, period, EV, Vadjust, Vrent, Vnoadjust, achoice, dadjchoice, dchoice, cadjchoice, cchoice)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Legacy module from BGLV code that plots all policy function choices
        !    over certain time periods. Outputs fixed-format spreadsheet.
        !
        !    Modified: 07/06/2019
        !
        !    PARAMETERS
        !
        !    period: Integer. Indicates the policy time period to look at, with
        !    the corresponding aggregate house price level. Corresponds to "p"
        !    in other functions.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IMPLICIT NONE
        INTEGER, dimension(:), INTENT(IN) :: agelist
        INTEGER, INTENT(IN) :: period
        REAL(8), dimension(:,:,:,:,:,:), INTENT(IN) :: EV, Vadjust, Vrent, Vnoadjust, achoice, dadjchoice, dchoice, cadjchoice, cchoice
        INTEGER, dimension(1) :: ageIter
        INTEGER :: i, i2, j, k, t, s, l, t2, u, u2
        CHARACTER(LEN=8) :: fname

        !OPEN (UNIT=21, FILE="vfunc1.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        !OPEN (UNIT=22, FILE="vfunc2.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        !OPEN (UNIT=23, FILE="vfunc3.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        !OPEN (UNIT=24, FILE="vfunc4.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        !OPEN (UNIT=25, FILE="vfunc5.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        !OPEN (UNIT=26, FILE="vfunc6.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

        s=1
        ageIter = shape(agelist)
        do while (s <= ageIter)
        
        t = agelist(s)
        write(fname, '(A6, I0.3)') "vfunc_" // t
        OPEN (UNIT=21, FILE=fname, STATUS="NEW", ACTION="WRITE", POSITION="REWIND")
        do i=1, zgridsize, 2
            do i2=2, 2
              do j=1, Dgridsize/2, 4
                do k=1, agridsize
                    do l=period,period
                        write(21, '(14F40.12)') hpnodes(l), znodes(i), Dnodes(j), anodes(k), &
                               EV(i, i2, j, k, t, l), Vadjust(i, i2, j, k, t, l), &
                               Vrent(i, i2, j, k, t, l), Vnoadjust(i, i2, j, k, t, l), &
                               achoice(i, i2, j, k, t, l), dadjchoice(i, i2, j, k, t, l), &
                               dchoice(i, i2, j, k, t, l), cadjchoice(i, i2, j, k, t, l),  &
                               cchoice(i, i2, j, k, t, l),  income(i, t)
                    end do
                end do
              end do
            end do
         end do
         s = s + 1
         CLOSE(21)
         end do


      OPEN (UNIT=22, FILE="vfuncLC.txt", STATUS="NEW", ACTION="WRITE", POSITION="REWIND")
      do i=1, zgridsize, 2
          do i2=2, 2
              do k=6, agridsize, 6
                  do l=period,period
                  do t=1, Tdie
                  j = 1  ! Renters
                      write(22, '(14F40.12)') hpnodes(l), znodes(i), anodes(k), real(t), EV(i, i2, j, k, t, l), Vadjust(i, i2, j, k, t, l), Vrent(i, i2, j, k, t, l),Vnoadjust(i, i2, j, k, t, l),  achoice(i, i2, j, k, t, l), dadjchoice(i, i2, j, k, t, l), dchoice(i, i2, j, k, t, l), cadjchoice(i, i2, j, k, t, l),  cchoice(i, i2, j, k, t, l),  income(i, t) !, mpc_pol(i, j, k, t, l)
                      end do
                  end do
              end do
          end do
    end do
    CLOSE(22)




    end subroutine output_vfuncs ! %>

end module lifecycle_solveDP
