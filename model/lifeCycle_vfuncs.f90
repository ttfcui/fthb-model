module lifecycle_vfuncs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE VFUNCS contain the main objective functions and constraints in the
! dynamic programming problem, in a sufficient general form.
! It also contains a subroutine that compares the utility of adjusting without
! the policy or adjusting with the policy (e.g. in a 2-period policy, to adjust
! now or adjust later)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    USE share
    USE lifeCycle_algs
    IMPLICIT NONE

    CONTAINS

    REAL(8) FUNCTION asset_shift(aprime, atrans, shifts)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Small function that crunches the vector of future income
        !    possible with temporary, IID shocks, values that are embedded into
        !    the asset state (so no need for one more state).
        !
        !    Modified: 06/22/2017
        !
        !    PARAMETERS
        !
        !    aprime: Future asset policy considered in DP program.
        !
        !    atrans: Actual financial assets next period (since "a"
        !    is actually voluntary equity). Needed to determine which
        !    interest rate applies on the asset. 
        !
        !    shifts: The amount of perturbation to the chosen asset level
        !    next period, coming from income shocks or temp. policy receipt.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: aprime, atrans, shifts

        ! Embedding temporary shocks into asset state space
            if (atrans < 0 .AND. atrans + shifts < 0) then
                    asset_shift = min(max(aprime+(shifts)/(1+rborrow),0.0),&
                                      anodes(agridsize))
            else if (atrans < 0) then
                    asset_shift = min(max(aprime+(shifts)/(1+r)+&
                                      (rborrow-r)*atrans,0.0),anodes(agridsize))
            else if (atrans >= 0 .AND. atrans + shifts >= 0) then
                    asset_shift = min(max(aprime+(shifts)/(1+r),0.0),&
                                      anodes(agridsize))
            else
                    asset_shift = min(max(aprime+(shifts)/(1+rborrow)+&
                                      (r-rborrow)*atrans,0.0),anodes(agridsize))
            end if


    END FUNCTION ! %>

    REAL(8) FUNCTION incomeCall(state)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Function wraps income calculations for each
        !    household (net taxes, unemployment benefits, etc)
        ! 
        !    Modified: 06/07/2018
        !
        !    PARAMETERS
        !
        !    state: The entire state variable fed into each
        !    value function call (though only 4 used)
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), DIMENSION(9), INTENT(IN), TARGET :: state
        REAL(8), POINTER :: zindex, t, balholder, unemp

        zindex => state(1)
        unemp => state(2)
        t => state(7)
        balholder => state(8)

        incomeCall = (unemp - 1.0)*(income(zindex, t) - inctax(zindex, t)) +&
                     (unemp - 2.0)*(unemptrans) + balholder

    END FUNCTION ! %>

    REAL(8) FUNCTION wealthCall(state, currentincome, assetholder)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Function wraps income calculations for each
        !    household (net taxes, unemployment benefits, etc)
        ! 
        !    Modified: 06/07/2018
        !
        !    PARAMETERS
        !
        !    state: 3/8 states encoded in state variable
        !    (housing, assets, price)
        !
        !    currentincome: HH flow income calculated by
        !    incomeCall.
        !
        !    assetholder: Calculation of liquid assets
        !    held based off of vol. equity and housing
        !    state variables.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), DIMENSION(3), INTENT(IN), TARGET :: state
        REAL(8), INTENT(IN) :: currentincome, assetholder
        REAL(8), POINTER :: a, D, hpholder

        D => state(1)
        a => state(2)
        hpholder => state(3)

        wealthCall = currentincome + D*(1-maint*delta-dtau)*exp(hpholder)&
            + assetholder

    END FUNCTION ! %>

    SUBROUTINE assetReturns(state, aprime, Dprime, presfact, dfact, &
                            thetaholder, inctax, assetholder, futliquid)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        ! 
        !    
 
        IMPLICIT NONE
        REAL(8), DIMENSION(3), INTENT(IN), TARGET :: state
        REAL(8), INTENT(IN) :: aprime, Dprime, presfact, dfact, thetaholder, inctax
        REAL(8), POINTER :: a, D, hpholder
        REAL(8), DIMENSION(2) :: finwealth
        REAL(8), INTENT(OUT) :: assetholder, futliquid

        D => state(1)
        a => state(2)
        hpholder => state(3)

        ! First element is current financial assets, 2nd is next period
        finwealth = (/ a + (-1.0+theta)*((1.0-scrapped)*D*presfact*exp(hpholder)&
                         + scrapped*-scrapvalue),&
                       aprime + (1.0-scrapped)*Dprime*(-1+theta+delta-maint*delta+dfact)&
                               &*exp(hpholder+hpdelta) + (-1.0+theta)*scrapped*-scrapvalue /)

        if (finwealth(1) < 0.0) then
                assetholder = (1+rborrow)*finwealth(1)&
                              + himdflag*MAX(rborrow*finwealth(1), inctax)
        else
                assetholder = (1+r)*finwealth(1)
        end if

        futliquid = finwealth(2)

    END SUBROUTINE  ! %>

    REAL(8) FUNCTION presCons(netwealth, Dprime, hpholder, downfact, othcosts)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        ! delete all the pricefact params that occur from the presCons calls 
        ! 
        !    

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: netwealth, Dprime, hpholder, downfact, othcosts

        presCons = netwealth - downfact*Dprime*&
                  &exp(hpholder) + othcosts


    END FUNCTION ! %>

    REAL(8) FUNCTION presVal(consumption, Dprime, dfact, utilfact)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        ! 
        !    

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: consumption, Dprime, dfact, utilfact

        ! TODO: what's the deal with the "+0"? is that an income
        ! effect/Stone-Geary kind of thing?
        presVal = (consumption+0)**elasticity2*(&
                  (1.0+dfact)*Dprime+0)**(1-elasticity2)
        if (elasticity .ne. 1) then
            presVal = utilfact*(presVal)**(1-elasticity)/(1-elasticity)
        else
            presVal = utilfact*log(presVal)
        end if

    END FUNCTION ! %>

    FUNCTION EV_linpol(awgt, Dwgt, zmin, zmax, EVinput)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Function linearly interpolates future
        !    expectation terms for (a, D) values off the grid.
        ! 
        !    Modified: 06/22/2017
        !
        !    PARAMETERS
        !
        !    *wgt: Weight on lower a/D term on the grid space.
        !
        !    zmin/zmax: Given the current AR(1) income state,
        !    marks the states into which agent can transfer
        !    with nontrivial probability (filtered out in
        !    module lifeCycle_calibrate)
        !
        !    EVinput: Expectation array fed into current
        !    DP problem (e.g. array where policy still active)
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: awgt, dwgt
        INTEGER, INTENT(IN) :: zmin, zmax
        REAL(8), DIMENSION(:, :, :), INTENT(IN) :: EVinput
        REAL(8), DIMENSION(zmax-zmin+1) :: EV_linpol

        EV_linpol = 0.0

         ! in the case of non-adjust function
        if (SIZE(EVinput, 2) < 2) then
            EV_linpol(:)=awgt*EVinput(zmin:zmax,1,1)+&
                         (1-awgt)*EVinput(zmin:zmax,1,2)
        else
            EV_linpol(:)=awgt*dwgt*EVinput(zmin:zmax,1,1)+&
                        (1-awgt)*dwgt*EVinput(zmin:zmax,1,2)+&
                        awgt*(1-dwgt)*EVinput(zmin:zmax,2,1)+&
                        (1-awgt)*(1-dwgt)*EVinput(zmin:zmax,2,2)
        end if
    END FUNCTION ! %>
        
    REAL(8) FUNCTION continuation(discount, prob, exval)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Small function that takes an array of
        !    probability and EV values and performs the
        !    dot product, while also discounting the future.
        ! 
        !    Modified: 06/22/2017
        !
        !    PARAMETERS
        !
        !    discount: Real number reflecting current discount factor.
        !
        !    Prob: Vector of probabilities for possible future
        !    states (must be less than 1, should sum to 1 too)
        !
        !    exval: Vector of expected values corresponding
        !    to possible future states.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: discount
        REAL(8), DIMENSION(:), INTENT(IN) :: prob, exval

        ! Tiny function that outputs discounted continuation
        ! value given two future states and outcomes.
        if (MAXVAL(prob) > 1 .AND. SIZE(prob) /= SIZE(exval)) call nrerror("Probability undefined")
        continuation = discount*DOT_PRODUCT(prob, exval)

    END FUNCTION ! %>
    REAL(8) FUNCTION temp_expect_loop(EVmat, aval, tmpnodes, probtmp, weightDprimel, &
                                      pwealth, zindex, zmin, zmax, transfers)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        ! 
        !    

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: zindex, zmin, zmax
        REAL(8), DIMENSION(:, :, :), INTENT(IN), TARGET :: EVmat
        REAL(8), DIMENSION(:), INTENT(IN) :: tmpnodes, Probtmp
        REAL(8), INTENT(IN) :: aval, weightDprimel, pwealth
        REAL(8), OPTIONAL :: transfers
        REAL(8), DIMENSION(:), ALLOCATABLE :: aprime_shift, temp_expect, weightaprimel
        REAL(8), DIMENSION(:,:,:), POINTER :: EVexpect
        INTEGER, DIMENSION(:), ALLOCATABLE :: aprimel_vec, aprimeh_vec
        INTEGER :: i, tmpgridsize

        tmpgridsize = size(tmpnodes)
        allocate(aprime_shift(tmpgridsize), temp_expect(tmpgridsize))

        if (present(transfers)) then
            do i=1, tmpgridsize
                aprime_shift(i) = asset_shift(aval, pwealth, &
                                              transfers + tmpnodes(i))
            end do
        else
            do i=1, tmpgridsize
                aprime_shift(i) = asset_shift(aval, pwealth,tmpnodes(i))
            end do
        end if
 
        call weightPrimevec(aprime_shift, anodes, weightaprimel, &
             aprimel_vec, aprimeh_vec)

        do i = 1, tmpgridsize
            EVexpect => EVmat(:, :, aprimel_vec(i):aprimeh_vec(i))
            if (zindex == zmin .AND. zindex == zmax) then
                temp_expect(i:i) = EV_linpol(weightaprimel(i), weightDprimel, zmin, zmax, EVexpect)
            else
            temp_expect(i) = DOT_PRODUCT(Probz(zindex,zmin:zmax), &
                EV_linpol(weightaprimel(i), weightDprimel, zmin, zmax, EVexpect))
            end if
        end do
           
        temp_expect_loop = DOT_PRODUCT(Probtmp, temp_expect)

    END FUNCTION ! %>
 
    REAL(8) FUNCTION expectation(zindex, Dval, aval, tval, hpval, price, uval, pwealth, Down,&
                                 expecting, rental, adjust, transfers)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Function calculates the expectation term after considering
        !    all possible shocks: AR(1) innovations, temporary income shocks,
        !    moving shock. Also takes into account if policy still exists
        !    in the future.
        ! 
        !    Modified: 08/09/2019
        !
        !    PARAMETERS
        !
        !    zindex/uval: Index of state space on the income process.
        !
        !    *val: Asset, durable and remaining state values in the
        !    FOLLOWING period.
        !
        !    pwealth: 2-vector of actual financial assets, net of household
        !    borrowing (so not voluntary equity)
        !
        !    expecting: Boolean. If TRUE, assumes next array stems from
        !    array where the policy is still active. 
        !
        !    rental: Boolean. If TRUE, understands function is used for rental
        !    problem: housing next period "Dprime" is guaranteed to be zero and
        !    irrelevant (hence interpolation is only over the assets dimension.)
        !    
        !    adjust: Boolean. If TRUE, understands function is used for
        !    adjustment to a new owned home. Important because this activates
        !    stochastic transaction costs.
        !
        !    transfers: Real. If nonzero, interpret as the temporary policy's face value.    
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        INTEGER, INTENT(IN), TARGET :: zindex
        REAL(8), INTENT(IN) :: aval, Dval, tval, hpval, price, uval, pwealth, Down
        LOGICAL, INTENT(IN) :: expecting, rental, adjust
        REAL(8), OPTIONAL :: transfers
        REAL(8), DIMENSION(:,:,:), POINTER :: EVmat, EVmatM, EVmatMR, EVmatU, EVmatD
        REAL(8), DIMENSION(:), POINTER :: tmpgrid, probtmp
        REAL(8), DIMENSION(:), ALLOCATABLE, TARGET:: costsgrid
        REAL, POINTER :: deathrisk, emprisk
        INTEGER, POINTER :: zmin, zmax
        REAL(8), DIMENSION(2), TARGET :: null_vec
        REAL(8) :: weightDprimel, weightaprimel, hpcurrent
        REAL(8) :: aprime_shift
        INTEGER :: i, ue_marker, Dprimel, Dprimeh, aprimel, aprimeh
        REAL(8), DIMENSION(4) :: EVholder, shockVecU, shockvecM, shockVec

        ! need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.
        ! Find nearest grid points and linearly interpolate.

        ! NOTE THE DIFFERENCE WITH THE ADJUST FUNCTION. We only interpolate
        ! over a because rental housing does not last for longer than one
        ! period. See indexing on EV as well.
        call weightPrimeL(Dval,Dnodes,weightDprimel,Dprimel,Dprimeh)
        ! TODO: are there transaction costs to renting?
        ! Maybe not as that should be absorbed into quadratic costs.
        null_vec = (/ 0.0, 1.0 /)
        tmpgrid = (/ sigma_temp, -sigma_temp /)  ! => null_vec(1:1)
        probtmp = (/ 0.5, 0.5 /)  ! => null_vec(2:2)
        if (rental) then
            weightDprimel = 0.0
            Dprimel = 1
            Dprimeh = Dprimel
        else if (adjust) then
            allocate(costsgrid(transgridsize))
            costsgrid = 0.0! transcst*(Dval*price) 
            !tmpgrid = 0.0! => costsgrid 
            probtmp => transcstProbs
        end if
 
        if (tval - 1 <= Tretire .AND. uval == 2) then
            zmin => minexpectationz(zindex)
            zmax => maxexpectationz(zindex)
            emprisk => empprob(tval)
            ue_marker = 1
            if (tval-1 == Tretire) ue_marker = 2

            shockVecU = (/ emprisk, 1.0-emprisk, emprisk, 1.0-emprisk /)
            if (movprob > 0.0 .AND. Down > 0.0) then
                shockVecM = (/ 1.0-movprob, 1.0-movprob, movprob, movprob/)
            else if (movprobR > 0.0 .AND. Down == 0.0) then
                shockVecM = (/ 1.0-movprobR, 1.0-movprobR, movprobR, movprobR/)
            else
                shockVecM = (/ 1.0, 1.0, 0.0, 0.0 /)
            end if
            shockVec = shockVecU*shockVecM
            EVholder= 0.0  ! Unaccessed value if mov/unemp shocks inapplicable
 
            ! Use pointers to simplify which EV array to reference
            ! in the case of anticipated policy
            ! This is the case with employment
            if (expecting) then
                EVmat => EVpol(:, 2, Dprimel:Dprimeh, :, tval, hpval)
                EVmatM => EVpolMov(:, 2, Dprimel:Dprimeh, :, tval, hpval)
                EVmatMR => EVpolMovR(:, 2, Dprimel:Dprimeh, :, tval, hpval)
            else
                EVmat => EV(:, 2, Dprimel:Dprimeh, :, tval, hpval)
                EVmatM => EVMov(:, 2, Dprimel:Dprimeh, :, tval, hpval)
                EVmatMR => EVMovR(:, 2, Dprimel:Dprimeh, :, tval, hpval)
            end if
 
 
            EVholder(1) = temp_expect_loop(EVmat, aval, tmpgrid, probtmp, weightDprimel,&
                                           pwealth, zindex, zmin, zmax, transfers)
            if (movprob > 0.0 .AND. Down > 0.0) then
                EVholder(3) = temp_expect_loop(EVmatM, aval, tmpgrid, probtmp, weightDprimel,&
                                               pwealth, zindex, zmin, zmax, transfers)
            else if (movprobR> 0.0 .AND. Down == 0.0) then
                EVholder(3) = temp_expect_loop(EVmatMR, aval, tmpgrid, probtmp, weightDprimel,&
                                               pwealth, zindex, zmin, zmax, transfers)
            end if
 
            if (empprob(tval) < 1.0) then
            if (expecting) then
                EVmatU => EVpol(:, ue_marker, Dprimel:Dprimeh, :, tval, hpval)
                EVmatM => EVpolMov(:, ue_marker, Dprimel:Dprimeh, :, tval, hpval)
                EVmatMR => EVpolMovR(:, ue_marker, Dprimel:Dprimeh, :, tval, hpval)
            else
                EVmat => EV(:, ue_marker, Dprimel:Dprimeh, :, tval, hpval)
                EVmatM => EVMov(:, ue_marker, Dprimel:Dprimeh, :, tval, hpval)
                EVmatMR => EVMovR(:, ue_marker, Dprimel:Dprimeh, :, tval, hpval)
            end if

            EVholder(2) = temp_expect_loop(EVmatU, aval, tmpgrid, probtmp, weightDprimel,&
                                           pwealth, zindex, zmin, zmax, transfers)
            if (movprob > 0.0 .AND. Down > 0.0) then
                EVholder(4) = temp_expect_loop(EVmatM, aval, tmpgrid, probtmp, weightDprimel,&
                                               pwealth, zindex, zmin, zmax, transfers)
            else if (movprobR> 0.0 .AND. Down == 0.0) then
                EVholder(4) = temp_expect_loop(EVmatMR, aval, tmpgrid, probtmp, weightDprimel,&
                                               pwealth, zindex, zmin, zmax, transfers)
            end if
            end if

            expectation = continuation(beta2, shockVec, EVholder(1:4))
 
        else
            ! TODO: do retired people just not move? Modelling decision
            if (expecting) then
                EVmatD => EVpol(:, 2, Dprimel:Dprimeh, :, Tdie+1, hpval)
                EVmatU => EVpol(:, 2, Dprimel:Dprimeh, :, tval, hpval)
                EVmat => EVpol(:, uval, Dprimel:Dprimeh, :, tval, hpval)
            else
                EVmatD => EV(:, 2, Dprimel:Dprimeh, :, Tdie+1, hpval)
                EVmatU => EV(:, 2, Dprimel:Dprimeh, :, tval, hpval)
                EVmat => EV(:, uval, Dprimel:Dprimeh, :, tval, hpval)
            end if

            EVholder(1) = temp_expect_loop(EVmat, aval, tmpgrid, probtmp, weightDprimel,&
                                           pwealth, zindex, zindex, zindex, transfers)
 
            if (tval-1 > Tretire) then
                deathrisk => deathprob(tval-Tretire-1)
                ! When retired, no risk in income but risk in death
                EVholder(2) = temp_expect_loop(EVmatD, aval, tmpgrid, probtmp, weightDprimel,&
                                               pwealth, zindex, zindex, zindex, transfers)
                expectation = continuation(beta2, (/ 1.0-deathrisk, deathrisk /),&
                                                  EVholder(1:2))
            else
                ! When unemployed, stochastic decision of staying in
                ! unemployment or rejoining at previous income level
                EVholder(2) = temp_expect_loop(EVmatU, aval, tmpgrid, probtmp, weightDprimel,&
                                               pwealth, zindex, zindex, zindex, transfers)
                expectation = continuation(beta2, (/ unempprob, 1.0-unempprob /),&
                                                  EVholder(1:2))
            end if

        end if

        if (expectation/= expectation) then
            write(*,*) "err0"
        end if
    END FUNCTION ! %>

    REAL(8) FUNCTION valfuncadjust(p, state, pol_opts)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    The objective function in the case where agent adjusts into
        !    owned housing, possibly paying transaction costs/receiving a policy/
        !    scrapping old durable. This function is general and can emulate
        !    many policies based on parameters set upon model initialization.
        ! 
        !    Modified: 08/09/2019
        !
        !    PARAMETERS
        !
        !    p: Vector of length 2, representing the choice variables
        !    (a' and D', C is residual from budget constraint)
        !
        !    state: Vector of state agent is currently in, both idiosyncratic
        !    as well as aggregate (price vector, transfer level)
        !
        !    pol_opts: Vector of length 2 specifying policy options not already
        !    covered by model parameters.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IMPLICIT NONE
        REAL(8), DIMENSION(2), TARGET :: p
        REAL(8), DIMENSION(9), INTENT(IN), TARGET :: state
        LOGICAL, DIMENSION(2), INTENT(IN) :: pol_opts
        REAL(8), POINTER :: a, D, hpholder, zindex, t, hpindex, balholder, unemp, rentimpute
        REAL(8), POINTER :: aprime, Dprime
        REAL(8) :: consumption, currentincome, hpcurrent, transcost
        REAL(8) :: laborsupply
        REAL(8) :: assetholder, futliquid, EVholder
        REAL(8) :: thetaholder, transfers, scrapholder, adjshifts, policyfactor, &
                   rebateholder, scrapvholder, transferfuture, transfer_down, dep

        !write(*,*) state

        zindex => state(1)
        unemp => state(2)
        D => state(3)
        a => state(4)
        hpholder => state(5)
        hpindex => state(6)
        t => state(7)
        balholder => state(8)
        if (hpindex==0) then
            hpcurrent = 1.0
        else
            hpcurrent = hpindex
        end if
        rentimpute => state(9)

        aprime => p(1)
        Dprime => p(2)
        

        if (polLevel(2) > 0.0) then  ! Fixed durable adjustment size (even if algorithm minimizes in 2 dimensions)
            Dprime => polLevel(2)
        end if
        if (aprime>anodes(agridsize)) then
            aprime => anodes(agridsize)
        end if
        if (Dprime>Dmax) then
            Dprime => Dnodes(dgridsize)
        end if

        thetaholder=(1-polLevel(3))*theta   ! Aggregate temp shock to DP level
        scrapholder=scrapped
        scrapvholder=scrapvalue
        if (pol_opts(1)) then
            transfers=1.0
            if (.NOT. pctageflag(zindex, t)) then  ! adjtransfer expressed in nominal terms
                transferfuture=(1.0-eta_transfer)*adjtransfer(zindex, t)
                transfer_down = eta_transfer*adjtransfer(zindex, t)/(exp(hpholder)*Dprime)
            else  ! adjtransfer means *percentage* of housing value, so the above times housing value
                transferfuture=(1.0-eta_transfer)*adjtransfer(zindex, t)*exp(hpholder)*Dprime
                transfer_down = eta_transfer*adjtransfer(zindex, t)
            end if
        else
            transfers=0.0
            transferfuture=transfers
            transfer_down=transfers
        end if
        ! Inclusive of unconditional transfer less extra cost of purchase
        ! (related to temp shock to DP)
        currentincome = incomeCall(state) + polLevel(1) - polLevel(4)
        dep = maint*delta

        if (aprime>=borrowconstraint .AND. Dprime>=Dmin) then
            ! liquidity const. today
            call assetReturns(state(3:5), aprime, Dprime, (1.0-delta)/(1.0+rborrow),&
                              discountflag*transfer_down, thetaholder, &
                              inctax(zindex, t), assetholder, futliquid)
            !transcost = F2*((Dprime-rentimpute)/((Dprime+rentimpute)))**(2.0)*((Dprime+rentimpute)) + F*D*(1-maint*delta-dtau)*exp(hpholder)
            transcost = F2*((Dprime-0.0)/((Dprime+rentimpute)))**(2.0)*((Dprime+rentimpute)) + F*D*(1-maint*delta-dtau)*exp(hpholder)
            ! Down payment w/ transaction costs (possibly rebated)
            rebateholder = 1.0! 1.0 + F/thetaholder
            ! Straight monetary transfer (for both CARS and FTHB designs)
            adjshifts = (1.0-downflag)*(1.0-discountflag)*(&
                        transfers*adjtransfer(zindex, t)*eta_transfer)
            ! CARS-specific scrappage policy parameters
            adjshifts = adjshifts - scrapholder*(-scrapvholder*(1.0-transfers)&
                        +D*(1-dep-dtau)*exp(hpholder)) 
            
            ! The value of policyfactor is face value of subsidy normalized by house value,
            ! Normalized again by the down payment factor thetaholder
            policyfactor = transfer_down/thetaholder
            ! Then, policyfactor is only turned on if we specified the sequence of
            ! policy dummies in share module such that the present policy is interpreted
            ! as one on the down payment.
            policyfactor = downflag*(1.0-discountflag)*policyfactor

            consumption = presCons(wealthCall(state(3:5), currentincome, assetholder) - aprime, &
                                   Dprime, hpholder, thetaholder*rebateholder*(1.0-policyfactor), &
                                   -transcost + adjshifts)
            !if (zindex == 10 .AND. D == 0.0 .AND. t == 2 .AND. 0.2 < a .AND. a < 0.3) then
            !write(*,*) consumption
            !end if
            ! Price level of down payment (possibly implied temp subsidy on down)

            if (consumption > 0) then
                valfuncadjust = presVal(consumption, Dprime, discountflag*transfer_down, utilscale(t))

                EVholder=expectation(INT(zindex), (1-delta+dep+discountflag*transfer_down)*Dprime,&
                                     aprime, t+1, hpindex+1, exp(hpholder), unemp, futliquid, D,&
                                     pol_opts(2), .FALSE.,.TRUE.,transferfuture- &
                                     downflag*(1.0-discountflag)*(1.0-exp(-transfer_down))*Dprime)

                valfuncadjust= valfuncadjust + EVholder
                if (valfuncadjust /= valfuncadjust) then  ! NA value
                    write(*,*) EVholder
                    valfuncadjust=-3.0e6
                    return
                end if
            else
                valfuncadjust=-2.0e6
            end if
        else
             valfuncadjust=-1.0e6
        end if
        valfuncadjust=-valfuncadjust  !powell minimizes

    END FUNCTION valfuncadjust ! %>

    REAL(8) FUNCTION valfuncrent(p, state, pol_opts)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    The objective function in the case where agent enters the rental
        !    market, consuming in the present period rental housing services
        !    that disappear in the following period.
        ! 
        !    Modified: 08/09/2019
        !
        !    PARAMETERS
        !
        !    p: Vector of length 2, representing the choice variables
        !    (a' and D', C is residual from budget constraint). Note that
        !    current housing consumption D' is chosen - it's just not carried
        !    into the next period.
        !
        !    state: Vector of state agent is currently in, both idiosyncratic
        !    as well as aggregate (price vector, transfer level)
        !
        !    pol_opts: Vector of length 2 specifying policy options not already
        !    covered by model parameters.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IMPLICIT NONE
        REAL(8), DIMENSION(2), TARGET :: p
        REAL(8), DIMENSION(9), INTENT(IN), TARGET :: state
        LOGICAL, DIMENSION(2), INTENT(IN) :: pol_opts
        REAL(8), POINTER :: a, D, hpholder, zindex, t, hpindex, balholder, unemp
        REAL(8), POINTER :: aprime, Dprime
        REAL(8) :: consumption, currentincome, hpcurrent
        REAL(8) :: laborsupply
        REAL(8) :: assetholder, futliquid, EVholder
        REAL(8) :: rentholder, rentpay, transferholder, dep

        !write(*,*) state

        zindex => state(1)
        unemp => state(2)
        D => state(3)
        a => state(4)
        hpholder => state(5)
        hpindex => state(6)
        t => state(7)
        balholder => state(8)
        if (hpindex==0) then
            hpcurrent = 1.0
        else
            hpcurrent = hpindex
        end if

        currentincome = incomeCall(state)+ polLevel(1)
        dep = maint*delta

        aprime => p(1)
        Dprime => p(2)
        if (aprime>anodes(agridsize)) then
            aprime => anodes(agridsize)
        end if
        if (Dprime>Dmax) then
            Dprime => Dnodes(dgridsize)
        end if

        if (pol_opts(1)) then
        transferholder=renttransfer
            if (pctageflag(zindex, t)) then
                transferholder=transferholder*Dprime*exp(rentelasticity*hpholder)   
            end if
        else
        transferholder=0
        end if

        ! First element is current financial assets, 2nd is next period
        ! finwealth = (/ a-D*(1-prF-dtau)*(1-theta)*exp(hpholder), aprime /)

        if (aprime>=borrowconstraint .AND. Dprime>=Dnodes(1)) then
            call assetReturns(state(3:5), aprime, Dprime, (1.0-delta)/(1.0+r),0.0, 1.0,&
                              inctax(zindex, t), assetholder, futliquid)
            if (t <= Tretire) then
                rentholder=r_rental(hpcurrent)
            else
                rentholder=r_rental_retire(hpcurrent)
            end if
            rentpay = rentholder*Dprime*exp(rentelasticity*hpholder)

            ! TODO: Other tax credits for renters?
            consumption = presCons(wealthCall(state(3:5), currentincome, assetholder) - aprime, &
                                   0.0, hpholder, 0.0, -rentpay+transferholder)
            if (consumption > 0) then
            ! %< need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.  Find nearest grid points and linearly interpolate.
                valfuncrent = presVal(consumption, Dprime, -1.0+rentUtil, utilscale(t))
                ! NOTE THE DIFFERENCE WITH THE ADJUST FUNCTION. We only interpolate
                ! over a because rental housing does not last for longer than one
                ! period. See indexing on EV as well.
                EVholder=expectation(INT(zindex), 0.0, aprime, t+1, hpindex+1,exp(rentelasticity*hpholder),&
                                     unemp, futliquid, 0.0, pol_opts(2), .TRUE.,.FALSE.)
                valfuncrent= valfuncrent + EVholder
                if (valfuncrent /= valfuncrent)  then ! NA value
                    valfuncrent=-6.0e6
                    return
                end if
            ! %>
            else
                valfuncrent=-5.0e6
            end if
        else
             valfuncrent=-4.0e6
        end if
        valfuncrent=-valfuncrent  !powell minimizes

    END FUNCTION valfuncrent ! %>

    REAL(8) FUNCTION valfuncnoadjust(aprime, state, pol_opts)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    The objective function in the case where agent does not adjust
        !    her housing good, therefore choosing only over assets a'.
        !    General enough to cover the case of both forced depreciation
        !    of durable and forced maintenance (though not endongenous choice
        !    of maintenace).
        ! 
        !    Modified: 08/09/2019
        !
        !    PARAMETERS
        !
        !    aprime: Real number representing asset policy a'.
        !
        !    state: Vector of state agent is currently in, both idiosyncratic
        !    as well as aggregate (price vector, transfer level)
        !
        !    pol_opts: Vector of length 2 specifying policy options not already
        !    covered by model parameters.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), INTENT(IN), TARGET :: aprime
        REAL(8), DIMENSION(9), INTENT(IN), TARGET :: state
        LOGICAL, INTENT(IN) :: pol_opts
        REAL(8), POINTER :: a, Dindex, hpholder, zindex, t, hpindex, balholder, unemp
        REAL(8), POINTER :: D
        REAL(8) :: consumption, aprime_adj, currentincome, hpcurrent
        REAL(8) :: laborsupply
        REAL(8) :: assetholder, futliquid, EVholder
        REAL(8) :: dep

        !write(*,*) state

        zindex => state(1)
        unemp => state(2)
        Dindex => state(3)
        a => state(4)
        hpholder => state(5)
        hpindex => state(6)
        t => state(7)
        balholder => state(8)
        if (hpindex==0) then
            hpcurrent = 1.0
        else
            hpcurrent = hpindex
        end if

        currentincome = incomeCall(state)+ polLevel(1)
        dep = maint*delta

        D => Dnodes(Dindex)
        if (D>Dmax) then
            D => Dnodes(dgridsize)
        end if
        aprime_adj = max(aprime+(1-theta)*D*hpdelta, borrowconstraint)
        if (aprime_adj>anodes(agridsize)) then
            aprime_adj = anodes(agridsize)
        end if

        if (aprime_adj>=borrowconstraint .AND. (D>=Dmin .OR. maint == 0.0 )) then
            call assetReturns((/ D, a, hpholder /), aprime, D, (1.0-delta)/(1.0+rborrow),&
                               0.0, theta, inctax(zindex, t), assetholder, futliquid)

            consumption = presCons(wealthCall((/ D, a, hpholder/), currentincome, assetholder) - aprime_adj, &
                                   D, hpholder, theta, 0.0)

            if (consumption>0) then
            ! %< need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.  Find nearest grid points and linearly interpolate.
                valfuncnoadjust = presVal(consumption, D, 0.0, utilscale(t))

                ! NOTE THE DIFFERENCE WITH THE ADJUST FUNCTION. We interpolate
                ! over a AND d because depreciation could occur even w/out adjustment.
                EVholder=expectation(INT(zindex), (1-delta+dep)*D,&
                                     aprime, t+1, hpindex+1, exp(hpholder), unemp, futliquid,&
                                     D, pol_opts, .FALSE.,.FALSE.)
                valfuncnoadjust= valfuncnoadjust + EVholder
                if (valfuncnoadjust /= valfuncnoadjust)  then ! NA value
                    valfuncnoadjust=-9.0e6
                    return
                end if
            ! %>
            else
                valfuncnoadjust=-8.0e6
            end if
        else
             valfuncnoadjust=-7.0e6
        end if
        valfuncnoadjust=-valfuncnoadjust  !powell minimizes

    END FUNCTION valfuncnoadjust ! %>

    SUBROUTINE adj_func_comp(p, state, aArray, Darray, rentalArray, cArray, &
                             choiceArray, EVarray, nocreditBool, price,&
                             consOut, rentOut, choiceOut, welfOut)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Subroutine evaluates utility from two separate optimization problems
        !    problems (e.g. taking up the policy now vs. waiting to take up
        !    policy later), then takes the max of the two.
        !    her housing good, therefore choosing only over assets a'.
        !    General enough to cover the case of both forced depreciation
        ! 
        !    Modified: 05/12/2019
        !
        IMPLICIT NONE
        REAL(8), DIMENSION(2), INTENT(INOUT) :: p
        REAL(8), DIMENSION(:), INTENT(IN) :: state
        LOGICAL, DIMENSION(:), INTENT(IN) :: nocreditBool
        REAL(8), INTENT(IN) :: price
        REAL(8), INTENT(INOUT) :: consOut, rentOut, choiceOut, welfOut
        REAL(8), dimension(:,:,:,:,:,:), INTENT(IN) :: aArray,DArray,rentalArray, cArray,choiceArray, EVarray
        REAL(8) :: rentchoice, constemp, renttemp, choicetemp, welftemp,&
            weightDl, weightal
        INTEGER :: Dl, Dh, al, ah

        REAL(8), DIMENSION(2) :: testutil, newteststate

        call pol_linworking(state, aArray, Darray, cArray, choiceArray, EVarray, &
            newteststate(1), newteststate(2),&
            constemp, renttemp, choicetemp, welftemp)

        call weightprimeL(state(3), Dnodes, weightDl, Dl, Dh)
        call weightprimeL(state(4), anodes, weightal, al, ah)

        rentchoice = rentalArray(state(1), state(2), Dl, al, state(5), state(6))

        ! Inputted function: adjustment with credit but with extra conditions
        testutil(1) = valfuncadjust((/p(2), p(1) /),&
            (/ state(1:4), price, state(6), state(5), balancer_internal, rentchoice /),&
             (/.TRUE., .FALSE./))
        ! Computed within subroutine: adjustment without credit, but possibly
        ! expecting one as well
        testutil(2) = valfuncadjust(newteststate,&
            (/ state(1:4), price, state(6), state(5), balancer_internal, rentchoice /),&
             nocreditBool)
        if (testutil(1) <= testutil(2)) then
            return
        else
            p(2) = newteststate(1)
            p(1) = newteststate(2)
            consOut = constemp
            rentOut = renttemp
            choiceOut = choicetemp
            welfOut = welftemp
        end if
    END SUBROUTINE ! %>

        REAL(8) FUNCTION get_aprime(Dprime, substate, downfact, currentincome, assetholder, adjshifts, policyfactor, thetaholder)
!            Returns the REAL(8) value of aprime along the BC governing the simplex given aprime. Can be used to obtain a (Q,D) pair
!            PARAMS:
!            - <REAL(8)> Dprime: Dprime value for the (Q,D) pair
!            - <REAL(8)>  downfact, currentincome, assetholder, adjshifts
!            - <REAL(8), DIMENSION(3)> substate: array containing information about current D, a, and hpholder
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: Dprime, downfact, currentincome, assetholder, adjshifts, policyfactor, thetaholder
            REAL(8), DIMENSION(4), INTENT(IN), TARGET :: substate
            REAL(8), POINTER :: D, a, hpholder, rentimpute
            REAL(8) :: myThet, myDelt, myKap, rebateholder
            ! Initialize pointers
            D => substate(1)
            a => substate(2)
            hpholder => substate(3)
            rentimpute => substate(4)

            rebateholder = 1  ! Could we generalize this factor?

            myDelt = downflag*(1-discountflag)*policyfactor/thetaholder
            myThet = exp(hpholder)*thetaholder*rebateholder
            myKap = wealthCall(substate(1:3), currentincome, assetholder) + adjshifts - (1.0-scrapped)*F*D*(1-maint*delta-dtau)*exp(hpholder)
            ! The -0.01 is a perturbation
            !get_aprime = myKap - myThet*(1.0 - myDelt)*Dprime - F2*((Dprime-rentimpute)**2.0)/(Dprime+rentimpute) - offset_consumption
            get_aprime = myKap - myThet*(1.0 - myDelt)*Dprime - F2*((Dprime-0.0)**2.0)/(Dprime+rentimpute) - offset_consumption
            get_aprime = min(get_aprime, amax)
            get_aprime = max(get_aprime, amin)
        END FUNCTION get_aprime


        REAL(8) FUNCTION get_Dprime(aprime, substate, downfact, currentincome, assetholder, adjshifts, policyfactor, thetaholder)
            !pass adjtransfer as param
!            Returns the REAL(8) value of Dprime along the BC governing the simplex given aprime. Can be used to obtain a (Q,D) pair
!            PARAMS:
!            - <REAL(8)> aprime: aprime (Q) value for the (Q,D) pair
!            - <REAL(8)> hpholder, downfact, currentincome, assetholder, adjshifts
!            - <REAL(8), DIMENSION(3)> substate: array containing information about current D, a, and hpholder
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: aprime, downfact, currentincome, assetholder, adjshifts, policyfactor, thetaholder
            REAL(8), DIMENSION(4), INTENT(IN), TARGET :: substate
            REAL(8), POINTER :: D, a, hpholder, rentimpute
            REAL(8) :: myThet, myKap, myDelt, myA, myB, myC, rebateholder
            ! Initialize pointers
            D => substate(1)
            a => substate(2)
            hpholder => substate(3)
            rentimpute => substate(4)

            rebateholder = 1  ! Could we generalize this factor?
            !Rewrite BC as a' + D'(myThet) + (F2*(D'-D)**2)/(D'+D) = myKap
            myDelt = downflag*(1-discountflag)*policyfactor/thetaholder 
            myThet = exp(hpholder)*thetaholder*rebateholder !exp(hpholder)*thetaholder*rebateholder*downflag*(1-discountflag)
            myKap = wealthCall(substate(1:3), currentincome, assetholder) + adjshifts -aprime - &
                (1.0-scrapped)*F*D*(1-maint*delta-dtau)*exp(hpholder) - offset_consumption
            myA = myThet + F2
            myB = myThet*rentimpute - myThet*myDelt - myKap
            myC = -myKap*rentimpute - myThet*myDelt*rentimpute
            !myB = -myKap + aprime + (myThet - myDelt + delta*exp(hpholder) - 2.0*F2)*0.0
            !myC = (-myKap+aprime - myDelt)*0.0 + (F2 - delta*exp(hpholder))*0.0**2 
            !possible for it to be <0?
            if (F2 <= 1e-7) then
                get_Dprime = myKap/myThet + myDelt
            else
            if (sqrt(myB**(2.0)-4*myA*myC) > 0) then
                get_Dprime = (-myB + sqrt(myB**(2.0)-4*myA*myC))/2/myA
            else
                get_Dprime = (-myB - sqrt(myB**(2.0)-4*myA*myC))/2/myA
            end if
            end if
            get_Dprime = min(get_Dprime, Dmax)
            get_Dprime = max(Dmin, get_Dprime)
        END FUNCTION get_Dprime

        SUBROUTINE get_simplex(amoebaGrid, state, currentincome, assetholder, myUsercost, pol_opts)
            !For now assume r=rborrow
            !            Populates amoebaGrid with 3 (Q,D) pairs which form the vertices of the simplex after accounting for subsidy.
            !            Intuition: a simplex is constructed from 3 constraints, and the vertices are the coordinates of interest:
            !            1) aprime >= borrowconstraint
            !            2) Dprime >= Dmin
            !            3) C >= 0
            !            PARAMS:
            !            - <REAL(8)> hpholder, downfact, currentincome, assetholder, othcost
            !            - <REAL(8), DIMENSION(8)> state
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: currentincome, assetholder, myUsercost
            LOGICAL, DIMENSION(2), INTENT(IN) :: pol_opts
            REAL(8), DIMENSION(9), INTENT(IN), TARGET :: state
            REAL(8), POINTER :: a, D, hpholder, zindex, t, hpindex, balholder, unemp, rentimpute
            REAL(8), DIMENSION(2,3), INTENT(OUT) :: amoebaGrid
            REAL(8), DIMENSION(4) :: substate
            REAL(8) :: thetaholder, transfers, adjshifts, myTheta, transcost,&
                    rebateholder, transfer_down, dep, hpcurrent, downfact
            ! Initialize pointers
            zindex => state(1); unemp => state(2); D => state(3); a => state(4)
            hpholder => state(5); hpindex => state(6); t => state(7); balholder => state(8); rentimpute => state(9)
            if (hpindex==0) then
                hpcurrent = 1.0
            else
                hpcurrent = hpindex
            end if
            ! Incorporate subsidy
            myTheta = 1-(1-thetamatlab)*myUsercost
            thetaholder=(1-polLevel(3))*myTheta
            if (pol_opts(1)) then
                transfers=1.0
                if (.NOT. pctageflag(zindex, t))  then ! adjtransfer expressed in nominal terms
                    transfer_down = eta_transfer*adjtransfer(zindex, t)/exp(hpholder)
                end if
            else
                transfers=0.0
                transfer_down=transfers
            end if
            dep = maint*delta
            ! Down payment w/ transaction costs (possibly rebated)
            rebateholder = 1.0!+(1.0-transfers*rebatedflag) !should just be 1.0
            ! Straight monetary transfer (for both CARS and FTHB designs)
            adjshifts = (1.0-downflag)*(1.0-discountflag)*(&
                    transfers*adjtransfer(zindex, t)*eta_transfer)
            ! CARS-specific scrappage policy parameters
            adjshifts = adjshifts - scrapped*(-scrapvalue*(1.0-transfers)&
                    +D*(1-dep-dtau)*exp(hpholder))
            ! Actual grid construction
            substate = (/D, a, hpholder, rentimpute /)
            downfact = thetaholder * rebateholder
            ! (borrowconstraint, Dmin)
            amoebaGrid(1, 1) = borrowconstraint
            amoebaGrid(2, 1) = Dmin
            ! (Q when BC @ Dmin, Dmin)
            amoebaGrid(1, 2) = get_aprime(Dmin, substate, downfact, currentincome, assetholder, adjshifts, &
                                          transfer_down/Dmin, thetaholder)
            amoebaGrid(2, 2) = Dmin
            ! (borrowconstraint, D when BC @ borrowconstraint)
            amoebaGrid(1, 3) = borrowconstraint
            ! 
            amoebaGrid(2, 3) = get_Dprime(borrowconstraint, substate, downfact, currentincome, assetholder, adjshifts,&
                                          transfer_down*exp(hpholder), thetaholder)
            if (D>=0.5 .AND. zindex == 10 .AND. a <= 0.25) then
                !write(*,*) amoebaGrid
            end if
            !if (t < 20) then
            !write(*,*) amoebaGrid
            !end if
            !if (D <= 0.1 .AND. a <= 0.1 .and. currentincome <= 1.0) then
            !    write(*,*) state(1:4), currentincome
            !    write(*,*) amoebaGrid
            !    transcost = F2*((amoebaGrid(2, 3)-D)/((amoebaGrid(2, 3)+D)/2.0))**(2.0)*((amoebaGrid(2, 3)+D)/2.0) + F*D*(1-maint*delta-dtau)*exp(hpholder)
            !    write(*,*) presCons(wealthCall(state(3:5), currentincome, assetholder) - amoebaGrid(1,3), &
            !                        amoebaGrid(2, 3), hpholder, thetaholder*rebateholder*(1.0-transfer_down/amoebagrid(2,3)), &
            !                        -transcost + adjshifts)
            !end if
            !write(*,*) amoebaGrid
        END SUBROUTINE

subroutine get_simplex_rent(amoebaGrid, state, currentincome, assetholder)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: currentincome, assetholder
    REAL(8), DIMENSION(9), INTENT(IN), TARGET :: state
    REAL(8), POINTER :: a, D, hpholder, zindex, t, hpindex, balholder, unemp
    REAL(8), DIMENSION(2,3), INTENT(OUT) :: amoebaGrid
    REAL(8), DIMENSION(3) :: substate
    REAL(8) :: hpcurrent, rentholder
    ! Initialize pointers
    zindex => state(1); unemp => state(2); D => state(3); a => state(4)
    hpholder => state(5); hpindex => state(6); t => state(7); balholder => state(8)
    if (hpindex==0) then
        hpcurrent = 1.0
    else
        hpcurrent = hpindex
    end if
    if (t <= Tretire) then
        rentholder=r_rental(hpcurrent)
    else
        rentholder=r_rental_retire(hpcurrent)
    end if
    ! Actual grid construction
    substate = (/D, a, hpholder/)
    ! (borrowconstraint, Dnodes(1))
    amoebaGrid(1, 1) = borrowconstraint
    amoebaGrid(2, 1) = Dnodes(2)
    ! (Q when BC @ Dnodes(1), Dnodes(1))
    amoebaGrid(1, 2) = get_aprime_rent(Dnodes(2), substate, currentincome, assetholder, rentholder)
    amoebaGrid(2, 2) = Dnodes(2)
    ! (borrowconstraint, D when BC @ borrowconstraint)
    amoebaGrid(1, 3) = borrowconstraint
    amoebaGrid(2, 3) = get_Dprime_rent(borrowconstraint, substate, currentincome, assetholder, rentholder)
!    write(*,*) amoebaGrid

end subroutine get_simplex_rent


REAL(8) FUNCTION get_aprime_rent(Dprime, substate, currentincome, assetholder, rentholder)
    IMPLICIT NONE
    REAL(8), DIMENSION(3), INTENT(IN), TARGET :: substate
    REAL(8), POINTER :: a, D, hpholder
    REAL(8), INTENT(IN) :: Dprime, currentincome, assetholder, rentholder
    D => substate(1)
    a => substate(2)
    hpholder => substate(3)
    get_aprime_rent = wealthCall(substate, currentincome, assetholder) - Dprime*rentholder*exp(rentelasticity*hpholder)- offset_consumption
    get_aprime_rent = min(get_aprime_rent, amax)
    get_aprime_rent = max(get_aprime_rent, amin)
END FUNCTION get_aprime_rent


REAL(8) FUNCTION get_Dprime_rent(aprime, substate, currentincome, assetholder, rentholder)
    IMPLICIT NONE
    REAL(8), DIMENSION(3), INTENT(IN), TARGET :: substate
    REAL(8), POINTER :: a, D, hpholder
    REAL(8), INTENT(IN) :: aprime, currentincome, assetholder, rentholder
    D => substate(1)
    a => substate(2)
    hpholder => substate(3)
    get_Dprime_rent = (wealthCall(substate, currentincome, assetholder) - aprime - offset_consumption)/(rentholder*exp(rentelasticity*hpholder)) 
    get_Dprime_rent = min(get_Dprime_rent, Dmax)
    get_Dprime_rent = max(0.0, get_Dprime_rent)
END FUNCTION get_Dprime_rent


end module lifecycle_vfuncs
