module lifecycle_vfuncs

    USE share
    USE lifeCycle_algs
    IMPLICIT NONE

    CONTAINS

    REAL(8) FUNCTION asset_shift(aprime, atrans, shifts) ! %<
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: aprime, atrans, shifts

        ! Embedding temporary shocks into asset state space
            if (atrans < 0 .AND. atrans + shifts < 0) then
                    asset_shift = min(max(aprime+(shifts)/(1+rborrow),amin),&
                                      anodes(agridsize))
            else if (atrans < 0) then
                    asset_shift = min(max(aprime+(shifts)/(1+r)+&
                                      (rborrow-r)*atrans,amin),anodes(agridsize))
            else if (atrans >= 0 .AND. atrans + shifts >= 0) then
                    asset_shift = min(max(aprime+(shifts)/(1+r),amin),&
                                      anodes(agridsize))
            else
                    asset_shift = min(max(aprime+(shifts)/(1+rborrow)+&
                                      (r-rborrow)*atrans,amin),anodes(agridsize))
            end if


    END FUNCTION ! %>


    FUNCTION EV_linpol(awgt, Dwgt, zmin, zmax, EVinput) ! %<
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

    REAL(8) FUNCTION continuation(discount, prob, e1, e2) ! %<
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: discount, prob, e1, e2

        ! Tiny function that outputs discounted continuation
        ! value given two future states and outcomes.
        if (prob > 1) call nrerror("Probability undefined")
        continuation = discount*(prob*e1 + (1-prob)*e2)

    END FUNCTION ! %>

    REAL(8) FUNCTION expectation(zindex, Dval, aval, tval, hpval, uval, pwealth,& ! %< 
                                 expecting, rental, transfers)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: zindex
        REAL(8), INTENT(IN) :: aval, Dval, tval, hpval, uval, pwealth
        LOGICAL, INTENT(IN) :: expecting, rental
        REAL(8), OPTIONAL :: transfers
        REAL(8), DIMENSION(:,:,:), POINTER :: EVmat, EVmatU, EVmatD
        INTEGER, POINTER :: zmin, zmax
        REAL(8) :: weightDprimel, weightaprimel, hpcurrent
        REAL(8), DIMENSION(tmpgridsize+1) :: aprime_shift
        REAL(8), DIMENSION(tmpgridsize) :: perm_expect
        INTEGER :: i, Dprimel, Dprimeh, aprimel, aprimeh
        REAL(8), DIMENSION(2) :: EVholder

        ! need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.
        ! Find nearest grid points and linearly interpolate.

        ! NOTE THE DIFFERENCE WITH THE ADJUST FUNCTION. We only interpolate
        ! over a because rental housing does not last for longer than one
        ! period. See indexing on EV as well.
        if (rental) then
            weightDprimel = 0.0
            Dprimel = 1
            Dprimeh = Dprimel
        else
            call weightPrimeL(Dval,Dnodes,weightDprimel,Dprimel,Dprimeh)
        end if
 
        if (present(transfers)) then
            aprime_shift(tmpgridsize+1) = asset_shift(aval, pwealth, transfers)
        else
            aprime_shift(tmpgridsize+1) = asset_shift(aval, pwealth, 0.0)
        end if
 
        if (tval-1 <= Tretire .AND. uval == 2) then
            zmin => minexpectationz(zindex)
            zmax => maxexpectationz(zindex)
            if (present(transfers)) then
                do i=1, tmpgridsize
                    aprime_shift(i) = asset_shift(aval, pwealth,&
                                                  transfers + tmpnodes(i))
                end do
            else
                do i=1, tmpgridsize
                    aprime_shift(i) = asset_shift(aval, pwealth,tmpnodes(i))
                end do
            end if
      
            ! TODO: Call future price as a state rather than a call to
            ! hpnodes
            do i=1, tmpgridsize
                call weightPrimeL(aprime_shift(i),anodes,weightaprimel,aprimel,aprimeh)
                ! TODO: Use pointers to simplify which EV array to reference
                ! in the case of anticipated policy
                if (expecting) then
                    EVmat => EVpol(:, 2, Dprimel:Dprimeh, aprimel:aprimeh, tval, hpval)
                    ! TODO: The unemployment indicator
                    EVmatU => EVpol(:, 1, Dprimel:Dprimeh, aprimel:aprimeh, tval, hpval)
                else
                    EVmat => EV(:, 2, Dprimel:Dprimeh, aprimel:aprimeh, tval, hpval)
                    ! TODO: The unemployment indicator
                    EVmatU => EV(:, 1, Dprimel:Dprimeh, aprimel:aprimeh, tval, hpval)
                end if
                perm_expect(i) = DOT_PRODUCT(Probz(zindex,zmin:zmax),&
                    EV_linpol(weightaprimel, weightDprimel, zmin, zmax, EVmat))
            end do
 
            EVholder(1)=DOT_PRODUCT(Probtmp, perm_expect)
            call weightPrimeL(aprime_shift(tmpgridsize+1),anodes,weightaprimel,&
                              aprimel,aprimeh)
            EVholder(2:2)=EV_linpol(weightaprimel, weightDprimel,zmin,zmax,EVmatU)
            expectation = continuation(beta2, empprob, EVholder(1), EVholder(2))
 
        else
            zmin => zindex
            zmax => zindex
            call weightPrimeL(aprime_shift(tmpgridsize+1),anodes,weightaprimel,aprimel,aprimeh)
            if (expecting) then
                EVmatD => EVpol(:, 2, Dprimel:Dprimeh, aprimel:aprimeh, Tdie+1, hpval)
                ! TODO: The unemployment indicator
                EVmatU => EVpol(:, 1, Dprimel:Dprimeh, aprimel:aprimeh, tval, hpval)
                EVmat => EVpol(:, 2, Dprimel:Dprimeh, aprimel:aprimeh, tval, hpval)
            else
                EVmatD => EV(:, 2, Dprimel:Dprimeh, aprimel:aprimeh, Tdie+1, hpval)
                ! TODO: The unemployment indicator
                EVmatU => EV(:, 1, Dprimel:Dprimeh, aprimel:aprimeh, tval, hpval)
                EVmat => EV(:, 2, Dprimel:Dprimeh, aprimel:aprimeh, tval, hpval)
            end if
            EVholder(1:1)=EV_linpol(weightaprimel, weightDprimel,zmin,zmax,EVmat)
 
            if (tval-1 > Tretire) then
                ! When retired, no risk in income but risk in death
                EVholder(2:2)=EV_linpol(weightaprimel, weightDprimel,zmin,zmax,EVmatD)
                expectation = continuation(beta2retire, 1.0-deathprob(tval-Tretire-1),&
                                           EVholder(1), EVholder(2))
            else
                ! When unemployed, stochastic decision of staying in
                ! unemployment or rejoining at previous income level
                EVholder(2:2)=EV_linpol(weightaprimel, weightDprimel,INT(zindex),EVmatU)
                expectation = continuation(beta2, unempprob, EVholder(1), EVholder(2))
            end if

        end if
    END FUNCTION ! %>


    REAL(8) FUNCTION valfuncadjust(p, state, pol_opts) ! %< p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
        IMPLICIT NONE
        REAL(8), DIMENSION(2), TARGET :: p
        REAL(8), DIMENSION(8), INTENT(IN), TARGET :: state
        LOGICAL, DIMENSION(2), INTENT(IN) :: pol_opts
        REAL(8), POINTER :: a, D, hpholder, zindex, t, hpindex, balholder, unemp
        REAL(8), POINTER :: aprime, Dprime
        REAL(8) :: consumption, currentincome, hpcurrent
        REAL(8) :: laborsupply
        REAL(8), DIMENSION(2) :: finwealth
        REAL(8) :: EVholder
        REAL(8) :: assetholder, thetaholder, transferholder, scrapholder,&
                   transferfuture, transfer_down, dep

        !write(*,*) state

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

        aprime => p(1)
        Dprime => p(2)
        if (aprime>anodes(agridsize)) then
            aprime => anodes(agridsize)
        end if

        thetaholder=theta
        if (pol_opts(1)) then
            transferholder=adjtransfer
            scrapholder=scrapped
            transferfuture=(1.0-eta_transfer)*transferholder
            transfer_down = eta_transfer*transferholder/(exp(hpholder)*Dprime)
        else
            transferholder=0.0
            scrapholder=0.0
            transferfuture=transferholder
            transfer_down=transferholder
        end if
        currentincome = (unemp - 1.0)*(income(zindex, t) - inctax(zindex, t)) + balholder
        dep = maint*delta
        ! First element is current financial assets, 2nd is next period
        finwealth = (/ a-D*(1-theta)*exp(hpholder),&
                       aprime-(1-delta+dep+discountflag*transfer_down)*Dprime&
                       &*(1-theta)*exp(hpholder+hpdelta) /)

        if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. &
            Dprime>=Dmin .AND. Dprime<=Dmax) then
            ! ((t <= Tretire .AND. Dprime>=Dmin) .OR. (t > Tretire .AND. Dprime > 0)) .AND. Dprime<=Dmax) then
            if (finwealth(1) < 0.0) then
                assetholder = (1+rborrow)*finwealth(1)
            else
                assetholder = (1+r)*finwealth(1)
            end if
            ! TODO: Maybe have a contemporaneous transfer too
            ! (like mortgage interest tax credit?)
            consumption=(currentincome+(1.0-scrapholder)*D*(1-prF)*(1-dep-dtau)*exp(hpholder)&
                         +(1.0-downflag)*(1.0-discountflag)*(transferholder-transferfuture)&
                         -thetaholder*Dprime*(1+(F-prF))*exp(hpholder-downflag*(1.0-discountflag)*&
                         &transfer_down/thetaholder)+assetholder-aprime)

            if (consumption > 0) then
            ! TODO: what's the deal with the "+0"? is that supposed to reflect
            ! a minimum house size?
                valfuncadjust=(consumption+0)**elasticity2*(&
                        (1+discountflag*transfer_down)*Dprime+0)**(1-elasticity2)
                if (elasticity .ne. 1) then
                    valfuncadjust=utilscale(t)*(valfuncadjust)**(1-elasticity)/(1-elasticity)
                else
                    valfuncadjust=utilscale(t)*log(valfuncadjust)
                end if

                EVholder=expectation(INT(zindex), (1-delta+dep+discountflag*transfer_down)*Dprime,&
                                     aprime, t+1, hpindex+1, unemp, finwealth(2),&
                                     pol_opts(2), .FALSE., transferfuture)
                valfuncadjust= valfuncadjust + EVholder
            else
                valfuncadjust=-1.0e13
            end if
        else
             valfuncadjust=-1.0e13
        end if
        valfuncadjust=-valfuncadjust  !powell minimizes

    END FUNCTION valfuncadjust ! %>

    REAL(8) FUNCTION valfuncrent(p, state, pol_opts) ! %< p is choice variables optimized over (a' and D', C is residual from budget constraint), state is current state
        IMPLICIT NONE
        REAL(8), DIMENSION(2), TARGET :: p
        REAL(8), DIMENSION(7), INTENT(IN), TARGET :: state
        LOGICAL, DIMENSION(2), INTENT(IN) :: pol_opts
        REAL(8), POINTER :: a, D, hpholder, zindex, t, hpindex, balholder, unemp
        REAL(8), POINTER :: aprime, Dprime
        REAL(8) :: consumption, currentincome, hpcurrent
        REAL(8) :: laborsupply
        REAL(8), DIMENSION(2) :: finwealth
        REAL(8) :: EVholder
        REAL(8) :: assetholder, rentholder, transferholder, dep

        !write(*,*) state

        zindex => state(1)
        D => state(2)
        a => state(3)
        hpholder => state(4)
        hpindex => state(5)
        t => state(6)
        balholder => state(7)
        if (hpindex==0) then
            hpcurrent = 1.0
        else
            hpcurrent = hpindex
        end if

        if (pol_opts(1)) then
        transferholder=renttransfer
        else
        transferholder=0
        end if
        currentincome = (unemp - 1.0)*(income(zindex, t) - inctax(zindex, t)) + balholder
        dep = maint*delta

        aprime => p(1)
        Dprime => p(2)
        if (aprime>anodes(agridsize)) then
            aprime => anodes(agridsize)
        end if
        ! First element is current financial assets, 2nd is next period
        finwealth = (/ a-D*(1-theta)*exp(hpholder), aprime /)

        if (aprime>=borrowconstraint .AND. aprime<=anodes(agridsize) .AND. Dprime>=Dnodes(1) .AND. Dprime<=Dmax_rent) then
            if ((a-D*(1-theta)*exp(hpholder))<0) then
                assetholder = (1+rborrow)*a - (1+rborrow)*(1-theta)*D*exp(hpholder)
            else
                assetholder = (1+r)*a - (1+r)*(1-theta)*D*exp(hpholder)
            end if
            if (t <= Tretire) then
                rentholder=r_rental(hpcurrent)
            else
                rentholder=r_rental_retire(hpcurrent)
            end if
            ! TODO: Other tax credits for renters?
                consumption=(currentincome+D*(1-prF)*(1-dep-dtau)*exp(hpholder)&
                    +transferholder-rentholder*Dprime*exp(rentelasticity*hpholder)&
                    +assetholder-aprime)
            if (consumption > 0) then
            ! %< need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.  Find nearest grid points and linearly interpolate.
                valfuncrent=(consumption+0)**elasticity2*(rentUtil*Dprime+0)**(1-elasticity2)
                if (elasticity .ne. 1) then
                    valfuncrent=utilscale(t)*(valfuncrent)**(1-elasticity)/(1-elasticity)
                else
                    valfuncrent=utilscale(t)*log(valfuncrent)
                end if
                ! NOTE THE DIFFERENCE WITH THE ADJUST FUNCTION. We only interpolate
                ! over a because rental housing does not last for longer than one
                ! period. See indexing on EV as well.

                EVholder=expectation(INT(zindex), 0.0, aprime, t+1, hpindex+1,&
                                     unemp, finwealth(2), pol_opts(2), .TRUE.)
                valfuncrent= valfuncrent + EVholder
            ! %>
            else
                valfuncrent=-1.0e13
            end if
        else
             valfuncrent=-1.0e13
        end if
        valfuncrent=-valfuncrent  !powell minimizes

    END FUNCTION valfuncrent ! %>

    REAL(8) FUNCTION valfuncnoadjust(aprime, state, pol_opts)  ! %< aprime is a', (given D'=D, C is residual from budget constraint), state is current state
        IMPLICIT NONE
        REAL(8), INTENT(IN), TARGET :: aprime
        REAL(8), DIMENSION(7), INTENT(IN), TARGET :: state
        LOGICAL, INTENT(IN) :: pol_opts
        REAL(8), POINTER :: a, Dindex, hpholder, zindex, t, hpindex, balholder, unemp
        REAL(8), POINTER :: D
        REAL(8) :: consumption, aprime_adj, currentincome, hpcurrent
        REAL(8) :: laborsupply
        REAL(8), DIMENSION(2) :: finwealth
        REAL(8) :: EVholder
        REAL(8) :: transferholder, assetholder, dep

        !write(*,*) state

        zindex => state(1)
        Dindex => state(2)
        a => state(3)
        hpholder => state(4)
        hpindex => state(5)
        t => state(6)
        balholder => state(7)
        if (hpindex==0) then
            hpcurrent = 1.0
        else
            hpcurrent = hpindex
        end if

        transferholder=noadjtransfer_internal
        currentincome = (unemp-1.0)*(income(zindex, t) - inctax(zindex, t)) + balholder
        dep = maint*delta

        D => Dnodes(Dindex)
        aprime_adj = max(aprime+(1-theta)*D*hpdelta, borrowconstraint)
        if (aprime_adj>anodes(agridsize)) then
            aprime_adj = anodes(agridsize)
        end if
        ! First element is current financial assets, 2nd is next period
        finwealth = (/ a-D*(1-theta)*exp(hpholder),&
                       aprime-(1-delta+dep)*D*(1-theta)*exp(hpholder+hpdelta) /)


        if (aprime_adj>=borrowconstraint .AND. aprime_adj<=anodes(agridsize)) then
            if ((a-D*(1-theta)*exp(hpholder))<0) then
                assetholder = (1+rborrow)*a - (1+rborrow)*(1-theta)*D*exp(hpholder)
            else
                assetholder = (1+r)*a - (1+r)*(1-theta)*D*exp(hpholder)
            end if
                consumption=(currentincome+D*(1-dep-dtau)*exp(hpholder)&
                    +transferholder-theta*D*exp(hpholder)&
                    +assetholder-aprime_adj)
            if (consumption>0) then
            ! %< need to figure out Vadjust and Vnoadjust at offgrid points tomorrow.  Find nearest grid points and linearly interpolate.
                valfuncnoadjust=(consumption+0)**elasticity2*(D+0)**(1-elasticity2)
                if (elasticity .ne. 1) then
                    valfuncnoadjust=utilscale(t)*(valfuncnoadjust)**(1-elasticity)/(1-elasticity)
                else
                    valfuncnoadjust=utilscale(t)*log(valfuncnoadjust)
                end if
                ! NOTE THE DIFFERENCE WITH THE ADJUST FUNCTION. We interpolate
                ! over a AND d because depreciation could occur even w/out adjustment.
                EVholder=expectation(INT(zindex), (1-delta+dep)*D,&
                                     aprime, t+1, hpindex+1, unemp, finwealth(2),&
                                     pol_opts, .FALSE.)
                valfuncnoadjust= valfuncnoadjust + EVholder
            ! %>
            else
                valfuncnoadjust=-1.0e13
            end if
        else
             valfuncnoadjust=-1.0e13
        end if
        valfuncnoadjust=-valfuncnoadjust  !powell minimizes

    END FUNCTION valfuncnoadjust ! %>

end module lifecycle_vfuncs
