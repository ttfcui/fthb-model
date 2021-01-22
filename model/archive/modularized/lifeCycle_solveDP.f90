module lifecycle_solveDP

    USE share
    USE lifecycle_algs
    USE lifecycle_vfuncs
    USE OMP_LIB
    IMPLICIT NONE

    CONTAINS

    subroutine solveworkingproblem(steady, price)
        IMPLICIT NONE
        INTEGER :: thread_count, actualHP
        LOGICAL, INTENT(IN) :: steady
        REAL(8), INTENT(IN) :: price
        REAL(8), dimension(:,:,:,:,:), ALLOCATABLE :: optpolicynoadjust, optpolicyadjust
        REAL(8), dimension(:,:,:,:,:,:), ALLOCATABLE :: pstartadjust
        REAL(8), dimension(:,:,:,:,:), ALLOCATABLE :: ystartadjust
        REAL(8), dimension(:,:,:,:), ALLOCATABLE :: ax, bx, cx, adjust
        REAL(8), dimension(:,:,:,:,:), ALLOCATABLE :: state

        REAL(8), DIMENSION(2) :: p2
        REAL(8) :: ftol

        ftol=5e-9

        if (steady) THEN
           actualHP = hpgridsize
           hpnodes(hpgridsize) = exp(price)
        ELSE
           actualHP = 1
        END IF


    ! Allocation of automatic arrays for OpenMP (skip this)
    ALLOCATE(ax(agridsize, Dgridsize, zgridsize, hpgridsize),&
             bx(agridsize, Dgridsize, zgridsize, hpgridsize),&
             cx(agridsize, Dgridsize, zgridsize, hpgridsize),&
             adjust(agridsize, Dgridsize, zgridsize, hpgridsize))
    ALLOCATE(optpolicynoadjust(agridsize, Dgridsize, zgridsize, hpgridsize, 2),&
             optpolicyadjust(agridsize, Dgridsize, zgridsize, hpgridsize, 2))
    ALLOCATE(pstartadjust(agridsize, Dgridsize, zgridsize, hpgridsize, 3, 2),&
             ystartadjust(agridsize, Dgridsize, zgridsize, hpgridsize, 3),&
             state(agridsize, Dgridsize, zgridsize, hpgridsize, 5))
            


        adjust=0

        do t=Tdie, 1,-1
            write(*,*) "Processing year state:", t
    !$OMP PARALLEL
    !$OMP DO
            do i=1, agridsize
                do j=1, Dgridsize
                    do k=1, zgridsize
                        do l=actualHP, hpgridsize
                            pstartadjust(i, j, k, l, 1, 1)=.01*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 1, 2)=.8*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 2, 1)=.25*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 2, 2)=.25*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 3, 1)=.5*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 3, 2)=.2*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))


                            state(i, j, k, l, 1)=anodes(i)
                            state(i, j, k, l, 2)=Dnodes(j)
                            state(i, j, k, l, 3)=k
                            state(i, j, k, l, 4)=l
                            state(i, j, k, l, 5)=t

                            ystartadjust(i, j, k, l, 1)=valfuncadjust2(pstartadjust(i, j, k, l, 1,:), state(i, j, k, l,:))
                            ystartadjust(i, j, k, l, 2)=valfuncadjust2(pstartadjust(i, j, k, l, 2,:), state(i, j, k, l,:))
                            ystartadjust(i, j, k, l, 3)=valfuncadjust2(pstartadjust(i, j, k, l, 3,:), state(i, j, k, l,:))


                            call amoeba(state(i, j, k, l,:), pstartadjust(i, j, k, l,:,:), ystartadjust(i, j, k, l,:), ftol, valfuncadjust)
                            Vadjust(i, j, k, l, t)=ystartadjust(i, j, k, l, 1)
                            achoiceadjust(i, j, k, l, t)=pstartadjust(i, j, k, l, 1, 1)
                            Dchoiceadjust(i, j, k, l, t)=pstartadjust(i, j, k, l, 1, 2)


                            pstartadjust(i, j, k, l, 1, 1)=0.0
                            pstartadjust(i, j, k, l, 1, 2)=.01*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 2, 1)=.05*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 2, 2)=.1*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 3, 1)=.4*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 3, 2)=.21*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))

                            ystartadjust(i, j, k, l, 1)=valfuncadjust2(pstartadjust(i, j, k, l, 1,:), state(i, j, k, l,:))
                            ystartadjust(i, j, k, l, 2)=valfuncadjust2(pstartadjust(i, j, k, l, 2,:), state(i, j, k, l,:))
                            ystartadjust(i, j, k, l, 3)=valfuncadjust2(pstartadjust(i, j, k, l, 3,:), state(i, j, k, l,:))

                            call amoeba(state(i, j, k, l,:), pstartadjust(i, j, k, l,:,:), ystartadjust(i, j, k, l,:), ftol, valfuncadjust)
                            if (ystartadjust(i, j, k, l, 1)<Vadjust(i, j, k, l, t)) then  !(again we're minimizing)
                                Vadjust(i, j, k, l, t)=ystartadjust(i, j, k, l, 1)
                                achoiceadjust(i, j, k, l, t)=pstartadjust(i, j, k, l, 1, 1)
                                Dchoiceadjust(i, j, k, l, t)=pstartadjust(i, j, k, l, 1, 2)
                            end if










                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            pstartadjust(i, j, k, l, 1, 1)=.01*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 1, 2)=.8*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 2, 1)=.25*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 2, 2)=.25*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 3, 1)=.5*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 3, 2)=.2*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))

                            state(i, j, k, l, 1)=anodes(i)
                            state(i, j, k, l, 2)=Dnodes(j)
                            state(i, j, k, l, 3)=k
                            state(i, j, k, l, 4)=l
                            state(i, j, k, l, 5)=t

                            ystartadjust(i, j, k, l, 1)=valfuncrent2(pstartadjust(i, j, k, l, 1,:), state(i, j, k, l,:))
                            ystartadjust(i, j, k, l, 2)=valfuncrent2(pstartadjust(i, j, k, l, 2,:), state(i, j, k, l,:))
                            ystartadjust(i, j, k, l, 3)=valfuncrent2(pstartadjust(i, j, k, l, 3,:), state(i, j, k, l,:))

                            call amoeba(state(i, j, k, l,:), pstartadjust(i, j, k, l,:,:), ystartadjust(i, j, k, l,:), ftol, valfuncrent)


                            Vrent(i, j, k, l, t)=ystartadjust(i, j, k, l, 1)
                            achoicerent(i, j, k, l, t)=pstartadjust(i, j, k, l, 1, 1)
                            Dchoicerent(i, j, k, l, t)=pstartadjust(i, j, k, l, 1, 2)


                            pstartadjust(i, j, k, l, 1, 1)=0.0
                            pstartadjust(i, j, k, l, 1, 2)=.01*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 2, 1)=.05*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 2, 2)=.1*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 3, 1)=.4*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))
                            pstartadjust(i, j, k, l, 3, 2)=.21*((1+r)*anodes(i)+income(k, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)+(1-F)*(1-delta)*Dnodes(j))


                            ystartadjust(i, j, k, l, 1)=valfuncrent2(pstartadjust(i, j, k, l, 1,:), state(i, j, k, l,:))
                            ystartadjust(i, j, k, l, 2)=valfuncrent2(pstartadjust(i, j, k, l, 2,:), state(i, j, k, l,:))
                            ystartadjust(i, j, k, l, 3)=valfuncrent2(pstartadjust(i, j, k, l, 3,:), state(i, j, k, l,:))

                            call amoeba(state(i, j, k, l,:), pstartadjust(i, j, k, l,:,:), ystartadjust(i, j, k, l,:), ftol, valfuncrent)


                            if (ystartadjust(i, j, k, l, 1)<Vrent(i, j, k, l, t)) then  !(again we're minimizing)
                                Vrent(i, j, k, l, t)=ystartadjust(i, j, k, l, 1)
                                achoicerent(i, j, k, l, t)=pstartadjust(i, j, k, l, 1, 1)
                                Dchoicerent(i, j, k, l, t)=pstartadjust(i, j, k, l, 1, 2)
                            end if


                            if (costlyequity==1) then

                                if (anodes(i)>=(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))) then
                                    ax(i, j, k, l)=(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                else
                                    ax(i, j, k, l)=anodes(i)
                                end if
                            bx(i, j, k, l)=ax(i, j, k, l)
                            else
                            ax(i, j, k, l)=0.0
                            bx(i, j, k, l)=borrowconstraint
                            end if
                            cx(i, j, k, l)=amax

                            state(i, j, k, l, 2)=j
                            Vnoadjust(i, j, k, l, t)=brentnew(ax(i, j, k, l), bx(i, j, k, l), cx(i, j, k, l), valfuncnoadjust, ftol, state(i, j, k, l, 1), state(i, j, k, l, 2), state(i, j, k, l, 3), state(i, j, k, l, 4), state(i, j, k, l, 5), achoicenoadjust(i, j, k, l, t))
                            Dchoicenoadjust(i, j, k, l, t)=Dnodes(j)

                            ! Gideon: computing cchoiceadjust, cchoicenoadjust and cchoicerent
                            if (anodes(i)-(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))<0) then
                                    cchoiceadjust(i, j, k, l, t)=income(k, t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoiceadjust(i, j, k, l, t)*exp(hpnodes(l))-achoiceadjust(i, j, k, l, t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                    cchoicenoadjust(i, j, k, l, t)=income(k, t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoicenoadjust(i, j, k, l, t)*exp(hpnodes(l))-achoicenoadjust(i, j, k, l, t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                    cchoicerent(i, j, k, l, t)=income(k, t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoicerent(i, j, k, l, t)*exp(rentelasticity*hpnodes(l))-achoicerent(i, j, k, l, t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                            else
                                    cchoiceadjust(i, j, k, l, t)=income(k, t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoicenoadjust(i, j, k, l, t)*exp(hpnodes(l))-achoicenoadjust(i, j, k, l, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                    cchoicenoadjust(i, j, k, l, t)=income(k, t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoicenoadjust(i, j, k, l, t)*exp(hpnodes(l))-achoicenoadjust(i, j, k, l, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                    cchoicerent(i, j, k, l, t)=income(k, t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoicerent(i, j, k, l, t)*exp(rentelasticity*hpnodes(l))-achoicerent(i, j, k, l, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                            end if

                            if (Vadjust(i, j, k, l, t)<Vnoadjust(i, j, k, l, t) .and. Vadjust(i, j, k, l, t)<Vrent(i, j, k, l, t)) then  ! since V = - V from minimization
                                achoice(i, j, k, l, t)=achoiceadjust(i, j, k, l, t)
                                Dchoice(i, j, k, l, t)=Dchoiceadjust(i, j, k, l, t)
                                if (anodes(i)-(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))<0) then
                                    cchoice(i, j, k, l, t)=income(k, t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i, j, k, l, t)*exp(hpnodes(l))-achoice(i, j, k, l, t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                    !write(*,*) "diff1", income(k, t)+(1+rborrow*exp(hpnodes(l))**)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i, j, k, l, t)*exp(hpnodes(l))-achoice(i, j, k, l, t)-(1+rborrow*exp(hpnodes(l))**)*(1-thetanodes(l))*(1-F)*Dnodes(j)*exp(hpnodes(l))-(income(k, t)+(1+rborrow*exp(hpnodes(l))**(-0.0))*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i, j, k, l, t)*exp(hpnodes(l))-achoice(i, j, k, l, t)-(1+rborrow*exp(hpnodes(l))**(-0.0))*(1-thetanodes(l))*(1-F)*Dnodes(j)*exp(hpnodes(l)))
                                else
                                    cchoice(i, j, k, l, t)=income(k, t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i, j, k, l, t)*exp(hpnodes(l))-achoice(i, j, k, l, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                end if
                                rentalindicator(i, j, k, l, t)=0
                                choiceindicator(i, j, k, l, t)=1
                            elseif (Vnoadjust(i, j, k, l, t)<Vadjust(i, j, k, l, t) .and. Vnoadjust(i, j, k, l, t)<Vrent(i, j, k, l, t)) then
                                achoice(i, j, k, l, t)=achoicenoadjust(i, j, k, l, t)
                                Dchoice(i, j, k, l, t)=Dchoicenoadjust(i, j, k, l, t)


                                if (anodes(i)-(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))<0) then
                                    cchoice(i, j, k, l, t)=income(k, t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i, j, k, l, t)*exp(hpnodes(l))-achoice(i, j, k, l, t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                    !write(*,*) "diff2", income(k, t)+(1+rborrow*exp(hpnodes(l))**)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i, j, k, l, t)*exp(hpnodes(l))-achoice(i, j, k, l, t)-(1+rborrow*exp(hpnodes(l))**)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))-(income(k, t)+(1+rborrow*exp(hpnodes(l))**(-0.0))*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i, j, k, l, t)*exp(hpnodes(l))-achoice(i, j, k, l, t)-(1+rborrow*exp(hpnodes(l))**(-0.0))*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l)))
                                else
                                    cchoice(i, j, k, l, t)=income(k, t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*exp(hpnodes(l))-thetanodes(l)*Dchoice(i, j, k, l, t)*exp(hpnodes(l))-achoice(i, j, k, l, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                end if
                                rentalindicator(i, j, k, l, t)=0
                                choiceindicator(i, j, k, l, t)=2
                            else
                                achoice(i, j, k, l, t)=achoicerent(i, j, k, l, t)
                                Dchoice(i, j, k, l, t)=Dchoicerent(i, j, k, l, t)
                                if (anodes(i)-(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))<0) then
                                    cchoice(i, j, k, l, t)=income(k, t)+(1+rborrow)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoice(i, j, k, l, t)*exp(rentelasticity*hpnodes(l))-achoice(i, j, k, l, t)-(1+rborrow)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                    !write(*,*) "diff3", income(k, t)+(1+rborrow*exp(hpnodes(l))**)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoice(i, j, k, l, t)*exp(hpnodes(l))-achoice(i, j, k, l, t)-(1+rborrow*exp(hpnodes(l))**)*(1-thetanodes(l))*(1-F)*Dnodes(j)*exp(hpnodes(l)) - (income(k, t)+(1+rborrow*exp(hpnodes(l))**(-0.0))*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoice(i, j, k, l, t)*exp(hpnodes(l))-achoice(i, j, k, l, t)-(1+rborrow*exp(hpnodes(l))**(-0.0))*(1-thetanodes(l))*(1-F)*Dnodes(j)*exp(hpnodes(l)))
                                else
                                    cchoice(i, j, k, l, t)=income(k, t)+(1+r)*anodes(i)+Dnodes(j)*(1-delta)*(1-F)*exp(hpnodes(l))-r_rental*Dchoice(i, j, k, l, t)*exp(rentelasticity*hpnodes(l))-achoice(i, j, k, l, t)-(1+r)*(1-thetanodes(l))*Dnodes(j)*exp(hpnodes(l))
                                end if
                                !write(*,*) l,(1+rborrow*exp(hpnodes(l))**)
                                rentalindicator(i, j, k, l, t)=1
                                choiceindicator(i, j, k, l, t)=3
                            end if
                           ! write(*,*) exp(hpnodes(l)), (1+rborrow*exp(hpnodes(l))**), income(6, l, 1)
                        end do
                    end do
                end do
            end do
    !$OMP END DO
    !$OMP END PARALLEL
            Vnoadjust(:,:,:,:, t)=-Vnoadjust(:,:,:,:, t)
            Vadjust(:,:,:,:, t)=-Vadjust(:,:,:,:, t)
            Vrent(:,:,:,:, t)=-Vrent(:,:,:,:, t)

            EV(:,:,:,:, t)=max(Vrent(:,:,:,:, t), max(Vnoadjust(:,:,:,:, t), Vadjust(:,:,:,:, t)))

           ! write(*,*) EV(1, 1, 2, 1, 2), EV(2, 1, 2, 1, 2)

            if (t==1) then
            write(*,*) "policy functions:"
            write(*,*) achoicerent(1, 1, 2, actualHP, 1), dchoicerent(1, 1, 2, actualHP, 1)
            write(*,*) achoiceadjust(1, 1, 2, actualHP, 1), dchoiceadjust(1, 1, 2, actualHP, 1)
            write(*,*) achoicenoadjust(1, 1, 2, actualHP, 1), dchoicenoadjust(1, 1, 2, actualHP, 1)
            write(*,*) achoice(1, 1, 2, actualHP, 1), dchoice(1, 1, 2, actualHP, 1), cchoice(1, 1, 2, actualHP, 1)

                             state(1, 1, 1, 1, 1)=anodes(1)
                            state(1, 1, 1, 1, 2)=Dnodes(1)
                            state(1, 1, 1, 1, 3)=2
                            state(1, 1, 1, 1, 4)=1
                            state(1, 1, 1, 1, 5)=1
                            pstartadjust(1, 1, 1, 1, 1, 1)=achoicerent(1, 1, 2, 1, 1)
                            pstartadjust(1, 1, 1, 1, 1, 2)=dchoicerent(1, 1, 2, 1, 1)

                            p2(1)=pstartadjust(1, 1, 2, 1, 1, 1)
                            p2(2)=pstartadjust(1, 1, 2, 1, 1, 2)
            ! write(*,*) valfuncrent3(p2, state(1, 1, 1, 1,:))
            ! p2(1)=achoicerent(1, 1, 2, 1, 1)+.000001
            ! write(*,*) valfuncrent3(p2, state(1, 1, 1, 1,:))
            end if

            if (priceshock==1) then
            cchoiceshock=cchoice
            Dchoiceshock=Dchoice
            end if
    end do

    DEALLOCATE(ax, bx, cx, adjust)

    end subroutine solveworkingproblem

    subroutine plotpolicy_constantwealth(p)
        USE OMP_LIB
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: p
        REAL(8) :: wealth, h, q, achoicelin, Dchoicelin, cchoicelin, rentallin, choicelin, welfare

        OPEN (UNIT=266, FILE="policy_constantwealth.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        266 format (I6.2, 99(F16.6))

        t=3
        do k=1, zgridsize
        do i=1, 50
            wealth=0+.1*i
            do j=1, 50
                h=0+(j-1)*.1
                q=wealth-theta*h
                call pol_linworking(q, h, REAL(k), p, REAL(t), achoicelin, Dchoicelin, cchoicelin, rentallin, choicelin, welfare)
                write(266, 266) k, wealth, h, Dchoicelin, cchoicelin, achoicelin
            end do
        end do
        end do

    end subroutine plotpolicy_constantwealth

    subroutine output_vfuncs(steady, price)
        IMPLICIT NONE
        LOGICAL, INTENT(IN) :: steady
        REAL(8), INTENT(IN) :: price
        INTEGER :: k2, t2, l2, actualHP

        OPEN (UNIT=1, FILE="vfunc1.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        1 format (99(F40.12))

        OPEN (UNIT=17, FILE="vfunc2.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        17 format (99(F40.12))

        OPEN (UNIT=18, FILE="vfunc3.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        18 format (99(F40.12))

        OPEN (UNIT=19, FILE="vfunc4.txt", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
        19 format (99(F40.12))

        if (steady) THEN
           actualHP = hpgridsize
           hpnodes(1) = exp(price)
        ELSE
           actualHP = 1
        END IF


        write(*,*) agridsize, Dgridsize, zgridsize, actualHP
        t=0
        do i=1, agridsize
            do j=1, Dgridsize
                do k=1, zgridsize
                    do l=actualHP,hpgridsize
                    t=t+1
                    t2=33
                    k2=k
                    l2=l
                        !call computempc(anodes(i), Dnodes(j), k2, l2, t2, mpc_pol(i, j, k, l, 33))
                        write(1, 1) EV(i, j, k, l, 33), achoice(i, j, k, l, 33), dchoice(i, j, k, l, 33), cchoice(i, j, k, l, 33), rentalindicator(i, j, k, l, 33), anodes(i), Dnodes(j), znodes(k), hpnodes(l), income(k, 1) !, mpc_pol(i, j, k, l, 33) , Vadjust(i, j, k, l, 1), Vnoadjust(i, j, k, l, 1)
                    end do
                end do
            end do
      end do

      t=0
        do i=1, agridsize
            do j=1, Dgridsize
                do k=1, zgridsize
                    do l=actualHP,hpgridsize
                    t=t+1
                        write(17, 17) EV(i, j, k, l, 34), achoice(i, j, k, l, 34), dchoice(i, j, k, l, 34), cchoice(i, j, k, l, 34), rentalindicator(i, j, k, l, 34), anodes(i), Dnodes(j), znodes(k), hpnodes(l), income(k, 1) !, Vadjust(i, j, k, l, 1), Vnoadjust(i, j, k, l, 1)
                    end do
                end do
            end do
      end do

      t=0
        do i=1, agridsize
            do j=1, Dgridsize
                do k=1, zgridsize
                    do l=actualHP,hpgridsize
                    t=t+1
                        write(18, 18) EV(i, j, k, l, 35), achoice(i, j, k, l, 35), dchoice(i, j, k, l, 35), cchoice(i, j, k, l, 35), rentalindicator(i, j, k, l, 35), anodes(i), Dnodes(j), znodes(k), hpnodes(l), income(k, 1) !, Vadjust(i, j, k, l, 1), Vnoadjust(i, j, k, l, 1)
                    end do
                end do
            end do
      end do

      t=0
        do i=1, agridsize
            do j=1, Dgridsize
                do k=1, zgridsize
                    do l=actualHP,hpgridsize
                    t=t+1
                        write(19, 19) EV(i, j, k, l, 36), achoice(i, j, k, l, 36), dchoice(i, j, k, l, 36), cchoice(i, j, k, l, 36), rentalindicator(i, j, k, l, 36), anodes(i), Dnodes(j), znodes(k), hpnodes(l), income(k, 1) !, Vadjust(i, j, k, l, 1), Vnoadjust(i, j, k, l, 1)
                    end do
                end do
            end do
      end do



      write(*,*) t
    end subroutine output_vfuncs

end module lifecycle_solveDP
