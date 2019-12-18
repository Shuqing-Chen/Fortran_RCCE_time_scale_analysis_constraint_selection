subroutine time_scale(ifPrint,htr,nSpecs,nrxn,phi,Temp,press,&
                      time_scale_rxn_f,time_scale_rxn_r, time_scale_rxn,&
                      stoichcoeff_f,stoichcoeff_b)

  use mod_hectars
    
    logical, intent(in) :: ifPrint
    integer, intent(in) :: nSpecs,nrxn
    type(hectars), intent(inout) :: htr
    real(8), intent(inout) :: phi(nSpecs)
    real(8), intent(inout) :: Temp
    real(8), intent(in) :: press
    real(8), intent(out) :: time_scale_rxn_f(nrxn)
    real(8), intent(out) :: time_scale_rxn_r(nrxn)
    real(8), intent(out) :: time_scale_rxn(nrxn)
    integer, intent(out) :: stoichcoeff_f(nSpecs,nrxn)
    integer, intent(out) :: stoichcoeff_b(nSpecs,nrxn)
    integer :: stoichcoeff(nSpecs,nrxn)
    real(8) :: FWDK(nrxn) 
    real(8) :: REVK(nrxn)
    real(8) :: NETK(nrxn)
    real(8) :: k_f(nrxn) 
    real(8) :: k_b(nrxn)    
    real(8) :: c(nSpecs) 
    real(8) :: Ru
    real(8) :: mlrMassAvg
    real(8) :: MWs(nSpecs)
    real(8) :: Jacobian_f(nSpecs,nrxn)
    real(8) :: Jacobian_b(nSpecs,nrxn)
    real(8) :: Jacobian(nSpecs,nrxn)
    real(8) :: wdot(nSpecs+1)
    real(8) :: wdot2(nSpecs+1)
    real(8) :: q_f(nrxn) 
    real(8) :: q_b(nrxn)    
    integer :: i,j,k
    
    stoichcoeff_f=0.d0
    stoichcoeff_b=0.d0
    stoichcoeff=0.d0
    do i=1,nrxn
    do j=1,nSpecs
       stoichcoeff_f(j,i) = reactantStoichCoeff(htr,j,i)
       stoichcoeff_b(j,i) = productStoichCoeff(htr,j,i)
       stoichcoeff(j,i) = stoichcoeff_b(j,i) - stoichcoeff_f(j,i)
!        if ( (stoichcoeff_f(j,i)/=0) .AND. (stoichcoeff_b(j,i)/=0) ) then
!        if (stoichcoeff_f(j,i)==stoichcoeff_b(j,i)) then
!          print*,"stoichcoeff_f(j,i)=stoichcoeff_b(j,i)",i,j
!        endif
!        endif
    enddo
    enddo
    
    call getFwdRatesOfProgress(htr,FWDK)
    call getRevRatesOfProgress(htr,REVK)
    call getNetRatesOfProgress(htr,NETK)

!    print*,"FWDK",FWDK
!    print*,"REVK",REVK
!    print*,"NETK1",FWDK-REVK
!    print*,"NETK2",NETK
        
    call setState_TPY(htr, Temp, press, phi(1:nSpecs))

    Ru = Runiv(htr)
    mlrMassAvg = meanMolecularWeight(htr)
    call getMolecularWeights(htr,MWs)
            
    do j=1,nSpecs
       c(j) = phi(j)*mlrMassAvg*press/(MWs(j)*Ru*Temp)
    end do  

    k_f=FWDK
    k_b=REVK
    do i=1,nrxn
    do j=1,nSpecs
       if (stoichcoeff_f(j,i)==1) then 
          k_f(i) = k_f(i)/c(j)
       elseif (stoichcoeff_f(j,i)==2) then
          k_f(i) = k_f(i)/(c(j)*c(j))
       endif
           
       if (stoichcoeff_b(j,i)==1) then 
          k_b(i) = k_b(i)/c(j)
       elseif (stoichcoeff_b(j,i)==2) then
          k_b(i) = k_b(i)/(c(j)*c(j))
       endif       
    end do
    end do 
    
!    print*,"k_f",k_f
!    print*,"k_b",k_b   
    
    Jacobian_f=0.d0
    Jacobian_b=0.d0
    
    do i=1,nrxn
    do j=1,nSpecs
        !if (stoichcoeff(j,i)==0) then
        !   Jacobian_f(j,i) = 0.d0
        !   Jacobian_b(j,i) = 0.d0
        !elseif (stoichcoeff(j,i)<0) then
        if (stoichcoeff_f(j,i)/=0) then
           Jacobian_f(j,i) = k_f(i)
           if (stoichcoeff_f(j,i)==1) then 
              Jacobian_f(j,i) = Jacobian_f(j,i)*stoichcoeff_f(j,i)
           elseif (stoichcoeff_f(j,i)==2) then
              Jacobian_f(j,i) = Jacobian_f(j,i)*stoichcoeff_f(j,i)*c(j)
           endif
           do kk=1,nSpecs
              if (stoichcoeff_f(kk,i)/=0) then
                 if (kk/=j) then 
                    if (stoichcoeff_f(kk,i)==1) then 
                       Jacobian_f(j,i) = Jacobian_f(j,i)*c(kk)
                    elseif (stoichcoeff_f(kk,i)==2) then
                       Jacobian_f(j,i) = Jacobian_f(j,i)*c(kk)*c(kk)
                    endif
                 endif
              endif
           enddo
        endif
        !else
        if (stoichcoeff_b(j,i)/=0) then
           Jacobian_b(j,i) = k_b(i)
           if (stoichcoeff_b(j,i)==1) then 
              Jacobian_b(j,i) = Jacobian_b(j,i)*stoichcoeff_b(j,i)
           elseif (stoichcoeff_b(j,i)==2) then
              Jacobian_b(j,i) = Jacobian_b(j,i)*stoichcoeff_b(j,i)*c(j)
           endif
           do kk=1,nSpecs
              if (stoichcoeff_b(kk,i)/=0) then
                 if (kk/=j) then 
                    if (stoichcoeff_b(kk,i)==1) then 
                       Jacobian_b(j,i) = Jacobian_b(j,i)*c(kk)
                    elseif (stoichcoeff_b(kk,i)==2) then
                       Jacobian_b(j,i) = Jacobian_b(j,i)*c(kk)*c(kk)
                    endif
                 endif
              endif
           enddo
        endif        
    end do
    end do 

    time_scale_rxn_f=0.d0
    time_scale_rxn_r=0.d0
    time_scale_rxn=0.d0  
      
    do i=1,nrxn
       do j=1,nSpecs
           time_scale_rxn_f(i)=time_scale_rxn_f(i)+Jacobian_f(j,i)*stoichcoeff_f(j,i)
           time_scale_rxn_r(i)=time_scale_rxn_r(i)+Jacobian_b(j,i)*stoichcoeff_b(j,i)
        enddo
    end do

    do i=1,nrxn
       do j=1,nSpecs
           Jacobian(j,i)=Jacobian_f(j,i)-Jacobian_b(j,i)
       enddo       
       do j=1,nSpecs
           time_scale_rxn(i)=time_scale_rxn(i)+Jacobian(j,i)*stoichcoeff(j,i)
       enddo
    end do
    
    do i=1,nrxn
        time_scale_rxn_f(i)=1.d0/time_scale_rxn_f(i)
        time_scale_rxn_r(i)=1.d0/time_scale_rxn_r(i)
        time_scale_rxn(i)=1.d0/time_scale_rxn(i)
        if (time_scale_rxn(i)<0) then
           time_scale_rxn(i)=-1.0d0*time_scale_rxn(i)
        endif
    end do
    
!    write (*,*) "c",c
!    write (*,*)
!    write (*,*) "Jacobian",Jacobian_b(:,3) 
    
    if (ifPrint .eqv. .true.) then
       do i=nrxn-20,nrxn
          print*,time_scale_rxn_f(i),time_scale_rxn_r(i)
       enddo
       open(unit=4000, file = 'stoichcoeff.txt')
       open(unit=5000, file = 'reationK.txt')
       open(unit=6000, file = 'Jacobian.txt')
       do i=1,nSpecs
          write(4000,'(133(i5)/)') stoichcoeff_f(i,:)
       enddo 
       write(4000,*)
       do i=1,nSpecs
          write(4000,'(133(i5)/)') stoichcoeff_b(i,:)
       enddo
        
       write(5000,'(133(d25.17)/)')FWDK
       write(5000,*)
       write(5000,'(133(d25.17)/)')REVK
       
       do i=1,nSpecs
          write(6000,'(133(d25.17)/)') Jacobian_f(i,:)
       enddo 
       write(6000,*)      
       do i=1,nSpecs
          write(6000,'(133(d25.17)/)') Jacobian_b(i,:)
       enddo
       close(4000)
       close(5000)
       close(6000)
    endif  

    call get_wdot_hectars(htr,phi,press,wdot) 
    
    wdot2=0.d0
    do j=1,nSpecs
       do i=1,nrxn
          wdot2(j)= wdot2(j)+stoichcoeff(j,i)*NETK(i)
       enddo
    enddo
    
    if (ifPrint .eqv. .true.) then
       do j=1,nSpecs
          write (*,*) "wdot",wdot(j),wdot2(j)
       enddo
       write (*,*)
       do i=1,nrxn
          write (*,*) "NETK",NETK(i),FWDK(i)-REVK(i) 
       enddo       
       !stop
    endif          
end subroutine time_scale
