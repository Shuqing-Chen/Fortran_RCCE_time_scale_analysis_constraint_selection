module constraints_mod
   use mod_hectars
   use Matrix_Mod
   implicit none
   private
   real(8) :: time
            
   public :: constraints

contains
!-----------------------------------------------
subroutine constraints(htr,nfast_rxns,fast_rxns)
  
  type(hectars), intent(inout) :: htr 
  integer, intent(inout) :: nfast_rxns
  integer, intent(inout) :: fast_rxns(nfast_rxns)
  integer :: i,j
  integer :: nelmnt, nspecs, nrxn
  character(5) :: elmntstr
  character(20) :: specsstr
  character(80) :: rxnstr
  integer, allocatable :: n_dot(:,:)
  integer, allocatable  :: new_species(:)
  integer :: nnew_species,nnew_species_local
  integer, allocatable  :: new_elements(:)
  integer, allocatable  :: j_new_elements(:)
  integer :: nnew_elements,nnew_elements_local
  integer, allocatable :: C_standard(:,:)
  integer :: nC_standard
  integer, allocatable :: dependent_C_standard(:)
  integer :: ndependent_C_standard
  integer, allocatable :: C_basic(:,:)
  integer :: nC_basic
  integer, allocatable :: C_basic_mole(:)
  integer :: is_C_basic_mole
  integer :: n_required_constraints
  integer :: nindependent_fast_rxns
  integer :: independent_fast_rxns(nfast_rxns)
    
  nelmnt = nElements(htr)
  nspecs = nSpecies(htr)
  nrxn = nReactions(htr)
  
  print*, nrxn
  do i=1,nelmnt
     elmntstr=getElementName(htr,i)
     print*,elmntstr
  enddo  

  do j=1,nSpecs
     specsstr=getSpeciesName(htr,j)
     print*,j,specsstr
  enddo
    
  do i=1,nrxn
     call getReactionString(htr,i,rxnstr)
!     print*,rxnstr
  enddo 
   
  allocate(n_dot(nspecs,nrxn))
  allocate(new_species(nSpecs))
  allocate(C_standard(nrxn,nSpecs))
  allocate(dependent_C_standard(nrxn))
  allocate(C_basic(nelmnt,nSpecs))
  allocate(C_basic_mole(nSpecs))
  allocate(new_elements(nelmnt))
  allocate(j_new_elements(nelmnt))
  
  do j=1,nSpecs
     do i=1,nrxn
        n_dot(j,i)=productStoichCoeff(htr,j,i)-reactantStoichCoeff(htr,j,i)
     enddo
  enddo

  ! nfast_rxns = 13
  ! fast_rxns = (/ 19, 23, 34, 50, 51, 53, 55, 61, 65, 77, 96, 104, 106 /)
  
  nC_standard=0
  C_standard=0
  nnew_species=0
  nC_basic=0
  C_basic=0
  is_C_basic_mole=0
  nnew_elements=0
  ndependent_C_standard=0
  nindependent_fast_rxns=nfast_rxns
  
  do i=1,nfast_rxns

     write(*,'(/,"Fast reaction index = ",i5)') fast_rxns(i)
     print*,"i",i
     call identify_new_species(htr,i,nfast_rxns,fast_rxns,nSpecs,&
     nnew_species,nnew_species_local,new_species)
     if (nnew_species_local>0) then
        call identify_new_elements(htr,nelmnt,nSpecs,nnew_species,&
             nnew_species_local,new_species,nnew_elements,nnew_elements_local,&
             new_elements,j_new_elements) 
     endif  
     if (nC_standard>0 ) then        
        call identify_dependent_C_standard(nrxn,nSpecs,n_dot,nC_standard,&
             C_standard,ndependent_C_standard,dependent_C_standard,fast_rxns(i))
     endif
     n_required_constraints = nnew_species_local - 1
     print*,"n_required_constraints",n_required_constraints
     
     
     if ( (nnew_species_local>0) ) then 
        independent_fast_rxns(i-(nfast_rxns-nindependent_fast_rxns)) = &
        fast_rxns(i)  
        call update_C_basic(htr,nelmnt,nSpecs,nnew_elements,new_elements,&
             nC_basic,C_basic,is_C_basic_mole,C_basic_mole,nnew_species,&
             nnew_species_local,new_species)       
        if (nnew_elements_local>0) then
           call creat_C_basic(htr,nelmnt,nSpecs,nnew_elements_local,&
                nnew_elements,new_elements,j_new_elements,nnew_species,&
                nnew_species_local,new_species,C_basic,nC_basic)
                n_required_constraints=n_required_constraints-nnew_elements_local
           print*,"n_required_constraints",n_required_constraints
        endif
        if (n_required_constraints > 0) then
        if (is_C_basic_mole==0) then
           call creat_C_basic_mole(nSpecs,is_C_basic_mole,C_basic_mole,&
                nnew_species,nnew_species_local,new_species)
           n_required_constraints = n_required_constraints - 1
           print*,"n_required_constraints",n_required_constraints
        endif
        endif
        if ( n_required_constraints > 0 ) then
        if (nnew_species_local>1) then
           call create_C_standard(i,nfast_rxns,fast_rxns,nSpecs,nrxn,n_dot,&
           nnew_species,nnew_species_local,new_species,nC_standard,C_standard)
           n_required_constraints=n_required_constraints-(nnew_species_local-1)
           print*,"n_required_constraints",n_required_constraints
           if (n_required_constraints<0) then
              C_standard(nC_standard,:)=0
              nC_standard = nC_standard - 1
              print*,"nC_standard",nC_standard
           endif           
        endif
        endif
        if ( (ndependent_C_standard>0) ) then
            call update_C_standard(i,nfast_rxns,fast_rxns,nSpecs,nrxn,n_dot,&
                 nnew_species,nnew_species_local,new_species,nC_standard,&
                 C_standard,ndependent_C_standard,dependent_C_standard)
        endif


     elseif ( (nnew_species_local==0) .AND. (ndependent_C_standard>0) ) then
        independent_fast_rxns(i-(nfast_rxns-nindependent_fast_rxns))=fast_rxns(i)        
        call update_C_standard_2(i,nfast_rxns,fast_rxns,nSpecs,nrxn,n_dot,&
             nnew_species,nnew_species_local,new_species,nC_standard,&
             C_standard,ndependent_C_standard,dependent_C_standard)

     elseif  ( (nnew_species_local==0) .AND. (ndependent_C_standard==0) ) then     
        nindependent_fast_rxns = nindependent_fast_rxns - 1 
        print*,"nindependent_fast_rxns",nindependent_fast_rxns

        
     endif
     print*,"C_standard = "
     write(*,'(22(i4,1X))') (C_standard(j,:),j=1,nC_standard)
 
  enddo

  print*,"nC_standard",nC_standard
  print*,"nnew_species",nnew_species
  print*,"new_species = "
  !write(*,'(i5,1x,A)')((new_species(j),getSpeciesName(htr,new_species(j))),j=1,nnew_species)
  print*,"C_basic = "
  write(*,'(22(i4,1X))') (C_basic(i,:),i=1,nC_basic)
  print*,"C_basic_mold = "
  write(*,'(22(i4,1X))') C_basic_mole(:)
  print*,"C_standard = "
  write(*,'(22(i4,1X))') (C_standard(i,:),i=1,nC_standard)
  print*,"nindependent_fast_rxns",nindependent_fast_rxns 
  do i=1,nindependent_fast_rxns
     call getReactionString(htr,independent_fast_rxns(i),rxnstr) 
     write(*,'(i5,1x,A)') independent_fast_rxns(i),rxnstr
  enddo
  
  open(unit=7000, file = 'constraints.txt')
  write(7000,'(22(i3,1X))') (C_basic(i,:),i=1,nC_basic)
  write(7000,'(22(i3,1X))') C_basic_mole(:)
  write(7000,'(22(i3,1X))') (C_standard(i,:),i=1,nC_standard)
  close(7000)

end subroutine constraints
!----------------------------------------------  
subroutine identify_new_species(htr,i,nfast_rxns,fast_rxns,nSpecs,nnew_species,&
nnew_species_local,new_species)
  
  type(hectars), intent(inout) :: htr 
  integer, intent(in) :: i
  integer, intent(in) :: nfast_rxns
  integer, intent(in) :: fast_rxns(nfast_rxns)
  integer, intent(in) :: nSpecs  
  integer, intent(inout) :: nnew_species
  integer, intent(out) :: nnew_species_local
  integer, intent(inout) :: new_species(nSpecs)  !!!!! see
  integer :: j,k 
  integer :: is_new
  
  nnew_species_local=0
  do j=1,nspecs
     if (reactantStoichCoeff(htr,j,fast_rxns(i))/=0) then
        is_new=1
        if (nnew_species>0) then
           do k=1,nnew_species
              if (new_species(k)==j) then
                 is_new=0
              endif
           enddo
        endif
        if (is_new==1) then
           nnew_species=nnew_species+1
           nnew_species_local=nnew_species_local+1
           new_species(nnew_species)=j
           print*,"new_Specs",j
        endif
     endif
     if (productStoichCoeff(htr,j,fast_rxns(i))/=0) then
        is_new=1
        if (nnew_species>0) then
           do k=1,nnew_species
              if (new_species(k)==j) then
                 is_new=0
              endif
           enddo
        endif
        if (is_new==1) then
           nnew_species=nnew_species+1
           nnew_species_local=nnew_species_local+1
           new_species(nnew_species)=j
           print*,"new_Specs",j
        endif
     endif        
  enddo
  print*,"nnew_species_local",nnew_species_local
end subroutine identify_new_species
!---------------------------------------------- 
subroutine identify_new_elements(htr,nelmnt,nSpecs,nnew_species,&
nnew_species_local,new_species,nnew_elements,nnew_elements_local,&
new_elements,j_new_elements)
type(hectars), intent(inout) :: htr
integer, intent(in) :: nelmnt,nSpecs
integer, intent(in) :: nnew_species
integer, intent(in) :: nnew_species_local
integer, intent(inout) :: new_species(nSpecs) 
integer, intent(inout) :: nnew_elements
integer, intent(out) :: nnew_elements_local
integer, intent(inout) :: new_elements(nelmnt)
integer, intent(out) :: j_new_elements(nelmnt)
integer :: i,j,k,jj,j_new_species
integer :: is_new
  
nnew_elements_local=0
do j=1,nnew_species_local
   do i=1,nelmnt
      is_new=1
      if (nAtoms(htr,new_species(nnew_species-nnew_species_local+j),i)/=0) then
         if (nnew_elements>0) then
            do k=1,nnew_elements
               if (new_elements(k)==i) then
                  is_new=0
               endif
            enddo
         endif
      else
         is_new=0
      endif
      if (is_new==1) then
         nnew_elements=nnew_elements+1
         nnew_elements_local=nnew_elements_local+1
         new_elements(nnew_elements)=i
         j_new_elements(nnew_elements)=&
         new_species(nnew_species-nnew_species_local+j)
         print*,"Specs",new_species(nnew_species-nnew_species_local+j),"elmnt",i    
!         j_new_species=new_species(nnew_species-nnew_species_local+j)
!         do  jj=nnew_species-nnew_species_local+j+1,nnew_species
!            new_species(jj-1)=new_species(jj)
!         enddo
!         new_species(nnew_species)=j_new_species
!         print*,"new_species"
!         do jj=1,nnew_species
!            print*,new_species(jj)
!         enddo
      endif
   enddo
enddo
print*,"nnew_elements_local",nnew_elements_local   
end subroutine identify_new_elements
!----------------------------------------------
subroutine identify_dependent_C_standard(nrxn,nSpecs,n_dot,nC_standard,&
C_standard,ndependent_C_standard,dependent_C_standard,fast_rxns_index)

integer, intent(in) :: nrxn,nSpecs
integer, intent(in) :: n_dot(nSpecs,nrxn)
integer, intent(in) :: nC_standard
integer, intent(in) :: C_standard(nrxn,nSpecs)
integer, intent(in) :: fast_rxns_index
integer, intent(out) :: ndependent_C_standard
integer, intent(out) :: dependent_C_standard(nC_standard)
integer :: C_standard_omega(nrxn,nrxn)
integer :: i,j,k
integer :: count

call C_standard_to_omega(nSpecs,nrxn,nC_standard,n_dot,C_standard,&
C_standard_omega)

count = 0
do i=1,nC_standard
   if (C_standard_omega(i,fast_rxns_index)/=0) then
      count = count +1 
      !dependent_C_standard(count)=C_standard_omega(i,fast_rxns_index)
      dependent_C_standard(count)=i
      print*,"dependent_C_standard",i
   endif
enddo
ndependent_C_standard=count
print*,"ndependent_C_standard",ndependent_C_standard
end subroutine identify_dependent_C_standard
!---------------------------------------------
subroutine update_C_basic(htr,nelmnt,nSpecs,nnew_elements,new_elements,&
nC_basic,C_basic,is_C_basic_mole,&
           C_basic_mole,nnew_species,nnew_species_local,new_species)
type(hectars), intent(inout) :: htr 
integer, intent(in) :: nelmnt,nSpecs,nnew_elements
integer, intent(in) :: new_elements(nelmnt)
integer, intent(inout) :: nC_basic
integer, intent(inout) :: C_basic(nelmnt,nSpecs)
integer, intent(in) :: is_C_basic_mole
integer, intent(inout) :: C_basic_mole(nSpecs)
integer, intent(in) :: nnew_species,nnew_species_local
integer, intent(in) :: new_species(nSpecs)
integer :: nAtoms_Specs(nSpecs,nelmnt)
integer :: i,j

do j=1,nSpecs
   do i=1,nelmnt
      nAtoms_Specs(j,i)=nAtoms(htr,j,i)
!      print*,j,i,nAtoms_Specs(j,i)
   enddo
enddo

if (nnew_species_local>0) then
   do j=1,nnew_species_local
      do i=1,nnew_elements
        C_basic(new_elements(i),new_species(nnew_species-nnew_species_local+j))&
        =C_basic(new_elements(i),new_species(nnew_species-nnew_species_local+j))&
        +nAtoms_Specs(new_species(nnew_species-nnew_species_local+j),new_elements(i))
      enddo
   enddo
   
   if (is_C_basic_mole==1) then
      do j=1,nnew_species_local
         C_basic_mole(new_species(nnew_species-nnew_species_local+j))=1
      enddo 
      print*,"update_C_basic_mole"
      write(*,'(26(i3,1X))'),C_basic_mole(:)  
   endif
endif

do i=1,nnew_elements
   print*,"update_C_basic"
   write(*,'(26(i3,1X))') C_basic(new_elements(i),:)
enddo

end subroutine update_C_basic
!----------------------------------------------
subroutine creat_C_basic(htr,nelmnt,nSpecs,nnew_elements_local,nnew_elements,&
new_elements,j_new_elements,nnew_species,nnew_species_local,new_species,&
C_basic,nC_basic)
type(hectars), intent(inout) :: htr 
integer, intent(in) :: nelmnt,nSpecs,nnew_elements_local,nnew_elements
integer, intent(in) :: new_elements(nelmnt)
integer, intent(in) :: j_new_elements(nelmnt)
integer, intent(in) :: nnew_species
integer, intent(in) :: nnew_species_local
integer, intent(in) :: new_species(nSpecs) 
integer, intent(inout) :: C_basic(nelmnt,nSpecs)
integer, intent(inout) :: nC_basic
integer :: i,j

do i=nnew_elements-nnew_elements_local+1,nnew_elements
   do j=nnew_species-nnew_species_local+1,nnew_species
      C_basic(new_elements(i),new_species(j))=nAtoms(htr,new_species(j),&
      new_elements(i))
   enddo
   nC_basic = nC_basic + 1
   print*,"create_C_basic"
   write(*,'(26(i3,1X))') C_basic(new_elements(i),:)
enddo

end subroutine creat_C_basic
!----------------------------------------------
subroutine creat_C_basic_mole(nSpecs,is_C_basic_mole,C_basic_mole,&
nnew_species,nnew_species_local,new_species)
integer, intent(in) :: nSpecs
integer, intent(out) :: is_C_basic_mole
integer, intent(out) :: C_basic_mole(nSpecs)
integer, intent(in) :: nnew_species,nnew_species_local
integer, intent(in) :: new_species(nSpecs)
integer :: j

is_C_basic_mole = 1
C_basic_mole = 0   
do j=1,nnew_species
   C_basic_mole(new_species(j))=1
enddo 
print*,"create_C_basic_mole"
write(*,'(26(i3,1X))') C_basic_mole(:)  

end subroutine creat_C_basic_mole
!----------------------------------------------
subroutine create_C_standard(i,nfast_rxns,fast_rxns,nSpecs,nrxn,n_dot,&
nnew_species,nnew_species_local,new_species,nC_standard,C_standard)
integer, intent(in) :: i
integer, intent(in) :: nfast_rxns
integer, intent(in) :: fast_rxns(nfast_rxns)
integer, intent(in) :: nSpecs,nrxn
integer, intent(in) :: n_dot(nSpecs,nrxn)
integer, intent(in) :: nnew_species
integer, intent(in) :: nnew_species_local
integer, intent(in) :: new_species(nSpecs) 
integer, intent(inout) :: nC_standard
integer, intent(inout) :: C_standard(nrxn,nSpecs)
integer :: coefficient1,coefficient2
integer :: j

do j=nnew_species-nnew_species_local+2,nnew_species
   coefficient1=n_dot(new_species(j),fast_rxns(i))
   coefficient2=n_dot(new_species(nnew_species-nnew_species_local+1),fast_rxns(i))   
   print*,"Specs being removed",new_species(nnew_species-nnew_species_local+1)
   print*,"Specs being updated",new_species(j)
   print*,"omega of Specs being removed",coefficient2
   print*,"omega of Specs being updated",coefficient1
   
   coefficient1=-1*coefficient1
   
   nC_standard=nC_standard+1

   C_standard(nC_standard,new_species(nnew_species-nnew_species_local+1))=&
   coefficient1
   C_standard(nC_standard,new_species(j))=coefficient2
   
   print*,"nC_standard",nC_standard
   print*,"create_C_standard"
   write(*,'(26(i3,1X))') C_standard(nC_standard,:)
enddo

end subroutine create_C_standard
!----------------------------------------------
subroutine update_C_standard(i,nfast_rxns,fast_rxns,nSpecs,nrxn,n_dot,&
nnew_species,nnew_species_local,new_species,nC_standard,C_standard,&
ndependent_C_standard,dependent_C_standard)
integer, intent(in) :: i
integer, intent(in) :: nfast_rxns
integer, intent(in) :: fast_rxns(nfast_rxns)
integer, intent(in) :: nSpecs,nrxn
integer, intent(in) :: n_dot(nSpecs,nrxn)
integer, intent(in) :: nnew_species
integer, intent(in) :: nnew_species_local
integer, intent(in) :: new_species(nSpecs) 
integer, intent(in) :: nC_standard
integer, intent(inout) :: C_standard(nrxn,nSpecs)
integer, intent(in) :: ndependent_C_standard
integer, intent(in) :: dependent_C_standard(nC_standard)
integer :: coefficient1,coefficient2
integer :: C_standard_omega(nrxn,nrxn)
integer :: ii,j,k

call C_standard_to_omega(nSpecs,nrxn,nC_standard,n_dot,C_standard,C_standard_omega)

do ii=1,ndependent_C_standard
   print*,"C_standard_omega(dependent_C_standard(ii),fast_rxns(i))",&
   C_standard_omega(dependent_C_standard(ii),fast_rxns(i))
   
   coefficient1=C_standard_omega(dependent_C_standard(ii),fast_rxns(i))
   coefficient2=n_dot(new_species(nnew_species-nnew_species_local+1),fast_rxns(i))   
   print*,"Specs being removed",new_species(nnew_species-nnew_species_local+1)
   print*,"C_standard being updated",dependent_C_standard(ii)
   print*,"omega of Specs being removed",coefficient2
   print*,"omega of C_standard being updated",coefficient1
   
   coefficient1=-1*coefficient1
   
   do j=1,nSpecs
      C_standard(dependent_C_standard(ii),j) = &
      coefficient2*C_standard(dependent_C_standard(ii),j)
   enddo
   C_standard(dependent_C_standard(ii),&
   new_species(nnew_species-nnew_species_local+1))=coefficient1
   print*,"update_C_standard"
   write(*,'(26(i3,1X))') C_standard(dependent_C_standard(ii),:)
   
   call C_standard_to_omega(nSpecs,nrxn,nC_standard,n_dot,C_standard,C_standard_omega)
   print*,"C_standard_omega(dependent_C_standard(ii),fast_rxns(i))",&
   C_standard_omega(dependent_C_standard(ii),fast_rxns(i))
enddo

end subroutine update_C_standard
!----------------------------------------------
subroutine update_C_standard_2(i,nfast_rxns,fast_rxns,nSpecs,nrxn,n_dot,&
nnew_species,nnew_species_local,new_species,nC_standard,C_standard,&
ndependent_C_standard,dependent_C_standard)
integer, intent(in) :: i
integer, intent(in) :: nfast_rxns
integer, intent(in) :: fast_rxns(nfast_rxns)
integer, intent(in) :: nSpecs,nrxn
integer, intent(in) :: n_dot(nSpecs,nrxn)
integer, intent(in) :: nnew_species
integer, intent(in) :: nnew_species_local
integer, intent(in) :: new_species(nSpecs) 
integer, intent(inout) :: nC_standard
integer, intent(inout) :: C_standard(nrxn,nSpecs)
integer, intent(in) :: ndependent_C_standard
integer, intent(in) :: dependent_C_standard(nC_standard)
integer :: coefficient1,coefficient2
integer :: C_standard_omega(nrxn,nrxn)
integer :: ii,j,k

call C_standard_to_omega(nSpecs,nrxn,nC_standard,n_dot,C_standard,C_standard_omega)

if (ndependent_C_standard>1) then
do ii=2,ndependent_C_standard
   print*,"C_standard_omega(dependent_C_standard(ii),fast_rxns(i))",&
   C_standard_omega(dependent_C_standard(ii),fast_rxns(i))
   
   coefficient1=C_standard_omega(dependent_C_standard(ii),fast_rxns(i))
   coefficient2=C_standard_omega(dependent_C_standard(1),fast_rxns(i))   
   print*,"C_standard being removed",dependent_C_standard(1)
   print*,"C_standard being updated",dependent_C_standard(ii)
   print*,"omega of C_standard being removed",coefficient2
   print*,"omega of C_standard being updated",coefficient1
   
  coefficient1=-1*coefficient1
   
   do j=1,nSpecs
      C_standard(dependent_C_standard(ii),j) = &
      coefficient2*C_standard(dependent_C_standard(ii),j)+&
      coefficient1*C_standard(dependent_C_standard(1),j)
   enddo
   
   print*,"update_C_standard_2"
   write(*,'(26(i3,1X))') C_standard(dependent_C_standard(ii),:)
   
   call C_standard_to_omega(nSpecs,nrxn,nC_standard,n_dot,C_standard,C_standard_omega)
   print*,"C_standard_omega(dependent_C_standard(ii),fast_rxns(i))",&
   C_standard_omega(dependent_C_standard(ii),fast_rxns(i))
enddo
endif

do ii=dependent_C_standard(1)+1,nC_standard
      C_standard(ii-1,1:nSpecs) = C_standard(ii,1:nSpecs)
enddo

nC_standard = nC_standard - 1
print*,"nC_standard",nC_standard
end subroutine update_C_standard_2
!----------------------------------------------
subroutine C_standard_to_omega(nSpecs,nrxn,nC_standard,n_dot,C_standard,&
C_standard_omega)
integer, intent(in) :: nSpecs,nrxn,nC_standard
integer, intent(in) :: n_dot(nSpecs,nrxn)
integer, intent(in) :: C_standard(nrxn,nSpecs)
integer, intent(out) :: C_standard_omega(nrxn,nrxn)
integer :: i,j,k

open(55,file='C_standard_omega.txt')

C_standard_omega=0

do i=1,nC_standard
   do k=1,nrxn
      do j=1,nSpecs
         C_standard_omega(i,k)=C_standard_omega(i,k)+C_standard(i,j)*n_dot(j,k)
      enddo
   enddo
   write(55,'(130(i4,1X))') C_standard_omega(i,1:nrxn)   
enddo 

close(55)
end subroutine C_standard_to_omega
!----------------------------------------------
subroutine omega_to_C_standard(nSpecs,nrxn,nC_standard,n_dot,C_standard,C_standard_omega)
integer, intent(in) :: nSpecs,nrxn,nC_standard
integer, intent(in) :: n_dot(nSpecs,nrxn)
integer, intent(out) :: C_standard(nrxn,nSpecs)
integer, intent(in) :: C_standard_omega(nrxn,nrxn)
integer :: n_dot_t(nrxn,nSpecs)
integer :: B(nSpecs,nSpecs)
integer :: B_inv(nSpecs,nSpecs)
integer :: C(nC_standard,nSpecs)
integer :: rowsC,colsC,ErrCode,indx(nSpecs)
integer :: j

call transpose_mtrx(n_dot,nSpecs,nrxn,n_dot_t)
call matrixproduct(n_dot,nSpecs,nrxn,n_dot_t,nrxn,nSpecs,B,rowsC,colsC,ErrCode)
do j=1,nSpecs
   write(*,'(26(i3,1X))'),B(j,:)
enddo
open(55,file='n_dot.txt')
do j=1,nSpecs
   write(55,'(130(i3,1X))'),n_dot(j,1:nrxn)
enddo
close(55)
call migs(B,nSpecs,B_inv,indx)
call matrixproduct(C_standard_omega,nrxn,nrxn,n_dot_t,nrxn,nSpecs,C,rowsC,&
     colsC,ErrCode)
call matrixproduct(C,nC_standard,nSpecs,B_inv,nSpecs,nSpecs,C_standard,rowsC,&
     colsC,ErrCode)
     
end subroutine omega_to_C_standard
!----------------------------------------------
subroutine transpose_mtrx(a,n_r,n_c,b)
integer, intent(in) :: n_r,n_c
integer, intent(in) :: a(n_r,n_c)
integer, intent(out) :: b(n_c,n_r)
integer :: i,j

do i=1,n_r
   do j=1,n_c
      b(j,i)=a(i,j)
   enddo  
enddo 

end subroutine transpose_mtrx       
end module constraints_mod
