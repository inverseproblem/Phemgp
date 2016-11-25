!------------------------------------------------------------------------
!
!    Copyright 2014,2015  Andrea Zunino 
!
!
!    Phemgp is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Phemgp is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Phemgp.  If not, see <http://www.gnu.org/licenses/>.
!
!------------------------------------------------------------------------
 
!#######################################################################
!#######################################################################

module realprec

  integer,parameter  :: pdigits=15  !15 
  integer,parameter  :: rexprange=307  !307 
    
  !! RESULT = SELECTED_REAL_KIND([P, R, RADIX]) 
  integer,parameter  :: dp = selected_real_kind(pdigits, rexprange)

  public :: pdigits,rexprange,dp

end module realprec

!#######################################################################

module  proptab

  use realprec
  implicit  none

  integer, parameter :: maxnumtab = 20
  integer, parameter :: lenstrnameprop = 24
  integer, parameter :: numnameprops = 50

  type tabledims
     character(len=1024) :: inptable(maxnumtab)
     integer :: nxpts,nypts,nzfiles,nprop,totnumphases
     real(dp) :: xmin,ymin,dx,dy
     real(dp), allocatable :: ztabs(:)
     character(len=128) :: title
     character(len=18)   :: xname,yname  
     character(len=lenstrnameprop) :: nameprop(numnameprops)
  end type tabledims
          
  type indexes
     integer :: vol_index
  end type indexes

  
  real(dp), allocatable, target :: propar(:,:,:,:)
  character(12), allocatable :: phasar(:,:,:,:)
  integer  :: proparlinelength 
  integer  :: phasarlinelength 
  integer, parameter :: rowlength = 2048
  character(len=12), allocatable :: lstphas(:)


  public :: tabledims,indexes,maxnumtab,propar,phasar,proparlinelength,phasarlinelength,lstphas,&
       redtab,mixprop,extractminpr,extractsingleminpr
  private :: rowlength,average,bilinear,xy2ij,checknodes,calcz


contains
  
  !!=====================================================================

  subroutine redtab(propar, phasar, tbd,lstphas)   
    
    !    propar: 3d array of real(dp) containing the bulk prop. and
    !             the properties of each phase. the first 3 indices contain
    !             p, t, followed by the no. of the phases present at 
    !             given conditions. the remaining elements are the bulks 
    !             properties followed by the prop. for each phase.
    !    phasar: 3d array of character containing  names of the phases 
    !             present at that conditions
    !    title: title of the table as appear inside the file
    !    nxpts, nypts: no. of points along i [x] and j [y] axes
    !    xmin, ymin: starting point of x and y
    !    dx, dy: increments in x and y direction
    !    xname, yname: label of x and y axis 
    !    nprop: no. of properties 
    !
    !    ------------------------------------------------------
    !     xpts = xpoint, ypts = ypoint
    !

    implicit none 

    character(len=12), allocatable :: lstphas(:)
    character(len=100) :: tfile,filenumstr
    character(1)  :: row(rowlength),substring(lenstrnameprop)
    character(24)  :: formatversion
    character(12)  :: pres,temp,listofphases(rowlength)
    real(dp), allocatable :: propar(:,:,:,:)
    character(12), allocatable :: phasar(:,:,:,:)
    character(12) :: curphase,emptystring,strfmt
    integer :: nprop,nphas,iph,kount,p, &
         nelem,ibulk,i,j,k,fle,maxlinelength,maxnphas, &
         noindipvar,istrs,istre,startstr(rowlength),endstr(rowlength)&
         ,mah,k2
    type(tabledims) :: tbd


    maxlinelength=0
    maxnphas=0

    if (allocated(tbd%ztabs) ) deallocate(tbd%ztabs)
    allocate(tbd%ztabs(tbd%nzfiles))

    !! starts reading tables
    !! ###############################################################
    !! ##### Tables MUST have the same number of properties !!!! #####
    !! ###############################################################
          
    !!------------ scan for phases present in the table ----------
    write(*,*)
    write(*,*) "Scanning for phase names..."
    write(*,*)

    !! remember, in redtab:  phasar=" "
    emptystring=" "  !!!phasar(1,1,1,phasarlinelength)  
    listofphases(1) = " "
    p=1
    maxnphas=0
    row(:)=' '
    substring(:)=' '

    do fle=1,tbd%nzfiles
       
       tfile = trim(adjustl(tbd%inptable(fle)))
       open(unit=10, file=tfile, status='old')
       
       read(10,*) formatversion
       read(10,*) tbd%title
       read(10,*) tbd%ztabs(fle)  !!! <<< NEW
       read(10,*) noindipvar
       read(10,*) tbd%xname
       read(10,*) tbd%xmin
       read(10,*) tbd%dx
       read(10,*) tbd%nxpts
       read(10,*) tbd%yname 
       read(10,*) tbd%ymin
       read(10,*) tbd%dy
       read(10,*) tbd%nypts
       read(10,*) tbd%nprop
       tbd%nprop=tbd%nprop-4 ! skipping first 4 columns |this header == all fields (# columns)
       read(10,*)   !! skipping properties etc strings

       ! loop on y nodes
       do j=1,tbd%nypts
          ! loop on x nodes
          do i=1,tbd%nxpts

             read(10,*) curphase, nphas
             !print*,curphase,"nphas", nphas
             if ( nphas > maxnphas ) then
                maxnphas=nphas
             endif

             LOOP_IDX4: do iph=1,nphas
                
                read(10,fmt='(a12)') curphase
                !print*,'curphase: ',curphase
                curphase=trim(adjustl(curphase))
                if (curphase==emptystring) then 
                   !! if curphase is an empty string jump to next element
                   cycle LOOP_IDX4 
                else if ( all(listofphases(1:p) /= curphase) )  then  
                   !!! ( all(listofphases(1:(p-1)) /= curphase) ) 
                   print*,'curphase: ',curphase
                   listofphases(p)=curphase 
                   p=p+1
                endif

             enddo LOOP_IDX4
             ! enddo LOOP_IDX4
          enddo
       enddo 
       
       close(10)
    enddo  !! ends z loop
     
    if (allocated(lstphas) ) deallocate(lstphas)
    allocate( lstphas(p-1) )  !! p-1 = total number of phases
    lstphas = listofphases(1:p-1)
    tbd%totnumphases = p-1
    print*,'total number of phases', p-1
    !print*,'max num. phases per one line',maxnphas 
    !!------- end scan -----------

    phasarlinelength= maxnphas
    proparlinelength= (phasarlinelength+1)*(tbd%nprop+1)-1
    
    print*,"############################################################"
    print*, "proparlinelength:", proparlinelength
    print*, "phasarlinelength:", phasarlinelength
    print*,"############################################################"

  
    ! allocate variable contained in module sizetab
    if (allocated(propar) ) deallocate(propar) !! in case one asks to reload tables from main menu
    if (allocated(phasar) ) deallocate(phasar) !! in case one asks to reload tables from main menu 
    allocate( propar(tbd%nxpts,tbd%nypts,tbd%nzfiles,proparlinelength) )
    allocate( phasar(tbd%nxpts,tbd%nypts,tbd%nzfiles,phasarlinelength) )
    propar=0.0  !! initialize arrays
    phasar=" "

    do fle=1,tbd%nzfiles
  
       print*, 'tbd%nzfiles',tbd%nzfiles
       tfile = trim(adjustl(tbd%inptable(fle)))
       print*,"--------------------------------------------"
       write (*,*) 'Reading file: ',tfile
       open(unit=10, file=tfile, status='old')
       read(10,*) formatversion
       read(10,*) tbd%title
       read(10,*) tbd%ztabs(fle)  !!! <<< NEW
       read(10,*) noindipvar
       read(10,*) tbd%xname
       read(10,*) tbd%xmin
       read(10,*) tbd%dx
       read(10,*) tbd%nxpts
       read(10,*) tbd%yname 
       read(10,*) tbd%ymin
       read(10,*) tbd%dy
       read(10,*) tbd%nypts
       read(10,*) tbd%nprop 
       tbd%nprop=tbd%nprop-4 ! this header == all fields (# columns): 
       !+ Name        Counter T(K)           P(bar)         phase,vo%  ......
       write(strfmt,'(i12)') rowlength
       strfmt="(a"//trim(adjustl(strfmt))//")"
       !print*,'^^^^ strfmt:',strfmt
       write(strfmt,'(i12)') rowlength
       strfmt="("//trim(adjustl(strfmt))//"(a1))"
       read(10,strfmt) (row(k),k=1,rowlength)

       !print *, "row:",(row(k),k=1,rowlength),"|END"
       istrs=1
       istre=1
       do k=1,(rowlength-1)
          if ( (row(k)==" ") .and. (row(k+1)/=" ") ) then
             startstr(istrs)=k+1
             !print *,istrs,startstr(istrs),row(k+1)
             istrs=istrs+1
          else if ( (row(k)/=" ") .and. (row(k+1)==" ") ) then
             endstr(istre)=k
             !print *,istre,endstr(istre),row(k)
             istre=istre+1
          end if
       end do

       k2=1
       ! skip first 2+noindipvar substrings
       do k=(2+noindipvar+1),istrs-1
          !print *,'k:',k,startstr(k),endstr(k)
          substring(1:(-startstr(k)+endstr(k)+1)) = &
               row(startstr(k):endstr(k))
          !print *, "substring:",substring," | ",row(startstr(k):endstr(k))
          do mah=1,(-startstr(k)+endstr(k)+1)
             tbd%nameprop(k2)(mah:mah)=substring(mah)
          enddo
          do mah=(-startstr(k)+endstr(k)+1)+1,lenstrnameprop !! len=24  !12
             tbd%nameprop(k2)(mah:mah)=" "
          enddo         
          ! print *, "substring:",substring
          ! print *,"tbd%nameprop:",tbd%nameprop(k2) 
          k2 = k2 + 1
       enddo

       write(*,*) "Format version: ",formatversion
       write(*,*) "Number of independent variables: ",noindipvar
       write(*,*) 'Title inside file: ',trim(tbd%title)
       write(*,*)  'z value (could be e.g. basalt fraction or water content): ',tbd%ztabs(fle)
       write(*,*) "x: (npts,xmin,x_incr) ",tbd%nxpts,tbd%xmin,tbd%dx
       write(*,*) "y: (npts,ymin,y_incr) ",tbd%nypts,tbd%ymin,tbd%dy
       write(*,*) "num. of tables:",tbd%nzfiles
       write(*,*) 'x name, y name: ',tbd%xname,tbd%yname
       write(*,*) "Number of properties: ",tbd%nprop
       write(*,*) 'Properties: ',(tbd%nameprop(k)//" ",k=1,tbd%nprop)
       write(*,*) "---------------------------------------------"
       write(*,*)
       
       ! tbd%nxpts*tbd%nypts lines correspond to the no. of points in x and y 
       ! loop on y nodes
       do j=1,tbd%nypts
          ! loop on x nodes
          do i=1,tbd%nxpts

             ! propar(i,j,fle,*ibulk+1*): ibulk+1 because first element (bulk) 
             ! contains the no. of phases
             read(10,fmt='(a12)',ADVANCE="no") curphase
             read(10,*) nphas, temp, pres, (propar(i,j,fle,ibulk+1),ibulk=1,tbd%nprop)

             propar(i,j,fle,1)= real(nphas)
          
             do iph=1,nphas
                kount=iph*tbd%nprop
                ! propar(i,j,fle,*nelem+2*): *nelem+2* because starts at tbd%nprop+2 
                read(10,fmt='(a12)',ADVANCE="no") phasar(i,j,fle,iph) !! to read '/' and other keywords in a string
                read(10,*) nphas, temp, pres, &
                    (propar(i,j,fle,nelem+2),nelem=kount,kount+tbd%nprop-1 )
             enddo

          end do !  end do for i      
       end do !  end do for j

       close(10)
       
     enddo !!### enddo on tbd%nzfiles ################

    return
  end subroutine redtab
     
  !!=====================================================================

  subroutine mixprop(propar,x,y,z,tbd,inp_req,indx,bulk_prop,phas_prop)

    implicit none
      
    integer :: i,j,k,num_preq_bulk,num_preq_phas
    real(dp),intent(in) :: x,y,z
    real(dp) ::zdec
    logical :: left,down
    logical :: onanode
    integer :: i_int,j_int,z1,z2,ii,jj
    real(dp) :: line_a(proparlinelength),line_b(proparlinelength) &
         ,line_c(proparlinelength),line_d(proparlinelength)
    real(dp),allocatable :: bulk_prop_int1(:),bulk_prop_int2(:)&
         ,phas_prop_int1(:),phas_prop_int2(:)
    real(dp),allocatable,intent(out) :: bulk_prop(:), phas_prop(:)
    real(dp), target,intent(in) :: propar(:,:,:,:)
    real(dp),allocatable :: av_phas_prop_a(:),av_phas_prop_b(:),&
         av_phas_prop_c(:), av_phas_prop_d(:),bulk_prop_a(:) &
         ,bulk_prop_b(:),bulk_prop_c(:),bulk_prop_d(:)
    real(dp), pointer :: propar1(:,:,:),propar2(:,:,:) 
    integer, allocatable :: preq_bulk(:),preq_phas(:,:)
    integer,intent(in) :: inp_req(:,:)
    type(tabledims),intent(in) :: tbd
    type(indexes),intent(in) :: indx ! which property to average on!

        
    num_preq_bulk = count(inp_req(1,:)>0)
    num_preq_phas = count(inp_req(2,:)>0) 
    allocate( preq_bulk(num_preq_bulk),preq_phas(2,num_preq_phas) )
    preq_bulk=0 ! initializing to zero all elements
    preq_phas=0
    
    ! requested bulk properties
    do j=1,num_preq_bulk
       preq_bulk(j) =inp_req(1,j) ! bulk props from perplex, not used
    enddo
    
    ! requested phase properties
    do j=1,num_preq_phas
       preq_phas(1,j) = inp_req(2,j)  !! which properties
       preq_phas(2,j) = inp_req(3,j)  !! average scheme
    enddo
       
    if (.not. allocated(bulk_prop) ) then
       allocate(bulk_prop(num_preq_bulk))
    endif
    allocate(bulk_prop_a(num_preq_bulk), &
         bulk_prop_b(num_preq_bulk),bulk_prop_c(num_preq_bulk),&
         bulk_prop_d(num_preq_bulk), bulk_prop_int1(num_preq_bulk))

    if (.not. allocated(phas_prop) ) then
       allocate(phas_prop(num_preq_phas))
    endif
    allocate(av_phas_prop_a(num_preq_phas),&
         av_phas_prop_b(num_preq_phas),av_phas_prop_c(num_preq_phas),&
         av_phas_prop_d(num_preq_phas),phas_prop_int1(num_preq_phas))


    if (tbd%nzfiles > 1) then
       !! compute the z fraction
       call calcz(z,tbd,z1,z2,zdec)
    else
       z1=1
    end if
   
    
    !----------------------------------------------------------------------
    !######### TABLE 1 ###################################
    
    propar1 => propar(:,:,z1,:) !! pointer to table 1

    call xy2ij(x,y,tbd,ii,jj,left,down)
    call checknodes(x,y,ii,jj,tbd,left,down,onanode,i_int,j_int)

    do k=1,proparlinelength
       line_a(k) = propar1(i_int,  j_int,  k)
       line_b(k) = propar1(i_int+1,j_int,  k)
       line_c(k) = propar1(i_int+1,j_int+1,k)
       line_d(k) = propar1(i_int,  j_int+1,k)
    enddo

    call average(line_a,tbd%nprop,preq_bulk,preq_phas,bulk_prop_a,av_phas_prop_a,indx)
    call average(line_b,tbd%nprop,preq_bulk,preq_phas,bulk_prop_b,av_phas_prop_b,indx)
    call average(line_c,tbd%nprop,preq_bulk,preq_phas,bulk_prop_c,av_phas_prop_c,indx)
    call average(line_d,tbd%nprop,preq_bulk,preq_phas,bulk_prop_d,av_phas_prop_d,indx)
    call bilinear(bulk_prop_a,bulk_prop_b,bulk_prop_c,bulk_prop_d, &
         tbd,x,y,i_int,j_int,bulk_prop_int1)
    call bilinear(av_phas_prop_a,av_phas_prop_b,av_phas_prop_c,av_phas_prop_d, &
         tbd,x,y,i_int,j_int,phas_prop_int1)
    
    !----------------------------------------------------------------------
    !######### TABLE 2 ###################################
    if (tbd%nzfiles > 1) then 

       allocate(bulk_prop_int2(num_preq_bulk))
       allocate(phas_prop_int2(num_preq_phas))

       propar2 => propar(:,:,z2,:) !! pointer to table 2
 
       call xy2ij(x,y,tbd,ii,jj,left,down)
       call checknodes(x,y,ii,jj,tbd,left,down,onanode,i_int,j_int)

       do k=1,proparlinelength
          line_a(k) = propar2(i_int,  j_int,  k)
          line_b(k) = propar2(i_int+1,j_int,  k)
          line_c(k) = propar2(i_int+1,j_int+1,k)
          line_d(k) = propar2(i_int,  j_int+1,k)
       enddo
       call average(line_a,tbd%nprop,preq_bulk,preq_phas,bulk_prop_a,av_phas_prop_a,indx)
       call average(line_b,tbd%nprop,preq_bulk,preq_phas,bulk_prop_b,av_phas_prop_b,indx)
       call average(line_c,tbd%nprop,preq_bulk,preq_phas,bulk_prop_c,av_phas_prop_c,indx)
       call average(line_d,tbd%nprop,preq_bulk,preq_phas,bulk_prop_d,av_phas_prop_d,indx) 
       call bilinear(bulk_prop_a,bulk_prop_b,bulk_prop_c,bulk_prop_d, &
            tbd,x,y,i_int,j_int,bulk_prop_int2)
       call bilinear(av_phas_prop_a,av_phas_prop_b,av_phas_prop_c,av_phas_prop_d, &
            tbd,x,y,i_int,j_int,phas_prop_int2)
 
       !----------------------------------------------------------------------
       !#####  arithmetic mean between tables 1 & 2 #####
       
       do i=1,num_preq_bulk
          !! arithmetic mean in z direction  
          bulk_prop(i) = (1.-zdec)*bulk_prop_int1(i)+zdec*bulk_prop_int2(i)
       enddo
       do i=1,num_preq_phas
          !! arithmetic mean in z direction 
          phas_prop(i) = (1.-zdec)*phas_prop_int1(i)+zdec*phas_prop_int2(i)
       enddo

       deallocate(bulk_prop_int2,phas_prop_int2 )    

    else

       !!!### only 1 table loaded
       do i=1,num_preq_bulk
          bulk_prop(i) = bulk_prop_int1(i)
       enddo
       do i=1,num_preq_phas
          phas_prop(i) = phas_prop_int1(i)
       enddo
    
    endif

    deallocate(preq_bulk,preq_phas)
    deallocate(bulk_prop_a, &
         bulk_prop_b,bulk_prop_c,&
         bulk_prop_d, bulk_prop_int1)
 
         
    deallocate(av_phas_prop_a,&
         av_phas_prop_b,av_phas_prop_c,&
         av_phas_prop_d,phas_prop_int1)
   
    return
  end subroutine mixprop
  
  !______________________________________________________________
  !______________________________________________________________
  
  subroutine average(line,nprop,preq_bulk,preq_phas,bulk_prop,av_phas_prop,indx)
    
    implicit none
    
    real(dp) :: line(proparlinelength) ,cur_prop_av,cur_prop_reuss,cur_prop_voigt,sumvol 
    integer :: nphas,nprop,l,k,i,j
    real(dp) :: bulk_prop(:),av_phas_prop(:)
    real(dp), allocatable :: cur_prop(:),vol(:)
    integer :: preq_bulk(:),preq_phas(:,:)
    type(indexes) :: indx

    nphas=int(line(1))
    allocate(cur_prop(nphas),vol(nphas))

    !! #########  bulk properties  ############
    bulk_prop=0.0
    if ( (maxval(preq_bulk)) > 0 ) then       
       do i=1,size(preq_bulk)
          !! i+1 because first is the no. of phases
          bulk_prop(i) = line(preq_bulk(i)+1)  !!CHECK CHECK CHECK
       enddo
    endif
    
    !! #########  phase properties  ############
    if ( (maxval(preq_phas(1,:))) > 0 ) then 
          
       l=1
       sumvol=0.0
       do k = (nprop+1+indx%vol_index),(nprop*nphas+1+nprop),nprop
          vol(l)=line(k) !!!!! <<<<<<<<<<<<<<<<<<<<<<<<<-----------<<<< 
          sumvol = sumvol + vol(l)
          l = l+1                 
       enddo
       !!!!!######### normalization ########
       vol=vol/sumvol  !! for all elements

       ! loop only on requested properties
       do j=1,size(preq_phas,2)
          i=preq_phas(1,j)  !! <<<<<<---<<<<<<

          ! loop on phases for a single property
          l=1
          do k = (nprop+1+i),(nprop*(nphas+1)+1),nprop
             cur_prop(l)=line(k)
             l = l+1
          enddo
 
          if (preq_phas(2,j)==1) then  
             
             !! ####  arithmetic mean #####
             cur_prop_av=0.0
             do k = 1,nphas 
                cur_prop_av=cur_prop_av+vol(k)*cur_prop(k) 
             enddo
             av_phas_prop(j)=cur_prop_av
             
          else if (preq_phas(2,j)==2) then
             
             !!####  VRH averaging  ####
             cur_prop_voigt=0.0d0
             cur_prop_reuss=0.0d0
             do k = 1,nphas
                cur_prop_voigt = cur_prop_voigt+vol(k)*cur_prop(k)
                cur_prop_reuss = cur_prop_reuss+vol(k)/cur_prop(k)
             enddo
             cur_prop_reuss = 1. / cur_prop_reuss
             cur_prop_av = 0.5d0*(cur_prop_reuss + cur_prop_voigt)
             av_phas_prop(j)=cur_prop_av

          else if (preq_phas(2,j)==3) then
             
             !!####  Voigt averaging  ####
             cur_prop_voigt=0.0d0
             do k = 1,nphas
                cur_prop_voigt = cur_prop_voigt+vol(k)*cur_prop(k)
             enddo
             av_phas_prop(j)=cur_prop_voigt

          else if (preq_phas(2,j)==4) then
             
             !!####  Reuss averaging  ####
             cur_prop_reuss=0.0d0
             do k = 1,nphas
                cur_prop_reuss = cur_prop_reuss+vol(k)/cur_prop(k)
             enddo
             cur_prop_reuss = 1. / cur_prop_reuss
             av_phas_prop(j)= cur_prop_reuss
             
          endif
       enddo  ! for nprop...
    endif !! ends: if ( (maxval(preq_phas(1,:))) > 0 ) then 

    deallocate(cur_prop,vol)

    return
  end subroutine average

  !!=====================================================================
  
  subroutine bilinear(prop_node_a,prop_node_b,prop_node_c,prop_node_d, &
       tbd,x,y,i_int,j_int,prop_int)
    !----------------------------------------------------------------------
    !     cfr. numerical recipes in fortran 77 vol.1 pag 116
    !
    !----------------------------------------------------------------------
 
    implicit none
    integer :: k,i_int,j_int
    real(dp) :: xlow,xhigh,ylow,yhigh,t,u,no1,no2,no3,no4,x,y
    real(dp) :: prop_int(:),prop_node_a(:),&
         prop_node_b(:),prop_node_c(:),prop_node_d(:)
    type(tabledims) :: tbd
    !
    ! the nodes are anticlockwise numbered
    !
    ! no4    no3
    !
    ! no1    no2
    !
    ! no1 : the reference point; the origin for u and t.
    ! t and u go from 0.0 to 1.0
    !
    
    ! bilinear interpolation for propar
    !     i + nx*(j-1)
    xlow  = tbd%xmin + dble(i_int-1)  * tbd%dx
    xhigh = tbd%xmin + dble(i_int+1-1)* tbd%dx
    ylow  = tbd%ymin + dble(j_int-1)  * tbd%dy 
    yhigh = tbd%ymin + dble(j_int+1-1)* tbd%dy
    t = (x-xlow) / (xhigh-xlow)
    u = (y-ylow) / (yhigh-ylow)
    do k=1,size(prop_node_a)
       no1=prop_node_a(k)
       no2=prop_node_b(k)
       no3=prop_node_c(k)
       no4=prop_node_d(k)
       prop_int(k) = no1*(1.-t)*(1.-u)+t*(1.-u)*no2+t*u*no3+(1.-t)*u*no4
    enddo  
  
    return
  end subroutine bilinear
  
  !!!!##############################################################

  subroutine extractminpr(xvec,yvec,zvec,propar,phasar,nprop,inp_minpr,tbd,outfilemin1,outfilemin2) 
    
    implicit none
    
    integer :: i,j,z,itab,jtab,nprop,z1,z2,iph,&
         ipr,sngpr,nphases,l,num_minpr_req,pp
    real(dp) :: xvec,yvec,zdec,zvec
    integer, allocatable :: minpr_req(:)
    character(len=12), allocatable :: phasnamelst(:)
    character(len=12),allocatable :: phasar(:,:,:,:)
    real(dp), allocatable :: proplst(:)
    real(dp), allocatable :: propar(:,:,:,:)
    type(tabledims) :: tbd
    character(100) :: outfilemin1,outfilemin2
    logical :: left,down,onanode
    integer :: inp_minpr(:)
    
    call xy2ij(xvec,yvec,tbd,i,j,left,down)
    call checknodes(xvec,yvec,i,j,tbd,left,down,onanode,itab,jtab)
    call calcz(zvec,tbd,z1,z2,zdec)


    open(unit=31,file=outfilemin1,status='unknown',POSITION='APPEND')    
    open(unit=32,file=outfilemin2,status='unknown',POSITION='APPEND')    
    
    num_minpr_req = count(inp_minpr(:)>0)
    allocate(minpr_req(num_minpr_req))  
    

    ! requested bulk properties
    do j=1,num_minpr_req
       minpr_req(j) = inp_minpr(j) 
    enddo
    
    !  if (onanode.eqv..true.) then
    !! ###### IF ON A NODE ???? #####
    if (tbd%nzfiles == 1 ) z2=z1 !! no loop on z

    do z=z1,z2
       !print*, "z1,z2",z1,z2
       do j=jtab,jtab+1  !! #### how to keep "checknodes" routine order ??
          do i=itab,itab+1 
             
             nphases=int(propar(i,j,z,1))
             allocate(proplst(num_minpr_req*nphases))
             allocate(phasnamelst(nphases))
        
             pp=1 !! array index starts from 1!!              
             !! loop on phases
             do iph=1,nphases
                
                phasnamelst(iph) = phasar(i,j,z,iph)
                
                ! loop only on requested properties
                do l=1,num_minpr_req
                   
                   ipr=minpr_req(l)
                   sngpr = nprop*iph+1+ipr
                   proplst(pp) = propar(i,j,z,sngpr)
                   pp=pp+1
              
                enddo
             enddo
             
             !! ---------------------------
             if ( (i==itab).and.(j==jtab)) then 

                if (z==z1) then !! TABLE 1
                   write(31,*) "--- node A ---",xvec,yvec
                   write(31,*) phasnamelst
                   write(31,*) proplst
                else if((z==z2).and.(tbd%nzfiles > 1 ) ) then !! TABLE 2
                   write(32,*) "--- node A ---",xvec,yvec
                   write(32,*) phasnamelst
                   write(32,*) proplst
                endif

             else if ((i==itab+1).and.(j==jtab) ) then

                if (z==z1) then !! TABLE 1
                   write(31,*) "--- node B ---",xvec,yvec
                   write(31,*) phasnamelst
                   write(31,*) proplst
                else if((z==z2).and.(tbd%nzfiles > 1)) then !! TABLE 2
                   write(32,*) "--- node B ---",xvec,yvec
                   write(32,*) phasnamelst
                   write(32,*) proplst
                endif            

             else if ((i==itab+1).and.(j==jtab+1) ) then

                if (z==z1) then !! TABLE 1
                   write(31,*) "--- node C ---",xvec,yvec
                   write(31,*) phasnamelst
                   write(31,*) proplst
                else if((z==z2).and.(tbd%nzfiles > 1 )) then !! TABLE 2
                   write(32,*) "--- node C ---",xvec,yvec
                   write(32,*) phasnamelst
                   write(32,*) proplst
                endif            

             else if ((i==itab).and.(j==jtab+1) ) then

                if (z==z1) then !! TABLE 1
                   write(31,*) "--- node D ---",xvec,yvec
                   write(31,*) phasnamelst
                   write(31,*) proplst
                else if((z==z2).and.(tbd%nzfiles > 1 )) then !! TABLE 2
                   write(32,*) "--- node D ---",xvec,yvec
                   write(32,*) phasnamelst
                   write(32,*) proplst
                endif            

             endif

           deallocate(proplst,phasnamelst)
           
          enddo
       enddo
    enddo
    
    deallocate(minpr_req)

    close(31)
    close(32)

    return
  end subroutine extractminpr

  !!----------------------------------------------------------------------
  !!=====================================================================

  subroutine extractsingleminpr(xvec,yvec,zvec,propar,phasar,nprop,inp_singleminpr,&
       req_singlemin,tbd,propssinglemin,phaseispresent) 

    implicit none
    
    integer :: i,j,z,itab,jtab,nprop,z1,z2,iph,&
         ipr,sngpr,nphases,l,num_minpr_req,pp,idxreqphase
    real(dp),intent(in) :: xvec,yvec,zvec
    integer, allocatable :: minpr_req(:)
    character(len=12),intent(in) :: phasar(:,:,:,:)
    character(len=12) :: req_singlemin
    real(dp),intent(out) :: propssinglemin(:)
    real(dp), allocatable :: tab_1_nd_a(:),&
         tab_2_nd_a(:),tab_1_nd_b(:),tab_2_nd_b(:),&
         tab_1_nd_c(:),tab_2_nd_c(:),tab_1_nd_d(:),tab_2_nd_d(:),&
         prop_int_tab_1(:),prop_int_tab_2(:)
    real(dp),intent(in) :: propar(:,:,:,:)
    real(dp) :: zdec
    type(tabledims) :: tbd

    logical :: left,down,onanode
    integer :: inp_singleminpr(:),nphasepresent
    logical,intent(out) :: phaseispresent

    
    call xy2ij(xvec,yvec,tbd,i,j,left,down)
    call checknodes(xvec,yvec,i,j,tbd,left,down,onanode,itab,jtab)

    if (tbd%nzfiles == 1 ) then
       z1=1
       z2=z1 !! no loop on z
    else
       call calcz(zvec,tbd,z1,z2,zdec)
    end if
    
    num_minpr_req = count(inp_singleminpr(:)>0)
    allocate(minpr_req(num_minpr_req))  
    !!allocate(propssinglemin(num_minpr_req))
    allocate( tab_1_nd_a(num_minpr_req),tab_2_nd_a(num_minpr_req),&
         tab_1_nd_b(num_minpr_req),tab_2_nd_b(num_minpr_req),&
         tab_1_nd_c(num_minpr_req),tab_2_nd_c(num_minpr_req),&
         tab_1_nd_d(num_minpr_req),tab_2_nd_d(num_minpr_req),&
         prop_int_tab_1(num_minpr_req),prop_int_tab_2(num_minpr_req) )
    
        
    ! requested properties
    do j=1,num_minpr_req
       minpr_req(j) = inp_singleminpr(j) 
    enddo
    
    nphasepresent=0 !! initializing variable
    phaseispresent=.false.
       
 
    do z=z1,z2
       do j=jtab,jtab+1  !! #### how to keep "checknodes" routine order ??
          do i=itab,itab+1 
             
             nphases=int(propar(i,j,z,1))
             idxreqphase=0 !! means phase is not present

             pp=1 !! array index starts from 1!!              
             !! loop on phases
             do iph=1,nphases
                !print*,i,j,z,trim(phasar(i,j,z,iph)),trim(req_singlemin)
                if (trim(phasar(i,j,z,iph))==(trim(req_singlemin))) then
                   idxreqphase=iph
                   nphasepresent=nphasepresent+1 
                endif
             enddo
             
             if (idxreqphase==0) then
                !! set to zero the weight of the node if phase is not present
                propssinglemin=0 
             else 
                ! loop only on requested properties
                do l=1,num_minpr_req
                   ipr=minpr_req(l)
                   sngpr = nprop * idxreqphase +1+ipr
                   propssinglemin(pp) = propar(i,j,z,sngpr)
                   pp=pp+1
                enddo
             endif
             
             !! -----------------------------------
             if ( (i==itab).and.(j==jtab)) then 
                if (z==z1) then !! TABLE 1
                   tab_1_nd_a = propssinglemin
                else if((z==z2).and.(tbd%nzfiles > 1 )) then !! TABLE 2
                   tab_2_nd_a = propssinglemin
                endif
             else if ((i==itab+1).and.(j==jtab) ) then
                if (z==z1) then !! TABLE 1
                   tab_1_nd_b = propssinglemin
                else if((z==z2).and.(tbd%nzfiles > 1 )) then !! TABLE 2
                   tab_2_nd_b = propssinglemin
                endif            
             else if ((i==itab+1).and.(j==jtab+1) ) then
                if (z==z1) then !! TABLE 1
                   tab_1_nd_c = propssinglemin
                else if((z==z2).and.(tbd%nzfiles > 1 )) then !! TABLE 2
                   tab_2_nd_c = propssinglemin
                endif            
             else if ((i==itab).and.(j==jtab+1) ) then
                if (z==z1) then !! TABLE 1
                   tab_1_nd_d = propssinglemin
                else if((z==z2).and.(tbd%nzfiles > 1 )) then !! TABLE 2
                   tab_2_nd_d = propssinglemin
                endif            
             endif

          enddo
       enddo
    enddo

    !print*, 'nphasepresent',nphasepresent    

    if (  ((tbd%nzfiles==1).and.(nphasepresent==4)) .or. &
        ((tbd%nzfiles>1).and.(nphasepresent==8)) ) then
       
       call bilinear(tab_1_nd_a,tab_1_nd_b,tab_1_nd_c,tab_1_nd_d, &
            tbd,xvec,yvec,itab,jtab,prop_int_tab_1)

       if (tbd%nzfiles > 1 ) then
          call bilinear(tab_2_nd_a,tab_2_nd_b,tab_2_nd_c,tab_2_nd_d, &
               tbd,xvec,yvec,itab,jtab,prop_int_tab_2)
          do i=1,num_minpr_req
             !! arithmetic mean in z direction  
             propssinglemin(i) = (1.-zdec)*prop_int_tab_1(i)+zdec*prop_int_tab_2(i)
          enddo
       else
          
          propssinglemin = prop_int_tab_1
  
       endif

       phaseispresent = .true.
       
    else
       
       phaseispresent = .false.
       !print*,"### Phase not present, skipping x,y,z =",xvec,yvec,zvec
   
    endif
    
    
    deallocate(minpr_req)
    deallocate( tab_1_nd_a,tab_2_nd_a,&
         tab_1_nd_b,tab_2_nd_b,&
         tab_1_nd_c,tab_2_nd_c,&
         tab_1_nd_d,tab_2_nd_d,&
         prop_int_tab_1,prop_int_tab_2 )

    return
  end subroutine extractsingleminpr

  !!===========================================================
  
  subroutine xy2ij (x,y,tbd,i,j,left,down)
    !----------------------------------------------------------------------
    ! xy2ij - identifies the grid point associated with coordinate x-y 
    !      and the sector of the x-y coordinate relative to the nodal point        
    !----------------------------------------------------------------------

    implicit none
    integer :: i, j
    logical :: left, down
    real(dp) :: res, x, y
    type(tabledims) :: tbd

    res = (x-tbd%xmin)/tbd%dx + 1d0
    i = int(res)
    if( i.lt.1 ) then
       i=1
    elseif (i.ge.tbd%nxpts) then
       i=tbd%nxpts
    endif
    if (res-dble(i).gt.0.5d0) then
       i = i + 1
       left = .true.
    else 
       left = .false.
    end if

    res = (y-tbd%ymin)/tbd%dy + 1d0
    j = int(res) 
    if( j.lt.1 ) then
       j=1
    elseif (j.ge.tbd%nypts) then
       j=tbd%nypts
    endif
    if (res-dble(j).gt.0.5d0) then
       j = j + 1
       down = .true.
    else 
       down = .false.
    end if

    return
  end subroutine xy2ij

  !#######################################################################

  subroutine checknodes(x,y,i,j,tbd,left,down,onanode,i_int,j_int)
    
    implicit none
    logical ::  onanode,left,down
    real(dp) :: res_x,res_y,x,y
    integer :: i_int,j_int,i,j
    type(tabledims) :: tbd

    ! checking whether we're on a node, if so then you can just use 
    ! properties as kept in propar.
    res_x = (x-tbd%xmin)/tbd%dx
    res_y = (y-tbd%ymin)/tbd%dy 
    if (res_x.le.1d-6.and.res_y.le.1d-6) then
       onanode=.true.
       ! i_int=i !! for future modifications...
       ! j_int=j
       ! return
    else
       onanode=.false.
    endif

    !! determining the 4 nodes
    ! the origin i_int,j_int determines the lower left corner

    ! i axis
    if (left .eqv. .true.) then
       i_int = i-1
    else
       i_int = i
    endif

    ! check bounds
    if (i_int .lt. 1 ) then
       i_int = 1
    else if (i_int .ge. tbd%nxpts ) then
       i_int = tbd%nxpts-1
    endif

    ! j axis
    if (down .eqv. .true.) then
       j_int = j-1
    else
       j_int = j
    endif

    ! check bounds
    if (j_int .le. 1 ) then
       j_int = 1
    else if (j_int .ge. tbd%nypts ) then
       j_int = tbd%nypts-1
    endif

    return
  end subroutine checknodes

  !#######################################################################

  subroutine calcz(z,tbd,z1,z2,zdec)
    
    implicit none
    type(tabledims) :: tbd
    real(dp) :: z,zdec,ztab1,ztab2,zmin,zmax
    integer :: z1,z2
    integer :: loc1(1),loc2(1)

    zmin = minval(tbd%ztabs)
    zmax = maxval(tbd%ztabs)

    !!### Mapping of z into indexes ###
    if ((z < zmin).or.(z>zmax) ) then
       print*,"### Error z<zmin or z>zmax ###"
       print*,'z:',z,zmin,zmax
       stop
    endif

    if (tbd%nzfiles==1) then
       z1=1
       z2=z1
       zdec=0.0d0
    else
       loc1=minloc(abs(tbd%ztabs - z))
       loc2=minloc(abs(tbd%ztabs - z), MASK = (tbd%ztabs /= tbd%ztabs(loc1(1)) ) )

       !! 1 MUST BE < 2
       ztab1=min( tbd%ztabs(loc1(1)),tbd%ztabs(loc2(1)) ) 
       ztab2=max( tbd%ztabs(loc1(1)),tbd%ztabs(loc2(1)) ) 
       z1=min(loc1(1),loc2(1))
       z2=max(loc1(1),loc2(1))

       zdec=(z-ztab1)/(ztab2-ztab1)

    endif

    ! print*,'z:',z
    ! print*,'zdec:',zdec
    ! print*,'z1,z2',z1,z2

    return
  end subroutine calcz
  !!============================================

end module proptab
  

!#######################################################################
!#######################################################################

module interaction

  use realprec
  
contains

  subroutine interactive()

    use proptab
    implicit none

    integer :: i,j,pt,stat,choiceend,choicemain,nreqminpr,&
         w,nreqsinglephas,ipha
    logical :: okchk
    integer :: nreqphas,nreqbulk
    real(dp) :: xvec,yvec,zvec,xlowb,xupb,ylowb,yupb,zlowb,zupb, &
         tol_bound_x,tol_bound_y
    real(dp), allocatable :: bulk_prop(:),phas_prop(:),propssinglemin(:)
    character(8)   :: astring
    character(len=12),allocatable :: req_phases(:)
    character(100) :: phaselstfile,outfile,xypts,outfilemin1,outfilemin2,xyzstring
    character(100),allocatable :: phaoutfile(:)
    character(8) :: phanu
    type(indexes) :: indx
    type(tabledims) :: tbd
    integer :: inp_req(3,50),inp_minpr(50),ienttable,nmaxpt,itab
    logical :: entiretable,phaseispresent
    character :: char

    
    tbd%nzfiles=-1
    inp_minpr=-1
    phaselstfile = 'list_req_phases.txt'

    
80  write(*,*) "_________________________________________________________"  
    write(*,*)


85  do itab=1,size(tbd%inptable) 
       write(*,*) "Enter the name of the input table file:"
       write(*,*) " * Leave blank to end *"
       read(*,fmt='(a1024)')  tbd%inptable(itab)(:)
       write(*,*) " "
       if (tbd%inptable(itab)=='') exit
       tbd%inptable(itab)=adjustl(tbd%inptable(itab))
       print *, "You requested: ",trim(tbd%inptable(itab))
       print *,""
    enddo
    tbd%nzfiles=itab-1
    print*, "You requested ",tbd%nzfiles," tables"
    if (tbd%nzfiles <= 0) then
       write(*,*) " Error, number of tables must be >= 1 "
       goto 85 
    endif


    write(*,*) " "
    write(*,*) "-----------------------------------------------------"


    !!-------------  LOAD TABLES --------------
    call redtab(propar, phasar, tbd, lstphas)   

    !! allocate vector of strings for requested phase names
    if ( allocated(req_phases) ) deallocate(req_phases)
    allocate(req_phases(tbd%totnumphases))
    req_phases = ""
    

90  write(*,*) "_________________________________________________________"  
    write(*,*)
    write(*,*) "What would you like to compute?"
    write(*,*)
    write(*,*) " [1] - Compute bulk (aggregate) properties"
    write(*,*) " [2] - Compute single phase(s) properties"
    write(*,*) ""
    write(*,*) "Only for debug:"
    write(*,*) " [3] - Extract phase properties from nearest nodes of the tables"
    read(*,*) choicemain 
    write(*,*) 
95  write(*,*) "Compute entire table OR input specific points from an x-y-z file?"
    write(*,*) " [1] - Compute entire table "
    write(*,*) " [2] - Input specific points from an ASCII file"
    read(*,*) ienttable
    if (ienttable == 1) then
       write(*,*) "You chose to compute the entire table."
       write(*,*)
       if (tbd%nzfiles == 1) then
          zvec=tbd%ztabs(1)
          write(*,*) "zvec is set to ",zvec
          write(*,*) " because there is only 1 input table (NO interpolation in the z direction)."
       else 
          write(*,*) "Enter the requested z (e.g., basalt fraction &
               &or other depending on input tables)"
          read(*,*) zvec
          write(*,*)
       endif
       entiretable=.true.
    else if (ienttable == 2) then 
       write(*,*) "Enter the filename containing the input x,y,z points:" 
       read(*,*)  xypts 
       write(*,*) 
       entiretable=.false.
    else  
       write(*,*) " Error, please choose 1 or 2."
       write(*,*) 
       goto 95
    endif

!!!################################################

    inp_req=0 !! zeroes all elements
    nreqbulk=0
    nreqphas=0
    nreqminpr=0

    select case(choicemain)
       
    case (1) !!############ Bulk phase properties  ############

       write(*,*) "Enter the name of the output file"
       read(*,*)  outfile
       write(*,*) 

       write(*,*) "Enter the requested properties"
       do i=1,tbd%nprop
          write(astring,'(i8)') i
          write(*,*) "["//trim(adjustl(astring))//"]-> ",trim(adjustl(tbd%nameprop(i)))
       enddo
       write(*,*) "0 - CONTINUE"

       nreqphas=0
       do i=1,size(inp_req,2)
110       read(*,*) inp_req(2,i)
          !! check for a valid entry
          if ( (inp_req(2,i) .lt. 0) .or. (inp_req(2,i) .gt. tbd%nprop) ) then
             print*,"Wrong entry number, retry:"
             go to 110
          endif
          !!-----end check-----------------
          if (inp_req(2,i)==0) then
             nreqphas=i-1
             exit
          endif
       enddo

       do i=1,nreqphas
          write(astring,'(i8)') inp_req(2,i)
          write(*,*) 'Which averaging scheme for "'//trim(adjustl(tbd%nameprop(inp_req(2,i) )))//'" ?'
          write(*,*)"[1=arithmetic, 2=VRH, 3=Voigt, 4=Reuss]"
120       read(*,*) inp_req(3,i)
          !! check for a valid entry
          if ( (inp_req(3,i) <1) .or. ( inp_req(3,i)>4 ) ) then
             print*,"Wrong entry number, retry:"
             go to 120
          end if
       enddo
       !!#### very important to initialize indx%vol_index,           ####
       !!#### otherwise problems 'cause in "mixprop" we may not know ####
       !!#### what to average on                                      ####
       indx%vol_index =-100

       !!if (inp_req(3,1) .ne. 0) then
       write(*,*) "Which property do you want to average on?"
       write(*,*) " (usually  phase,vo%) "
       do i=1,tbd%nprop
          write(astring,'(i8)') i
          write(*,*) "["//trim(adjustl(astring))//"]-> ",trim(adjustl(tbd%nameprop(i)))
       enddo
130    read(*,*) indx%vol_index
       !! check for a valid entry
       if ( (indx%vol_index .lt. 1) .or. (indx%vol_index .gt. tbd%nprop) ) then
          print*,"Wrong entry number, retry:"
          go to 130
       endif
       !!-----end check-----------------
       !!endif

       write(*,*) "You requested the following properties: "
       if (nreqphas/=0) then 
          write(*,*) ( (trim(adjustl(tbd%nameprop(inp_req(2,i))))) //"    ", i=1,nreqphas )
       endif
       write(*,*)


    case (2) !!############ Single phase(s) properties  ############

       write(*,*) "Enter the name of the output file"
       read(*,*)  outfile
       write(*,*) 

       write(*,*) "Enter the number of the phase(s) you are interested in:"
       write(*,*)
       write(*,*) "[a] ALL phases (present at requested conditions)"
       write(*,*) "    or choose one or more amongst: "
       do i=1,size(lstphas)
          write(astring,'(i8)') i
          write(*,*) "["//trim(adjustl(astring))//"]-> ",trim(adjustl(lstphas(i)))
       enddo
       write(*,*) "0 - CONTINUE"

       read(*,*) char
       if (char=='a') then
          do i=1,tbd%totnumphases
             req_phases(i) = trim(adjustl(lstphas(i)))
          end do
          nreqsinglephas=tbd%totnumphases
       else 
          write(ipha,'(a1)') char
          nreqsinglephas=1
          i=0 ! read before
          do while (nreqsinglephas<=tbd%totnumphases)
             i=i+1
             if (i>1) read(*,*) ipha
             !! check for a valid entry
             if ( (ipha .lt. 0) .or. (ipha .gt. tbd%totnumphases) ) then
                print*,"Wrong entry number, retry:"
                cycle
             else if ( (ipha .gt. 0) .and. (ipha .le. tbd%totnumphases) ) then 
                req_phases(i) = trim(adjustl(lstphas(ipha)))
                nreqsinglephas=nreqsinglephas+1
             endif
             !!-----end check-----------------
             if (ipha==0) then
                exit
             endif
          end do
       end if

       write(*,*)
       write(*,*) "You requested the following phases: "
       !!print*,'nreqsinglephas',nreqsinglephas
       write(*,*) ( (trim(req_phases(i))//"  "),i=1,nreqsinglephas)
       write(*,*)
       
       !! allocate file name array for several phases
       if (allocated(phaoutfile)) deallocate(phaoutfile)
       allocate(phaoutfile(nreqsinglephas))
           
       write(*,*) "Enter the requested properties"
       do i=1,tbd%nprop
          write(astring,'(i8)') i
          write(*,*) "["//trim(adjustl(astring))//"]-> ",trim(adjustl(tbd%nameprop(i)))
       enddo
       write(*,*) "0 - CONTINUE"
       nreqminpr=0
       do i=1,size(inp_minpr)
150       read(*,*) inp_minpr(i)
          !! check for a valid entry
          if ( (inp_minpr(i) .lt. 0) .or. (inp_minpr(i) .gt. tbd%nprop) ) then
             print*,"Wrong entry number, retry:"
             go to 150
          endif
          !!-----end check-----------------
          if (inp_minpr(i)==0) then
             nreqminpr=i-1
             exit
          endif
       enddo
       write(*,*)
       

    case (3)!!############ Phase properties from nearest nodes of the tables ############

       write(*,*) "Enter the **root** of the name of the output file"
       read(*,*)  outfile
       write(*,*) 

       write(*,*) "Enter the requested properties"
       do i=1,tbd%nprop
          write(astring,'(i8)') i
          write(*,*) "["//trim(adjustl(astring))//"]-> ",trim(adjustl(tbd%nameprop(i)))
       enddo
       write(*,*) "0 - CONTINUE"
       nreqminpr=0
       do i=1,size(inp_minpr)
105       read(*,*) inp_minpr(i)
          !! check for a valid entry
          if ( (inp_minpr(i) .lt. 0) .or. (inp_minpr(i) .gt. tbd%nprop) ) then
             print*,"Wrong entry number, retry:"
             go to 105
          endif
          !!-----end check-----------------
          if (inp_minpr(i)==0) then
             nreqminpr=i-1
             exit
          endif
       enddo  

    end select  !!!----------------------------

!!!############ END SELECT  ##############################

    if (tbd%nzfiles > 1) then
       xyzstring=tbd%xname//"  "//tbd%yname//"  z "
    else if (tbd%nzfiles.eq.1) then
       xyzstring=tbd%xname//"  "//tbd%yname
    endif

    !!  -------- open files, etc. on the ground of previous choice

    if ( choicemain.eq.1 ) then

       open(unit=22,file=outfile,status='replace')
       write(22,*) trim(adjustl(tbd%title))
       write(22,*) trim(adjustl(xyzstring)),"    ",( trim(adjustl(tbd%nameprop(inp_req(2,i)))) // "     ", i=1,nreqphas )


    else if (choicemain.eq.2) then

       if (allocated(propssinglemin)) deallocate(propssinglemin)
       allocate(propssinglemin(nreqminpr))

       write(*,*) "A list of requested phases and corresponding file number "
       write(*,*) " will be written to file: ",phaselstfile
   
       open(unit=60,file=phaselstfile,status='replace')
       
       do ipha = 1,nreqsinglephas
          write(phanu,"(i3)") ipha
          phaoutfile(ipha) = trim(trim(adjustl(outfile))//'_'//trim(adjustl(phanu)))
          write(60,*) phanu//'  '//req_phases(ipha)
          
          open(unit=40,file=phaoutfile(ipha),status='replace')
          write(40,*) trim(adjustl(tbd%title))
          write(40,*) req_phases(ipha) 
          write(40,*)trim(adjustl(xyzstring)),"    ",( trim(adjustl(tbd%nameprop(inp_minpr(i)))) // "     ", i=1,nreqminpr ) 
          close(40)
       end do
       close(60)
       
    else if (choicemain.eq.3)  then

       outfilemin1 = trim(adjustl(outfile)) // "_tab1.ascii" 
       outfilemin2 = trim(adjustl(outfile)) // "_tab2.ascii"
       print*,outfilemin1
       print*,outfilemin2
       !! to overwrite old files:
       open(unit=31,file=outfilemin1,status='replace')    
       open(unit=32,file=outfilemin2,status='replace')    
       close(31)
       close(32)
       
    endif

    ! !!######################## compute properties ###################################

    if ( entiretable.eqv..true. ) then
       ! !!############ compute properties for entire table  ######################

       nmaxpt=tbd%nxpts*tbd%nypts

       do j=0,tbd%nypts-1
          yvec=tbd%ymin+(j)*tbd%dy
          do i=0,tbd%nxpts-1
             xvec=tbd%xmin+(i)*tbd%dx

             if ( choicemain.eq.1 ) then
                call mixprop(propar,xvec,yvec,zvec,tbd,inp_req,indx &
                     ,bulk_prop,phas_prop ) 
                if (tbd%nzfiles == 1)  then
                   write(22,*) xvec,yvec,phas_prop
                else
                   write(22,*) xvec,yvec,zvec,phas_prop
                endif

             else if (choicemain.eq.2) then


                logical firsttspha(nreqsinglephas)
                
                do ipha = 1,nreqsinglephas
                   call extractsingleminpr(xvec,yvec,zvec,propar,phasar,tbd%nprop,inp_minpr,&
                        req_phases(ipha),tbd,propssinglemin,phaseispresent)
                   if (phaseispresent) then
                      
                      if (tbd%nzfiles == 1)  then
                         open(unit=40,file=phaoutfile(ipha),status='old',position='append')
                         write(40,*) xvec,yvec,propssinglemin
                         close(40)
                      else
                         open(unit=40,file=phaoutfile(ipha),status='old',position='append')
                         write(40,*) xvec,yvec,zvec,propssinglemin
                         close(40)
                      endif
                   end if
                end do
                
             else if (choicemain.eq.3) then
                call extractminpr(xvec,yvec,zvec,propar,phasar,tbd%nprop,inp_minpr,tbd,&
                     outfilemin1,outfilemin2) 

             endif

          enddo
       enddo

    else if ( entiretable.eqv..false. ) then
       ! !!############ compute properties for selected points  #####################

       open(unit=21,file=xypts,status='old')
       print*, "Reading  input x,y,z file..."

       tol_bound_x = 0.001*tbd%dx
       tol_bound_y  = 0.001*tbd%dy
       xlowb=tbd%xmin - tol_bound_x
       xupb=tbd%xmin+tbd%dx*(tbd%nxpts-1) + tol_bound_x
       ylowb=tbd%ymin - tol_bound_x
       yupb=tbd%ymin+tbd%dy*(tbd%nypts-1) + tol_bound_x
       zlowb=0.0
       zupb=tbd%nzfiles-1
       print*," Table bounds: "
       print*," x: ",xlowb,xupb
       print*," y: ",ylowb,yupb
       print*," z: ",zlowb,zupb
       print*, " "
       print*, "Calculating properties... "
       print*, " "

       !!#######################################################
       !! read until the end of file (->stat.lt.0)
       pt=0
       LOOP_xyzpts: do 
          read(21,*,iostat=stat) xvec,yvec,zvec
          pt=pt+1
          if (stat .lt. 0) exit 

          !! check for tables bounds
          if ( (xvec.lt.xlowb) .or. (xvec.gt.xupb) ) then 
             print*, "####  x  out of table's range.   ####"
             print*, " x = ",xvec
             print*, "####  Check x-y-z input file.    ####"
             print*, "####  Skipping this x-y-z point. ####"
             cycle LOOP_xyzpts 
          endif
          if ( (yvec.lt.ylowb) .or. (yvec.gt.yupb) ) then 
             print*, "####  y  out of table's range.  ####" 
             print*, " y = ",yvec
             print*, "####  Check x-y-z input file.   ####"
             print*, "####  Skipping this x-y-z point. ####"
             cycle LOOP_xyzpts 
          endif
          !!    mapping of z STILL WORK IN PROGRESS
          if ( (zvec.lt.zlowb) .or. (zvec.gt.zupb) ) then 
             print*, "####  z  out of table's range.  ####"
             print*, " z = ",zvec
             print*, "####  Check x-y-z input file.   ####"
             print*, "####  Skipping this x-y-z point. ####"
             cycle LOOP_xyzpts 
          endif

          
          if ( choicemain.eq.1 ) then
             !! now output phas_prop  
             call mixprop(propar,xvec,yvec,zvec,tbd,inp_req,indx &
                  ,bulk_prop,phas_prop ) 
             if (tbd%nzfiles == 1)  then
                write(22,*) xvec,yvec,phas_prop
             else
                write(22,*) xvec,yvec,zvec,phas_prop
             endif
             
          else if (choicemain.eq.2) then        
             do ipha = 1,nreqsinglephas
                call extractsingleminpr(xvec,yvec,zvec,propar,phasar,tbd%nprop,inp_minpr,&
                    req_phases(ipha),tbd,propssinglemin,phaseispresent)
                if (phaseispresent) then
                   if (tbd%nzfiles == 1)  then
                      open(unit=40,file=phaoutfile(ipha),status='old',position='append')
                      write(40,*) xvec,yvec,propssinglemin
                      close(40)
                   else
                      open(unit=40,file=phaoutfile(ipha),status='old',position='append')
                      write(40,*) xvec,yvec,zvec,propssinglemin
                      close(40)
                   endif
                end if
             end do
             
          else if (choicemain.eq.3) then
             call extractminpr(xvec,yvec,zvec,propar,phasar,tbd%nprop,inp_minpr,tbd,&
                  outfilemin1,outfilemin2)
          endif

       enddo LOOP_xyzpts  !! end loop on requested x-y-z points

       close(21) !! close x-y-z input file


    endif   !!! if ( entiretable.eqv..true. ) then


    !!#######################################################
    
    if ( choicemain.eq.1 ) then
       if (allocated(bulk_prop)) deallocate(bulk_prop)
       if (allocated(phas_prop)) deallocate(phas_prop)
       close(22) !! close output file
       print*,
       print*, "Output written to: ",outfile

    else if (choicemain.eq.2 ) then

       deallocate(propssinglemin)
       do ipha = 1,nreqsinglephas
          write(*,*) "Output written to: ",trim(phaoutfile(ipha))
       end do
       write(*,*)
       write(*,*) "Notice that the second line of each of the output files contains "
       write(*,*) " the name of the phase relative to that file number"
       
    else if (choicemain.eq.3 ) then
       write(*,*) 
       write(*,*) "Output written to: ",outfile

    endif


!!!############# back to main menu or exit  ##############
    choiceend=0
200 write(*,*) "______________________________________________"
    write(*,*)
    write(*,*) " Do you want to continue? "
    write(*,*) " [1] Go back to prompt to enter the input file "
    write(*,*) " [2] Go back to prompt for requested properties "
    write(*,*) " [3] Quit the program "
    read(*,*) choiceend

    if ( (choiceend.lt.1).or.(choiceend.gt.3)) then
       print*," Wrong number, please retry: "
       go to 200
    endif

    select case(choiceend)
    case (1)
       go to 80
    case (2)
       go to 90
    case (3)
       print*,"Bye bye!"
       stop
    end select

  end subroutine interactive

end module interaction


!#######################################################################
!#######################################################################

program interactivetab
  use interaction
  implicit none
  call interactive()
end program interactivetab

!#######################################################################
!#######################################################################
