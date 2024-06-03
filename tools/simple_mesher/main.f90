program meshof

   use ofgeom

   implicit none

   integer :: nbx,nby,nbz,nvx,nvy,nvz,ifac,im,ifacbc(6) 
   real :: lx,ly,lz,dx,dy,dz
   real :: x(NVTXMAX),y(NVTXMAX),z(NVTXMAX)

   type(vtxlist)  :: listv
   type(block), allocatable :: bmat(:,:,:)
   
   type(facelist) :: bcface,ibface
   type(blocklist) ::listb

   type(toposet), allocatable :: foam_toposet(:)

   logical :: constantmesh,check1,check2,mask, create_topoSet
   integer :: nmask=0,ntopoSet=0
   integer, allocatable :: vmask_start(:,:),vmask_end(:,:)
   character*6 :: auxname,auxname2

   integer :: i,j,k,ic,is,js,i_rack

   ! user vars
   real :: pos_hvac,pos_rack1,x_rack,x_sep,dbattery_x,dbattery_y
   real :: y0, mid_gap,side_gap,y_rack,pos_rack2,pos_rack4,pos_mod
   real :: zhvac,z_rack,xsep_rack_block,x_block,Lbattery,Lfan,Lbat_y
   real :: xtopo(100,3),pos_modulex,pos_moduley,pos_set,dx0,pos_rack(5)

   write(*,*) " Starting "
   write(*,*) " __________________________________________ "
   
   ! initialise soms vars
   x(:) = 0.0
   y(:) = 0.0
   z(:) = 0.0

   mask = .false.
   create_topoSet = .false.

   ! vtxmatrix  (mm)
   write(*,*) " Creating vtx ..."

   !========================  (problem dependent)
   ! X
   lx = 12192

   pos_hvac  =  985
   pos_rack1 =  1204
   x_rack    =  1368
   x_sep     =  410

   Lbattery  = x_rack/3.0
   Lfan = 100

   x(1) = 0.0                 !<<---------------
   
   x(2) =  pos_hvac           ! hvac 
   x(3) =  pos_rack1
   x(4) = x(3) + x_rack       ! rack 1 
   x(5) = x(4) + x_sep
   pos_rack2 = x(5)
   x(6) = x(5) + x_rack*2     ! rack 2 and 3
   x(7) = x(6) + x_sep
   pos_rack4 = x(7)
   x(8) = x(7) + x_rack*2     ! rack 4 and 5
   ! HVAC 2 TO DO
   xsep_rack_block = 146
   x_block = 1200!1000
   x(9)  = x(8) + xsep_rack_block ! block
   x(10) = x(9)  + x_block
   x(11) = lx  - 985            ! hvac 
   !
   x(12) = lx                  !<<---------------
   nvx  = 12

   write(*,*) " POINTS IN X"
   do i= 1,nvx
      write(*,*) i,x(i)
   end do 

   ! Y 
   ly = 2438.0

   y0   = 60  !origin of inside coordinates
   side_gap = 192
   y_rack = 725 
   mid_gap = 380
   y(1) = -0.5*ly             !<<---------------
   
   y(3) = -0.5*mid_gap
   y(2) = y(3) - y_rack
   
   y(4) = -y(3) 
   y(5) = -y(2) 

   y(6) = 0.5*ly               !<<---------------
   nvy = 6

   write(*,*) " POINTS IN Y"
   do j= 1,nvy
      write(*,*) j,y(j)
   end do   
   
   ! Z
   lz = 2896.0

   zhvac  = 2052
   z_rack = 2200

   z(1) = 0.0
   z(2) =  zhvac
   z(3) =  z_rack
   z(4) =  lz
   nvz = 4
   
   ! change to m
   x(:) = x(:)/1000.0
   y(:) = y(:)/1000.0
   z(:) = z(:)/1000.0

   ! spacing contant in each direction
   dx0 = 20.0 !  40 m


   dx = dx0/1000.0   
   dy = dx
   dz = dx
   constantmesh  = .true.
  
   ! masking
   mask=.true.

   if (mask) write(*,*) " Masking blocks ...   "
   ! number of masks regions
   ! hvac1 + racks(6) + 2 blocks + hvac2
   nmask = 1 + 6 + 2 +1
   ! allocate space for masks
   allocate (vmask_start(nmask,3),vmask_end(nmask,3))

   if (mask) write(*,*)  "                    ... DONE",nmask
    
   ! vertex index that define contiguous mask
   
   ! HVAC
   im = 1 
   vmask_start(im,1) = 1  !i
   vmask_end(im,1)   = 2

   vmask_start(im,2) = 2  !j
   vmask_end(im,2)   = 5

   vmask_start(im,3) = 1  !k
   vmask_end(im,3)   = 2
   
   ! RACKS 1 -3
   do im = 2,4   
      vmask_start(im,1) = 3 + 2*(im-2) ! im=2-> 3, im=3->5  im=4->7
      vmask_end(im,1)   = 4 + 2*(im-2) ! im=2-> 4        6        8
      vmask_start(im,2) = 2  
      vmask_end(im,2)   = 3 
      vmask_start(im,3) = 1
      vmask_end(im,3)   = 3
   end do

   ! RACKS 1 -3
   do im = 5,7   
      vmask_start(im,1) = 3 + 2*(im-5) 
      vmask_end(im,1)   = 4 + 2*(im-5) 
      vmask_start(im,2) = 4  
      vmask_end(im,2)   = 5 
      vmask_start(im,3) = 1
      vmask_end(im,3)   = 3
   end do

   ! BLOCKS  (masks 8-9)
   do im = 8,9
      ic = im - 8 ! 0 and 1
      vmask_start(im,1) = 9  !i
      vmask_end(im,1)   = 10
      vmask_start(im,2) = 2 + ic*2 
      vmask_end(im,2)   = 3 + ic*2               
      vmask_start(im,3) = 1  !k
      vmask_end(im,3)   = 3
   end do
   im = 10
   ! HVAC 2
   vmask_start(im,1) = 11  !i
   vmask_end(im,1)   = 12
   vmask_start(im,2) = 2 
   vmask_end(im,2)   = 5             
   vmask_start(im,3) = 1  !k
   vmask_end(im,3)   = 2

   create_topoSet = .true.

   ntopoSet  = 26*5*2

   allocate( foam_toposet(ntopoSet)  )

   ! racks 1 to 3
   Lfan = 100
   Lbat_y = 115

  

   write(*,*) "dxloc=",dx0," pos_set",pos_set

   auxname = "tfan"
   auxname2 = "fan"  !<---------- mpatch name

   dbattery_x = Lbattery
   dbattery_y = z_rack/9.0

   ! rack position
   pos_rack(1) = pos_rack1 
   pos_rack(2) = pos_rack2
   pos_rack(3) = pos_rack2 + x_rack
   pos_rack(4) = pos_rack4 
   pos_rack(5) = pos_rack4 + x_rack
   
   ! Racks 1-5
   im = 0
   do i_rack = 1,5 ! was 5
      pos_mod = pos_rack(i_rack)
      ! SOUTH
      pos_set =  -0.5*mid_gap - y_rack  
      do js = 1,9
      do is = 1,3      
         if ((is + js) > 2) then
            im = im + 1         
            pos_modulex =  pos_mod + 0.5*dbattery_x + (is -1)*dbattery_x
            pos_moduley =  0.5*dbattery_y + (js -1)*dbattery_y
            foam_toposet(im)%box1(1) = pos_modulex - Lfan/2.0
            foam_toposet(im)%box2(1) = pos_modulex + Lfan/2.0
            foam_toposet(im)%box1(2) = pos_set - dx0/3.0
            foam_toposet(im)%box2(2) = pos_set + dx0/3.0
            foam_toposet(im)%box1(3) = pos_moduley - Lfan/2.0
            foam_toposet(im)%box2(3) = pos_moduley + Lfan/2.0         
            ! Define topo and patch name
            foam_toposet(im)%name = setname(auxname,im) 
            foam_toposet(im)%patchname = setname(auxname2,im) 
         end if
      end do
      end do
      ! NORTH
      pos_set =  0.5*mid_gap + y_rack   
      do js = 1,9
      do is = 1,3      
         if ((is + js) > 2) then
            im = im + 1         
            pos_modulex =  pos_mod + 0.5*dbattery_x + (is -1)*dbattery_x
            pos_moduley =  0.5*dbattery_y + (js -1)*dbattery_y
            foam_toposet(im)%box1(1) = pos_modulex - Lfan/2.0
            foam_toposet(im)%box2(1) = pos_modulex + Lfan/2.0
            foam_toposet(im)%box1(2) = pos_set - dx0/3.0
            foam_toposet(im)%box2(2) = pos_set + dx0/3.0
            foam_toposet(im)%box1(3) = pos_moduley - Lfan/2.0
            foam_toposet(im)%box2(3) = pos_moduley + Lfan/2.0         
            ! Define topo and patch name
            foam_toposet(im)%name = setname(auxname,im) 
            foam_toposet(im)%patchname = setname(auxname2,im) 
         end if
      end do
      end do
   end do

   write(*,*) "TopSet=",ntopoSet," im=",im
   ntopoSet = im

   ! change from mm to m
   do im = 1,ntopoSet
      foam_toposet(im)%box1(:) = foam_toposet(im)%box1(:)/1000.0
      foam_toposet(im)%box2(:) = foam_toposet(im)%box2(:)/1000.0              
      ! print topoSets (two files)
      call export_faceSet(ioTopo,foam_toposet(im) )
      call export_facePatch(ioPatch,foam_toposet(im))
   end do   

   !stop









   


   !===============================================

   ! from vertices array create list
   call create_listvtx(listv,x,y,z,nvx,nvy,nvz)

   write(*,*) " number of points=",listv%nvtx

   write(*,*) "                    ... DONE"
   
   ! for debug
   !call print_listvtx(listv)
   
   call export_listvtx(iovtx,listv)

   write(*,*) " Creating block ..."
   ! blocks
   nbx = nvx  - 1
   nby = nvy  - 1
   nbz = nvz  - 1
   write(*,*) " number blocks x y z  ",nbx,nby,nbz
   allocate(bmat(nbx,nby,nbz))
   ! create blocks from vertex list and connect
   do i=1,nbx
   do j=1,nby
   do k=1,nbz
      call  create_blockmap(bmat(i,j,k),i,j,k)      
   end do
   end do
   end do
   write(*,*) "                    ... DONE"

   
   ! number of grid points
   if (constantmesh) then
      write(*,*) " constant spacing  delta ~",dx
   do i=1,nbx
   do j=1,nby
   do k=1,nbz      
      bmat(i,j,k)%nx = floor((x(i+1) - x(i))/dx)
      bmat(i,j,k)%ny = floor((y(j+1) - y(j))/dy)
      bmat(i,j,k)%nz = floor((z(k+1) - z(k))/dz)
   end do
   end do
   end do
   else
      write(*,*) " ERROR only constantmesh option valid"
      stop
   end if  

   ! check connectivity of grid points correct
   write(*,*) " checking for connectivty problems.."
   do i=1,nbx
   do j=1,nby
   do k=1,nbz   
      ! check x ------------
      if (i < nbx) then
         check1  = bmat(i,j,k)%ny .eq. bmat(i+1,j,k)%ny
         check1  = check1 .and. (bmat(i,j,k)%nz .eq. bmat(i+1,j,k)%nz)           
      else
         check1  = .true.
      end if      
      if (i > 1)  then  
         check2 = bmat(i,j,k)%ny .eq. bmat(i-1,j,k)%ny           
         check2 = check2 .and. (bmat(i,j,k)%nz .eq. bmat(i-1,j,k)%nz )
      else
         check2 = .true.
      end if   

      if (.not.(check1 .and. check2)) then 
         write(*,*) "ERROR ... x-dir",check1,check2
         write(*,*) "  in block i j k",i,j,k
         call print_block(bmat(i,j,k))
         stop
      end if   
      ! y
      if (j < nby) then
         check1  = bmat(i,j,k)%nx .eq. bmat(i,j+1,k)%nx
         check1  = check1 .and. (bmat(i,j,k)%nz .eq. bmat(i,j+1,k)%nz)           
      else
         check1  = .true.
      end if      
      if (j > 1)  then  
         check2 = bmat(i,j,k)%nx .eq. bmat(i,j-1,k)%nx           
         check2 = check2 .and. (bmat(i,j,k)%nz .eq. bmat(i,j-1,k)%nz )
      else
         check2 = .true.
      end if   

      if (.not.(check1 .and. check2)) then 
         write(*,*) "ERROR ... y-dir",check1,check2
         write(*,*) "  in block i j k",i,j,k
         call print_block(bmat(i,j,k))
         if (check1) call print_block(bmat(i,j-1,k))
         if (check2) call print_block(bmat(i,j+1,k))
         stop
      end if 
      ! z  
      if (k < nbz) then
         check1  = bmat(i,j,k)%nx .eq. bmat(i,j,k+1)%nx
         check1  = check1 .and. (bmat(i,j,k)%ny .eq. bmat(i,j,k+1)%ny)           
      else
         check1  = .true.
      end if      
      if (k > 1)  then  
         check2 = bmat(i,j,k)%nx .eq. bmat(i,j,k-1)%nx           
         check2 = check2 .and. (bmat(i,j,k)%ny .eq. bmat(i,j,k-1)%ny )
      else
         check2 = .true.
      end if   

      if (.not.(check1 .and. check2)) then 
         write(*,*) "ERROR ... z-dir",check1,check2
         write(*,*) "  in block i j k",i,j,k
         call print_block(bmat(i,j,k))
         stop
      end if 

      
    end do
    end do
   end do

   write(*,*) "                    ... NO ERRORS"
   

   if (mask) then
      write(*,*) " Masking blocks.."
      ! detect block to eliminate
      do im = 1,nmask
         do i = vmask_start(im,1),vmask_end(im,1)-1
         do j = vmask_start(im,2),vmask_end(im,2)-1  
         do k = vmask_start(im,3),vmask_end(im,3)-1 
            bmat(i,j,k)%active =.false.
         end do
         end do   
         end do
      end do 
      write(*,*) "                    ... DONE"
   end if   


   write(*,*) " Creating faces bc list ..."   

   ! collect external boundaries  (only if blocks are active)
   ! NORTH 
   ifac = 0   
   do i=1,nbx
   do k=1,nbz
      if (bmat(i,nby,k)%active) then
         ifac = ifac + 1      
         bcface%f(ifac,inorth) = bmat(i,nby,k)%north
      end if   
   end do
   end do 
   bcface%nfac(inorth) = ifac
   bcface%name(inorth) = cnorth
   ! SOUTH 
   ifac = 0
   do i=1,nbx
   do k=1,nbz      
      if (bmat(i,1,k)%active) then
         ifac = ifac + 1      
         bcface%f(ifac,isouth) = bmat(i,1,k)%south         
      end if   
   end do
   end do
   bcface%nfac(isouth) = ifac
   bcface%name(isouth)  = csouth
   ! EAST 
   ifac = 0
   do j=1,nby
   do k=1,nbz
      if (bmat(nbx,j,k)%active) then
         ifac = ifac + 1
         bcface%f(ifac,ieast) = bmat(nbx,j,k)%east
      end if
   end do
   end do  
   bcface%nfac(ieast) = ifac
   bcface%name(ieast)  = ceast
   ! WEST 
   ifac = 0
   do j=1,nby
   do k=1,nbz      
      if (bmat(1,j,k)%active) then
         ifac = ifac + 1      
         bcface%f(ifac,iwest) = bmat(1,j,k)%west         
      end if   
   end do
   end do
   bcface%nfac(iwest) = ifac
   bcface%name(iwest)  = cwest
   ! TOP
   ifac = 0
   do i=1,nbx
   do j=1,nby
      if (bmat(i,j,nbz)%active) then
         ifac = ifac + 1
         bcface%f(ifac,itop) = bmat(i,j,nbz)%top
      end if
   end do
   end do
   bcface%nfac(itop) = ifac
   bcface%name(itop)  = ctop
   ! BOT   
   ifac  = 0
   do i=1,nbx
   do j=1,nby   
      if (bmat(i,j,1)%active) then 
      ifac = ifac + 1      
      bcface%f(ifac,ibottom) = bmat(i,j,1)%bottom         
      end if
   end do
   end do
   bcface%nfac(ibottom) = ifac
   bcface%name(ibottom) = cbot

   write(*,*) " number of external faces="
   write(*,*) "  north=",bcface%nfac(inorth)," south =",bcface%nfac(isouth)
   write(*,*) "  east =",bcface%nfac(ieast) ," west  =",bcface%nfac(iwest)
   write(*,*) "  top  =",bcface%nfac(itop)  ," bottom=",bcface%nfac(ibottom)
   write(*,*) "------------------------------------------------------------"

   call export_listfaces(iofac,bcface)

         
   write(*,*) "                    ... DONE"
    

   ! for debug
   !call print_block(bmat(nbx,1,1))


   ! add new faces created after masking
   if (mask) then      
      write(*,*) " Creating new faces after masking..."   
      ifacbc(:) = 0
      do k= 1,nbz
      do j= 1,nby
      do i= 1,nbx
         ! detect block not.active
         if (.not. bmat(i,j,k)%active) then            
            if ((j < nby) .and.( bmat(i,j+1,k)%active)) then  ! add north
               ifacbc(inorth) = ifacbc(inorth) + 1
               ifac  =  ifacbc(inorth)
               ibface%f(ifac,inorth) = bmat(i,j,k)%north 
               ibface%nfac(inorth) = ifac               
            end if               
            
            if ((j > 1)   .and.( bmat(i,j-1,k)%active)) then  ! add south
               ifacbc(isouth) = ifacbc(isouth) + 1
               ifac = ifacbc(isouth)
               ibface%f(ifac,isouth) = bmat(i,j,k)%south
               ibface%nfac(isouth) = ifac               
            end if               
            if ((i < nbx) .and.( bmat(i+1,j,k)%active)) then  ! add east
               ifacbc(ieast) = ifacbc(ieast) + 1
               ifac = ifacbc(ieast)
               ibface%f(ifac,ieast) = bmat(i,j,k)%east
               ibface%nfac(ieast) = ifac               
            end if               
            if ((i > 1)   .and.( bmat(i-1,j,k)%active)) then  ! add west
               ifacbc(iwest) = ifacbc(iwest) + 1
               ifac = ifacbc(iwest)
               ibface%f(ifac,iwest) = bmat(i,j,k)%west
               ibface%nfac(iwest) = ifac               
            end if               
            if ((k < nbz) .and.( bmat(i,j,k+1)%active)) then  ! add top
               ifacbc(itop) = ifacbc(itop) + 1
               ifac = ifacbc(itop)
               ibface%f(ifac,itop) = bmat(i,j,k)%top
               ibface%nfac(itop) = ifac               
            end if               
            if ((k > 1)   .and.( bmat(i,j,k-1)%active)) then  ! add bottom
               ifacbc(ibottom) = ifacbc(ibottom) + 1
               ifac = ifacbc(ibottom)
               ibface%f(ifac,ibottom) = bmat(i,j,k)%bottom
               ibface%nfac(ibottom) = ifac               
            end if                           
         end if   
      end do
      end do
      end do        

      write(*,*) " number of internal faces after masking"
      write(*,*) "  north=",ibface%nfac(inorth)," south =",ibface%nfac(isouth)
      write(*,*) "  east =",ibface%nfac(ieast) ," west  =",ibface%nfac(iwest)
      write(*,*) "  top  =",ibface%nfac(itop)  ," bottom=",ibface%nfac(ibottom)
      write(*,*) "------------------------------------------------------------"
   

      call export_newfaces(iofac,ibface,"inner")

      write(*,*) "                    ... DONE"

   else
      
      write(iofac,*) ");" 

   end if   

   write(*,*) " Creating block list ..."

   ! from blockmap array create list block (removing inactive blocks)
   call create_listblock (listb,bmat,nbx,nby,nbz)

   write(*,*) " number of blocks=",listb%nblocks

   write(*,*) "                    ... DONE"

   call export_listblock(ioblck,listb)





   
 

   


   

   



end program meshof
