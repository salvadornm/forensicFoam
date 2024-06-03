module ofgeom

implicit none


  integer, parameter :: NVTXMAX=500
  integer, parameter :: NBLKMAX=200

  integer,parameter :: inorth=1,ieast=2,isouth=3,iwest=4,itop=5,ibottom=6

  character*8, parameter :: cnorth="northbc",csouth="southbc",ctop="topbc",cbot="bottombc"
  character*8, parameter :: ceast ="eastbc" ,cwest ="westbc"

  integer,parameter :: iovtx=77,iofac=79,ioblck=78,ioTopo=65,ioPatch=66


  type vtx
    real  :: x,y,z
    integer :: id
  end type vtx

  type face
    integer :: p(4)
    integer :: id 
  end type face
 
  type block
    integer   :: p(0:7)   ! index of points
    integer :: id
    integer :: nx,ny,nz
    logical :: active
    type(face) :: north,east,west,south,top,bottom
  end type block

  type facelist
    integer      :: nfac(6)
    character*8  :: name(6)
    type(face)   :: f(NBLKMAX,6)    
  end type facelist 

  type vtxlist
    integer     :: nvtx
    type(vtx)   :: v(NVTXMAX)    
  end type vtxlist 

  type blocklist
    integer     :: nblocks
    type(block) :: b(NBLKMAX)
    type(vtx)   :: v(NVTXMAX)
    integer :: nfaces,nboundaries       
  end type blocklist


  type toposet
    character*8 :: name,patchname
    real :: box1(3),box2(3)
    integer :: tipe,action,source
  end type toposet

  integer, allocatable, private :: connect(:,:,:)

  contains

  ! -----------------------------------------------------
  ! creates vertex list and creates mapping i,j,k --> ijk
  ! -----------------------------------------------------
  subroutine create_listvtx(lvtx,xv,yv,zv,nx,ny,nz)
    implicit none

    type(vtxlist) :: lvtx
    real, intent(in) :: xv(:),yv(:),zv(:)
    integer :: nx,ny,nz
    ! local
    integer :: i,j,k, count

    allocate ( connect(nx,ny,nz) )

    count = 0
    do k = 1,nz
    do j = 1,ny
    do i = 1,nx
      count = count + 1
      lvtx%v(count)%id = count
      lvtx%v(count)%x = xv(i)
      lvtx%v(count)%y = yv(j)
      lvtx%v(count)%z = zv(k)            
      connect(i,j,k) = count 
    end do
    end do
    end do      
    lvtx%nvtx = count

  end subroutine create_listvtx
  ! -----------------------------------------------------
  subroutine print_listvtx(lvtx)

    implicit none

    type(vtxlist) :: lvtx

    write(*,*) " number vertices =",lvtx%nvtx
    write(*,*) " first vertex",lvtx%v(1)%x,lvtx%v(1)%y,lvtx%v(1)%z
    write(*,*) " last  vertex" ,lvtx%v(lvtx%nvtx)%x,lvtx%v(lvtx%nvtx)%y,lvtx%v(lvtx%nvtx)%z

    ! write(*,*) " connect =",connect(1,1,1), connect(2,1,1),connect(1,2,1),connect(2,2,1)

  end subroutine print_listvtx

  ! -----------------------------------------------------
  ! export vertex list to unit OpenFOAM style
  ! -----------------------------------------------------
  subroutine export_listvtx(unit,lvtx)

    implicit none
 
    integer, intent(in) :: unit
    type(vtxlist) :: lvtx 

    integer :: i

    write(unit,*) "vertices"
    write(unit,*) "( "    
    do i=1,lvtx%nvtx 
    write(unit,'(a,3f12.6,a)') "   (",lvtx%v(i)%x,lvtx%v(i)%y,lvtx%v(i)%z," )"
    end do
    write(unit,*) " );"    
    
  end subroutine export_listvtx
  ! ----------------------------------------------------
  ! creates block from mapping
  ! ----------------------------------------------------
  subroutine create_blockmap(b,i,j,k)

    implicit none

    type(block) :: b
    integer, intent(in) :: i,j,k

    b%nx = 0
    b%ny = 0
    b%nz = 0
    b%id = 0

    b%active = .true.
    
    ! order  sorting
    !     3 <--  2         7 <--  6
    !     |      |         |      | 
    !     |      |         |      |
    !     0  --> 1         4  --> 5

    if (.not. allocated(connect)) then
      write(*,*) " ERROR ... mapping vertex not defined"
      write(*,*) " call create_listvtx before"
      stop       
    end if  

    b%p(0) = connect(i,j,k) 
    b%p(1) = connect(i+1,j,k) 
    b%p(2) = connect(i+1,j+1,k) 
    b%p(3) = connect(i,j+1,k) 
    b%p(4) = connect(i,j,k+1) 
    b%p(5) = connect(i+1,j,k+1) 
    b%p(6) = connect(i+1,j+1,k+1) 
    b%p(7) = connect(i,j+1,k+1) 

    ! faces
    ! west  (0,4,7,3)
    b%west%p(1) = b%p(0)
    b%west%p(2) = b%p(4)
    b%west%p(3) = b%p(7)
    b%west%p(4) = b%p(3)
    ! east (1,2,6,5)
    b%east%p(1) = b%p(1)
    b%east%p(2) = b%p(2)
    b%east%p(3) = b%p(6)
    b%east%p(4) = b%p(5)
    ! south  (0,1,5,4)
    b%south%p(1) = b%p(0)
    b%south%p(2) = b%p(1)
    b%south%p(3) = b%p(5)
    b%south%p(4) = b%p(4)
    ! north  (3,7,6,2)
    b%north%p(1) = b%p(3)
    b%north%p(2) = b%p(7)
    b%north%p(3) = b%p(6)
    b%north%p(4) = b%p(2)
    ! bottom  (0,3,2,1)
    b%bottom%p(1) = b%p(0)
    b%bottom%p(2) = b%p(3)
    b%bottom%p(3) = b%p(2)
    b%bottom%p(4) = b%p(1)
    ! top  (4,5,6,7)
    b%top%p(1) = b%p(4)
    b%top%p(2) = b%p(5)
    b%top%p(3) = b%p(6)
    b%top%p(4) = b%p(7)    

  end subroutine create_blockmap  
  ! ----------------------------------------------------
  ! create list of blocks
  ! ----------------------------------------------------
  subroutine create_listblock(lblck,b,nx,ny,nz)
    implicit none

    type(blocklist) :: lblck 
    type(block), intent(in) :: b(:,:,:)
    integer, intent(in) :: nx,ny,nz
    ! local
    integer :: i,j,k, count

    count = 0
    do k = 1,nz
    do j = 1,ny
    do i = 1,nx
      if (b(i,j,k)%active) then
        count = count + 1
        lblck%b(count) = b(i,j,k)
      end if
    end do
    end do
    end do      
    lblck%nblocks = count
  end subroutine create_listblock
  ! ----------------------------------------------------
  ! print information onf one block
  ! ----------------------------------------------------
  subroutine print_block(b)

    implicit none

    type(block) :: b

    write(*,*) "  ****************************"
    write(*,*) " block  id=",b%id, "ACTIVE=",b%active
    write(*,*) " block  nx ny nz=",b%nx,b%ny,b%nz 
    write(*,*) " block  points (fortran notation)"
    write(*,'(A,8i3,A)') " hex (",b%p(0:7),")"
    write(*,*) " faces (fortran notation)"
    write(*,*) " N",b%north%p(:)
    write(*,*) " S",b%south%p(:)
    write(*,*) " E",b%east%p(:)
    write(*,*) " W",b%west%p(:)
    write(*,*) " T",b%top%p(:)
    write(*,*) " B",b%bottom%p(:)
    write(*,*) " OpenFOAM Notation"
    write(*,'(A,8i4,A,3i4,A)') " hex (",b%p(0:7)-1,") (",b%nx,b%ny,b%nz,") simpleGrading (1 1 1)"

    write(*,*) "  *************************** "
    
  end subroutine print_block
  ! -----------------------------------------------------
  ! export vertex list to unit OpenFOAM style (hence the -1)
  ! -----------------------------------------------------
  subroutine export_listblock(unit,lblck)

    implicit none
 
    integer, intent(in) :: unit
    type(blocklist) :: lblck 

    integer :: i

    write(unit,*) "blocks"
    write(unit,*) "( "    
    do i=1,lblck%nblocks
    write(unit,'(A,8i4,A,3i4,A)') " hex (",lblck%b(i)%p(0:7)-1,") (", & 
      & lblck%b(i)%nx,lblck%b(i)%ny,lblck%b(i)%nz,") simpleGrading (1 1 1)"
    end do
    write(unit,*) " );"    
    
  end subroutine export_listblock
  
  ! -----------------------------------------------------
  ! export face list to unit OpenFOAM style (hence the -1)
  ! -----------------------------------------------------
  subroutine export_listfaces(unit,lfac)

    implicit none
 
    integer, intent(in) :: unit
    type(facelist) :: lfac
    
    integer :: i,bc

    write(unit,*) "boundary" 
    write(unit,*) "( "   
    do bc = 1,6
      write(unit,*) lfac%name(bc)
      write(unit,*) "  {   "
      write(unit,*) "     type patch;   "
      write(unit,*) "      faces        "
      write(unit,*) "      (            "
      
      do i=1,lfac%nfac(bc)
        write(unit,'(A,4i4,A)') "     (",lfac%f(i,bc)%p(:)-1,")"
      end do  
      write(unit,*) "      );            "
      write(unit,*) "  }   "
    end do

   ! write(unit,*) ");"    

  end subroutine  export_listfaces
  ! -----------------------------------------------------
  ! export newly created face list to unit OpenFOAM style (hence the -1)
  ! -----------------------------------------------------
  subroutine export_newfaces(unit,lfac,bcname)

    implicit none
 
    integer, intent(in) :: unit
    type(facelist) :: lfac
    character(*) :: bcname

    integer :: i,bc

    write(unit,*) bcname
    write(unit,*) "  {   "
    write(unit,*) "     type patch;   "
    write(unit,*) "      faces        "
    write(unit,*) "      (            "
    
    do bc = 1,6
      do i=1,lfac%nfac(bc)
        write(unit,'(A,4i4,A)') "     (",lfac%f(i,bc)%p(:)-1,")"
      end do  
    end do      
    write(unit,*) "      );            "
    write(unit,*) "  }   "
    
    write(unit,*) ");"    

  end subroutine  export_newfaces
  ! -----------------------------------------------------
  ! export topoSet list to unit OpenFOAM style 
  ! -----------------------------------------------------
  subroutine export_faceSet(unit,set)

    implicit none
 
    integer, intent(in) :: unit
    type(toposet),intent(in) :: set
        
    write(unit,*) "  {   "
    write(unit,*) "name    ",trim(set%name)," ;"
    write(unit,*) "type    faceSet;"
    write(unit,*) "action  new;"
    write(unit,*) "source  boxToFace;" 
    write(unit,'(A,3F10.4,A,3F10.4,A)') " box (",set%box1(1:3)," )(",set%box2(1:3)," );"
    write(unit,*) "  }   "
    write(unit,*) "     "
    
  end subroutine   export_faceSet
  ! -----------------------------------------------------
  ! export topoSet list to unit OpenFOAM style 
  ! -----------------------------------------------------
  subroutine export_facePatch(unit,set)

    implicit none
 
    integer, intent(in) :: unit
    type(toposet), intent(in) :: set
        
    write(unit,*) " {   "
    write(unit,*) "name    ",trim(set%patchname)," ;"
    write(unit,*) "patchInfo"
    write(unit,*) "  {   "
    write(unit,*) "type patch;"
    write(unit,*) "  }   "
    write(unit,*) "constructFrom set;" 
    write(unit,*) "set  ",trim(set%name),";"
    write(unit,*) " }   "
    write(unit,*) "     "
    
  end subroutine   export_facePatch


  ! -----------------------------------------------------
  !  name and integer return8
  ! -----------------------------------------------------
  function setname(name,aux) result(sname)

    implicit none

    integer, intent(in) :: aux
    character*6, intent(in) :: name
    character*8 :: sname
    !
    character*6 ::str1

    write(str1, '(i0)') aux
    sname = trim(name) // trim(adjustl(str1)) 
      
  end function  setname

  


end module ofgeom
