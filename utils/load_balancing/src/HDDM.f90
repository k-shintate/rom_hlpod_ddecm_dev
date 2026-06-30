!> このプログラムは、monolis_solid の実行プロセス数と同じプロセス数で実行することが前提
!> コマンドライン引数 : deflation基底取得のための領域分割数
!
! 対応する partition 結果:
!  - graph_nedelec_elem_test.dat.* : nodal graph
!  - graph_nedelec_elem.dat.*      : conn graph
!  - graph_elem_conn.dat.*         : additional conn graph
!  - node_coordinate_elem.dat.*    : graph_nedelec_elem.dat 基準の real distval
!  - elem_bool.dat.*               : graph_nedelec_elem.dat 基準の int distval
!  - nedelec_elem.dat.*            : graph_nedelec_elem.dat 基準の int distval
!  - nedelec_edge_sign.dat.*       : graph_nedelec_elem.dat 基準の int distval
!  - D_bc_ned.dat.*                : graph_nedelec_elem_test.dat 基準の BC
!
! 重要:
!  graph_nedelec_elem.dat 基準の distval ごとに ndof を別管理する。
!    node_coordinate_elem.dat  -> ndof_coord
!    elem_bool.dat             -> ndof_bool
!    nedelec_elem.dat          -> ndof_nedelec_elem
!    nedelec_edge_sign.dat     -> ndof_sign
!
! MPI 並列時:
!  - metagraph.dat.id.<rank> を読み、rank ごとの domain_id list に従ってマージする
!  - n_graphs は metagraph.dat.id.<rank> のヘッダから取得する
!
! graph_elem_conn について:
!  - GEDATSU 側は n_conn_vertex=0 に未対応なので、driver 側で回避する
!  - local rank 内で graph_elem_conn が全て空なら、GEDATSU merge を呼ばずに空 graph を作る
!  - local rank 内に非空 graph がある場合は、非空 graph だけを抽出して GEDATSU merge に渡す
!  - 出力は常に行う。空 graph の場合も 0 graph ファイルを出力する
!
program merge_test_cio_full
  use mod_gedatsu
  use mod_monolis_utils
  implicit none

  logical, parameter :: FILE_IS_C0 = .true.
  logical, parameter :: LOG_ENABLE = .true.

  character(monolis_charlen), parameter :: IN_PART_DIR  = './parted.1/'
  character(monolis_charlen), parameter :: OUT_PART_DIR = './parted.0/'

  character(monolis_charlen), parameter :: NODAL_GRAPH_BASE  = 'graph_nedelec_elem_test.dat'
  character(monolis_charlen), parameter :: NEDELEC_CONN_BASE = 'graph_nedelec_elem.dat'
  character(monolis_charlen), parameter :: ELEM_CONN_BASE    = 'graph_elem_conn.dat'

  character(monolis_charlen), parameter :: NODE_COORD_BASE   = 'node_coordinate_elem.dat'
  character(monolis_charlen), parameter :: BC_BASE           = 'D_bc_ned.dat'

  character(monolis_charlen), parameter :: ELEM_BOOL_BASE    = 'elem_bool.dat'
  character(monolis_charlen), parameter :: NEDELEC_ELEM_BASE = 'nedelec_elem.dat'
  character(monolis_charlen), parameter :: NEDELEC_SIGN_BASE = 'nedelec_edge_sign.dat'

  type(gedatsu_graph), allocatable :: nodal_graphs(:)
  type(gedatsu_graph), allocatable :: conn_graphs(:)
  type(gedatsu_graph), allocatable :: elem_conn_graphs(:)

  ! graph_elem_conn merge 用の非空 graph 抽出配列
  type(gedatsu_graph), allocatable :: elem_conn_graphs_nonzero(:)
  type(gedatsu_graph), allocatable :: elem_conn_nodal_graphs_nonzero(:)

  type(gedatsu_graph) :: merged_nodal_graph
  type(gedatsu_graph) :: merged_conn_graph
  type(gedatsu_graph) :: merged_elem_conn_graph

  type(monolis_COM), allocatable :: monoCOMs(:)
  type(monolis_COM) :: merged_monoCOM

  type(monolis_list_I), allocatable :: n_dof_list_coord(:)
  type(monolis_list_R), allocatable :: list_coord(:)

  type(monolis_list_I), allocatable :: n_dof_list_bc(:)
  type(monolis_list_I), allocatable :: n_dof_list_load(:)

  type(monolis_list_I), allocatable :: list_ibc(:)
  type(monolis_list_I), allocatable :: list_iload(:)

  type(monolis_list_R), allocatable :: list_rbc(:)
  type(monolis_list_R), allocatable :: list_rload(:)

  !------------------------------------------------------------
  ! graph_nedelec_elem.dat 基準の I 系 distval
  !------------------------------------------------------------
  type(monolis_list_I), allocatable :: n_dof_list_bool(:)
  type(monolis_list_I), allocatable :: list_bool(:)
  integer(kint), allocatable :: merged_n_dof_list_bool(:)
  integer(kint), allocatable :: merged_array_bool(:)
  character(monolis_charlen) :: label_bool

  type(monolis_list_I), allocatable :: n_dof_list_nedelec_elem(:)
  type(monolis_list_I), allocatable :: list_nedelec_elem(:)
  integer(kint), allocatable :: merged_n_dof_list_nedelec_elem(:)
  integer(kint), allocatable :: merged_array_nedelec_elem(:)
  character(monolis_charlen) :: label_nedelec_elem

  type(monolis_list_I), allocatable :: n_dof_list_sign(:)
  type(monolis_list_I), allocatable :: list_sign(:)
  integer(kint), allocatable :: merged_n_dof_list_sign(:)
  integer(kint), allocatable :: merged_array_sign(:)
  character(monolis_charlen) :: label_sign

  integer(kint), allocatable :: arrayI_bool(:)
  integer(kint), allocatable :: arrayI_bool_val(:)

  integer(kint), allocatable :: arrayI_nedelec_elem(:)
  integer(kint), allocatable :: arrayI_nedelec_elem_val(:)

  integer(kint), allocatable :: arrayI_sign(:)
  integer(kint), allocatable :: arrayI_sign_val(:)

  integer(kint), allocatable :: arrayI_coord(:)
  real(kdouble), allocatable :: arrayR_coord(:)

  integer(kint) :: i, j, k, idx
  integer(kint) :: nnode, ndof
  integer(kint) :: nbc, nload
  integer(kint) :: iS, iE
  integer(kint) :: nnode_nodal_local
  integer(kint) :: nedge_local
  integer(kint) :: nedge_merged
  integer(kint) :: ncoord_merged
  integer(kint) :: nnode_merged

  ! graph_elem_conn 用
  integer(kint) :: n_elem_conn_local
  integer(kint) :: n_elem_conn_graphs_nonzero
  integer(kint) :: icount

  ! distval ごとの ndof
  integer(kint) :: ndof_coord
  integer(kint) :: ndof_bool
  integer(kint) :: ndof_nedelec_elem
  integer(kint) :: ndof_sign

  character(monolis_charlen) :: fname, cnum, label

  integer(kint), allocatable :: ibc(:,:), iload(:,:)

  integer(kint), allocatable :: merged_n_dof_list_coord(:)
  real(kdouble), allocatable :: merged_array_coord(:)

  integer(kint), allocatable :: merged_n_dof_list_bc(:)
  integer(kint), allocatable :: merged_n_dof_list_load(:)
  integer(kint), allocatable :: merged_array_ibc(:)
  integer(kint), allocatable :: merged_array_iload(:)

  integer(kint), allocatable :: tmpI(:), perm(:)
  real(kdouble), allocatable :: tmpR(:)

  real(kdouble), allocatable :: coord(:,:)
  real(kdouble), allocatable :: rbc(:)
  real(kdouble), allocatable :: rload(:)
  real(kdouble), allocatable :: merged_array_rbc(:)
  real(kdouble), allocatable :: merged_array_rload(:)

  integer(kint) :: ndomain_basis
  integer(kint) :: n_graphs
  integer(kint) :: funit
  character(monolis_charlen) :: ndomain_basis_
  integer(kint), allocatable :: domain_id(:)

  integer(kint) :: myrank, nprocs
  logical :: ex

  integer(kint) :: sum_internal

  call monolis_mpi_initialize()

  myrank = monolis_mpi_get_global_my_rank()
  nprocs = monolis_mpi_get_global_comm_size()

  call logp(myrank, "START program")

  !============================================================
  ! argv
  !============================================================
  call logp(myrank, "Read argv")

  if(command_argument_count() < 1) stop "Error : argv is need."

  call get_command_argument(1, ndomain_basis_)
  read(ndomain_basis_,*) ndomain_basis

  call logp(myrank, "ndomain_basis = "//trim(ndomain_basis_))

  if(ndomain_basis < nprocs) then
    write(*,*) "Error: ndomain_basis < nprocs"
    write(*,*) "ndomain_basis=", ndomain_basis
    write(*,*) "nprocs=", nprocs
    call monolis_mpi_finalize()
    stop
  endif

  if(ndomain_basis == 1 .and. nprocs == 1) then
    call monolis_mpi_finalize()
    stop
  endif

  !============================================================
  ! n_graphs
  !
  ! MPI 並列時は metagraph.dat.id.<rank> のヘッダから、
  ! この rank が読む domain 数を取得する。
  ! 1 proc 時は従来どおり argv の ndomain_basis 個を読む。
  !============================================================
  call logp(myrank, "Read n_graphs")

  if(nprocs == 1) then
    n_graphs = ndomain_basis
  else
    write(cnum, '(i0)') myrank
    fname = trim(OUT_PART_DIR)//"metagraph.dat.id."//trim(cnum)

    inquire(file=fname, exist=ex)

    if(.not. ex) then
      write(*,*) "Error: metagraph.dat.id.* not found: ", trim(fname)
      call monolis_mpi_finalize()
      stop
    endif

    open(newunit=funit, file=fname, status='old')
      read(funit,*) n_graphs
    close(funit)
  endif

  call logp_i(myrank, "n_graphs", n_graphs)

  !============================================================
  ! initialize / allocate
  !============================================================
  call gedatsu_graph_initialize(merged_nodal_graph)
  call gedatsu_graph_initialize(merged_conn_graph)
  call gedatsu_graph_initialize(merged_elem_conn_graph)
  call monolis_com_initialize_by_self(merged_monoCOM)

  allocate(nodal_graphs(n_graphs))
  allocate(conn_graphs(n_graphs))
  allocate(elem_conn_graphs(n_graphs))
  allocate(monoCOMs(n_graphs))

  allocate(n_dof_list_coord(n_graphs))
  allocate(list_coord(n_graphs))

  allocate(n_dof_list_bc(n_graphs))
  allocate(n_dof_list_load(n_graphs))

  allocate(list_ibc(n_graphs))
  allocate(list_iload(n_graphs))
  allocate(list_rbc(n_graphs))
  allocate(list_rload(n_graphs))

  call monolis_list_initialize_I(n_dof_list_coord, n_graphs)
  call monolis_list_initialize_R(list_coord,       n_graphs)

  call monolis_list_initialize_I(n_dof_list_bc,   n_graphs)
  call monolis_list_initialize_I(n_dof_list_load, n_graphs)

  call monolis_list_initialize_I(list_ibc,   n_graphs)
  call monolis_list_initialize_I(list_iload, n_graphs)
  call monolis_list_initialize_R(list_rbc,   n_graphs)
  call monolis_list_initialize_R(list_rload, n_graphs)

  allocate(n_dof_list_bool(n_graphs))
  allocate(list_bool(n_graphs))
  call monolis_list_initialize_I(n_dof_list_bool, n_graphs)
  call monolis_list_initialize_I(list_bool,       n_graphs)

  allocate(n_dof_list_nedelec_elem(n_graphs))
  allocate(list_nedelec_elem(n_graphs))
  call monolis_list_initialize_I(n_dof_list_nedelec_elem, n_graphs)
  call monolis_list_initialize_I(list_nedelec_elem,       n_graphs)

  allocate(n_dof_list_sign(n_graphs))
  allocate(list_sign(n_graphs))
  call monolis_list_initialize_I(n_dof_list_sign, n_graphs)
  call monolis_list_initialize_I(list_sign,       n_graphs)

  call monolis_alloc_I_1d(domain_id, n_graphs)

  !============================================================
  ! domain_id
  !
  ! domain_id は入力ファイル番号に使う。
  ! MPI 並列時は metagraph.dat.id.<rank> から読み込む。
  ! vertex_domain_id は GEDATSU 側の実装に合わせて myrank を入れる。
  !============================================================
  call logp(myrank, "Read domain_id list")

  if(nprocs == 1) then
    do i = 1, n_graphs
      domain_id(i) = i - 1
    enddo
  else
    write(cnum, '(i0)') myrank
    fname = trim(OUT_PART_DIR)//"metagraph.dat.id."//trim(cnum)

    inquire(file=fname, exist=ex)

    if(.not. ex) then
      write(*,*) "Error: metagraph.dat.id.* not found: ", trim(fname)
      call monolis_mpi_finalize()
      stop
    endif

    open(newunit=funit, file=fname, status='old')
      read(funit,*)
      read(funit,*)

      do i = 1, n_graphs
        read(funit,*) domain_id(i)
      enddo
    close(funit)
  endif

  call logp(myrank, "domain_id list read OK")

  !============================================================
  ! read loop
  !============================================================
  call logp(myrank, "Begin per-graph read loop")

  sum_internal = 0

  ndof_coord        = -1
  ndof_bool         = -1
  ndof_nedelec_elem = -1
  ndof_sign         = -1

  do i = 1, n_graphs
    call logp_i(myrank, "READ GRAPH i", i)
    call logp_i(myrank, "domain_id(i)", domain_id(i))

    call gedatsu_graph_initialize(nodal_graphs(i))
    call gedatsu_graph_initialize(conn_graphs(i))
    call gedatsu_graph_initialize(elem_conn_graphs(i))
    call monolis_com_initialize_by_self(monoCOMs(i))

    write(cnum, '(i0)') domain_id(i)

    !============================================================
    ! nodal graph
    !============================================================
    fname = trim(IN_PART_DIR)//trim(NODAL_GRAPH_BASE)//"."//trim(cnum)

    call logp(myrank, "read nodal graph: "//trim(fname))

    call monolis_input_graph(fname, &
      & nodal_graphs(i)%n_vertex, &
      & nodal_graphs(i)%vertex_id, &
      & nodal_graphs(i)%index, &
      & nodal_graphs(i)%item)

    if(FILE_IS_C0) then
      nodal_graphs(i)%item = nodal_graphs(i)%item + 1
    endif

    fname = trim(IN_PART_DIR)//trim(NODAL_GRAPH_BASE)//".n_internal."//trim(cnum)
    call monolis_input_internal_vertex_number(fname, nodal_graphs(i)%n_internal_vertex)

    call logp_i(myrank, "nodal n_internal(read)", nodal_graphs(i)%n_internal_vertex)

    sum_internal = sum_internal + nodal_graphs(i)%n_internal_vertex

    nnode_nodal_local = nodal_graphs(i)%n_vertex
    call logp_i(myrank, "nodal n_vertex", nnode_nodal_local)

    call monolis_dealloc_I_1d(nodal_graphs(i)%vertex_id)

    fname = trim(IN_PART_DIR)//trim(NODAL_GRAPH_BASE)//".id."//trim(cnum)
    call monolis_input_global_id(fname, nodal_graphs(i)%n_vertex, nodal_graphs(i)%vertex_id)

    call monolis_alloc_I_1d(nodal_graphs(i)%vertex_domain_id, nodal_graphs(i)%n_vertex)
    nodal_graphs(i)%vertex_domain_id(:) = myrank

    !============================================================
    ! com table
    !============================================================
    call logp(myrank, "read com tables for domain "//trim(cnum))

    call read_com_domain_suffix( &
      & monoCOMs(i), &
      & trim(IN_PART_DIR), &
      & trim(NODAL_GRAPH_BASE), &
      & trim(cnum), &
      & FILE_IS_C0)

    !============================================================
    ! graph_nedelec_elem.dat : conn graph
    !============================================================
    fname = trim(IN_PART_DIR)//trim(NEDELEC_CONN_BASE)//"."//trim(cnum)

    call logp(myrank, "read nedelec conn graph: "//trim(fname))

    call monolis_input_graph(fname, &
      & conn_graphs(i)%n_vertex, &
      & conn_graphs(i)%vertex_id, &
      & conn_graphs(i)%index, &
      & conn_graphs(i)%item)

    if(FILE_IS_C0) then
      conn_graphs(i)%item = conn_graphs(i)%item + 1
    endif

    fname = trim(IN_PART_DIR)//trim(NEDELEC_CONN_BASE)//".n_internal."//trim(cnum)
    call monolis_input_internal_vertex_number(fname, conn_graphs(i)%n_internal_vertex)

    call monolis_dealloc_I_1d(conn_graphs(i)%vertex_id)

    fname = trim(IN_PART_DIR)//trim(NEDELEC_CONN_BASE)//".id."//trim(cnum)
    call monolis_input_global_id(fname, conn_graphs(i)%n_vertex, conn_graphs(i)%vertex_id)

    call monolis_alloc_I_1d(conn_graphs(i)%vertex_domain_id, conn_graphs(i)%n_vertex)
    conn_graphs(i)%vertex_domain_id(:) = myrank

    nedge_local = conn_graphs(i)%n_vertex
    call logp_i(myrank, "nedelec conn n_vertex", nedge_local)

    !============================================================
    ! graph_elem_conn.dat : additional conn graph
    !============================================================
    fname = trim(IN_PART_DIR)//trim(ELEM_CONN_BASE)//"."//trim(cnum)

    inquire(file=fname, exist=ex)

    if(.not. ex) then
      write(*,*) "Error: graph_elem_conn file not found: ", trim(fname)
      call monolis_mpi_finalize()
      stop
    endif

    call logp(myrank, "read elem conn graph: "//trim(fname))

    call monolis_input_graph(fname, &
      & elem_conn_graphs(i)%n_vertex, &
      & elem_conn_graphs(i)%vertex_id, &
      & elem_conn_graphs(i)%index, &
      & elem_conn_graphs(i)%item)

    if(FILE_IS_C0) then
      elem_conn_graphs(i)%item = elem_conn_graphs(i)%item + 1
    endif

    fname = trim(IN_PART_DIR)//trim(ELEM_CONN_BASE)//".n_internal."//trim(cnum)
    call monolis_input_internal_vertex_number(fname, elem_conn_graphs(i)%n_internal_vertex)

    call monolis_dealloc_I_1d(elem_conn_graphs(i)%vertex_id)

    fname = trim(IN_PART_DIR)//trim(ELEM_CONN_BASE)//".id."//trim(cnum)
    call monolis_input_global_id(fname, elem_conn_graphs(i)%n_vertex, elem_conn_graphs(i)%vertex_id)

    call monolis_alloc_I_1d(elem_conn_graphs(i)%vertex_domain_id, elem_conn_graphs(i)%n_vertex)
    elem_conn_graphs(i)%vertex_domain_id(:) = myrank

    call logp_i(myrank, "elem conn n_vertex", elem_conn_graphs(i)%n_vertex)

    !============================================================
    ! node_coordinate_elem.dat : conn graph 基準, R distval
    !============================================================
    fname = trim(IN_PART_DIR)//trim(NODE_COORD_BASE)//"."//trim(cnum)

    call logp(myrank, "read node_coordinate_elem distval: "//trim(fname))

    call monolis_input_distval_R(fname, label, nnode, ndof, coord)

    if(ndof <= 0) then
      write(*,*) "Error: invalid ndof(node_coordinate_elem) read=", ndof, &
        & " file=", trim(fname), " rank=", myrank
      call monolis_mpi_finalize()
      stop "invalid ndof(node_coordinate_elem)"
    endif

    call set_or_check_ndof( &
      & myrank, &
      & trim(NODE_COORD_BASE), &
      & ndof, &
      & ndof_coord, &
      & fname)

    if(nnode /= nedge_local) then
      write(*,*) "Error: node_coordinate_elem nvertex mismatch: coord=", nnode, &
        & " conn=", nedge_local, " rank=", myrank, " file=", trim(fname)
      call monolis_mpi_finalize()
      stop
    endif

    call logp_i(myrank, "node_coordinate_elem nnode", nnode)
    call logp_i(myrank, "node_coordinate_elem ndof", ndof_coord)

    call monolis_dealloc_I_1d(arrayI_coord)
    call monolis_alloc_I_1d(arrayI_coord, nedge_local)

    arrayI_coord(:) = ndof_coord

    call monolis_list_set_I(n_dof_list_coord, i, nedge_local, arrayI_coord)

    call monolis_dealloc_R_1d(arrayR_coord)
    call monolis_alloc_R_1d(arrayR_coord, nedge_local*ndof_coord)

    do j = 1, nedge_local
      do k = 1, ndof_coord
        idx = ndof_coord*(j-1) + k
        arrayR_coord(idx) = coord(k,j)
      enddo
    enddo

    call monolis_list_set_R(list_coord, i, nedge_local*ndof_coord, arrayR_coord)

    !============================================================
    ! elem_bool.dat : conn graph 基準, I distval
    !============================================================
    fname = trim(IN_PART_DIR)//trim(ELEM_BOOL_BASE)//"."//trim(cnum)

    call logp(myrank, "read elem_bool (I, conn-based): "//trim(fname))

    call read_distval_fixndof_I_header2( &
      & fname, label_bool, nnode, arrayI_bool, arrayI_bool_val, ex)

    if(ex) then
      if(nnode /= nedge_local) then
        write(*,*) "Error: elem_bool nvertex mismatch: bool=", nnode, &
          & " conn=", nedge_local, " rank=", myrank, " file=", trim(fname)
        call monolis_mpi_finalize()
        stop
      endif

      call set_or_check_ndof_from_list( &
        & myrank, &
        & trim(ELEM_BOOL_BASE), &
        & nnode, &
        & arrayI_bool, &
        & ndof_bool, &
        & fname)

      call monolis_list_set_I(n_dof_list_bool, i, nedge_local, arrayI_bool)
      call monolis_list_set_I(list_bool, i, int(size(arrayI_bool_val), kint), arrayI_bool_val)

      call logp_i(myrank, "elem_bool ndof", ndof_bool)
      call logp(myrank, "elem_bool -> list OK")
    else
      call logp(myrank, "elem_bool not found: set empty")

      call monolis_dealloc_I_1d(arrayI_bool)
      call monolis_alloc_I_1d(arrayI_bool, nedge_local)
      arrayI_bool(:) = 0

      call monolis_dealloc_I_1d(arrayI_bool_val)
      call monolis_alloc_I_1d(arrayI_bool_val, 0)

      call monolis_list_set_I(n_dof_list_bool, i, nedge_local, arrayI_bool)
      call monolis_list_set_I(list_bool,       i, 0, arrayI_bool_val)
    endif

    !============================================================
    ! nedelec_elem.dat : conn graph 基準, I distval
    !============================================================
    fname = trim(IN_PART_DIR)//trim(NEDELEC_ELEM_BASE)//"."//trim(cnum)

    call logp(myrank, "read nedelec_elem (I, conn-based): "//trim(fname))

    call read_distval_fixndof_I_header2( &
      & fname, label_nedelec_elem, nnode, &
      & arrayI_nedelec_elem, arrayI_nedelec_elem_val, ex)

    if(ex) then
      if(nnode /= nedge_local) then
        write(*,*) "Error: nedelec_elem nvertex mismatch: nedelec_elem=", nnode, &
          & " conn=", nedge_local, " rank=", myrank, " file=", trim(fname)
        call monolis_mpi_finalize()
        stop
      endif

      call set_or_check_ndof_from_list( &
        & myrank, &
        & trim(NEDELEC_ELEM_BASE), &
        & nnode, &
        & arrayI_nedelec_elem, &
        & ndof_nedelec_elem, &
        & fname)

      call monolis_list_set_I( &
        & n_dof_list_nedelec_elem, i, nedge_local, arrayI_nedelec_elem)

      call monolis_list_set_I( &
        & list_nedelec_elem, i, int(size(arrayI_nedelec_elem_val), kint), arrayI_nedelec_elem_val)

      call logp_i(myrank, "nedelec_elem ndof", ndof_nedelec_elem)
      call logp(myrank, "nedelec_elem -> list OK")
    else
      call logp(myrank, "nedelec_elem not found: set empty")

      call monolis_dealloc_I_1d(arrayI_nedelec_elem)
      call monolis_alloc_I_1d(arrayI_nedelec_elem, nedge_local)
      arrayI_nedelec_elem(:) = 0

      call monolis_dealloc_I_1d(arrayI_nedelec_elem_val)
      call monolis_alloc_I_1d(arrayI_nedelec_elem_val, 0)

      call monolis_list_set_I( &
        & n_dof_list_nedelec_elem, i, nedge_local, arrayI_nedelec_elem)

      call monolis_list_set_I( &
        & list_nedelec_elem, i, 0, arrayI_nedelec_elem_val)
    endif

    !============================================================
    ! nedelec_edge_sign.dat : conn graph 基準, I distval
    !============================================================
    fname = trim(IN_PART_DIR)//trim(NEDELEC_SIGN_BASE)//"."//trim(cnum)

    call logp(myrank, "read nedelec_edge_sign (I, conn-based): "//trim(fname))

    call read_distval_fixndof_I_header2( &
      & fname, label_sign, nnode, arrayI_sign, arrayI_sign_val, ex)

    if(ex) then
      if(nnode /= nedge_local) then
        write(*,*) "Error: sign nvertex mismatch: sign=", nnode, &
          & " conn=", nedge_local, " rank=", myrank, " file=", trim(fname)
        call monolis_mpi_finalize()
        stop
      endif

      call set_or_check_ndof_from_list( &
        & myrank, &
        & trim(NEDELEC_SIGN_BASE), &
        & nnode, &
        & arrayI_sign, &
        & ndof_sign, &
        & fname)

      call monolis_list_set_I(n_dof_list_sign, i, nedge_local, arrayI_sign)
      call monolis_list_set_I(list_sign, i, int(size(arrayI_sign_val), kint), arrayI_sign_val)

      call logp_i(myrank, "nedelec_edge_sign ndof", ndof_sign)
      call logp(myrank, "nedelec_edge_sign -> list OK")
    else
      call logp(myrank, "nedelec_edge_sign not found: set empty")

      call monolis_dealloc_I_1d(arrayI_sign)
      call monolis_alloc_I_1d(arrayI_sign, nedge_local)
      arrayI_sign(:) = 0

      call monolis_dealloc_I_1d(arrayI_sign_val)
      call monolis_alloc_I_1d(arrayI_sign_val, 0)

      call monolis_list_set_I(n_dof_list_sign, i, nedge_local, arrayI_sign)
      call monolis_list_set_I(list_sign,       i, 0, arrayI_sign_val)
    endif

    !============================================================
    ! D_bc_ned.dat : nodal graph 基準, R BC
    !============================================================
    fname = trim(IN_PART_DIR)//trim(BC_BASE)//"."//trim(cnum)

    call logp(myrank, "read D_bc_ned: "//trim(fname))

    call monolis_input_bc_R(fname, nbc, ndof, ibc, rbc)

    if(ndof /= 1) then
      write(*,*) "Error: ndof(bc) read=", ndof, &
        & " file=", trim(fname), " rank=", myrank
      call monolis_mpi_finalize()
      stop "ndof isn't 1 (bc)"
    endif

    call logp_i(myrank, "nbc", nbc)

    if(FILE_IS_C0) then
      if(nbc > 0) ibc(1,:) = ibc(1,:) + 1
      if(nbc > 0) ibc(2,:) = ibc(2,:) + 1
    endif

    !------------------------------------------------------------
    ! BC -> list
    !------------------------------------------------------------
    call logp(myrank, "bc -> list (nodal-based)")

    call monolis_dealloc_I_1d(tmpI)
    call monolis_alloc_I_1d(tmpI, nnode_nodal_local)
    tmpI(:) = 0

    do j = 1, nbc
      if(ibc(1,j) < 1 .or. ibc(1,j) > nnode_nodal_local) then
        write(*,*) "Error: bc node id out of range: node=", ibc(1,j), &
          & " nnode=", nnode_nodal_local, " rank=", myrank, " file=", trim(fname)
        call monolis_mpi_finalize()
        stop
      endif

      tmpI(ibc(1,j)) = tmpI(ibc(1,j)) + 1
    enddo

    call monolis_list_set_I(n_dof_list_bc, i, nnode_nodal_local, tmpI)

    call monolis_dealloc_I_1d(tmpI)
    call monolis_dealloc_I_1d(perm)

    call monolis_alloc_I_1d(tmpI, nbc)
    call monolis_alloc_I_1d(perm, nbc)

    if(nbc > 0) then
      tmpI(:) = ibc(1,:)

      call monolis_get_sequence_array_I(perm, nbc, 1, 1)
      call monolis_qsort_I_2d(tmpI, perm, 1, nbc)

      do j = 1, nbc
        idx = perm(j)
        tmpI(j) = ibc(2,idx)
      enddo

      call monolis_list_set_I(list_ibc, i, nbc, tmpI)

      call monolis_dealloc_R_1d(tmpR)
      call monolis_alloc_R_1d(tmpR, nbc)

      do j = 1, nbc
        idx = perm(j)
        tmpR(j) = rbc(idx)
      enddo

      call monolis_list_set_R(list_rbc, i, nbc, tmpR)
    else
      call monolis_list_set_I(list_ibc, i, 0, tmpI)

      call monolis_dealloc_R_1d(tmpR)
      call monolis_alloc_R_1d(tmpR, 0)
      call monolis_list_set_R(list_rbc, i, 0, tmpR)
    endif

    !============================================================
    ! load は今回使わないので空扱い
    !============================================================
    call logp(myrank, "load is not used: set empty load list")

    nload = 0

    call monolis_dealloc_I_1d(tmpI)
    call monolis_alloc_I_1d(tmpI, nnode_nodal_local)
    tmpI(:) = 0

    call monolis_list_set_I(n_dof_list_load, i, nnode_nodal_local, tmpI)

    call monolis_dealloc_I_1d(tmpI)
    call monolis_alloc_I_1d(tmpI, 0)

    call monolis_list_set_I(list_iload, i, 0, tmpI)

    call monolis_dealloc_R_1d(tmpR)
    call monolis_alloc_R_1d(tmpR, 0)

    call monolis_list_set_R(list_rload, i, 0, tmpR)

    call logp_i(myrank, "DONE GRAPH i", i)
  enddo

  call logp(myrank, "End per-graph read loop")
  call logp_i(myrank, "sum n_internal(targets)", sum_internal)

  !============================================================
  ! merge
  !============================================================
  call logp(myrank, "Merge nodal subgraphs")

  call gedatsu_merge_nodal_subgraphs( &
    & n_graphs, nodal_graphs, monoCOMs, &
    & merged_nodal_graph, merged_monoCOM, ORDER_DOMAIN_ID)

  call logp_i(myrank, "sum_internal final", sum_internal)

  call logp(myrank, "Merge graph_nedelec_elem connectivity subgraphs")

  call gedatsu_merge_connectivity_subgraphs( &
    & n_graphs, nodal_graphs, merged_nodal_graph, merged_monoCOM, &
    & n_graphs, conn_graphs, merged_conn_graph)

  !============================================================
  ! graph_elem_conn
  !
  ! GEDATSU 側は n_conn_vertex=0 に未対応なので、
  ! driver 側で以下を行う:
  !
  ! 1. local rank 内の graph_elem_conn の n_vertex 合計を見る
  ! 2. 合計 0 なら GEDATSU merge を呼ばず、空 graph を作る
  ! 3. 合計 > 0 なら、n_vertex > 0 の graph だけを抽出して GEDATSU merge に渡す
  !============================================================
  n_elem_conn_local = 0
  n_elem_conn_graphs_nonzero = 0

  do i = 1, n_graphs
    n_elem_conn_local = n_elem_conn_local + elem_conn_graphs(i)%n_vertex

    if(elem_conn_graphs(i)%n_vertex > 0) then
      n_elem_conn_graphs_nonzero = n_elem_conn_graphs_nonzero + 1
    endif
  enddo

  call logp_i(myrank, "n_elem_conn_local", n_elem_conn_local)
  call logp_i(myrank, "n_elem_conn_graphs_nonzero", n_elem_conn_graphs_nonzero)

  if(n_elem_conn_local == 0) then
    call logp(myrank, "Create empty graph_elem_conn: local n_vertex = 0")

    call make_empty_graph(merged_elem_conn_graph)
  else
    call logp(myrank, "Merge graph_elem_conn connectivity subgraphs using nonzero graphs only")

    allocate(elem_conn_graphs_nonzero(n_elem_conn_graphs_nonzero))
    allocate(elem_conn_nodal_graphs_nonzero(n_elem_conn_graphs_nonzero))

    icount = 0

    do i = 1, n_graphs
      if(elem_conn_graphs(i)%n_vertex <= 0) cycle

      icount = icount + 1

      elem_conn_graphs_nonzero(icount) = elem_conn_graphs(i)
      elem_conn_nodal_graphs_nonzero(icount) = nodal_graphs(i)
    enddo

    call gedatsu_merge_connectivity_subgraphs( &
      & n_elem_conn_graphs_nonzero, elem_conn_nodal_graphs_nonzero, merged_nodal_graph, merged_monoCOM, &
      & n_elem_conn_graphs_nonzero, elem_conn_graphs_nonzero, merged_elem_conn_graph)

    deallocate(elem_conn_graphs_nonzero)
    deallocate(elem_conn_nodal_graphs_nonzero)
  endif

  ! node_coordinate_elem は conn graph 基準
  call logp(myrank, "Merge node_coordinate_elem distval R (conn-based)")

  call gedatsu_merge_distval_R( &
    & n_graphs, conn_graphs, merged_conn_graph, &
    & n_dof_list_coord, list_coord, &
    & merged_n_dof_list_coord, merged_array_coord)

  ! BC/load は nodal graph 基準
  call logp(myrank, "Merge bc/load distvals (nodal-based)")

  call gedatsu_merge_distval_I( &
    & n_graphs, nodal_graphs, merged_nodal_graph, &
    & n_dof_list_bc, list_ibc, &
    & merged_n_dof_list_bc, merged_array_ibc)

  call gedatsu_merge_distval_R( &
    & n_graphs, nodal_graphs, merged_nodal_graph, &
    & n_dof_list_bc, list_rbc, &
    & merged_n_dof_list_bc, merged_array_rbc)

  call gedatsu_merge_distval_I( &
    & n_graphs, nodal_graphs, merged_nodal_graph, &
    & n_dof_list_load, list_iload, &
    & merged_n_dof_list_load, merged_array_iload)

  call gedatsu_merge_distval_R( &
    & n_graphs, nodal_graphs, merged_nodal_graph, &
    & n_dof_list_load, list_rload, &
    & merged_n_dof_list_load, merged_array_rload)

  ! elem_bool / nedelec_elem / sign は conn graph 基準
  call logp(myrank, "Merge elem_bool distval I (conn-based)")

  call gedatsu_merge_distval_I( &
    & n_graphs, conn_graphs, merged_conn_graph, &
    & n_dof_list_bool, list_bool, &
    & merged_n_dof_list_bool, merged_array_bool)

  call logp(myrank, "Merge nedelec_elem distval I (conn-based)")

  call gedatsu_merge_distval_I( &
    & n_graphs, conn_graphs, merged_conn_graph, &
    & n_dof_list_nedelec_elem, list_nedelec_elem, &
    & merged_n_dof_list_nedelec_elem, merged_array_nedelec_elem)

  call logp(myrank, "Merge nedelec_edge_sign distval I (conn-based)")

  call gedatsu_merge_distval_I( &
    & n_graphs, conn_graphs, merged_conn_graph, &
    & n_dof_list_sign, list_sign, &
    & merged_n_dof_list_sign, merged_array_sign)

  !============================================================
  ! list -> node_coordinate_elem
  !============================================================
  call logp(myrank, "Rebuild node_coordinate_elem array")

  ncoord_merged = merged_conn_graph%n_vertex

  if(ndof_coord <= 0) then
    write(*,*) "Error: ndof_coord is invalid:", ndof_coord, " rank=", myrank
    call monolis_mpi_finalize()
    stop "invalid ndof_coord"
  endif

  call check_merged_ndof_list( &
    & myrank, &
    & trim(NODE_COORD_BASE), &
    & ncoord_merged, &
    & merged_n_dof_list_coord, &
    & ndof_coord)

  call monolis_dealloc_R_2d(coord)
  call monolis_alloc_R_2d(coord, ndof_coord, ncoord_merged)

  do i = 1, ncoord_merged
    do j = 1, ndof_coord
      idx = ndof_coord*(i-1) + j
      coord(j,i) = merged_array_coord(idx)
    enddo
  enddo

  !============================================================
  ! list -> BC
  !============================================================
  call logp(myrank, "Rebuild bc arrays")

  nnode_merged = merged_nodal_graph%n_vertex
  nbc = size(merged_array_rbc)

  call monolis_dealloc_I_2d(ibc)
  call monolis_dealloc_R_1d(rbc)

  call monolis_alloc_I_2d(ibc, 2, nbc)
  call monolis_alloc_R_1d(rbc, nbc)

  iS = 1

  do i = 1, nnode_merged
    j = merged_n_dof_list_bc(i)

    if(j /= 0) then
      iE = iS + j - 1

      ibc(1,iS:iE) = i
      ibc(2,iS:iE) = merged_array_ibc(iS:iE)
      rbc(iS:iE)   = merged_array_rbc(iS:iE)

      iS = iE + 1
    endif
  enddo

  !============================================================
  ! list -> load
  ! 今回は空のまま
  !============================================================
  call logp(myrank, "Rebuild empty load arrays")

  nload = size(merged_array_rload)

  call monolis_dealloc_I_2d(iload)
  call monolis_dealloc_R_1d(rload)

  call monolis_alloc_I_2d(iload, 2, nload)
  call monolis_alloc_R_1d(rload, nload)

  iS = 1

  do i = 1, nnode_merged
    j = merged_n_dof_list_load(i)

    if(j /= 0) then
      iE = iS + j - 1

      iload(1,iS:iE) = i
      iload(2,iS:iE) = merged_array_iload(iS:iE)
      rload(iS:iE)   = merged_array_rload(iS:iE)

      iS = iE + 1
    endif
  enddo

  !============================================================
  ! output
  !============================================================
  call logp(myrank, "Output merged common files")

  call logp_i(myrank, "ndof_coord",        ndof_coord)
  call logp_i(myrank, "ndof_bool",         ndof_bool)
  call logp_i(myrank, "ndof_nedelec_elem", ndof_nedelec_elem)
  call logp_i(myrank, "ndof_sign",         ndof_sign)

  call output_all_c( &
    & trim(OUT_PART_DIR), &
    & merged_nodal_graph, &
    & merged_conn_graph, &
    & merged_monoCOM, &
    & ncoord_merged, ndof_coord, coord, &
    & nbc, ibc, rbc, &
    & nload, iload, rload, &
    & FILE_IS_C0, sum_internal)

  !============================================================
  ! graph_elem_conn output
  !
  ! 空 graph の場合も output_graph_c 側で 0 graph を出力する。
  !============================================================
  call logp_i(myrank, "merged_elem_conn_graph n_vertex", merged_elem_conn_graph%n_vertex)
  call logp(myrank, "Output graph_elem_conn")

  call output_graph_c( &
    & trim(OUT_PART_DIR), &
    & trim(ELEM_CONN_BASE), &
    & merged_elem_conn_graph, &
    & FILE_IS_C0)

  nedge_merged = merged_conn_graph%n_vertex

  write(cnum, '(i0)') myrank

  call logp(myrank, "Output elem_bool / nedelec_elem / nedelec_edge_sign per-rank")

  call write_distval_fixndof_I_header2( &
    & trim(OUT_PART_DIR)//trim(ELEM_BOOL_BASE)//"."//trim(cnum), &
    & label_bool, nedge_merged, ndof_bool, &
    & merged_n_dof_list_bool, merged_array_bool)

  call write_distval_fixndof_I_header2( &
    & trim(OUT_PART_DIR)//trim(NEDELEC_ELEM_BASE)//"."//trim(cnum), &
    & label_nedelec_elem, nedge_merged, ndof_nedelec_elem, &
    & merged_n_dof_list_nedelec_elem, merged_array_nedelec_elem)

  call write_distval_fixndof_I_header2( &
    & trim(OUT_PART_DIR)//trim(NEDELEC_SIGN_BASE)//"."//trim(cnum), &
    & label_sign, nedge_merged, ndof_sign, &
    & merged_n_dof_list_sign, merged_array_sign)

  call logp(myrank, "Finalize MPI")

  call monolis_mpi_finalize()
  stop

contains

  !============================================================
  ! rank付きログ
  !============================================================
  subroutine logp(rank, msg)
    integer(kint), intent(in) :: rank
    character(*), intent(in) :: msg

    if(LOG_ENABLE) then
      write(*,'(A,I0,A,1X,A)') '[rank ', rank, ']', trim(msg)
      call flush(6)
    endif
  end subroutine logp

  subroutine logp_i(rank, key, val)
    integer(kint), intent(in) :: rank
    character(*), intent(in) :: key
    integer(kint), intent(in) :: val
    character(64) :: s

    write(s, '(I0)') val
    call logp(rank, trim(key)//" = "//trim(s))
  end subroutine logp_i

  !============================================================
  ! 空 graph を作る
  !
  ! n_vertex = 0
  ! n_internal_vertex = 0
  ! index は CSR 形式なので 1 要素だけ持たせて index(1)=0
  !============================================================
  subroutine make_empty_graph(g)
    type(gedatsu_graph), intent(inout) :: g

    call gedatsu_graph_initialize(g)

    g%n_vertex = 0
    g%n_internal_vertex = 0

    call monolis_alloc_I_1d(g%vertex_id, 0)
    call monolis_alloc_I_1d(g%vertex_domain_id, 0)

    call monolis_alloc_I_1d(g%index, 1)
    g%index(1) = 0

    call monolis_alloc_I_1d(g%item, 0)
  end subroutine make_empty_graph

  !============================================================
  ! distval ごとの ndof を保存・チェックする
  !============================================================
  subroutine set_or_check_ndof(rank, name, ndof_read, ndof_store, fname)
    integer(kint), intent(in) :: rank
    character(*), intent(in) :: name
    integer(kint), intent(in) :: ndof_read
    integer(kint), intent(inout) :: ndof_store
    character(*), intent(in) :: fname

    if(ndof_read < 0) then
      write(*,*) "Error: invalid ndof"
      write(*,*) "name=", trim(name)
      write(*,*) "ndof_read=", ndof_read
      write(*,*) "file=", trim(fname)
      write(*,*) "rank=", rank
      call monolis_mpi_finalize()
      stop
    endif

    if(ndof_store < 0) then
      ndof_store = ndof_read
    else
      if(ndof_read /= ndof_store) then
        write(*,*) "Error: ndof mismatch"
        write(*,*) "name=", trim(name)
        write(*,*) "read=", ndof_read
        write(*,*) "expected=", ndof_store
        write(*,*) "file=", trim(fname)
        write(*,*) "rank=", rank
        call monolis_mpi_finalize()
        stop
      endif
    endif
  end subroutine set_or_check_ndof

  !============================================================
  ! I系 distval の n_dof_node から固定 ndof を取得してチェック
  !============================================================
  subroutine set_or_check_ndof_from_list(rank, name, nnode, n_dof_node, ndof_store, fname)
    integer(kint), intent(in) :: rank
    character(*), intent(in) :: name
    integer(kint), intent(in) :: nnode
    integer(kint), intent(in) :: n_dof_node(:)
    integer(kint), intent(inout) :: ndof_store
    character(*), intent(in) :: fname

    integer(kint) :: i
    integer(kint) :: ndof_read

    if(nnode == 0) return

    ndof_read = n_dof_node(1)

    do i = 1, nnode
      if(n_dof_node(i) /= ndof_read) then
        write(*,*) "Error: variable ndof inside one fixed distval is not supported"
        write(*,*) "name=", trim(name)
        write(*,*) "file=", trim(fname)
        write(*,*) "i=", i
        write(*,*) "n_dof_node(i)=", n_dof_node(i)
        write(*,*) "expected=", ndof_read
        write(*,*) "rank=", rank
        call monolis_mpi_finalize()
        stop
      endif
    enddo

    call set_or_check_ndof(rank, name, ndof_read, ndof_store, fname)
  end subroutine set_or_check_ndof_from_list

  !============================================================
  ! merge 後の n_dof_list が固定 ndof になっているか確認
  !============================================================
  subroutine check_merged_ndof_list(rank, name, nnode, n_dof_node, ndof_expected)
    integer(kint), intent(in) :: rank
    character(*), intent(in) :: name
    integer(kint), intent(in) :: nnode
    integer(kint), intent(in) :: n_dof_node(:)
    integer(kint), intent(in) :: ndof_expected

    integer(kint) :: i

    do i = 1, nnode
      if(n_dof_node(i) /= ndof_expected) then
        write(*,*) "Error: merged ndof mismatch"
        write(*,*) "name=", trim(name)
        write(*,*) "i=", i
        write(*,*) "n_dof_node(i)=", n_dof_node(i)
        write(*,*) "expected=", ndof_expected
        write(*,*) "rank=", rank
        call monolis_mpi_finalize()
        stop
      endif
    enddo
  end subroutine check_merged_ndof_list

  !============================================================
  ! com read
  !============================================================
  subroutine read_com_domain_suffix(com, dir, base, suf, is_c0)
    type(monolis_COM), intent(inout) :: com
    character(*), intent(in) :: dir
    character(*), intent(in) :: base
    character(*), intent(in) :: suf
    logical, intent(in) :: is_c0

    character(monolis_charlen) :: fs, fr
    logical :: exs, exr

    fs = trim(dir)//trim(base)//".send."//trim(suf)
    fr = trim(dir)//trim(base)//".recv."//trim(suf)

    inquire(file=fs, exist=exs)
    inquire(file=fr, exist=exr)

    if(.not. (exs .and. exr)) then
      write(*,*) "Error: com file not found"
      write(*,*) "send:", trim(fs)
      write(*,*) "recv:", trim(fr)
      call monolis_mpi_finalize()
      stop
    endif

    call monolis_input_send_com_table(fs, com)
    call monolis_input_recv_com_table(fr, com)

    if(is_c0) then
      if(associated(com%send_item)) com%send_item = com%send_item + 1
      if(associated(com%recv_item)) com%recv_item = com%recv_item + 1
    endif
  end subroutine read_com_domain_suffix

  !============================================================
  ! 固定 ndof reader : I 版
  !
  ! 想定形式:
  !   label
  !   nnode ndof
  !   values...
  !============================================================
  subroutine read_distval_fixndof_I_header2(fname, label, nnode, n_dof_node, vals, exist)
    character(*), intent(in) :: fname
    character(monolis_charlen), intent(out) :: label
    integer(kint), intent(out) :: nnode
    integer(kint), allocatable, intent(out) :: n_dof_node(:)
    integer(kint), allocatable, intent(out) :: vals(:)
    logical, intent(out) :: exist

    integer :: unit, ios
    integer(kint) :: nn, nd, i
    logical :: ex_local

    inquire(file=fname, exist=ex_local)

    if(.not. ex_local) then
      exist = .false.
      label = ''
      nnode = 0

      call monolis_dealloc_I_1d(n_dof_node)
      call monolis_dealloc_I_1d(vals)

      call monolis_alloc_I_1d(n_dof_node, 0)
      call monolis_alloc_I_1d(vals, 0)

      return
    endif

    open(newunit=unit, file=fname, status='old')

    read(unit, '(A)', iostat=ios) label

    if(ios /= 0) then
      exist = .false.
      close(unit)
      return
    endif

    read(unit, *, iostat=ios) nn, nd

    if(ios /= 0 .or. nn < 0 .or. nd < 0) then
      exist = .false.
      close(unit)
      return
    endif

    nnode = nn

    call monolis_dealloc_I_1d(n_dof_node)
    call monolis_alloc_I_1d(n_dof_node, nnode)

    n_dof_node(:) = nd

    call monolis_dealloc_I_1d(vals)
    call monolis_alloc_I_1d(vals, nnode*nd)

    do i = 1, nnode
      if(nd > 0) then
        read(unit, *, iostat=ios) vals((i-1)*nd+1:i*nd)
      else
        read(unit, *, iostat=ios)
      endif

      if(ios /= 0) then
        exist = .false.
        close(unit)
        return
      endif
    enddo

    close(unit)

    exist = .true.
  end subroutine read_distval_fixndof_I_header2

  !============================================================
  ! 固定 ndof writer : I 版
  !
  ! n_dof_node(1) から ndof を決めず、distvalごとの ndof_out を使う。
  !============================================================
  subroutine write_distval_fixndof_I_header2(fname, label, nnode, ndof_out, n_dof_node, vals)
    character(*), intent(in) :: fname
    character(monolis_charlen), intent(in) :: label
    integer(kint), intent(in) :: nnode
    integer(kint), intent(in) :: ndof_out
    integer(kint), intent(in) :: n_dof_node(:)
    integer(kint), intent(in) :: vals(:)

    integer :: unit
    integer(kint) :: i

    if(nnode < 0) return

    if(ndof_out < 0) then
      write(*,*) "Error: ndof_out is invalid in writer"
      write(*,*) "file=", trim(fname)
      write(*,*) "ndof_out=", ndof_out
      stop
    endif

    do i = 1, nnode
      if(n_dof_node(i) /= ndof_out) then
        write(*,*) "Error: n_dof_node mismatch in writer"
        write(*,*) "file=", trim(fname)
        write(*,*) "i=", i
        write(*,*) "n_dof_node=", n_dof_node(i)
        write(*,*) "ndof_out=", ndof_out
        stop
      endif
    enddo

    if(int(size(vals), kint) /= nnode*ndof_out) then
      write(*,*) "Error: vals size mismatch in writer"
      write(*,*) "file=", trim(fname)
      write(*,*) "size(vals)=", size(vals)
      write(*,*) "expected=", nnode*ndof_out
      stop
    endif

    open(newunit=unit, file=fname, status='replace')

      write(unit, '(A)') trim(label)
      write(unit, *) nnode, ndof_out

      if(ndof_out > 0) then
        do i = 1, nnode
          write(unit, *) vals((i-1)*ndof_out+1:i*ndof_out)
        enddo
      else
        do i = 1, nnode
          write(unit, *)
        enddo
      endif

    close(unit)
  end subroutine write_distval_fixndof_I_header2

  !============================================================
  ! graph 出力 C形式
  !
  ! n_vertex = 0 の場合も、
  !   basename
  !   basename.n_internal
  !   basename.id
  ! を出力する。
  !============================================================
  subroutine output_graph_c(outdir, basename, g, is_c0)
    character(*), intent(in) :: outdir
    character(*), intent(in) :: basename
    type(gedatsu_graph), intent(in) :: g
    logical, intent(in) :: is_c0

    character(monolis_charlen) :: fname
    integer(kint), allocatable :: vertex_id_tmp(:)
    integer(kint), allocatable :: index_tmp(:)
    integer(kint), allocatable :: item_tmp(:)

    if(g%n_vertex == 0) then
      call monolis_alloc_I_1d(vertex_id_tmp, 0)
      call monolis_alloc_I_1d(index_tmp, 1)
      call monolis_alloc_I_1d(item_tmp, 0)

      index_tmp(1) = 0

      fname = monolis_get_global_output_file_name(outdir, "", trim(basename))
      call monolis_output_graph(fname, 0, vertex_id_tmp, index_tmp, item_tmp)

      fname = monolis_get_global_output_file_name(outdir, "", trim(basename)//".n_internal")
      call monolis_output_internal_vertex_number(fname, 0)

      fname = monolis_get_global_output_file_name(outdir, "", trim(basename)//".id")
      call monolis_output_global_id(fname, 0, vertex_id_tmp)

      call monolis_dealloc_I_1d(vertex_id_tmp)
      call monolis_dealloc_I_1d(index_tmp)
      call monolis_dealloc_I_1d(item_tmp)

      return
    endif

    call monolis_alloc_I_1d(vertex_id_tmp, size(g%vertex_id))
    call monolis_alloc_I_1d(index_tmp,    size(g%index))
    call monolis_alloc_I_1d(item_tmp,     size(g%item))

    vertex_id_tmp = g%vertex_id
    index_tmp     = g%index
    item_tmp      = g%item

    if(is_c0) item_tmp = item_tmp - 1

    fname = monolis_get_global_output_file_name(outdir, "", trim(basename))
    call monolis_output_graph(fname, g%n_vertex, vertex_id_tmp, index_tmp, item_tmp)

    fname = monolis_get_global_output_file_name(outdir, "", trim(basename)//".n_internal")
    call monolis_output_internal_vertex_number(fname, g%n_internal_vertex)

    fname = monolis_get_global_output_file_name(outdir, "", trim(basename)//".id")
    call monolis_output_global_id(fname, g%n_vertex, vertex_id_tmp)

    call monolis_dealloc_I_1d(vertex_id_tmp)
    call monolis_dealloc_I_1d(index_tmp)
    call monolis_dealloc_I_1d(item_tmp)
  end subroutine output_graph_c

  !============================================================
  ! 共通ファイル出力
  !
  ! graph_nedelec_elem_test.dat.n_internal は
  ! merged_nodal_graph%n_internal_vertex ではなく、
  ! 読み込んだ internal 数の合計 sum_internal をそのまま出力する。
  !============================================================
  subroutine output_all_c( &
      & outdir, gnod, gcon, com, &
      & ncoord, ndof_coord_out, coord_out, &
      & nbc, ibc, rbc, &
      & nload, iload, rload, &
      & is_c0, n_internal_out)

    character(*), intent(in) :: outdir
    type(gedatsu_graph), intent(in) :: gnod
    type(gedatsu_graph), intent(in) :: gcon
    type(monolis_COM), intent(in) :: com

    integer(kint), intent(in) :: ncoord
    integer(kint), intent(in) :: ndof_coord_out
    real(kdouble), intent(in) :: coord_out(:,:)

    integer(kint), intent(in) :: nbc
    integer(kint), intent(in) :: nload
    integer(kint), intent(in) :: ibc(:,:)
    integer(kint), intent(in) :: iload(:,:)
    real(kdouble), intent(in) :: rbc(:)
    real(kdouble), intent(in) :: rload(:)

    logical, intent(in) :: is_c0
    integer(kint), intent(in) :: n_internal_out

    character(monolis_charlen) :: fname, label
    integer(kint), allocatable :: item_tmp(:)
    integer(kint), allocatable :: index_tmp(:)
    integer(kint), allocatable :: ibc_tmp(:,:)

    type(monolis_COM) :: com_out

    !------------------------------------------------------------
    ! nodal graph
    !------------------------------------------------------------
    call monolis_alloc_I_1d(index_tmp, size(gnod%index))
    call monolis_alloc_I_1d(item_tmp,  size(gnod%item))

    index_tmp = gnod%index
    item_tmp  = gnod%item

    if(is_c0) item_tmp = item_tmp - 1

    fname = monolis_get_global_output_file_name(outdir, "", trim(NODAL_GRAPH_BASE))
    call monolis_output_graph(fname, gnod%n_vertex, gnod%vertex_id, index_tmp, item_tmp)

    fname = monolis_get_global_output_file_name(outdir, "", trim(NODAL_GRAPH_BASE)//".n_internal")
    call monolis_output_internal_vertex_number(fname, n_internal_out)

    fname = monolis_get_global_output_file_name(outdir, "", trim(NODAL_GRAPH_BASE)//".id")
    call monolis_output_global_id(fname, gnod%n_vertex, gnod%vertex_id)

    call monolis_dealloc_I_1d(index_tmp)
    call monolis_dealloc_I_1d(item_tmp)

    !------------------------------------------------------------
    ! com
    !------------------------------------------------------------
    com_out = com

    if(is_c0) then
      if(associated(com_out%send_item)) com_out%send_item = com_out%send_item - 1
      if(associated(com_out%recv_item)) com_out%recv_item = com_out%recv_item - 1
    endif

    fname = monolis_get_global_output_file_name(outdir, "", trim(NODAL_GRAPH_BASE)//".send")
    call monolis_output_send_com_table(fname, com_out)

    fname = monolis_get_global_output_file_name(outdir, "", trim(NODAL_GRAPH_BASE)//".recv")
    call monolis_output_recv_com_table(fname, com_out)

    !------------------------------------------------------------
    ! graph_nedelec_elem.dat
    !------------------------------------------------------------
    call monolis_alloc_I_1d(index_tmp, size(gcon%index))
    call monolis_alloc_I_1d(item_tmp,  size(gcon%item))

    index_tmp = gcon%index
    item_tmp  = gcon%item

    if(is_c0) item_tmp = item_tmp - 1

    fname = monolis_get_global_output_file_name(outdir, "", trim(NEDELEC_CONN_BASE))
    call monolis_output_graph(fname, gcon%n_vertex, gcon%vertex_id, index_tmp, item_tmp)

    fname = monolis_get_global_output_file_name(outdir, "", trim(NEDELEC_CONN_BASE)//".n_internal")
    call monolis_output_internal_vertex_number(fname, gcon%n_internal_vertex)

    fname = monolis_get_global_output_file_name(outdir, "", trim(NEDELEC_CONN_BASE)//".id")
    call monolis_output_global_id(fname, gcon%n_vertex, gcon%vertex_id)

    call monolis_dealloc_I_1d(index_tmp)
    call monolis_dealloc_I_1d(item_tmp)

    !------------------------------------------------------------
    ! node_coordinate_elem.dat
    !------------------------------------------------------------
    fname = monolis_get_global_output_file_name(outdir, "", trim(NODE_COORD_BASE))
    label = '#node_coordinate_elem'

    call monolis_output_distval_R(fname, label, ncoord, ndof_coord_out, coord_out)

    !------------------------------------------------------------
    ! D_bc_ned.dat
    !------------------------------------------------------------
    call monolis_alloc_I_2d(ibc_tmp, 2, nbc)

    if(nbc > 0) then
      ibc_tmp = ibc

      if(is_c0) then
        ibc_tmp(1,:) = ibc_tmp(1,:) - 1
        ibc_tmp(2,:) = ibc_tmp(2,:) - 1
      endif
    endif

    fname = monolis_get_global_output_file_name(outdir, "", trim(BC_BASE))
    call monolis_output_bc_R(fname, nbc, 1, ibc_tmp, rbc)

    call monolis_dealloc_I_2d(ibc_tmp)
  end subroutine output_all_c

end program merge_test_cio_full
