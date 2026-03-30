!> このプログラムは、monolis_solid の実行プロセス数と同じプロセス数で実行することが前提
!> コマンドライン引数 : deflation基底取得のための領域分割数
!
! 追加要件:
!  - elem_bool.dat / nedelec_edge_sign.dat は graph_nedelec_elem の id 基準で扱う（conn graph基準）
!  - elem_bool.dat / nedelec_edge_sign.dat は int の生データ（distval形式）なので、I系で読む・マージする
!  - どこまで実行したか分かるように、rank付きログを随所に出す
!  - 出力は rank毎に  ./parted.0/elem_bool.dat.<rank> / ./parted.0/nedelec_edge_sign.dat.<rank> へ出す
!
program merge_test_cio_full
  use mod_gedatsu
  use mod_monolis_utils
  implicit none

  logical, parameter :: FILE_IS_C0 = .true.   ! C形式で固定（自動判定なし）
  logical, parameter :: LOG_ENABLE = .true.   ! ログ出力ON/OFF

  type(gedatsu_graph), allocatable :: nodal_graphs(:), conn_graphs(:)
  type(gedatsu_graph) :: merged_nodal_graph, merged_conn_graph
  type(monolis_COM), allocatable :: monoCOMs(:)
  type(monolis_COM) :: merged_monoCOM

  type(monolis_list_I), allocatable :: n_dof_list_node(:), n_dof_list_bc(:), n_dof_list_load(:)
  type(monolis_list_I), allocatable :: list_ibc(:), list_iload(:)
  type(monolis_list_R), allocatable :: list_node(:), list_rbc(:), list_rload(:)

  !------------------------------------------------------------
  ! elem_bool / nedelec_edge_sign は int（I系）で扱う：conn_graph基準
  !------------------------------------------------------------
  type(monolis_list_I), allocatable :: n_dof_list_bool(:)
  type(monolis_list_I), allocatable :: list_bool(:)
  integer(kint), allocatable :: merged_n_dof_list_bool(:)
  integer(kint), allocatable :: merged_array_bool(:)
  character(monolis_charlen) :: label_bool

  type(monolis_list_I), allocatable :: n_dof_list_sign(:)
  type(monolis_list_I), allocatable :: list_sign(:)
  integer(kint), allocatable :: merged_n_dof_list_sign(:)
  integer(kint), allocatable :: merged_array_sign(:)
  character(monolis_charlen) :: label_sign

  integer(kint) :: i, j, k, idx, nnode, ndof, nbc, nload, iS, iE
  integer(kint) :: nnode_node, nedge_local, nedge_merged
  character(monolis_charlen) :: fname, cnum, label

  integer(kint), allocatable :: ibc(:,:), iload(:,:)
  integer(kint), allocatable :: merged_n_dof_list_node(:), merged_n_dof_list_bc(:), merged_n_dof_list_load(:)
  integer(kint), allocatable :: merged_array_ibc(:), merged_array_iload(:)

  ! bc/load の一時配列
  integer(kint), allocatable :: tmpI(:), perm(:)
  real(kdouble), allocatable :: tmpR(:)

  ! distvalごとの arrayI/arrayR（node）
  integer(kint), allocatable :: arrayI_node(:)
  real(kdouble), allocatable :: arrayR_node(:)

  ! distvalごとの arrayI/arrayVal（bool/sign）
  integer(kint), allocatable :: arrayI_bool(:), arrayI_sign(:)
  integer(kint), allocatable :: arrayI_bool_val(:), arrayI_sign_val(:)

  real(kdouble), allocatable :: node(:,:), rbc(:), rload(:)
  real(kdouble), allocatable :: merged_array_node(:), merged_array_rbc(:), merged_array_rload(:)

  integer(kint) :: ndomain_basis, n_graphs, funit
  character(monolis_charlen) :: ndomain_basis_
  integer(kint), allocatable :: domain_id(:)

  integer(kint) :: myrank, nprocs
  logical :: ex

  ! 追加：マージ対象の n_internal 合計を保持
  integer(kint) :: sum_internal

  call monolis_mpi_initialize()
  myrank = monolis_mpi_get_global_my_rank()
  nprocs = monolis_mpi_get_global_comm_size()

  call logp(myrank, "START program")

  !> 基底取得領域数の取得
  call logp(myrank, "Read argv")
  if(command_argument_count() < 1) stop "Error : argv is need."
  call get_command_argument(1, ndomain_basis_)
  read(ndomain_basis_,*) ndomain_basis
  call logp(myrank, "ndomain_basis = "//trim(ndomain_basis_))

  if(ndomain_basis < nprocs) stop "Error : ndomain_basis is invalid."
  if(ndomain_basis == 1 .and. nprocs == 1) then
    call monolis_mpi_finalize()
    stop
  endif

  !> n_graphs の取得
  call logp(myrank, "Read n_graphs")
  if(nprocs==1)then
    fname = "./parted.0/parted.1/metagraph.dat"
    inquire(file=fname, exist=ex); if(.not.ex) stop "Error: metagraph.dat not found"
    open(newunit=funit, file=fname, status='old')
      read(funit,*) n_graphs
    close(funit)
  else
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "metagraph.dat.n_internal")
    inquire(file=fname, exist=ex); if(.not.ex) stop "Error: metagraph.dat.n_internal not found"
    call monolis_input_internal_vertex_number(fname, n_graphs)
  endif
  call logp_i(myrank, "n_graphs", n_graphs)

  call gedatsu_graph_initialize(merged_nodal_graph)
  call gedatsu_graph_initialize(merged_conn_graph)
  call monolis_com_initialize_by_self(merged_monoCOM)

  allocate(nodal_graphs(n_graphs))
  allocate(conn_graphs(n_graphs))
  allocate(monoCOMs(n_graphs))

  allocate(n_dof_list_node(n_graphs))
  allocate(n_dof_list_bc(n_graphs))
  allocate(n_dof_list_load(n_graphs))
  allocate(list_node(n_graphs))
  allocate(list_ibc(n_graphs))
  allocate(list_iload(n_graphs))
  allocate(list_rbc(n_graphs))
  allocate(list_rload(n_graphs))

  call monolis_list_initialize_I(n_dof_list_node, n_graphs)
  call monolis_list_initialize_I(n_dof_list_bc,   n_graphs)
  call monolis_list_initialize_I(n_dof_list_load, n_graphs)
  call monolis_list_initialize_I(list_ibc,        n_graphs)
  call monolis_list_initialize_I(list_iload,      n_graphs)
  call monolis_list_initialize_R(list_node,       n_graphs)
  call monolis_list_initialize_R(list_rbc,        n_graphs)
  call monolis_list_initialize_R(list_rload,      n_graphs)

  allocate(n_dof_list_bool(n_graphs))
  allocate(list_bool(n_graphs))
  call monolis_list_initialize_I(n_dof_list_bool, n_graphs)
  call monolis_list_initialize_I(list_bool,       n_graphs)

  allocate(n_dof_list_sign(n_graphs))
  allocate(list_sign(n_graphs))
  call monolis_list_initialize_I(n_dof_list_sign, n_graphs)
  call monolis_list_initialize_I(list_sign,       n_graphs)

  call monolis_alloc_I_1d(domain_id, n_graphs)

  !> domain_id の取得
  call logp(myrank, "Read domain_id list")
  if(nprocs==1)then
    do i = 1, n_graphs
      domain_id(i) = i - 1
    enddo
  else
    write(cnum, '(i0)') myrank
    fname = './parted.0/metagraph.dat.id.' // trim(cnum)
    inquire(file=fname, exist=ex); if(.not.ex) stop "Error: metagraph.dat.id.* not found"
    open(newunit=funit, file=fname, status='old')
      read(funit,*)
      read(funit,*)
      do i = 1, n_graphs
        read(funit,*) domain_id(i)
      enddo
    close(funit)
  endif
  call logp(myrank, "domain_id list read OK")

  !> ---------------
  !> グラフ読み込み
  !> ---------------
  call logp(myrank, "Begin per-graph read loop")
  sum_internal = 0

  do i = 1, n_graphs
    call logp_i(myrank, "READ GRAPH i", i)

    call gedatsu_graph_initialize(nodal_graphs(i))
    call gedatsu_graph_initialize(conn_graphs(i))
    call monolis_com_initialize_by_self(monoCOMs(i))

    write(cnum,'(i0)') domain_id(i)

    !> ===== nodal graph =====
    fname = './parted.0/parted.1/graph.dat.' // trim(cnum)
    call logp(myrank, "read nodal graph: "//trim(fname))
    call monolis_input_graph(fname, &
      & nodal_graphs(i)%n_vertex, nodal_graphs(i)%vertex_id, nodal_graphs(i)%index, nodal_graphs(i)%item)

    if(FILE_IS_C0) then
      nodal_graphs(i)%item = nodal_graphs(i)%item + 1
    endif

    fname = './parted.0/parted.1/graph.dat.n_internal.' // trim(cnum)
    call monolis_input_internal_vertex_number(fname, nodal_graphs(i)%n_internal_vertex)
    call logp_i(myrank, "nodal n_internal(read)", nodal_graphs(i)%n_internal_vertex)

    ! 合計（要求：マージ対象の合計）
    sum_internal = sum_internal + nodal_graphs(i)%n_internal_vertex

    call monolis_dealloc_I_1d(nodal_graphs(i)%vertex_id)
    fname = './parted.0/parted.1/graph.dat.id.' // trim(cnum)
    call monolis_input_global_id(fname, nodal_graphs(i)%n_vertex, nodal_graphs(i)%vertex_id)

    call monolis_alloc_I_1d(nodal_graphs(i)%vertex_domain_id, nodal_graphs(i)%n_vertex)
    nodal_graphs(i)%vertex_domain_id(:) = myrank

    !> ===== com =====
    call logp(myrank, "read com tables for domain "//trim(cnum))
    call read_com_domain_suffix(monoCOMs(i), './parted.0/parted.1/', 'graph.dat', trim(cnum), FILE_IS_C0)

    !> ===== connectivity (graph_nedelec_elem) =====
    fname = './parted.0/parted.1/graph_nedelec_elem.dat.' // trim(cnum)
    call logp(myrank, "read conn graph: "//trim(fname))
    call monolis_input_graph(fname, &
      & conn_graphs(i)%n_vertex, conn_graphs(i)%vertex_id, conn_graphs(i)%index, conn_graphs(i)%item)

    if(FILE_IS_C0) then
      conn_graphs(i)%item = conn_graphs(i)%item + 1
    endif

    fname = './parted.0/parted.1/graph_nedelec_elem.dat.n_internal.' // trim(cnum)
    call monolis_input_internal_vertex_number(fname, conn_graphs(i)%n_internal_vertex)

    call monolis_dealloc_I_1d(conn_graphs(i)%vertex_id)
    fname = './parted.0/parted.1/graph_nedelec_elem.dat.id.' // trim(cnum)
    call monolis_input_global_id(fname, conn_graphs(i)%n_vertex, conn_graphs(i)%vertex_id)

    call monolis_alloc_I_1d(conn_graphs(i)%vertex_domain_id, conn_graphs(i)%n_vertex)
    conn_graphs(i)%vertex_domain_id(:) = myrank

    nedge_local = conn_graphs(i)%n_vertex
    call logp_i(myrank, "conn n_vertex", nedge_local)

    !> ===== distval (node) =====
    fname = './parted.0/parted.1/node_distval.dat.' // trim(cnum)
    call logp(myrank, "read node distval: "//trim(fname))
    call monolis_input_distval_R(fname, label, nnode, ndof, node)
    if(ndof /= 3) then
      write(*,*) "Error: ndof(node) read=", ndof, " file=", trim(fname), " rank=", myrank
      call monolis_mpi_finalize()
      stop "ndof isn't 3 (node)"
    endif
    nnode_node = nnode
    call logp_i(myrank, "node nnode", nnode_node)

    ! node distval -> list
    call logp(myrank, "node distval -> list")
    call monolis_dealloc_I_1d(arrayI_node)
    call monolis_alloc_I_1d(arrayI_node, nnode_node)
    arrayI_node(:) = ndof
    call monolis_list_set_I(n_dof_list_node, i, nnode_node, arrayI_node)

    call monolis_dealloc_R_1d(arrayR_node)
    call monolis_alloc_R_1d(arrayR_node, nnode_node*ndof)
    do j = 1, nnode_node
      do k = 1, ndof
        idx = ndof*(j-1) + k
        arrayR_node(idx) = node(k,j)
      enddo
    enddo
    call monolis_list_set_R(list_node, i, nnode_node*ndof, arrayR_node)

    !============================================================
    ! elem_bool.dat.<domain> : conn graph基準, int distval
    !============================================================
    fname = './parted.0/parted.1/elem_bool.dat.' // trim(cnum)
    call logp(myrank, "read elem_bool (I, conn-based): "//trim(fname))
    call read_distval_fixndof_I_header2(fname, label_bool, nnode, arrayI_bool, arrayI_bool_val, ex)

    if(ex) then
      if(nnode /= nedge_local) then
        write(*,*) "Error: elem_bool nvertex mismatch: bool=", nnode, " conn=", nedge_local, " rank=", myrank, " file=", trim(fname)
        call monolis_mpi_finalize()
        stop
      endif
      call monolis_list_set_I(n_dof_list_bool, i, nedge_local, arrayI_bool)
      call monolis_list_set_I(list_bool,       i, size(arrayI_bool_val), arrayI_bool_val)
      call logp(myrank, "elem_bool -> list OK")
    else
      call logp(myrank, "elem_bool not found (set empty)")
      call monolis_dealloc_I_1d(arrayI_bool)
      call monolis_alloc_I_1d(arrayI_bool, nedge_local)
      arrayI_bool(:)=0
      call monolis_dealloc_I_1d(arrayI_bool_val)
      call monolis_alloc_I_1d(arrayI_bool_val, 0)
      call monolis_list_set_I(n_dof_list_bool, i, nedge_local, arrayI_bool)
      call monolis_list_set_I(list_bool,       i, 0, arrayI_bool_val)
    endif

    !============================================================
    ! nedelec_edge_sign.dat.<domain> : conn graph基準, int distval
    !============================================================
    fname = './parted.0/parted.1/nedelec_edge_sign.dat.' // trim(cnum)
    call logp(myrank, "read sign (I, conn-based): "//trim(fname))
    call read_distval_fixndof_I_header2(fname, label_sign, nnode, arrayI_sign, arrayI_sign_val, ex)

    if(ex) then
      if(nnode /= nedge_local) then
        write(*,*) "Error: sign nvertex mismatch: sign=", nnode, " conn=", nedge_local, " rank=", myrank, " file=", trim(fname)
        call monolis_mpi_finalize()
        stop
      endif
      call monolis_list_set_I(n_dof_list_sign, i, nedge_local, arrayI_sign)
      call monolis_list_set_I(list_sign,       i, size(arrayI_sign_val), arrayI_sign_val)
      call logp(myrank, "sign -> list OK")
    else
      call logp(myrank, "sign not found (set empty)")
      call monolis_dealloc_I_1d(arrayI_sign)
      call monolis_alloc_I_1d(arrayI_sign, nedge_local)
      arrayI_sign(:)=0
      call monolis_dealloc_I_1d(arrayI_sign_val)
      call monolis_alloc_I_1d(arrayI_sign_val, 0)
      call monolis_list_set_I(n_dof_list_sign, i, nedge_local, arrayI_sign)
      call monolis_list_set_I(list_sign,       i, 0, arrayI_sign_val)
    endif

    !============================================================
    ! bc/load（node基準でマージ）
    !============================================================
    fname = './parted.0/parted.1/D_bc.dat.' // trim(cnum)
    call logp(myrank, "read bc: "//trim(fname))
    call monolis_input_bc_R(fname, nbc, ndof, ibc, rbc)
    if(ndof /= 1) then
      write(*,*) "Error: ndof(bc) read=", ndof, " file=", trim(fname), " rank=", myrank
      call monolis_mpi_finalize()
      stop "ndof isn't 1 (bc)"
    endif
    call logp_i(myrank, "nbc", nbc)

    fname = './parted.0/parted.1/D_bc.dat.' // trim(cnum)
    call logp(myrank, "read load: "//trim(fname))
    call monolis_input_bc_R(fname, nload, ndof, iload, rload)
    if(ndof /= 1) then
      write(*,*) "Error: ndof(load) read=", ndof, " file=", trim(fname), " rank=", myrank
      call monolis_mpi_finalize()
      stop "ndof isn't 1 (load)"
    endif
    call logp_i(myrank, "nload", nload)

    if(FILE_IS_C0) then
      if(nbc  > 0) ibc(1,:)   = ibc(1,:)   + 1
      if(nload> 0) iload(1,:) = iload(1,:) + 1
      if(nbc  > 0) ibc(2,:)   = ibc(2,:)   + 1
      if(nload> 0) iload(2,:) = iload(2,:) + 1
    endif

    ! bc -> list（node基準）
    call logp(myrank, "bc -> list (node-based)")
    call monolis_dealloc_I_1d(tmpI)
    call monolis_alloc_I_1d(tmpI, nnode_node)
    tmpI(:) = 0
    do j = 1, nbc
      tmpI(ibc(1,j)) = tmpI(ibc(1,j)) + 1
    enddo
    call monolis_list_set_I(n_dof_list_bc, i, nnode_node, tmpI)

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

    ! load -> list（node基準）
    call logp(myrank, "load -> list (node-based)")
    call monolis_dealloc_I_1d(tmpI)
    call monolis_alloc_I_1d(tmpI, nnode_node)
    tmpI(:) = 0
    do j = 1, nload
      tmpI(iload(1,j)) = tmpI(iload(1,j)) + 1
    enddo
    call monolis_list_set_I(n_dof_list_load, i, nnode_node, tmpI)

    call monolis_dealloc_I_1d(tmpI)
    call monolis_dealloc_I_1d(perm)
    call monolis_alloc_I_1d(tmpI, nload)
    call monolis_alloc_I_1d(perm, nload)
    if(nload > 0) then
      tmpI(:) = iload(1,:)
      call monolis_get_sequence_array_I(perm, nload, 1, 1)
      call monolis_qsort_I_2d(tmpI, perm, 1, nload)
      do j = 1, nload
        idx = perm(j)
        tmpI(j) = iload(2,idx)
      enddo
      call monolis_list_set_I(list_iload, i, nload, tmpI)

      call monolis_dealloc_R_1d(tmpR)
      call monolis_alloc_R_1d(tmpR, nload)
      do j = 1, nload
        idx = perm(j)
        tmpR(j) = rload(idx)
      enddo
      call monolis_list_set_R(list_rload, i, nload, tmpR)
    else
      call monolis_list_set_I(list_iload, i, 0, tmpI)
      call monolis_dealloc_R_1d(tmpR)
      call monolis_alloc_R_1d(tmpR, 0)
      call monolis_list_set_R(list_rload, i, 0, tmpR)
    endif

    call logp_i(myrank, "DONE GRAPH i", i)
  enddo
  call logp(myrank, "End per-graph read loop")
  call logp_i(myrank, "sum n_internal(targets)", sum_internal)

  ! ---------------
  ! 結合
  ! ---------------
  call logp(myrank, "Merge nodal subgraphs")
  call gedatsu_merge_nodal_subgraphs( &
    & n_graphs, nodal_graphs, monoCOMs, merged_nodal_graph, merged_monoCOM, ORDER_DOMAIN_ID)

  ! ★合計をログで確認
  call logp_i(myrank, "sum_internal (final)", sum_internal)

  call logp(myrank, "Merge connectivity subgraphs")
  call gedatsu_merge_connectivity_subgraphs( &
    & n_graphs, nodal_graphs, merged_nodal_graph, merged_monoCOM, &
    & n_graphs, conn_graphs,  merged_conn_graph)

  call logp(myrank, "Merge node distval (nodal-based)")
  call gedatsu_merge_distval_R(n_graphs, nodal_graphs, merged_nodal_graph, &
    & n_dof_list_node, list_node, merged_n_dof_list_node, merged_array_node)

  call logp(myrank, "Merge bc/load distvals (nodal-based)")
  call gedatsu_merge_distval_I(n_graphs, nodal_graphs, merged_nodal_graph, &
    & n_dof_list_bc, list_ibc, merged_n_dof_list_bc, merged_array_ibc)
  call gedatsu_merge_distval_R(n_graphs, nodal_graphs, merged_nodal_graph, &
    & n_dof_list_bc, list_rbc, merged_n_dof_list_bc, merged_array_rbc)

  call gedatsu_merge_distval_I(n_graphs, nodal_graphs, merged_nodal_graph, &
    & n_dof_list_load, list_iload, merged_n_dof_list_load, merged_array_iload)
  call gedatsu_merge_distval_R(n_graphs, nodal_graphs, merged_nodal_graph, &
    & n_dof_list_load, list_rload, merged_n_dof_list_load, merged_array_rload)

  call logp(myrank, "Merge elem_bool/sign (I) distvals (conn-based)")
  call gedatsu_merge_distval_I(n_graphs, conn_graphs, merged_conn_graph, &
    & n_dof_list_bool, list_bool, merged_n_dof_list_bool, merged_array_bool)

  call gedatsu_merge_distval_I(n_graphs, conn_graphs, merged_conn_graph, &
    & n_dof_list_sign, list_sign, merged_n_dof_list_sign, merged_array_sign)

  ! -------------------------
  ! list -> distval (node)
  ! -------------------------
  call logp(myrank, "Rebuild node distval array")
  nnode = merged_nodal_graph%n_vertex
  ndof  = 3

  call monolis_dealloc_R_2d(node)
  call monolis_alloc_R_2d(node, ndof, nnode)
  do i = 1, nnode
    do j = 1, ndof
      idx = ndof*(i-1) + j
      node(j,i) = merged_array_node(idx)
    enddo
  enddo

  ! -------------------------
  ! list -> bc
  ! -------------------------
  call logp(myrank, "Rebuild bc arrays")
  nbc = size(merged_array_rbc)

  call monolis_dealloc_I_2d(ibc)
  call monolis_dealloc_R_1d(rbc)
  call monolis_alloc_I_2d(ibc, 2, nbc)
  call monolis_alloc_R_1d(rbc, nbc)

  iS = 1
  do i = 1, nnode
    j = merged_n_dof_list_bc(i)
    if(j /= 0)then
      iE = iS + j - 1
      ibc(1,iS:iE) = i
      ibc(2,iS:iE) = merged_array_ibc(iS:iE)
      rbc(iS:iE)   = merged_array_rbc(iS:iE)
      iS = iE + 1
    endif
  enddo

  ! -------------------------
  ! list -> load
  ! -------------------------
  call logp(myrank, "Rebuild load arrays")
  nload = size(merged_array_rload)

  call monolis_dealloc_I_2d(iload)
  call monolis_dealloc_R_1d(rload)
  call monolis_alloc_I_2d(iload, 2, nload)
  call monolis_alloc_R_1d(rload, nload)

  iS = 1
  do i = 1, nnode
    j = merged_n_dof_list_load(i)
    if(j /= 0)then
      iE = iS + j - 1
      iload(1,iS:iE) = i
      iload(2,iS:iE) = merged_array_iload(iS:iE)
      rload(iS:iE)   = merged_array_rload(iS:iE)
      iS = iE + 1
    endif
  enddo

  !> ---------------
  !> 出力（C形式で出す）
  !> ---------------
  call logp(myrank, "Output merged (common) files")

  ! ★修正ポイント：
  !   graph.dat.n_internal は merged_nodal_graph%n_internal_vertex を信用せず
  !   sum_internal を「そのまま」出す
  call output_all_c("./parted.0/", merged_nodal_graph, merged_conn_graph, merged_monoCOM, &
    & nnode, ndof, node, nbc, ibc, rbc, nload, iload, rload, FILE_IS_C0, sum_internal)

  ! conn graph 由来の vertex数
  nedge_merged = merged_conn_graph%n_vertex
  write(cnum,'(i0)') myrank

  call logp(myrank, "Output elem_bool/sign per-rank")
  call write_distval_fixndof_I_header2( &
    & "./parted.0/elem_bool.dat."//trim(cnum), label_bool, nedge_merged, &
    & merged_n_dof_list_bool, merged_array_bool )

  call write_distval_fixndof_I_header2( &
    & "./parted.0/nedelec_edge_sign.dat."//trim(cnum), label_sign, nedge_merged, &
    & merged_n_dof_list_sign, merged_array_sign )

  call logp(myrank, "Finalize MPI")
  call monolis_mpi_finalize()
  stop

contains

  ! ------------------------
  ! rank付きログ
  ! ------------------------
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
    write(s,'(I0)') val
    call logp(rank, trim(key)//" = "//trim(s))
  end subroutine logp_i

  ! ---- com read (domain suffix send/recv) ----
  subroutine read_com_domain_suffix(com, dir, base, suf, is_c0)
    type(monolis_COM), intent(inout) :: com
    character(*), intent(in) :: dir, base, suf
    logical, intent(in) :: is_c0
    character(monolis_charlen) :: fs, fr
    logical :: exs, exr

    fs = trim(dir)//trim(base)//'.send.'//trim(suf)
    fr = trim(dir)//trim(base)//'.recv.'//trim(suf)
    inquire(file=fs, exist=exs)
    inquire(file=fr, exist=exr)
    if(.not.(exs .and. exr)) stop "Error: com file not found"

    call monolis_input_send_com_table(fs, com)
    call monolis_input_recv_com_table(fr, com)

    if(is_c0) then
      if(associated(com%send_item)) com%send_item = com%send_item + 1
      if(associated(com%recv_item)) com%recv_item = com%recv_item + 1
    endif
  end subroutine read_com_domain_suffix

  ! ============================================================
  ! 固定ndof reader（I版）
  ! ============================================================
  subroutine read_distval_fixndof_I_header2(fname, label, nnode, n_dof_node, vals, exist)
    character(*), intent(in) :: fname
    character(monolis_charlen), intent(out) :: label
    integer(kint), intent(out) :: nnode
    integer(kint), allocatable, intent(out) :: n_dof_node(:)
    integer(kint), allocatable, intent(out) :: vals(:)
    logical, intent(out) :: exist

    integer :: unit, ios
    integer(kint) :: nn, nd, i
    logical :: ex

    inquire(file=fname, exist=ex)
    if(.not.ex) then
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

    read(unit,'(A)', iostat=ios) label
    if(ios /= 0) then
      exist = .false.; close(unit); return
    endif

    read(unit,*, iostat=ios) nn, nd
    if(ios /= 0 .or. nn < 0 .or. nd < 0) then
      exist = .false.; close(unit); return
    endif
    nnode = nn

    call monolis_dealloc_I_1d(n_dof_node)
    call monolis_alloc_I_1d(n_dof_node, nnode)
    n_dof_node(:) = nd

    call monolis_dealloc_I_1d(vals)
    call monolis_alloc_I_1d(vals, nnode*nd)

    do i = 1, nnode
      read(unit,*, iostat=ios) vals((i-1)*nd+1:i*nd)
      if(ios /= 0) then
        exist = .false.; close(unit); return
      endif
    enddo

    close(unit)
    exist = .true.
  end subroutine read_distval_fixndof_I_header2

  ! ============================================================
  ! 固定ndof writer（I版）
  ! ============================================================
  subroutine write_distval_fixndof_I_header2(fname, label, nnode, n_dof_node, vals)
    character(*), intent(in) :: fname
    character(monolis_charlen), intent(in) :: label
    integer(kint), intent(in) :: nnode
    integer(kint), intent(in) :: n_dof_node(:)
    integer(kint), intent(in) :: vals(:)

    integer :: unit
    integer(kint) :: nd, i

    if(nnode == 0) return
    nd = n_dof_node(1)

    open(newunit=unit, file=fname, status='replace')
      write(unit,'(A)') trim(label)
      write(unit,*) nnode, nd
      do i = 1, nnode
        write(unit,*) vals((i-1)*nd+1:i*nd)
      enddo
    close(unit)
  end subroutine write_distval_fixndof_I_header2

  ! ---- 出力（C形式へ戻す） ----
  !
  ! ★修正:
  !   - graph.dat.n_internal の出力値を引数 n_internal_out で固定
  !   - com%send_item/recv_item を intent(in) から書き換えない（ローカルコピーで 0-start 化）
  !
  subroutine output_all_c(outdir, gnod, gcon, com, nnode, ndof, node, nbc, ibc, rbc, nload, iload, rload, is_c0, n_internal_out)
    character(*), intent(in) :: outdir
    type(gedatsu_graph), intent(in) :: gnod, gcon
    type(monolis_COM), intent(in) :: com
    integer(kint), intent(in) :: nnode, ndof, nbc, nload
    real(kdouble), intent(in) :: node(:,:)
    integer(kint), intent(in) :: ibc(:,:), iload(:,:)
    real(kdouble), intent(in) :: rbc(:), rload(:)
    logical, intent(in) :: is_c0
    integer(kint), intent(in) :: n_internal_out

    character(monolis_charlen) :: fname, label
    integer(kint), allocatable :: item_tmp(:)
    integer(kint), allocatable :: index_tmp(:)
    integer(kint), allocatable :: ibc_tmp(:,:), iload_tmp(:,:)

    ! com 出力用のローカルコピー（0-start化したいのは item だけ）
    type(monolis_COM) :: com_out

    call monolis_alloc_I_1d(index_tmp, size(gnod%index))
    call monolis_alloc_I_1d(item_tmp,  size(gnod%item))
    index_tmp = gnod%index
    item_tmp  = gnod%item
    if(is_c0) item_tmp = item_tmp - 1

    fname = monolis_get_global_output_file_name(outdir, "", "graph.dat")
    call monolis_output_graph(fname, gnod%n_vertex, gnod%vertex_id, index_tmp, item_tmp)

    fname = monolis_get_global_output_file_name(outdir, "", "graph.dat.n_internal")
    call monolis_output_internal_vertex_number(fname, n_internal_out)

    fname = monolis_get_global_output_file_name(outdir, "", "graph.dat.id")
    call monolis_output_global_id(fname, gnod%n_vertex, gnod%vertex_id)

    !> com（C形式へ戻す：itemだけ 0-start へ） ※ com を直接触らない
    com_out = com
    if(is_c0) then
      if(associated(com_out%send_item)) com_out%send_item = com_out%send_item - 1
      if(associated(com_out%recv_item)) com_out%recv_item = com_out%recv_item - 1
    endif

    fname = monolis_get_global_output_file_name(outdir, "", "graph.dat.send")
    call monolis_output_send_com_table(fname, com_out)
    fname = monolis_get_global_output_file_name(outdir, "", "graph.dat.recv")
    call monolis_output_recv_com_table(fname, com_out)

    ! ---- conn graph ----
    index_tmp = gcon%index
    call monolis_dealloc_I_1d(item_tmp)
    call monolis_alloc_I_1d(item_tmp, size(gcon%item))
    item_tmp  = gcon%item
    if(is_c0) item_tmp = item_tmp - 1

    fname = monolis_get_global_output_file_name(outdir, "", "graph_nedelec_elem.dat")
    call monolis_output_graph(fname, gcon%n_vertex, gcon%vertex_id, index_tmp, item_tmp)

    fname = monolis_get_global_output_file_name(outdir, "", "graph_nedelec_elem.dat.n_internal")
    call monolis_output_internal_vertex_number(fname, gcon%n_internal_vertex)

    fname = monolis_get_global_output_file_name(outdir, "", "graph_nedelec_elem.dat.id")
    call monolis_output_global_id(fname, gcon%n_vertex, gcon%vertex_id)

    fname = monolis_get_global_output_file_name(outdir, "", "node_distval.dat")
    label = '#node'
    call monolis_output_distval_R(fname, label, nnode, ndof, node)

    call monolis_alloc_I_2d(ibc_tmp, 2, nbc)
    call monolis_alloc_I_2d(iload_tmp, 2, nload)
    if(nbc>0) then
      ibc_tmp = ibc
      if(is_c0) then
        ibc_tmp(1,:) = ibc_tmp(1,:) - 1
        ibc_tmp(2,:) = ibc_tmp(2,:) - 1
      endif
    endif
    if(nload>0) then
      iload_tmp = iload
      if(is_c0) then
        iload_tmp(1,:) = iload_tmp(1,:) - 1
        iload_tmp(2,:) = iload_tmp(2,:) - 1
      endif
    endif

    fname = monolis_get_global_output_file_name(outdir, "", "D_bc.dat")
    call monolis_output_bc_R(fname, nbc, 1, ibc_tmp, rbc)

    fname = monolis_get_global_output_file_name(outdir, "", "D_bc.dat")
    call monolis_output_bc_R(fname, nload, 1, iload_tmp, rload)

    call monolis_dealloc_I_1d(index_tmp)
    call monolis_dealloc_I_1d(item_tmp)
    call monolis_dealloc_I_2d(ibc_tmp)
    call monolis_dealloc_I_2d(iload_tmp)
  end subroutine output_all_c

end program merge_test_cio_full
