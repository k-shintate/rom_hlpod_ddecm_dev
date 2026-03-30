!======================================================================
! merge_test_param_optional.f90
!
! 目的:
!   - graph/conn/node/bc/load の「ファイル名(データ名)」と
!     「0-start(C0)/1-start(F1)フォーマット」を引数で与えられる merge プログラム
!   - ただし、--graph/--conn/--node/--bc/--load が "引数で与えられない" 場合は
!     その対象は「読み込みもしない」「マージもしない」「出力もしない」
!
! 前提:
!   - Fortranプログラムとして動作
!   - CSR item は Fortran配列参照に使うため、内部では必ず 1-start に正規化する
!     (入力が C0 の場合は +1 して内部保持)
!   - CSR index は offset 配列として扱う（入力 base を維持して良い想定）
!   - COM は send/recv を読む（関数名は環境に合わせて read_com_with_spec 内を調整）
!
! 例:
!   (graph + com + conn + node + bc + load 全部)
!   mpirun -np 4 ./merge_test_param_optional \
!     --ndomain 16 --top ./parted.0 --part parted.1 --out ./merged \
!     --graph graph.dat --conn connectivity.dat --node node.dat --bc bc.dat --load load.dat \
!     --nodal-csr C0,C0 --conn-csr C0,C0 --com C0,C0 --bcfmt C0,C0
!
!   (graph + com だけ)
!   mpirun -np 4 ./merge_test_param_optional \
!     --ndomain 16 --top ./parted.0 --part parted.1 --out ./merged \
!     --graph graph.dat --nodal-csr C0,C0 --com C0,C0
!
!   (node だけ) ※graph無しなら単独マージはできないので通常は意味が薄いが、
!              本コードでは "node 単体は読まずにスキップ" ではなく
!              nodal_graph が無いと node の結合先が無いのでエラーにする
!              (必要ならここも silent-skip に変えられます)
!======================================================================

module mod_merge_io_spec
  use mod_monolis_utils
  implicit none

  integer(kint), parameter :: FMT_C0 = 0
  integer(kint), parameter :: FMT_F1 = 1

  type :: t_csr_format
    integer(kint) :: index_base = FMT_F1
    integer(kint) :: item_base  = FMT_F1
  end type

  type :: t_bc_format
    integer(kint) :: node_base = FMT_F1
    integer(kint) :: dir_base  = FMT_F1
  end type

  type :: t_com_format
    integer(kint) :: index_base = FMT_F1
    integer(kint) :: item_base  = FMT_F1
    logical       :: has_domain_suffix = .true.
  end type

  type :: t_io_names
    character(monolis_charlen) :: top_dir   = "./parted.0"
    character(monolis_charlen) :: part_dir  = "parted.1"
    character(monolis_charlen) :: graph     = "graph.dat"
    character(monolis_charlen) :: conn      = "connectivity.dat"
    character(monolis_charlen) :: node      = "node.dat"
    character(monolis_charlen) :: bc        = "bc.dat"
    character(monolis_charlen) :: load      = "load.dat"
    character(monolis_charlen) :: meta_nint = "metagraph.dat.n_internal"
    character(monolis_charlen) :: meta_id_prefix = "metagraph.dat.id."
    character(monolis_charlen) :: out_dir   = "./parted.0"
  end type

  type :: t_merge_spec
    type(t_csr_format) :: nodal_csr
    type(t_csr_format) :: conn_csr
    type(t_com_format) :: com
    type(t_bc_format)  :: bc
    type(t_io_names)   :: name

    logical :: do_graph = .false.
    logical :: do_conn  = .false.
    logical :: do_node  = .false.
    logical :: do_bc    = .false.
    logical :: do_load  = .false.
    logical :: do_com   = .false.
  end type

contains

  pure function parse_base(s) result(base)
    character(*), intent(in) :: s
    integer(kint) :: base
    character(:), allocatable :: t
    t = adjustl(trim(s))
    select case(t)
    case("C0","c0","0")
      base = FMT_C0
    case("F1","f1","1")
      base = FMT_F1
    case default
      base = FMT_F1
    end select
  end function parse_base

  pure function join2(a, b) result(p)
    character(*), intent(in) :: a, b
    character(monolis_charlen) :: p
    character(monolis_charlen) :: aa, bb
    aa = trim(a); bb = trim(b)
    if(len_trim(aa) == 0) then
      p = bb
      return
    endif
    if(aa(len_trim(aa):len_trim(aa)) == "/") then
      p = trim(aa)//trim(bb)
    else
      p = trim(aa)//"/"//trim(bb)
    endif
  end function join2

  pure function join3(a, b, c) result(p)
    character(*), intent(in) :: a, b, c
    character(monolis_charlen) :: p
    p = join2(join2(a,b), c)
  end function join3

  subroutine get_arg_value(key, val, found)
    character(*), intent(in) :: key
    character(monolis_charlen), intent(out) :: val
    logical, intent(out) :: found
    integer :: n, i
    character(monolis_charlen) :: a
    found = .false.
    val = ""
    n = command_argument_count()
    if(n < 1) return
    do i = 1, n-1
      call get_command_argument(i, a)
      if(trim(a) == trim(key)) then
        call get_command_argument(i+1, val)
        found = .true.
        return
      endif
    enddo
  end subroutine get_arg_value

  subroutine split2_csv(s, a, b)
    character(*), intent(in) :: s
    character(monolis_charlen), intent(out) :: a, b
    integer :: p
    a = ""; b = ""
    p = index(trim(s), ",")
    if(p <= 0) then
      a = trim(s)
      b = ""
    else
      a = trim(s(1:p-1))
      b = trim(s(p+1:))
    endif
  end subroutine split2_csv

  subroutine split3_csv(s, a, b, c)
    character(*), intent(in) :: s
    character(monolis_charlen), intent(out) :: a, b, c
    character(monolis_charlen) :: rest
    call split2_csv(s, a, rest)
    call split2_csv(rest, b, c)
  end subroutine split3_csv

end module mod_merge_io_spec

!======================================================================

program merge_test_param_optional
  use mod_gedatsu
  use mod_monolis_utils
  use mod_merge_io_spec
  implicit none

  type(t_merge_spec) :: spec

  type(gedatsu_graph), allocatable :: nodal_graphs(:), conn_graphs(:)
  type(gedatsu_graph) :: merged_nodal_graph, merged_conn_graph
  type(monolis_COM), allocatable :: monoCOMs(:)
  type(monolis_COM) :: merged_monoCOM

  type(monolis_list_I), allocatable :: n_dof_list_node(:), n_dof_list_bc(:), n_dof_list_load(:)
  type(monolis_list_I), allocatable :: list_ibc(:), list_iload(:)
  type(monolis_list_R), allocatable :: list_node(:), list_rbc(:), list_rload(:)

  integer(kint) :: i, j, k, idx, nnode, ndof, nbc, nload
  character(monolis_charlen) :: fname, cnum, label

  integer(kint), allocatable :: ibc(:,:), iload(:,:)
  integer(kint), allocatable :: merged_n_dof_list_node(:), merged_n_dof_list_bc(:), merged_n_dof_list_load(:)
  integer(kint), allocatable :: merged_array_ibc(:), merged_array_iload(:)
  real(kdouble), allocatable :: node(:,:), rbc(:), rload(:)
  real(kdouble), allocatable :: merged_array_node(:), merged_array_rbc(:), merged_array_rload(:)

  integer(kint), allocatable :: array_I(:), perm(:)
  real(kdouble), allocatable :: array_R(:)

  integer(kint) :: ndomain_basis, n_graphs, funit
  integer(kint), allocatable :: domain_id(:)
  integer(kint) :: myrank, nprocs
  logical :: ex

  call monolis_mpi_initialize()
  myrank = monolis_mpi_get_global_my_rank()
  nprocs = monolis_mpi_get_global_comm_size()

  call setup_defaults(spec)
  call parse_args(spec, ndomain_basis)

  ! ---- sanity: "do_conn/node/bc/load/com" needs graph (merged_nodal_graph) ----
  if((spec%do_conn .or. spec%do_node .or. spec%do_bc .or. spec%do_load .or. spec%do_com) &
     .and. .not. spec%do_graph) then
    if(myrank == 0) then
      write(*,*) "Error: --conn/--node/--bc/--load/--com requires --graph (nodal graph)."
    endif
    call monolis_mpi_finalize()
    stop
  endif

  if(ndomain_basis < nprocs) then
    if(myrank == 0) write(*,*) "Error: ndomain_basis < nprocs"
    call monolis_mpi_finalize()
    stop
  endif

  ! ---- get n_graphs ----
  if(nprocs == 1)then
    fname = join3(spec%name%top_dir, spec%name%part_dir, "metagraph.dat")
    inquire(file=fname, exist=ex)
    if(.not.ex) stop "Error: metagraph.dat not found"
    open(newunit=funit, file=fname, status='old')
      read(funit,*) n_graphs
    close(funit)
  else
    fname = join3(spec%name%top_dir, spec%name%part_dir, "metagraph.dat.n_internal.<rank>")
    inquire(file=fname, exist=ex)
    if(.not.ex) stop "Error: metagraph.dat.n_internal not found"
    call monolis_input_internal_vertex_number(fname, n_graphs)
  endif

  if(myrank == 0) then
    write(*,*) "n_graphs=", n_graphs
    write(*,*) "do_graph/do_conn/do_node/do_bc/do_load/do_com=", &
      spec%do_graph, spec%do_conn, spec%do_node, spec%do_bc, spec%do_load, spec%do_com
  endif

  ! ---- init merged objects only if needed ----
  if(spec%do_graph) call gedatsu_graph_initialize(merged_nodal_graph)
  if(spec%do_conn ) call gedatsu_graph_initialize(merged_conn_graph)
  if(spec%do_com  ) call monolis_com_initialize_by_self(merged_monoCOM)

  ! ---- allocate arrays ----
  allocate(nodal_graphs(n_graphs))
  allocate(conn_graphs(n_graphs))
  allocate(monoCOMs(n_graphs))
  call monolis_alloc_I_1d(domain_id, n_graphs)

  if(spec%do_node) then
    allocate(n_dof_list_node(n_graphs))
    allocate(list_node(n_graphs))
    call monolis_list_initialize_I(n_dof_list_node, n_graphs)
    call monolis_list_initialize_R(list_node, n_graphs)
  endif

  if(spec%do_bc) then
    allocate(n_dof_list_bc(n_graphs))
    allocate(list_ibc(n_graphs))
    allocate(list_rbc(n_graphs))
    call monolis_list_initialize_I(n_dof_list_bc, n_graphs)
    call monolis_list_initialize_I(list_ibc, n_graphs)
    call monolis_list_initialize_R(list_rbc, n_graphs)
  endif

  if(spec%do_load) then
    allocate(n_dof_list_load(n_graphs))
    allocate(list_iload(n_graphs))
    allocate(list_rload(n_graphs))
    call monolis_list_initialize_I(n_dof_list_load, n_graphs)
    call monolis_list_initialize_I(list_iload, n_graphs)
    call monolis_list_initialize_R(list_rload, n_graphs)
  endif

  ! ---- domain_id ----
  if(nprocs == 1)then
    do i = 1, n_graphs
      domain_id(i) = i - 1
    enddo
  else
    write(cnum,'(i0)') myrank
    fname = join2(spec%name%top_dir, trim(spec%name%meta_id_prefix)//trim(cnum))
    inquire(file=fname, exist=ex)
    if(.not.ex) stop "Error: metagraph.dat.id.<rank> not found"
    open(newunit=funit, file=fname, status='old')
      read(funit,*)
      read(funit,*)
      do i = 1, n_graphs
        read(funit,*) domain_id(i)
      enddo
    close(funit)
  endif

  !> ---------------
  !> read each subgraph
  !> ---------------
  do i = 1, n_graphs
    call gedatsu_graph_initialize(nodal_graphs(i))
    call gedatsu_graph_initialize(conn_graphs(i))
    call monolis_com_initialize_by_self(monoCOMs(i))

    write(cnum,'(i0)') domain_id(i)

    ! ---- nodal graph ----
    if(spec%do_graph) then
      fname = join3(spec%name%top_dir, spec%name%part_dir, &
                    trim(spec%name%graph)//"."//trim(cnum))
      call monolis_input_graph(fname, &
        & nodal_graphs(i)%n_vertex, nodal_graphs(i)%vertex_id, &
        & nodal_graphs(i)%index, nodal_graphs(i)%item)
      call normalize_csr_item_fortran(nodal_graphs(i)%item, spec%nodal_csr%item_base)

      fname = join3(spec%name%top_dir, spec%name%part_dir, &
                    trim(spec%name%graph)//".n_internal."//trim(cnum))
      call monolis_input_internal_vertex_number(fname, nodal_graphs(i)%n_internal_vertex)

      call monolis_dealloc_I_1d(nodal_graphs(i)%vertex_id)
      fname = join3(spec%name%top_dir, spec%name%part_dir, &
                    trim(spec%name%graph)//".id."//trim(cnum))
      call monolis_input_global_id(fname, nodal_graphs(i)%n_vertex, nodal_graphs(i)%vertex_id)

      call monolis_alloc_I_1d(nodal_graphs(i)%vertex_domain_id, nodal_graphs(i)%n_vertex)
      nodal_graphs(i)%vertex_domain_id(:) = myrank
    endif

    ! ---- com ----
    if(spec%do_com) then
      call read_com_with_spec(monoCOMs(i), &
        join3(spec%name%top_dir, spec%name%part_dir, ""), &
        spec%name%graph, trim(cnum), spec%com)
    endif

    ! ---- connectivity graph ----
    if(spec%do_conn) then
      fname = join3(spec%name%top_dir, spec%name%part_dir, &
                    trim(spec%name%conn)//"."//trim(cnum))
      call monolis_input_graph(fname, &
        & conn_graphs(i)%n_vertex, conn_graphs(i)%vertex_id, &
        & conn_graphs(i)%index, conn_graphs(i)%item)
      call normalize_csr_item_fortran(conn_graphs(i)%item, spec%conn_csr%item_base)

      fname = join3(spec%name%top_dir, spec%name%part_dir, &
                    trim(spec%name%conn)//".n_internal."//trim(cnum))
      call monolis_input_internal_vertex_number(fname, conn_graphs(i)%n_internal_vertex)

      call monolis_dealloc_I_1d(conn_graphs(i)%vertex_id)
      fname = join3(spec%name%top_dir, spec%name%part_dir, &
                    trim(spec%name%conn)//".id."//trim(cnum))
      call monolis_input_global_id(fname, conn_graphs(i)%n_vertex, conn_graphs(i)%vertex_id)

      call monolis_alloc_I_1d(conn_graphs(i)%vertex_domain_id, conn_graphs(i)%n_vertex)
      conn_graphs(i)%vertex_domain_id(:) = myrank
    endif

    ! ---- node(distval) ----
    if(spec%do_node) then
      fname = join3(spec%name%top_dir, spec%name%part_dir, &
                    trim(spec%name%node)//"."//trim(cnum))
      call monolis_input_distval_R(fname, label, nnode, ndof, node)
      if(ndof /= 3) then
        write(*,*) "Error: ndof(node)=", ndof, " file=", trim(fname), " rank=", myrank
        stop "ndof isn't 3 (node)"
      endif

      call monolis_dealloc_I_1d(array_I)
      call monolis_alloc_I_1d(array_I, nnode)
      array_I(:) = ndof
      call monolis_list_set_I(n_dof_list_node, i, nnode, array_I)

      call monolis_dealloc_R_1d(array_R)
      call monolis_alloc_R_1d(array_R, nnode*ndof)
      do j = 1, nnode
        do k = 1, ndof
          idx = ndof*(j-1) + k
          array_R(idx) = node(k,j)
        enddo
      enddo
      call monolis_list_set_R(list_node, i, nnode*ndof, array_R)
    endif

    ! ---- bc ----
    if(spec%do_bc) then
      fname = join3(spec%name%top_dir, spec%name%part_dir, &
                    trim(spec%name%bc)//"."//trim(cnum))
      call monolis_input_bc_R(fname, nbc, ndof, ibc, rbc)
      if(ndof /= 3) then
        write(*,*) "Error: ndof(bc)=", ndof, " file=", trim(fname), " rank=", myrank
        stop "ndof isn't 3 (bc)"
      endif
      call normalize_bc_fortran(ibc, spec%bc%node_base, spec%bc%dir_base)
      call bc_to_list(nnode, nbc, ibc, rbc, i, n_dof_list_bc, list_ibc, list_rbc)
    endif

    ! ---- load ----
    if(spec%do_load) then
      fname = join3(spec%name%top_dir, spec%name%part_dir, &
                    trim(spec%name%load)//"."//trim(cnum))
      call monolis_input_bc_R(fname, nload, ndof, iload, rload)
      if(ndof /= 3) then
        write(*,*) "Error: ndof(load)=", ndof, " file=", trim(fname), " rank=", myrank
        stop "ndof isn't 3 (load)"
      endif
      call normalize_bc_fortran(iload, spec%bc%node_base, spec%bc%dir_base)
      call bc_to_list(nnode, nload, iload, rload, i, n_dof_list_load, list_iload, list_rload)
    endif
  enddo

  !> ---------------
  !> merge
  !> ---------------
  if(spec%do_graph) then
    call gedatsu_merge_nodal_subgraphs(n_graphs, nodal_graphs, monoCOMs, &
      & merged_nodal_graph, merged_monoCOM, ORDER_NODAL_ID)
  endif

  if(spec%do_conn) then
    call gedatsu_merge_connectivity_subgraphs(n_graphs, nodal_graphs, &
      & merged_nodal_graph, merged_monoCOM, n_graphs, conn_graphs, merged_conn_graph)
  endif

  if(spec%do_node) then
    call gedatsu_merge_distval_R(n_graphs, nodal_graphs, merged_nodal_graph, &
      & n_dof_list_node, list_node, merged_n_dof_list_node, merged_array_node)
  endif

  if(spec%do_bc) then
    call gedatsu_merge_distval_I(n_graphs, nodal_graphs, merged_nodal_graph, &
      & n_dof_list_bc, list_ibc, merged_n_dof_list_bc, merged_array_ibc)
    call gedatsu_merge_distval_R(n_graphs, nodal_graphs, merged_nodal_graph, &
      & n_dof_list_bc, list_rbc, merged_n_dof_list_bc, merged_array_rbc)
  endif

  if(spec%do_load) then
    call gedatsu_merge_distval_I(n_graphs, nodal_graphs, merged_nodal_graph, &
      & n_dof_list_load, list_iload, merged_n_dof_list_load, merged_array_iload)
    call gedatsu_merge_distval_R(n_graphs, nodal_graphs, merged_nodal_graph, &
      & n_dof_list_load, list_rload, merged_n_dof_list_load, merged_array_rload)
  endif

  !> ---------------
  !> output (only enabled targets)
  !> ---------------
  call output_enabled(spec, merged_nodal_graph, merged_conn_graph, merged_monoCOM, &
    & merged_n_dof_list_node, merged_array_node, &
    & merged_n_dof_list_bc, merged_array_ibc, merged_array_rbc, &
    & merged_n_dof_list_load, merged_array_iload, merged_array_rload)

  call monolis_mpi_finalize()
  if(myrank == 0) write(*,*) "DONE"
contains

  subroutine setup_defaults(spec)
    type(t_merge_spec), intent(inout) :: spec
    ! base デフォルトは C0（あなたの要件に合わせる）
    spec%nodal_csr%index_base = FMT_C0
    spec%nodal_csr%item_base  = FMT_C0
    spec%conn_csr%index_base  = FMT_C0
    spec%conn_csr%item_base   = FMT_C0
    spec%com%index_base       = FMT_C0
    spec%com%item_base        = FMT_C0
    spec%com%has_domain_suffix = .true.
    spec%bc%node_base         = FMT_C0
    spec%bc%dir_base          = FMT_C0

    spec%name%top_dir         = "./parted.0"
    spec%name%part_dir        = "parted.1"
    spec%name%graph           = "graph.dat"
    spec%name%conn            = "connectivity.dat"
    spec%name%node            = "node.dat"
    spec%name%bc              = "bc.dat"
    spec%name%load            = "load.dat"
    spec%name%out_dir         = "./parted.0"

    ! 重要: 何も指定されない限りマージしない（あなたの要望）
    spec%do_graph = .false.
    spec%do_conn  = .false.
    spec%do_node  = .false.
    spec%do_bc    = .false.
    spec%do_load  = .false.
    spec%do_com   = .false.
  end subroutine setup_defaults

  subroutine parse_args(spec, ndomain_basis)
    type(t_merge_spec), intent(inout) :: spec
    integer(kint), intent(out) :: ndomain_basis
    character(monolis_charlen) :: v, a1, a2, a3
    logical :: ok

    ! 必須: --ndomain
    call get_arg_value("--ndomain", v, ok)
    if(.not.ok) stop "Error: --ndomain is required"
    read(v,*) ndomain_basis

    call get_arg_value("--top", v, ok);   if(ok) spec%name%top_dir  = trim(v)
    call get_arg_value("--part", v, ok);  if(ok) spec%name%part_dir = trim(v)
    call get_arg_value("--out", v, ok);   if(ok) spec%name%out_dir  = trim(v)

    ! ---- 対象は「引数があるときだけ do_* = true」 ----
    call get_arg_value("--graph", v, ok)
    if(ok) then
      spec%name%graph = trim(v)
      spec%do_graph = .true.
      spec%do_com   = .true.  ! graph を読むなら com も通常は必要なのでON（不要なら --nocom を作る等）
    endif

    call get_arg_value("--conn", v, ok)
    if(ok) then
      spec%name%conn = trim(v)
      spec%do_conn = .true.
    endif

    call get_arg_value("--node", v, ok)
    if(ok) then
      spec%name%node = trim(v)
      spec%do_node = .true.
    endif

    call get_arg_value("--bc", v, ok)
    if(ok) then
      spec%name%bc = trim(v)
      spec%do_bc = .true.
    endif

    call get_arg_value("--load", v, ok)
    if(ok) then
      spec%name%load = trim(v)
      spec%do_load = .true.
    endif

    call get_arg_value("--nodal-csr", v, ok)
    if(ok) then
      call split2_csv(v, a1, a2)
      spec%nodal_csr%index_base = parse_base(a1)
      spec%nodal_csr%item_base  = parse_base(a2)
    endif

    call get_arg_value("--conn-csr", v, ok)
    if(ok) then
      call split2_csv(v, a1, a2)
      spec%conn_csr%index_base = parse_base(a1)
      spec%conn_csr%item_base  = parse_base(a2)
    endif

    call get_arg_value("--bcfmt", v, ok)
    if(ok) then
      call split2_csv(v, a1, a2)
      spec%bc%node_base = parse_base(a1)
      spec%bc%dir_base  = parse_base(a2)
    endif

    call get_arg_value("--com", v, ok)
    if(ok) then
      call split3_csv(v, a1, a2, a3)
      spec%com%index_base = parse_base(a1)
      spec%com%item_base  = parse_base(a2)
      if(len_trim(a3) == 0) then
        spec%com%has_domain_suffix = .true.
      else
        spec%com%has_domain_suffix = (trim(a3) /= "nosuffix")
      endif
      spec%do_com = .true.
    endif
  end subroutine parse_args

  subroutine normalize_csr_item_fortran(item, item_base)
    integer(kint), intent(inout) :: item(:)
    integer(kint), intent(in)    :: item_base
    if(item_base == FMT_C0) item = item + 1
  end subroutine normalize_csr_item_fortran

  subroutine normalize_bc_fortran(ibc, node_base, dir_base)
    integer(kint), intent(inout) :: ibc(:,:)
    integer(kint), intent(in) :: node_base, dir_base
    if(size(ibc,2) == 0) return
    if(node_base == FMT_C0) ibc(1,:) = ibc(1,:) + 1
    if(dir_base  == FMT_C0) ibc(2,:) = ibc(2,:) + 1
  end subroutine normalize_bc_fortran

  subroutine read_com_with_spec(com, dir, base, suf, fmt)
    type(monolis_COM), intent(inout) :: com
    character(*), intent(in) :: dir, base, suf
    type(t_com_format), intent(in) :: fmt
    character(monolis_charlen) :: fs, fr
    logical :: exs, exr

    if(fmt%has_domain_suffix) then
      fs = trim(dir)//trim(base)//".send."//trim(suf)
      fr = trim(dir)//trim(base)//".recv."//trim(suf)
    else
      fs = trim(dir)//trim(base)//".send"
      fr = trim(dir)//trim(base)//".recv"
    endif

    inquire(file=fs, exist=exs)
    inquire(file=fr, exist=exr)
    if(.not.(exs .and. exr)) then
      write(*,*) "Error: com file not found:", trim(fs), trim(fr)
      stop
    endif

    !------------------------------------------------------------
    ! ★ここはあなたの環境の関数名に合わせてください★
    !------------------------------------------------------------
    call monolis_input_send_com_table(fs, com)
    call monolis_input_recv_com_table(fr, com)

    if(fmt%item_base == FMT_C0) then
      if(associated(com%send_item)) com%send_item = com%send_item + 1
      if(associated(com%recv_item)) com%recv_item = com%recv_item + 1
    endif
  end subroutine read_com_with_spec

  subroutine bc_to_list(nnode, nbc, ibc, rbc, ig, n_dof_list, list_i, list_r)
    integer(kint), intent(in) :: nnode, nbc, ig
    integer(kint), intent(in) :: ibc(:,:)
    real(kdouble), intent(in) :: rbc(:)
    type(monolis_list_I), intent(inout) :: n_dof_list(:)
    type(monolis_list_I), intent(inout) :: list_i(:)
    type(monolis_list_R), intent(inout) :: list_r(:)

    integer(kint), allocatable :: cnt(:), tmp_i(:), p(:)
    real(kdouble), allocatable :: tmp_r(:)
    integer(kint) :: j, id

    call monolis_alloc_I_1d(cnt, nnode)
    cnt(:) = 0
    do j = 1, nbc
      cnt(ibc(1,j)) = cnt(ibc(1,j)) + 1
    enddo
    call monolis_list_set_I(n_dof_list, ig, nnode, cnt)

    if(nbc > 0) then
      call monolis_alloc_I_1d(tmp_i, nbc)
      call monolis_alloc_I_1d(p, nbc)
      tmp_i(:) = ibc(1,:)
      call monolis_get_sequence_array_I(p, nbc, 1, 1)
      call monolis_qsort_I_2d(tmp_i, p, 1, nbc)
      do j = 1, nbc
        id = p(j)
        tmp_i(j) = ibc(2,id)
      enddo
      call monolis_list_set_I(list_i, ig, nbc, tmp_i)

      call monolis_alloc_R_1d(tmp_r, nbc)
      do j = 1, nbc
        id = p(j)
        tmp_r(j) = rbc(id)
      enddo
      call monolis_list_set_R(list_r, ig, nbc, tmp_r)

      call monolis_dealloc_I_1d(tmp_i)
      call monolis_dealloc_I_1d(p)
      call monolis_dealloc_R_1d(tmp_r)
    else
      call monolis_alloc_I_1d(tmp_i, 0)
      call monolis_list_set_I(list_i, ig, 0, tmp_i)
      call monolis_dealloc_I_1d(tmp_i)

      call monolis_alloc_R_1d(tmp_r, 0)
      call monolis_list_set_R(list_r, ig, 0, tmp_r)
      call monolis_dealloc_R_1d(tmp_r)
    endif

    call monolis_dealloc_I_1d(cnt)
  end subroutine bc_to_list

  subroutine list_to_bc(nnode, n_dof_list, arr_i, arr_r, bc, br)
    integer(kint), intent(in) :: nnode
    integer(kint), intent(in) :: n_dof_list(:)
    integer(kint), intent(in) :: arr_i(:)
    real(kdouble), intent(in) :: arr_r(:)
    integer(kint), allocatable, intent(out) :: bc(:,:)
    real(kdouble), allocatable, intent(out) :: br(:)

    integer(kint) :: n, iS, iE, i, m

    n = size(arr_r)
    call monolis_alloc_I_2d(bc, 2, n)
    call monolis_alloc_R_1d(br, n)

    iS = 1
    do i = 1, nnode
      m = n_dof_list(i)
      if(m /= 0) then
        iE = iS + m - 1
        bc(1, iS:iE) = i
        bc(2, iS:iE) = arr_i(iS:iE)
        br(iS:iE)    = arr_r(iS:iE)
        iS = iE + 1
      endif
    enddo
  end subroutine list_to_bc

  subroutine output_enabled(spec, gnod, gcon, com, &
      nlist_node, arr_node, &
      nlist_bc, arr_ibc, arr_rbc, &
      nlist_load, arr_iload, arr_rload)
    type(t_merge_spec), intent(in) :: spec
    type(gedatsu_graph), intent(in) :: gnod, gcon
    type(monolis_COM), intent(in) :: com
    integer(kint), allocatable, intent(in) :: nlist_node(:), nlist_bc(:), nlist_load(:)
    real(kdouble),    allocatable, intent(in) :: arr_node(:), arr_rbc(:), arr_rload(:)
    integer(kint),    allocatable, intent(in) :: arr_ibc(:), arr_iload(:)

    character(monolis_charlen) :: outbase, fname, lab
    integer(kint), allocatable :: item_tmp(:), index_tmp(:)
    integer(kint), allocatable :: ibc(:,:), iload(:,:)
    real(kdouble), allocatable :: node(:,:), rbc(:), rload(:)
    integer(kint) :: nnode, ndof, nbc, nload, i, j, id

    outbase = trim(spec%name%out_dir)//"/"

    ! --- graph ---
    if(spec%do_graph) then
      call monolis_alloc_I_1d(index_tmp, size(gnod%index))
      call monolis_alloc_I_1d(item_tmp,  size(gnod%item))
      index_tmp = gnod%index
      item_tmp  = gnod%item
      if(spec%nodal_csr%item_base == FMT_C0) item_tmp = item_tmp - 1

      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%graph))
      call monolis_output_graph(fname, gnod%n_vertex, gnod%vertex_id, index_tmp, item_tmp)

      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%graph)//".n_internal")
      call monolis_output_internal_vertex_number(fname, gnod%n_internal_vertex)

      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%graph)//".id")
      call monolis_output_global_id(fname, gnod%n_vertex, gnod%vertex_id)

      call monolis_dealloc_I_1d(index_tmp)
      call monolis_dealloc_I_1d(item_tmp)
    endif

    ! --- com ---
    if(spec%do_com) then
      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%graph)//".send")
      call monolis_output_send_com_table(fname, com)
      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%graph)//".recv")
      call monolis_output_recv_com_table(fname, com)
    endif

    ! --- connectivity ---
    if(spec%do_conn) then
      call monolis_alloc_I_1d(index_tmp, size(gcon%index))
      call monolis_alloc_I_1d(item_tmp,  size(gcon%item))
      index_tmp = gcon%index
      item_tmp  = gcon%item
      if(spec%conn_csr%item_base == FMT_C0) item_tmp = item_tmp - 1

      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%conn))
      call monolis_output_graph(fname, gcon%n_vertex, gcon%vertex_id, index_tmp, item_tmp)

      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%conn)//".n_internal")
      call monolis_output_internal_vertex_number(fname, gcon%n_internal_vertex)

      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%conn)//".id")
      call monolis_output_global_id(fname, gcon%n_vertex, gcon%vertex_id)

      call monolis_dealloc_I_1d(index_tmp)
      call monolis_dealloc_I_1d(item_tmp)
    endif

    ! --- node(distval) ---
    if(spec%do_node) then
      nnode = gnod%n_vertex
      ndof  = 3
      call monolis_alloc_R_2d(node, ndof, nnode)
      do i = 1, nnode
        do j = 1, ndof
          id = ndof*(i-1) + j
          node(j,i) = arr_node(id)
        enddo
      enddo

      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%node))
      lab = "#node"
      call monolis_output_distval_R(fname, lab, nnode, ndof, node)
      call monolis_dealloc_R_2d(node)
    endif

    ! --- bc ---
    if(spec%do_bc) then
      nnode = gnod%n_vertex
      call list_to_bc(nnode, nlist_bc, arr_ibc, arr_rbc, ibc, rbc)
      nbc = size(rbc)

      if(nbc > 0) then
        if(spec%bc%node_base == FMT_C0) ibc(1,:) = ibc(1,:) - 1
        if(spec%bc%dir_base  == FMT_C0) ibc(2,:) = ibc(2,:) - 1
      endif

      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%bc))
      call monolis_output_bc_R(fname, nbc, 3, ibc, rbc)
      call monolis_dealloc_I_2d(ibc)
      call monolis_dealloc_R_1d(rbc)
    endif

    ! --- load ---
    if(spec%do_load) then
      nnode = gnod%n_vertex
      call list_to_bc(nnode, nlist_load, arr_iload, arr_rload, iload, rload)
      nload = size(rload)

      if(nload > 0) then
        if(spec%bc%node_base == FMT_C0) iload(1,:) = iload(1,:) - 1
        if(spec%bc%dir_base  == FMT_C0) iload(2,:) = iload(2,:) - 1
      endif

      fname = monolis_get_global_output_file_name(outbase, "", trim(spec%name%load))
      call monolis_output_bc_R(fname, nload, 3, iload, rload)
      call monolis_dealloc_I_2d(iload)
      call monolis_dealloc_R_1d(rload)
    endif
  end subroutine output_enabled

end program merge_test_param_optional
