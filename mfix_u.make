.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
mfix.exe : \
    AMBM.mod \
    BC.mod \
    BOUNDFUNIJK3.mod \
    BOUNDFUNIJK.mod \
    CDIST.mod \
    CHECK.mod \
    CHISCHEME.mod \
    COEFF.mod \
    CONSTANT.mod \
    CONT.mod \
    CORNER.mod \
    DRAG.mod \
    ENERGY.mod \
    FLDVAR.mod \
    FUNCTION.mod \
    FUNITS.mod \
    GEOMETRY.mod \
    IC.mod \
    INDICES.mod \
    IS.mod \
    KINTHEORY.mod \
    LEQSOL.mod \
    MACHINE.mod \
    MATRIX.mod \
    MFIX_NETCDF.mod \
    MFLUX.mod \
    OUTPUT.mod \
    PARALLEL.mod \
    PARAM1.mod \
    PARAM.mod \
    PARSE.mod \
    PGCOR.mod \
    PHYSPROP.mod \
    PSCOR.mod \
    RESIDUAL.mod \
    RUN.mod \
    RXNS.mod \
    SCALARS.mod \
    SCALES.mod \
    TAU_G.mod \
    TAU_S.mod \
    TIME_CPU.mod \
    TMP_ARRAY1.mod \
    TMP_ARRAY.mod \
    TOLERANC.mod \
    TRACE.mod \
    TURB.mod \
    UR_FACS.mod \
    USR.mod \
    VISC_G.mod \
    VISC_S.mod \
    VSHEAR.mod \
    XSI_ARRAY.mod \
    CUTCELL.mod \
    DASHBOARD.mod \
    POLYGON.mod \
    PROGRESS_BAR.mod \
    QUADRIC.mod \
    STL.mod \
    VTK.mod \
    MCHEM.mod \
    DES_BC.mod \
    DESGRID.mod \
    DES_IC.mod \
    DESMPI.mod \
    DESMPI_WRAPPER.mod \
    DES_RXNS.mod \
    DES_THERMO.mod \
    DISCRETELEMENT.mod \
    INTERPOLATION.mod \
    MFIX_PIC.mod \
    MPPIC_WALLBC.mod \
    RANDOMNO.mod \
    SENDRECVNODE.mod \
    COMPAR.mod \
    DBG_UTIL.mod \
    DEBUG.mod \
    GRIDMAP.mod \
    MPI.mod \
    MPI_UTILITY.mod \
    PARALLEL_MPI.mod \
    SENDRECV3.mod \
    SENDRECV.mod \
    GHDTHEORY.mod \
    QMOMK_BC.mod \
    QMOMK_COLLISION.mod \
    QMOMK_FLUXES.mod \
    QMOM_KINETIC_EQUATION.mod \
    QMOMK_PARAMETERS.mod \
    QMOMK_QUADRATURE.mod \
    QMOMK_TOOLS.mod \
    accum_resid.$(OBJ_EXT) \
    adjust_a_u_g.$(OBJ_EXT) \
    adjust_a_u_s.$(OBJ_EXT) \
    adjust_a_v_g.$(OBJ_EXT) \
    adjust_a_v_s.$(OBJ_EXT) \
    adjust_a_w_g.$(OBJ_EXT) \
    adjust_a_w_s.$(OBJ_EXT) \
    adjust_dt.$(OBJ_EXT) \
    adjust_eps.$(OBJ_EXT) \
    adjust_leq.$(OBJ_EXT) \
    adjust_rop.$(OBJ_EXT) \
    adjust_theta.$(OBJ_EXT) \
    allocate_arrays.$(OBJ_EXT) \
    bc_phi.$(OBJ_EXT) \
    bc_theta.$(OBJ_EXT) \
    b_m_p_star.$(OBJ_EXT) \
    bound_x.$(OBJ_EXT) \
    calc_cell.$(OBJ_EXT) \
    calc_coeff.$(OBJ_EXT) \
    calc_d.$(OBJ_EXT) \
    calc_dif_g.$(OBJ_EXT) \
    calc_dif_s.$(OBJ_EXT) \
    calc_drag.$(OBJ_EXT) \
    calc_e.$(OBJ_EXT) \
    calc_gama.$(OBJ_EXT) \
    calc_grbdry.$(OBJ_EXT) \
    calc_h.$(OBJ_EXT) \
    calc_k_cp.$(OBJ_EXT) \
    calc_k_g.$(OBJ_EXT) \
    calc_k_s.$(OBJ_EXT) \
    calc_mflux.$(OBJ_EXT) \
    calc_mu_g.$(OBJ_EXT) \
    calc_mu_s.$(OBJ_EXT) \
    calc_mw.$(OBJ_EXT) \
    calc_outflow.$(OBJ_EXT) \
    calc_p_star.$(OBJ_EXT) \
    calc_resid.$(OBJ_EXT) \
    calc_s_ddot_s.$(OBJ_EXT) \
    calc_trd_g.$(OBJ_EXT) \
    calc_trd_s.$(OBJ_EXT) \
    calc_u_friction.$(OBJ_EXT) \
    calc_vol_fr.$(OBJ_EXT) \
    calc_xsi.$(OBJ_EXT) \
    cal_d.$(OBJ_EXT) \
    check_ab_m.$(OBJ_EXT) \
    check_convergence.$(OBJ_EXT) \
    check_data_01.$(OBJ_EXT) \
    check_data_02.$(OBJ_EXT) \
    check_data_03.$(OBJ_EXT) \
    check_data_04.$(OBJ_EXT) \
    check_data_05.$(OBJ_EXT) \
    check_data_06.$(OBJ_EXT) \
    check_data_07.$(OBJ_EXT) \
    check_data_08.$(OBJ_EXT) \
    check_data_09.$(OBJ_EXT) \
    check_data_20.$(OBJ_EXT) \
    check_data_30.$(OBJ_EXT) \
    check_mass_balance.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    check_plane.$(OBJ_EXT) \
    cn_extrapol.$(OBJ_EXT) \
    compare.$(OBJ_EXT) \
    conv_dif_phi.$(OBJ_EXT) \
    conv_dif_u_g.$(OBJ_EXT) \
    conv_dif_u_s.$(OBJ_EXT) \
    conv_dif_v_g.$(OBJ_EXT) \
    conv_dif_v_s.$(OBJ_EXT) \
    conv_dif_w_g.$(OBJ_EXT) \
    conv_dif_w_s.$(OBJ_EXT) \
    conv_pp_g.$(OBJ_EXT) \
    conv_rop.$(OBJ_EXT) \
    conv_rop_g.$(OBJ_EXT) \
    conv_rop_s.$(OBJ_EXT) \
    conv_source_epp.$(OBJ_EXT) \
    copy_a.$(OBJ_EXT) \
    corner.$(OBJ_EXT) \
    correct_0.$(OBJ_EXT) \
    correct_1.$(OBJ_EXT) \
    dgtsl.$(OBJ_EXT) \
    dif_u_is.$(OBJ_EXT) \
    dif_v_is.$(OBJ_EXT) \
    dif_w_is.$(OBJ_EXT) \
    discretize.$(OBJ_EXT) \
    display_resid.$(OBJ_EXT) \
    drag_gs.$(OBJ_EXT) \
    drag_ss.$(OBJ_EXT) \
    eosg.$(OBJ_EXT) \
    equal.$(OBJ_EXT) \
    error_routine.$(OBJ_EXT) \
    exchange.$(OBJ_EXT) \
    exit.$(OBJ_EXT) \
    flow_to_vel.$(OBJ_EXT) \
    g_0.$(OBJ_EXT) \
    get_bc_area.$(OBJ_EXT) \
    get_data.$(OBJ_EXT) \
    get_eq.$(OBJ_EXT) \
    get_flow_bc.$(OBJ_EXT) \
    get_hloss.$(OBJ_EXT) \
    get_is.$(OBJ_EXT) \
    get_philoss.$(OBJ_EXT) \
    get_smass.$(OBJ_EXT) \
    get_stats.$(OBJ_EXT) \
    get_walls_bc.$(OBJ_EXT) \
    in_bin_512.$(OBJ_EXT) \
    in_bin_512i.$(OBJ_EXT) \
    init_ab_m.$(OBJ_EXT) \
    init_fvars.$(OBJ_EXT) \
    init_namelist.$(OBJ_EXT) \
    init_resid.$(OBJ_EXT) \
    iterate.$(OBJ_EXT) \
    k_epsilon_prop.$(OBJ_EXT) \
    kintheory_drag_ss.$(OBJ_EXT) \
    kintheory_energy_dissipation_ss.$(OBJ_EXT) \
    kintheory_u_s.$(OBJ_EXT) \
    kintheory_v_s.$(OBJ_EXT) \
    kintheory_w_s.$(OBJ_EXT) \
    leq_bicgs.$(OBJ_EXT) \
    leq_bicgst.$(OBJ_EXT) \
    leq_cg.$(OBJ_EXT) \
    leq_gmres.$(OBJ_EXT) \
    leq_sor.$(OBJ_EXT) \
    line_too_big.$(OBJ_EXT) \
    location_check.$(OBJ_EXT) \
    location.$(OBJ_EXT) \
    machine.$(OBJ_EXT) \
    make_upper_case.$(OBJ_EXT) \
    mark_phase_4_cor.$(OBJ_EXT) \
    mfix.$(OBJ_EXT) \
    mod_bc_i.$(OBJ_EXT) \
    mod_bc_j.$(OBJ_EXT) \
    mod_bc_k.$(OBJ_EXT) \
    open_file.$(OBJ_EXT) \
    open_files.$(OBJ_EXT) \
    out_array_c.$(OBJ_EXT) \
    out_array.$(OBJ_EXT) \
    out_array_kc.$(OBJ_EXT) \
    out_array_k.$(OBJ_EXT) \
    out_bin_512.$(OBJ_EXT) \
    out_bin_512i.$(OBJ_EXT) \
    out_bin_512r.$(OBJ_EXT) \
    out_bin_r.$(OBJ_EXT) \
    parse_line.$(OBJ_EXT) \
    parse_resid_string.$(OBJ_EXT) \
    parse_rxn.$(OBJ_EXT) \
    partial_elim.$(OBJ_EXT) \
    physical_prop.$(OBJ_EXT) \
    read_database.$(OBJ_EXT) \
    read_namelist.$(OBJ_EXT) \
    read_res0.$(OBJ_EXT) \
    read_res1.$(OBJ_EXT) \
    remove_comment.$(OBJ_EXT) \
    reset_new.$(OBJ_EXT) \
    rrates0.$(OBJ_EXT) \
    rrates.$(OBJ_EXT) \
    rrates_init.$(OBJ_EXT) \
    scalar_prop.$(OBJ_EXT) \
    seek_comment.$(OBJ_EXT) \
    seek_end.$(OBJ_EXT) \
    set_bc0.$(OBJ_EXT) \
    set_bc1.$(OBJ_EXT) \
    set_constants.$(OBJ_EXT) \
    set_constprop.$(OBJ_EXT) \
    set_flags.$(OBJ_EXT) \
    set_fluidbed_p.$(OBJ_EXT) \
    set_geometry1.$(OBJ_EXT) \
    set_geometry.$(OBJ_EXT) \
    set_ic.$(OBJ_EXT) \
    set_increments3.$(OBJ_EXT) \
    set_increments.$(OBJ_EXT) \
    set_index1a3.$(OBJ_EXT) \
    set_index1a.$(OBJ_EXT) \
    set_index1.$(OBJ_EXT) \
    set_l_scale.$(OBJ_EXT) \
    set_max2.$(OBJ_EXT) \
    set_mw_mix_g.$(OBJ_EXT) \
    set_outflow.$(OBJ_EXT) \
    set_ro_g.$(OBJ_EXT) \
    set_wall_bc.$(OBJ_EXT) \
    shift_dxyz.$(OBJ_EXT) \
    solve_continuity.$(OBJ_EXT) \
    solve_energy_eq.$(OBJ_EXT) \
    solve_epp.$(OBJ_EXT) \
    solve_granular_energy.$(OBJ_EXT) \
    solve_k_epsilon_eq.$(OBJ_EXT) \
    solve_lin_eq.$(OBJ_EXT) \
    solve_pp_g.$(OBJ_EXT) \
    solve_scalar_eq.$(OBJ_EXT) \
    solve_species_eq.$(OBJ_EXT) \
    solve_vel_star.$(OBJ_EXT) \
    source_granular_energy.$(OBJ_EXT) \
    source_phi.$(OBJ_EXT) \
    source_pp_g.$(OBJ_EXT) \
    source_rop_g.$(OBJ_EXT) \
    source_rop_s.$(OBJ_EXT) \
    source_u_g.$(OBJ_EXT) \
    source_u_s.$(OBJ_EXT) \
    source_v_g.$(OBJ_EXT) \
    source_v_s.$(OBJ_EXT) \
    source_w_g.$(OBJ_EXT) \
    source_w_s.$(OBJ_EXT) \
    tau_u_g.$(OBJ_EXT) \
    tau_u_s.$(OBJ_EXT) \
    tau_v_g.$(OBJ_EXT) \
    tau_v_s.$(OBJ_EXT) \
    tau_w_g.$(OBJ_EXT) \
    tau_w_s.$(OBJ_EXT) \
    test_lin_eq.$(OBJ_EXT) \
    time_march.$(OBJ_EXT) \
    transfer.$(OBJ_EXT) \
    transport_prop.$(OBJ_EXT) \
    undef_2_0.$(OBJ_EXT) \
    under_relax.$(OBJ_EXT) \
    update_old.$(OBJ_EXT) \
    usr0.$(OBJ_EXT) \
    usr1.$(OBJ_EXT) \
    usr2.$(OBJ_EXT) \
    usr3.$(OBJ_EXT) \
    usr_init_namelist.$(OBJ_EXT) \
    usr_write_out0.$(OBJ_EXT) \
    usr_write_out1.$(OBJ_EXT) \
    utilities.$(OBJ_EXT) \
    vavg_u_g.$(OBJ_EXT) \
    vavg_u_s.$(OBJ_EXT) \
    vavg_v_g.$(OBJ_EXT) \
    vavg_v_s.$(OBJ_EXT) \
    vavg_w_g.$(OBJ_EXT) \
    vavg_w_s.$(OBJ_EXT) \
    vf_gs_x.$(OBJ_EXT) \
    vf_gs_y.$(OBJ_EXT) \
    vf_gs_z.$(OBJ_EXT) \
    vtc_scalar.$(OBJ_EXT) \
    write_ab_m.$(OBJ_EXT) \
    write_ab_m_var.$(OBJ_EXT) \
    write_error.$(OBJ_EXT) \
    write_header.$(OBJ_EXT) \
    write_out0.$(OBJ_EXT) \
    write_out1.$(OBJ_EXT) \
    write_out3.$(OBJ_EXT) \
    write_res0.$(OBJ_EXT) \
    write_res1.$(OBJ_EXT) \
    write_spx0.$(OBJ_EXT) \
    write_spx1.$(OBJ_EXT) \
    write_table.$(OBJ_EXT) \
    write_usr0.$(OBJ_EXT) \
    write_usr1.$(OBJ_EXT) \
    xerbla.$(OBJ_EXT) \
    zero_array.$(OBJ_EXT) \
    zero_norm_vel.$(OBJ_EXT) \
    allocate_cut_cell_arrays.$(OBJ_EXT) \
    allocate_dummy_cut_cell_arrays.$(OBJ_EXT) \
    calc_vort_out.$(OBJ_EXT) \
    cartesian_grid_init_namelist.$(OBJ_EXT) \
    CG_set_bc0.$(OBJ_EXT) \
    CG_set_outflow.$(OBJ_EXT) \
    CG_source_u_g.$(OBJ_EXT) \
    CG_source_u_s.$(OBJ_EXT) \
    CG_source_v_g.$(OBJ_EXT) \
    CG_source_v_s.$(OBJ_EXT) \
    CG_source_w_g.$(OBJ_EXT) \
    CG_source_w_s.$(OBJ_EXT) \
    check_data_cartesian.$(OBJ_EXT) \
    cut_cell_preprocessing.$(OBJ_EXT) \
    deallocate_cut_cell_arrays.$(OBJ_EXT) \
    define_quadrics.$(OBJ_EXT) \
    dmp_cartesian.$(OBJ_EXT) \
    eval_usr_fct.$(OBJ_EXT) \
    get_alpha.$(OBJ_EXT) \
    get_connectivity.$(OBJ_EXT) \
    get_cut_cell_flags.$(OBJ_EXT) \
    get_cut_cell_volume_area.$(OBJ_EXT) \
    get_delh.$(OBJ_EXT) \
    get_master.$(OBJ_EXT) \
    get_poly_data.$(OBJ_EXT) \
    get_stl_data.$(OBJ_EXT) \
    set_Odxyz.$(OBJ_EXT) \
    update_dashboard.$(OBJ_EXT) \
    vtk_out.$(OBJ_EXT) \
    write_progress_bar.$(OBJ_EXT) \
    calc_jacobian.$(OBJ_EXT) \
    check_data_chem.$(OBJ_EXT) \
    dgpadm.$(OBJ_EXT) \
    exponential.$(OBJ_EXT) \
    fex.$(OBJ_EXT) \
    g_derivs.$(OBJ_EXT) \
    jac.$(OBJ_EXT) \
    mchem_init.$(OBJ_EXT) \
    mchem_odepack_init.$(OBJ_EXT) \
    mchem_time_march.$(OBJ_EXT) \
    misat_table_init.$(OBJ_EXT) \
    react.$(OBJ_EXT) \
    usrfg.$(OBJ_EXT) \
    add_part_to_link_list.$(OBJ_EXT) \
    calc_app_coh_force.$(OBJ_EXT) \
    calc_cap_coh_force.$(OBJ_EXT) \
    calc_cohesive_forces.$(OBJ_EXT) \
    calc_esc_coh_force.$(OBJ_EXT) \
    calc_square_well.$(OBJ_EXT) \
    calc_van_der_waals.$(OBJ_EXT) \
    check_link.$(OBJ_EXT) \
    check_sw_wall_interaction.$(OBJ_EXT) \
    check_vdw_wall_interaction.$(OBJ_EXT) \
    initialize_cohesion_parameters.$(OBJ_EXT) \
    initialize_coh_int_search.$(OBJ_EXT) \
    linked_interaction_eval.$(OBJ_EXT) \
    remove_part_from_link_list.$(OBJ_EXT) \
    unlinked_interaction_eval.$(OBJ_EXT) \
    update_search_grids.$(OBJ_EXT) \
    calc_attrition_des.$(OBJ_EXT) \
    calc_force_des.$(OBJ_EXT) \
    calc_rrate_des.$(OBJ_EXT) \
    calc_thermo_des.$(OBJ_EXT) \
    cell_near_wall.$(OBJ_EXT) \
    cfassign.$(OBJ_EXT) \
    cffctowall.$(OBJ_EXT) \
    cffctow.$(OBJ_EXT) \
    cfnewvalues.$(OBJ_EXT) \
    cfrelvel.$(OBJ_EXT) \
    cfslide.$(OBJ_EXT) \
    cfslidewall.$(OBJ_EXT) \
    cfupdateold.$(OBJ_EXT) \
    cfwallcontact.$(OBJ_EXT) \
    cfwallposvel.$(OBJ_EXT) \
    check_des_bc.$(OBJ_EXT) \
    check_des_data.$(OBJ_EXT) \
    check_des_ic.$(OBJ_EXT) \
    check_des_rxns.$(OBJ_EXT) \
    check_des_thermo.$(OBJ_EXT) \
    des_allocate_arrays.$(OBJ_EXT) \
    des_check_particle.$(OBJ_EXT) \
    des_functions.$(OBJ_EXT) \
    des_granular_temperature.$(OBJ_EXT) \
    des_init_arrays.$(OBJ_EXT) \
    des_init_bc.$(OBJ_EXT) \
    des_init_namelist.$(OBJ_EXT) \
    des_mass_inlet.$(OBJ_EXT) \
    des_physical_prop.$(OBJ_EXT) \
    des_reaction_model.$(OBJ_EXT) \
    des_rrates.$(OBJ_EXT) \
    des_set_ic.$(OBJ_EXT) \
    des_thermo_cond.$(OBJ_EXT) \
    des_thermo_conv.$(OBJ_EXT) \
    des_thermo_newvalues.$(OBJ_EXT) \
    des_thermo_rad.$(OBJ_EXT) \
    des_time_march.$(OBJ_EXT) \
    des_wallbc_preprocessing.$(OBJ_EXT) \
    drag_fgs.$(OBJ_EXT) \
    gas_drag.$(OBJ_EXT) \
    generate_particle_config.$(OBJ_EXT) \
    grid_based_neighbor_search.$(OBJ_EXT) \
    make_arrays_des.$(OBJ_EXT) \
    mppic_routines.$(OBJ_EXT) \
    neighbour.$(OBJ_EXT) \
    nsquare.$(OBJ_EXT) \
    octree.$(OBJ_EXT) \
    particles_in_cell.$(OBJ_EXT) \
    quadtree.$(OBJ_EXT) \
    read_des_restart.$(OBJ_EXT) \
    thermo_nbr.$(OBJ_EXT) \
    walledgecontact.$(OBJ_EXT) \
    wallfacecontact.$(OBJ_EXT) \
    wallnodecontact.$(OBJ_EXT) \
    write_des_data.$(OBJ_EXT) \
    write_des_restart.$(OBJ_EXT) \
    gaussj.$(OBJ_EXT) \
    odeint.$(OBJ_EXT) \
    rkck.$(OBJ_EXT) \
    rkqs.$(OBJ_EXT) \
    source_population_eq.$(OBJ_EXT) \
    usr_dqmom.$(OBJ_EXT) \
    adjust_eps_ghd.$(OBJ_EXT) \
    bulk_viscosity.$(OBJ_EXT) \
    calc_d_ghd.$(OBJ_EXT) \
    calc_external_forces.$(OBJ_EXT) \
    calc_nflux.$(OBJ_EXT) \
    chi_ij_GHD.$(OBJ_EXT) \
    cooling_rate.$(OBJ_EXT) \
    cooling_rate_tc.$(OBJ_EXT) \
    dufour_coeff.$(OBJ_EXT) \
    ghd.$(OBJ_EXT) \
    ghdmassflux.$(OBJ_EXT) \
    mass_mobility.$(OBJ_EXT) \
    ordinary_diff.$(OBJ_EXT) \
    partial_elim_ghd.$(OBJ_EXT) \
    pressure.$(OBJ_EXT) \
    shear_viscosity.$(OBJ_EXT) \
    source_ghd_granular_energy.$(OBJ_EXT) \
    thermal_conductivity.$(OBJ_EXT) \
    thermal_diffusivity.$(OBJ_EXT) \
    thermal_mobility.$(OBJ_EXT) \
    transport_coeff_ghd.$(OBJ_EXT) \
    qmomk_allocate_arrays.$(OBJ_EXT) \
    qmomk_gas_drag.$(OBJ_EXT) \
    qmomk_init_bc.$(OBJ_EXT) \
    qmomk_initial_conditions.$(OBJ_EXT) \
    qmomk_init_namelist.$(OBJ_EXT) \
    qmomk_make_arrays.$(OBJ_EXT) \
    qmomk_read_restart.$(OBJ_EXT) \
    qmomk_set_bc.$(OBJ_EXT) \
    qmomk_time_march.$(OBJ_EXT) \
    qmomk_write_restart.$(OBJ_EXT) \
    get_values.$(OBJ_EXT) \
    readTherm.$(OBJ_EXT) \
    blas90.a odepack.a dgtsv90.a
	$(LINK_CMD) $(LINK_FLAGS) \
    accum_resid.$(OBJ_EXT) \
    adjust_a_u_g.$(OBJ_EXT) \
    adjust_a_u_s.$(OBJ_EXT) \
    adjust_a_v_g.$(OBJ_EXT) \
    adjust_a_v_s.$(OBJ_EXT) \
    adjust_a_w_g.$(OBJ_EXT) \
    adjust_a_w_s.$(OBJ_EXT) \
    adjust_dt.$(OBJ_EXT) \
    adjust_eps.$(OBJ_EXT) \
    adjust_leq.$(OBJ_EXT) \
    adjust_rop.$(OBJ_EXT) \
    adjust_theta.$(OBJ_EXT) \
    allocate_arrays.$(OBJ_EXT) \
    ambm_mod.$(OBJ_EXT) \
    bc_mod.$(OBJ_EXT) \
    bc_phi.$(OBJ_EXT) \
    bc_theta.$(OBJ_EXT) \
    b_m_p_star.$(OBJ_EXT) \
    boundfunijk3_mod.$(OBJ_EXT) \
    boundfunijk_mod.$(OBJ_EXT) \
    bound_x.$(OBJ_EXT) \
    calc_cell.$(OBJ_EXT) \
    calc_coeff.$(OBJ_EXT) \
    calc_d.$(OBJ_EXT) \
    calc_dif_g.$(OBJ_EXT) \
    calc_dif_s.$(OBJ_EXT) \
    calc_drag.$(OBJ_EXT) \
    calc_e.$(OBJ_EXT) \
    calc_gama.$(OBJ_EXT) \
    calc_grbdry.$(OBJ_EXT) \
    calc_h.$(OBJ_EXT) \
    calc_k_cp.$(OBJ_EXT) \
    calc_k_g.$(OBJ_EXT) \
    calc_k_s.$(OBJ_EXT) \
    calc_mflux.$(OBJ_EXT) \
    calc_mu_g.$(OBJ_EXT) \
    calc_mu_s.$(OBJ_EXT) \
    calc_mw.$(OBJ_EXT) \
    calc_outflow.$(OBJ_EXT) \
    calc_p_star.$(OBJ_EXT) \
    calc_resid.$(OBJ_EXT) \
    calc_s_ddot_s.$(OBJ_EXT) \
    calc_trd_g.$(OBJ_EXT) \
    calc_trd_s.$(OBJ_EXT) \
    calc_u_friction.$(OBJ_EXT) \
    calc_vol_fr.$(OBJ_EXT) \
    calc_xsi.$(OBJ_EXT) \
    cal_d.$(OBJ_EXT) \
    cdist_mod.$(OBJ_EXT) \
    check_ab_m.$(OBJ_EXT) \
    check_convergence.$(OBJ_EXT) \
    check_data_01.$(OBJ_EXT) \
    check_data_02.$(OBJ_EXT) \
    check_data_03.$(OBJ_EXT) \
    check_data_04.$(OBJ_EXT) \
    check_data_05.$(OBJ_EXT) \
    check_data_06.$(OBJ_EXT) \
    check_data_07.$(OBJ_EXT) \
    check_data_08.$(OBJ_EXT) \
    check_data_09.$(OBJ_EXT) \
    check_data_20.$(OBJ_EXT) \
    check_data_30.$(OBJ_EXT) \
    check_mass_balance.$(OBJ_EXT) \
    check_mod.$(OBJ_EXT) \
    check_one_axis.$(OBJ_EXT) \
    check_plane.$(OBJ_EXT) \
    chischeme_mod.$(OBJ_EXT) \
    cn_extrapol.$(OBJ_EXT) \
    coeff_mod.$(OBJ_EXT) \
    compare.$(OBJ_EXT) \
    constant_mod.$(OBJ_EXT) \
    cont_mod.$(OBJ_EXT) \
    conv_dif_phi.$(OBJ_EXT) \
    conv_dif_u_g.$(OBJ_EXT) \
    conv_dif_u_s.$(OBJ_EXT) \
    conv_dif_v_g.$(OBJ_EXT) \
    conv_dif_v_s.$(OBJ_EXT) \
    conv_dif_w_g.$(OBJ_EXT) \
    conv_dif_w_s.$(OBJ_EXT) \
    conv_pp_g.$(OBJ_EXT) \
    conv_rop.$(OBJ_EXT) \
    conv_rop_g.$(OBJ_EXT) \
    conv_rop_s.$(OBJ_EXT) \
    conv_source_epp.$(OBJ_EXT) \
    copy_a.$(OBJ_EXT) \
    corner.$(OBJ_EXT) \
    corner_mod.$(OBJ_EXT) \
    correct_0.$(OBJ_EXT) \
    correct_1.$(OBJ_EXT) \
    dgtsl.$(OBJ_EXT) \
    dif_u_is.$(OBJ_EXT) \
    dif_v_is.$(OBJ_EXT) \
    dif_w_is.$(OBJ_EXT) \
    discretize.$(OBJ_EXT) \
    display_resid.$(OBJ_EXT) \
    drag_gs.$(OBJ_EXT) \
    drag_mod.$(OBJ_EXT) \
    drag_ss.$(OBJ_EXT) \
    energy_mod.$(OBJ_EXT) \
    eosg.$(OBJ_EXT) \
    equal.$(OBJ_EXT) \
    error_routine.$(OBJ_EXT) \
    exchange.$(OBJ_EXT) \
    exit.$(OBJ_EXT) \
    fldvar_mod.$(OBJ_EXT) \
    flow_to_vel.$(OBJ_EXT) \
    function_mod.$(OBJ_EXT) \
    funits_mod.$(OBJ_EXT) \
    g_0.$(OBJ_EXT) \
    geometry_mod.$(OBJ_EXT) \
    get_bc_area.$(OBJ_EXT) \
    get_data.$(OBJ_EXT) \
    get_eq.$(OBJ_EXT) \
    get_flow_bc.$(OBJ_EXT) \
    get_hloss.$(OBJ_EXT) \
    get_is.$(OBJ_EXT) \
    get_philoss.$(OBJ_EXT) \
    get_smass.$(OBJ_EXT) \
    get_stats.$(OBJ_EXT) \
    get_walls_bc.$(OBJ_EXT) \
    ic_mod.$(OBJ_EXT) \
    in_bin_512.$(OBJ_EXT) \
    in_bin_512i.$(OBJ_EXT) \
    indices_mod.$(OBJ_EXT) \
    init_ab_m.$(OBJ_EXT) \
    init_fvars.$(OBJ_EXT) \
    init_namelist.$(OBJ_EXT) \
    init_resid.$(OBJ_EXT) \
    is_mod.$(OBJ_EXT) \
    iterate.$(OBJ_EXT) \
    k_epsilon_prop.$(OBJ_EXT) \
    kintheory_drag_ss.$(OBJ_EXT) \
    kintheory_energy_dissipation_ss.$(OBJ_EXT) \
    kintheory_mod.$(OBJ_EXT) \
    kintheory_u_s.$(OBJ_EXT) \
    kintheory_v_s.$(OBJ_EXT) \
    kintheory_w_s.$(OBJ_EXT) \
    leq_bicgs.$(OBJ_EXT) \
    leq_bicgst.$(OBJ_EXT) \
    leq_cg.$(OBJ_EXT) \
    leq_gmres.$(OBJ_EXT) \
    leqsol_mod.$(OBJ_EXT) \
    leq_sor.$(OBJ_EXT) \
    line_too_big.$(OBJ_EXT) \
    location_check.$(OBJ_EXT) \
    location.$(OBJ_EXT) \
    machine.$(OBJ_EXT) \
    machine_mod.$(OBJ_EXT) \
    make_upper_case.$(OBJ_EXT) \
    mark_phase_4_cor.$(OBJ_EXT) \
    matrix_mod.$(OBJ_EXT) \
    mfix.$(OBJ_EXT) \
    mfix_netcdf_mod.$(OBJ_EXT) \
    mflux_mod.$(OBJ_EXT) \
    mod_bc_i.$(OBJ_EXT) \
    mod_bc_j.$(OBJ_EXT) \
    mod_bc_k.$(OBJ_EXT) \
    open_file.$(OBJ_EXT) \
    open_files.$(OBJ_EXT) \
    out_array_c.$(OBJ_EXT) \
    out_array.$(OBJ_EXT) \
    out_array_kc.$(OBJ_EXT) \
    out_array_k.$(OBJ_EXT) \
    out_bin_512.$(OBJ_EXT) \
    out_bin_512i.$(OBJ_EXT) \
    out_bin_512r.$(OBJ_EXT) \
    out_bin_r.$(OBJ_EXT) \
    output_mod.$(OBJ_EXT) \
    parallel_mod.$(OBJ_EXT) \
    param1_mod.$(OBJ_EXT) \
    param_mod.$(OBJ_EXT) \
    parse_line.$(OBJ_EXT) \
    parse_mod.$(OBJ_EXT) \
    parse_resid_string.$(OBJ_EXT) \
    parse_rxn.$(OBJ_EXT) \
    partial_elim.$(OBJ_EXT) \
    pgcor_mod.$(OBJ_EXT) \
    physical_prop.$(OBJ_EXT) \
    physprop_mod.$(OBJ_EXT) \
    pscor_mod.$(OBJ_EXT) \
    read_database.$(OBJ_EXT) \
    read_namelist.$(OBJ_EXT) \
    read_res0.$(OBJ_EXT) \
    read_res1.$(OBJ_EXT) \
    remove_comment.$(OBJ_EXT) \
    reset_new.$(OBJ_EXT) \
    residual_mod.$(OBJ_EXT) \
    rrates0.$(OBJ_EXT) \
    rrates.$(OBJ_EXT) \
    rrates_init.$(OBJ_EXT) \
    run_mod.$(OBJ_EXT) \
    rxns_mod.$(OBJ_EXT) \
    scalar_prop.$(OBJ_EXT) \
    scalars_mod.$(OBJ_EXT) \
    scales_mod.$(OBJ_EXT) \
    seek_comment.$(OBJ_EXT) \
    seek_end.$(OBJ_EXT) \
    set_bc0.$(OBJ_EXT) \
    set_bc1.$(OBJ_EXT) \
    set_constants.$(OBJ_EXT) \
    set_constprop.$(OBJ_EXT) \
    set_flags.$(OBJ_EXT) \
    set_fluidbed_p.$(OBJ_EXT) \
    set_geometry1.$(OBJ_EXT) \
    set_geometry.$(OBJ_EXT) \
    set_ic.$(OBJ_EXT) \
    set_increments3.$(OBJ_EXT) \
    set_increments.$(OBJ_EXT) \
    set_index1a3.$(OBJ_EXT) \
    set_index1a.$(OBJ_EXT) \
    set_index1.$(OBJ_EXT) \
    set_l_scale.$(OBJ_EXT) \
    set_max2.$(OBJ_EXT) \
    set_mw_mix_g.$(OBJ_EXT) \
    set_outflow.$(OBJ_EXT) \
    set_ro_g.$(OBJ_EXT) \
    set_wall_bc.$(OBJ_EXT) \
    shift_dxyz.$(OBJ_EXT) \
    solve_continuity.$(OBJ_EXT) \
    solve_energy_eq.$(OBJ_EXT) \
    solve_epp.$(OBJ_EXT) \
    solve_granular_energy.$(OBJ_EXT) \
    solve_k_epsilon_eq.$(OBJ_EXT) \
    solve_lin_eq.$(OBJ_EXT) \
    solve_pp_g.$(OBJ_EXT) \
    solve_scalar_eq.$(OBJ_EXT) \
    solve_species_eq.$(OBJ_EXT) \
    solve_vel_star.$(OBJ_EXT) \
    source_granular_energy.$(OBJ_EXT) \
    source_phi.$(OBJ_EXT) \
    source_pp_g.$(OBJ_EXT) \
    source_rop_g.$(OBJ_EXT) \
    source_rop_s.$(OBJ_EXT) \
    source_u_g.$(OBJ_EXT) \
    source_u_s.$(OBJ_EXT) \
    source_v_g.$(OBJ_EXT) \
    source_v_s.$(OBJ_EXT) \
    source_w_g.$(OBJ_EXT) \
    source_w_s.$(OBJ_EXT) \
    tau_g_mod.$(OBJ_EXT) \
    tau_s_mod.$(OBJ_EXT) \
    tau_u_g.$(OBJ_EXT) \
    tau_u_s.$(OBJ_EXT) \
    tau_v_g.$(OBJ_EXT) \
    tau_v_s.$(OBJ_EXT) \
    tau_w_g.$(OBJ_EXT) \
    tau_w_s.$(OBJ_EXT) \
    test_lin_eq.$(OBJ_EXT) \
    time_cpu_mod.$(OBJ_EXT) \
    time_march.$(OBJ_EXT) \
    tmp_array1_mod.$(OBJ_EXT) \
    tmp_array_mod.$(OBJ_EXT) \
    toleranc_mod.$(OBJ_EXT) \
    trace_mod.$(OBJ_EXT) \
    transfer.$(OBJ_EXT) \
    transport_prop.$(OBJ_EXT) \
    turb_mod.$(OBJ_EXT) \
    undef_2_0.$(OBJ_EXT) \
    under_relax.$(OBJ_EXT) \
    update_old.$(OBJ_EXT) \
    ur_facs_mod.$(OBJ_EXT) \
    usr0.$(OBJ_EXT) \
    usr1.$(OBJ_EXT) \
    usr2.$(OBJ_EXT) \
    usr3.$(OBJ_EXT) \
    usr_init_namelist.$(OBJ_EXT) \
    usr_mod.$(OBJ_EXT) \
    usr_write_out0.$(OBJ_EXT) \
    usr_write_out1.$(OBJ_EXT) \
    utilities.$(OBJ_EXT) \
    vavg_u_g.$(OBJ_EXT) \
    vavg_u_s.$(OBJ_EXT) \
    vavg_v_g.$(OBJ_EXT) \
    vavg_v_s.$(OBJ_EXT) \
    vavg_w_g.$(OBJ_EXT) \
    vavg_w_s.$(OBJ_EXT) \
    vf_gs_x.$(OBJ_EXT) \
    vf_gs_y.$(OBJ_EXT) \
    vf_gs_z.$(OBJ_EXT) \
    visc_g_mod.$(OBJ_EXT) \
    visc_s_mod.$(OBJ_EXT) \
    vshear_mod.$(OBJ_EXT) \
    vtc_scalar.$(OBJ_EXT) \
    write_ab_m.$(OBJ_EXT) \
    write_ab_m_var.$(OBJ_EXT) \
    write_error.$(OBJ_EXT) \
    write_header.$(OBJ_EXT) \
    write_out0.$(OBJ_EXT) \
    write_out1.$(OBJ_EXT) \
    write_out3.$(OBJ_EXT) \
    write_res0.$(OBJ_EXT) \
    write_res1.$(OBJ_EXT) \
    write_spx0.$(OBJ_EXT) \
    write_spx1.$(OBJ_EXT) \
    write_table.$(OBJ_EXT) \
    write_usr0.$(OBJ_EXT) \
    write_usr1.$(OBJ_EXT) \
    xerbla.$(OBJ_EXT) \
    xsi_array_mod.$(OBJ_EXT) \
    zero_array.$(OBJ_EXT) \
    zero_norm_vel.$(OBJ_EXT) \
    allocate_cut_cell_arrays.$(OBJ_EXT) \
    allocate_dummy_cut_cell_arrays.$(OBJ_EXT) \
    calc_vort_out.$(OBJ_EXT) \
    cartesian_grid_init_namelist.$(OBJ_EXT) \
    CG_set_bc0.$(OBJ_EXT) \
    CG_set_outflow.$(OBJ_EXT) \
    CG_source_u_g.$(OBJ_EXT) \
    CG_source_u_s.$(OBJ_EXT) \
    CG_source_v_g.$(OBJ_EXT) \
    CG_source_v_s.$(OBJ_EXT) \
    CG_source_w_g.$(OBJ_EXT) \
    CG_source_w_s.$(OBJ_EXT) \
    check_data_cartesian.$(OBJ_EXT) \
    cutcell_mod.$(OBJ_EXT) \
    cut_cell_preprocessing.$(OBJ_EXT) \
    dashboard_mod.$(OBJ_EXT) \
    deallocate_cut_cell_arrays.$(OBJ_EXT) \
    define_quadrics.$(OBJ_EXT) \
    dmp_cartesian.$(OBJ_EXT) \
    eval_usr_fct.$(OBJ_EXT) \
    get_alpha.$(OBJ_EXT) \
    get_connectivity.$(OBJ_EXT) \
    get_cut_cell_flags.$(OBJ_EXT) \
    get_cut_cell_volume_area.$(OBJ_EXT) \
    get_delh.$(OBJ_EXT) \
    get_master.$(OBJ_EXT) \
    get_poly_data.$(OBJ_EXT) \
    get_stl_data.$(OBJ_EXT) \
    polygon_mod.$(OBJ_EXT) \
    progress_bar_mod.$(OBJ_EXT) \
    quadric_mod.$(OBJ_EXT) \
    set_Odxyz.$(OBJ_EXT) \
    stl_mod.$(OBJ_EXT) \
    update_dashboard.$(OBJ_EXT) \
    vtk_mod.$(OBJ_EXT) \
    vtk_out.$(OBJ_EXT) \
    write_progress_bar.$(OBJ_EXT) \
    calc_jacobian.$(OBJ_EXT) \
    check_data_chem.$(OBJ_EXT) \
    dgpadm.$(OBJ_EXT) \
    exponential.$(OBJ_EXT) \
    fex.$(OBJ_EXT) \
    g_derivs.$(OBJ_EXT) \
    jac.$(OBJ_EXT) \
    mchem_init.$(OBJ_EXT) \
    mchem_mod.$(OBJ_EXT) \
    mchem_odepack_init.$(OBJ_EXT) \
    mchem_time_march.$(OBJ_EXT) \
    misat_table_init.$(OBJ_EXT) \
    react.$(OBJ_EXT) \
    usrfg.$(OBJ_EXT) \
    add_part_to_link_list.$(OBJ_EXT) \
    calc_app_coh_force.$(OBJ_EXT) \
    calc_cap_coh_force.$(OBJ_EXT) \
    calc_cohesive_forces.$(OBJ_EXT) \
    calc_esc_coh_force.$(OBJ_EXT) \
    calc_square_well.$(OBJ_EXT) \
    calc_van_der_waals.$(OBJ_EXT) \
    check_link.$(OBJ_EXT) \
    check_sw_wall_interaction.$(OBJ_EXT) \
    check_vdw_wall_interaction.$(OBJ_EXT) \
    initialize_cohesion_parameters.$(OBJ_EXT) \
    initialize_coh_int_search.$(OBJ_EXT) \
    linked_interaction_eval.$(OBJ_EXT) \
    remove_part_from_link_list.$(OBJ_EXT) \
    unlinked_interaction_eval.$(OBJ_EXT) \
    update_search_grids.$(OBJ_EXT) \
    calc_attrition_des.$(OBJ_EXT) \
    calc_force_des.$(OBJ_EXT) \
    calc_rrate_des.$(OBJ_EXT) \
    calc_thermo_des.$(OBJ_EXT) \
    cell_near_wall.$(OBJ_EXT) \
    cfassign.$(OBJ_EXT) \
    cffctowall.$(OBJ_EXT) \
    cffctow.$(OBJ_EXT) \
    cfnewvalues.$(OBJ_EXT) \
    cfrelvel.$(OBJ_EXT) \
    cfslide.$(OBJ_EXT) \
    cfslidewall.$(OBJ_EXT) \
    cfupdateold.$(OBJ_EXT) \
    cfwallcontact.$(OBJ_EXT) \
    cfwallposvel.$(OBJ_EXT) \
    check_des_bc.$(OBJ_EXT) \
    check_des_data.$(OBJ_EXT) \
    check_des_ic.$(OBJ_EXT) \
    check_des_rxns.$(OBJ_EXT) \
    check_des_thermo.$(OBJ_EXT) \
    des_allocate_arrays.$(OBJ_EXT) \
    des_bc_mod.$(OBJ_EXT) \
    des_check_particle.$(OBJ_EXT) \
    des_functions.$(OBJ_EXT) \
    des_granular_temperature.$(OBJ_EXT) \
    desgrid_mod.$(OBJ_EXT) \
    des_ic_mod.$(OBJ_EXT) \
    des_init_arrays.$(OBJ_EXT) \
    des_init_bc.$(OBJ_EXT) \
    des_init_namelist.$(OBJ_EXT) \
    des_mass_inlet.$(OBJ_EXT) \
    desmpi_mod.$(OBJ_EXT) \
    desmpi_wrapper_mod.$(OBJ_EXT) \
    des_physical_prop.$(OBJ_EXT) \
    des_reaction_model.$(OBJ_EXT) \
    des_rrates.$(OBJ_EXT) \
    des_rxns_mod.$(OBJ_EXT) \
    des_set_ic.$(OBJ_EXT) \
    des_thermo_cond.$(OBJ_EXT) \
    des_thermo_conv.$(OBJ_EXT) \
    des_thermo_mod.$(OBJ_EXT) \
    des_thermo_newvalues.$(OBJ_EXT) \
    des_thermo_rad.$(OBJ_EXT) \
    des_time_march.$(OBJ_EXT) \
    des_wallbc_preprocessing.$(OBJ_EXT) \
    discretelement_mod.$(OBJ_EXT) \
    drag_fgs.$(OBJ_EXT) \
    gas_drag.$(OBJ_EXT) \
    generate_particle_config.$(OBJ_EXT) \
    grid_based_neighbor_search.$(OBJ_EXT) \
    interpolation_mod.$(OBJ_EXT) \
    make_arrays_des.$(OBJ_EXT) \
    mfix_pic_mod.$(OBJ_EXT) \
    mppic_routines.$(OBJ_EXT) \
    mppic_wallbc_mod.$(OBJ_EXT) \
    neighbour.$(OBJ_EXT) \
    nsquare.$(OBJ_EXT) \
    octree.$(OBJ_EXT) \
    particles_in_cell.$(OBJ_EXT) \
    quadtree.$(OBJ_EXT) \
    randomno_mod.$(OBJ_EXT) \
    read_des_restart.$(OBJ_EXT) \
    sendrecvnode_mod.$(OBJ_EXT) \
    thermo_nbr.$(OBJ_EXT) \
    walledgecontact.$(OBJ_EXT) \
    wallfacecontact.$(OBJ_EXT) \
    wallnodecontact.$(OBJ_EXT) \
    write_des_data.$(OBJ_EXT) \
    write_des_restart.$(OBJ_EXT) \
    compar_mod.$(OBJ_EXT) \
    dbg_util_mod.$(OBJ_EXT) \
    debug_mod.$(OBJ_EXT) \
    gridmap_mod.$(OBJ_EXT) \
    mpi_mod.$(OBJ_EXT) \
    mpi_utility_mod.$(OBJ_EXT) \
    parallel_mpi_mod.$(OBJ_EXT) \
    sendrecv3_mod.$(OBJ_EXT) \
    sendrecv_mod.$(OBJ_EXT) \
    gaussj.$(OBJ_EXT) \
    odeint.$(OBJ_EXT) \
    rkck.$(OBJ_EXT) \
    rkqs.$(OBJ_EXT) \
    source_population_eq.$(OBJ_EXT) \
    usr_dqmom.$(OBJ_EXT) \
    adjust_eps_ghd.$(OBJ_EXT) \
    bulk_viscosity.$(OBJ_EXT) \
    calc_d_ghd.$(OBJ_EXT) \
    calc_external_forces.$(OBJ_EXT) \
    calc_nflux.$(OBJ_EXT) \
    chi_ij_GHD.$(OBJ_EXT) \
    cooling_rate.$(OBJ_EXT) \
    cooling_rate_tc.$(OBJ_EXT) \
    dufour_coeff.$(OBJ_EXT) \
    ghd.$(OBJ_EXT) \
    ghdmassflux.$(OBJ_EXT) \
    ghdtheory_mod.$(OBJ_EXT) \
    mass_mobility.$(OBJ_EXT) \
    ordinary_diff.$(OBJ_EXT) \
    partial_elim_ghd.$(OBJ_EXT) \
    pressure.$(OBJ_EXT) \
    shear_viscosity.$(OBJ_EXT) \
    source_ghd_granular_energy.$(OBJ_EXT) \
    thermal_conductivity.$(OBJ_EXT) \
    thermal_diffusivity.$(OBJ_EXT) \
    thermal_mobility.$(OBJ_EXT) \
    transport_coeff_ghd.$(OBJ_EXT) \
    qmomk_allocate_arrays.$(OBJ_EXT) \
    qmomk_bc_mod.$(OBJ_EXT) \
    qmomk_collision_mod.$(OBJ_EXT) \
    qmomk_fluxes_mod.$(OBJ_EXT) \
    qmomk_gas_drag.$(OBJ_EXT) \
    qmom_kinetic_equation_mod.$(OBJ_EXT) \
    qmomk_init_bc.$(OBJ_EXT) \
    qmomk_initial_conditions.$(OBJ_EXT) \
    qmomk_init_namelist.$(OBJ_EXT) \
    qmomk_make_arrays.$(OBJ_EXT) \
    qmomk_parameters_mod.$(OBJ_EXT) \
    qmomk_quadrature_mod.$(OBJ_EXT) \
    qmomk_read_restart.$(OBJ_EXT) \
    qmomk_set_bc.$(OBJ_EXT) \
    qmomk_time_march.$(OBJ_EXT) \
    qmomk_tools_mod.$(OBJ_EXT) \
    qmomk_write_restart.$(OBJ_EXT) \
    get_values.$(OBJ_EXT) \
    readTherm.$(OBJ_EXT) \
  -o mfix.exe $(LIB_FLAGS)
  
blas90.a : BLAS.o
	ar cr blas90.a BLAS.o
BLAS.o : BLAS.F
	$(FORTRAN_CMD) $(FORT_FLAGS) BLAS.F
dgtsv90.a : DGTSV.o
	ar cr dgtsv90.a DGTSV.o
DGTSV.o : DGTSV.F
	$(FORTRAN_CMD) $(FORT_FLAGS) DGTSV.F
odepack.a : ODEPACK.o
	ar cr odepack.a ODEPACK.o
ODEPACK.o : ODEPACK.F
	$(FORTRAN_CMD) $(FORT_FLAGS3) ODEPACK.F
AMBM.mod : ambm_mod.f \
            PARAM.mod \
            PARAM1.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ambm_mod.f 
BC.mod : bc_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) bc_mod.f 
BOUNDFUNIJK3.mod : boundfunijk3_mod.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk3_mod.f 
BOUNDFUNIJK.mod : boundfunijk_mod.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk_mod.f 
CDIST.mod : cdist_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) cdist_mod.f 
CHECK.mod : check_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_mod.f 
CHISCHEME.mod : chischeme_mod.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) chischeme_mod.f 
COEFF.mod : coeff_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) coeff_mod.f 
CONSTANT.mod : constant_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) constant_mod.f 
CONT.mod : cont_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) cont_mod.f 
CORNER.mod : corner_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) corner_mod.f 
DRAG.mod : drag_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) drag_mod.f 
ENERGY.mod : energy_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) energy_mod.f 
FLDVAR.mod : fldvar_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) fldvar_mod.f 
FUNCTION.mod : function_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) function_mod.f 
FUNITS.mod : funits_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) funits_mod.f 
GEOMETRY.mod : geometry_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) geometry_mod.f 
IC.mod : ic_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ic_mod.f 
INDICES.mod : indices_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) indices_mod.f 
IS.mod : is_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) is_mod.f 
KINTHEORY.mod : kintheory_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_mod.f 
LEQSOL.mod : leqsol_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) leqsol_mod.f 
MACHINE.mod : machine_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) machine_mod.f 
MATRIX.mod : matrix_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) matrix_mod.f 
MFIX_NETCDF.mod : mfix_netcdf_mod.f \
            MFIX_netcdf_constants.fi                                     \
            MFIX_netcdf_overloads.fi                                     \
            MFIX_netcdf_variables.fi                                     \
            MFIX_netcdf_misc.fi                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) mfix_netcdf_mod.f 
MFLUX.mod : mflux_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) mflux_mod.f 
OUTPUT.mod : output_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) output_mod.f 
PARALLEL.mod : parallel_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parallel_mod.f 
PARAM1.mod : param1_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) param1_mod.f 
PARAM.mod : param_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) param_mod.f 
PARSE.mod : parse_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_mod.f 
PGCOR.mod : pgcor_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pgcor_mod.f 
PHYSPROP.mod : physprop_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) physprop_mod.f 
PSCOR.mod : pscor_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pscor_mod.f 
RESIDUAL.mod : residual_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) residual_mod.f 
RUN.mod : run_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) run_mod.f 
RXNS.mod : rxns_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) rxns_mod.f 
SCALARS.mod : scalars_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scalars_mod.f 
SCALES.mod : scales_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scales_mod.f 
TAU_G.mod : tau_g_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_g_mod.f 
TAU_S.mod : tau_s_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_s_mod.f 
TIME_CPU.mod : time_cpu_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) time_cpu_mod.f 
TMP_ARRAY1.mod : tmp_array1_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array1_mod.f 
TMP_ARRAY.mod : tmp_array_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array_mod.f 
TOLERANC.mod : toleranc_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) toleranc_mod.f 
TRACE.mod : trace_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) trace_mod.f 
TURB.mod : turb_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) turb_mod.f 
UR_FACS.mod : ur_facs_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ur_facs_mod.f 
USR.mod : usr_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_mod.f 
VISC_G.mod : visc_g_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_g_mod.f 
VISC_S.mod : visc_s_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_s_mod.f 
VSHEAR.mod : vshear_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) vshear_mod.f 
XSI_ARRAY.mod : xsi_array_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) xsi_array_mod.f 
CUTCELL.mod : ./cartesian_grid/cutcell_mod.f \
            PARAM.mod \
            PARAM1.mod \
            PROGRESS_BAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cutcell_mod.f 
DASHBOARD.mod : ./cartesian_grid/dashboard_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/dashboard_mod.f 
POLYGON.mod : ./cartesian_grid/polygon_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/polygon_mod.f 
PROGRESS_BAR.mod : ./cartesian_grid/progress_bar_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/progress_bar_mod.f 
QUADRIC.mod : ./cartesian_grid/quadric_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/quadric_mod.f 
STL.mod : ./cartesian_grid/stl_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/stl_mod.f 
VTK.mod : ./cartesian_grid/vtk_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/vtk_mod.f 
MCHEM.mod : ./chem/mchem_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/mchem_mod.f 
DES_BC.mod : ./des/des_bc_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_bc_mod.f 
DESGRID.mod : ./des/desgrid_mod.f \
            PARAM1.mod \
            FUNITS.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod \
            CONSTANT.mod \
            DESMPI_WRAPPER.mod \
            DES_THERMO.mod \
            des/desgrid_functions.inc                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desgrid_mod.f 
DES_IC.mod : ./des/des_ic_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_ic_mod.f 
DESMPI.mod : ./des/desmpi_mod.f \
            PARALLEL_MPI.mod \
            MPI_UTILITY.mod \
            DISCRETELEMENT.mod \
            DESGRID.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            SENDRECV.mod \
            DES_BC.mod \
            DESMPI_WRAPPER.mod \
            SENDRECVNODE.mod \
            MFIX_PIC.mod \
            des/desgrid_functions.inc                                    \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desmpi_mod.f 
DESMPI_WRAPPER.mod : ./des/desmpi_wrapper_mod.f \
            PARALLEL_MPI.mod \
            MPI_UTILITY.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desmpi_wrapper_mod.f 
DES_RXNS.mod : ./des/des_rxns_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_rxns_mod.f 
DES_THERMO.mod : ./des/des_thermo_mod.f \
            PARAM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_mod.f 
DISCRETELEMENT.mod : ./des/discretelement_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/discretelement_mod.f 
INTERPOLATION.mod : ./des/interpolation_mod.f \
            CONSTANT.mod \
            DISCRETELEMENT.mod \
            GEOMETRY.mod \
            PARAM1.mod \
            COMPAR.mod \
            INDICES.mod \
            PARAM.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/interpolation_mod.f 
MFIX_PIC.mod : ./des/mfix_pic_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mfix_pic_mod.f 
MPPIC_WALLBC.mod : ./des/mppic_wallbc_mod.f \
            PARAM.mod \
            PARAM1.mod \
            DISCRETELEMENT.mod \
            BC.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            FUNITS.mod \
            MPI_UTILITY.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            RANDOMNO.mod \
            CUTCELL.mod \
            FLDVAR.mod \
            MFIX_PIC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mppic_wallbc_mod.f 
RANDOMNO.mod : ./des/randomno_mod.f \
            CONSTANT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/randomno_mod.f 
SENDRECVNODE.mod : ./des/sendrecvnode_mod.f \
            PARALLEL_MPI.mod \
            MPI_UTILITY.mod \
            DISCRETELEMENT.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            SENDRECV.mod \
            DESMPI_WRAPPER.mod \
            DESGRID.mod \
            function.inc                                                 \
            des/desgrid_functions.inc                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/sendrecvnode_mod.f 
COMPAR.mod : ./dmp_modules/compar_mod.f \
            MPI.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/compar_mod.f 
DBG_UTIL.mod : ./dmp_modules/dbg_util_mod.f \
            COMPAR.mod \
            GEOMETRY.mod \
            PARALLEL_MPI.mod \
            INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/dbg_util_mod.f 
DEBUG.mod : ./dmp_modules/debug_mod.f \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/debug_mod.f 
GRIDMAP.mod : ./dmp_modules/gridmap_mod.f \
            MPI_UTILITY.mod \
            PARALLEL_MPI.mod \
            GEOMETRY.mod \
            SENDRECV.mod \
            COMPAR.mod \
            RUN.mod \
            INDICES.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/gridmap_mod.f 
MPI.mod : ./dmp_modules/mpi_mod.f \
            mpif.h                                                      
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_mod.f 
MPI_UTILITY.mod : ./dmp_modules/mpi_utility_mod.f \
            GEOMETRY.mod \
            COMPAR.mod \
            PARALLEL_MPI.mod \
            DEBUG.mod \
            INDICES.mod \
            FUNITS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_utility_mod.f 
PARALLEL_MPI.mod : ./dmp_modules/parallel_mpi_mod.f \
            GEOMETRY.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/parallel_mpi_mod.f 
SENDRECV3.mod : ./dmp_modules/sendrecv3_mod.f \
            PARALLEL_MPI.mod \
            DEBUG.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            MPI.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv3_mod.f 
SENDRECV.mod : ./dmp_modules/sendrecv_mod.f \
            PARALLEL_MPI.mod \
            DEBUG.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            MPI.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv_mod.f 
GHDTHEORY.mod : ./GhdTheory/ghdtheory_mod.f \
            PARAM.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ghdtheory_mod.f 
QMOMK_BC.mod : ./qmomk/qmomk_bc_mod.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            BC.mod \
            QMOM_KINETIC_EQUATION.mod \
            QMOMK_QUADRATURE.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_bc_mod.f 
QMOMK_COLLISION.mod : ./qmomk/qmomk_collision_mod.f \
            CONSTANT.mod \
            QMOMK_PARAMETERS.mod \
            QMOMK_QUADRATURE.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_collision_mod.f 
QMOMK_FLUXES.mod : ./qmomk/qmomk_fluxes_mod.f \
            QMOMK_PARAMETERS.mod \
            QMOMK_COLLISION.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_fluxes_mod.f 
QMOM_KINETIC_EQUATION.mod : ./qmomk/qmom_kinetic_equation_mod.f \
            PARAM.mod \
            PARAM1.mod \
            QMOMK_PARAMETERS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmom_kinetic_equation_mod.f 
QMOMK_PARAMETERS.mod : ./qmomk/qmomk_parameters_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_parameters_mod.f 
QMOMK_QUADRATURE.mod : ./qmomk/qmomk_quadrature_mod.f \
            QMOMK_TOOLS.mod \
            QMOMK_PARAMETERS.mod \
            CONSTANT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_quadrature_mod.f 
QMOMK_TOOLS.mod : ./qmomk/qmomk_tools_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_tools_mod.f 
accum_resid.$(OBJ_EXT) : accum_resid.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            RESIDUAL.mod \
            RUN.mod 
adjust_a_u_g.$(OBJ_EXT) : adjust_a_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            USR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_u_s.$(OBJ_EXT) : adjust_a_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_v_g.$(OBJ_EXT) : adjust_a_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_v_s.$(OBJ_EXT) : adjust_a_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_w_g.$(OBJ_EXT) : adjust_a_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_w_s.$(OBJ_EXT) : adjust_a_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            SENDRECV.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_dt.$(OBJ_EXT) : adjust_dt.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
adjust_eps.$(OBJ_EXT) : adjust_eps.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
adjust_leq.$(OBJ_EXT) : adjust_leq.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            LEQSOL.mod 
adjust_rop.$(OBJ_EXT) : adjust_rop.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
adjust_theta.$(OBJ_EXT) : adjust_theta.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            COMPAR.mod \
            function.inc                                                
allocate_arrays.$(OBJ_EXT) : allocate_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            AMBM.mod \
            COEFF.mod \
            CONT.mod \
            DRAG.mod \
            ENERGY.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            PHYSPROP.mod \
            PSCOR.mod \
            RESIDUAL.mod \
            RXNS.mod \
            RUN.mod \
            SCALARS.mod \
            TURB.mod \
            TAU_G.mod \
            TAU_S.mod \
            TMP_ARRAY.mod \
            TMP_ARRAY1.mod \
            TRACE.mod \
            VISC_G.mod \
            VISC_S.mod \
            XSI_ARRAY.mod \
            VSHEAR.mod \
            MFLUX.mod \
            MCHEM.mod \
            GHDTHEORY.mod \
            KINTHEORY.mod \
            CDIST.mod 
bc_phi.$(OBJ_EXT) : bc_phi.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            GEOMETRY.mod \
            OUTPUT.mod \
            INDICES.mod \
            BC.mod \
            COMPAR.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
bc_theta.$(OBJ_EXT) : bc_theta.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            GEOMETRY.mod \
            OUTPUT.mod \
            INDICES.mod \
            BC.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            TURB.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
b_m_p_star.$(OBJ_EXT) : b_m_p_star.f \
            PARAM.mod \
            PARAM1.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            RXNS.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
bound_x.$(OBJ_EXT) : bound_x.f \
            PARAM.mod \
            PARAM1.mod 
calc_cell.$(OBJ_EXT) : calc_cell.f \
            PARAM.mod \
            PARAM1.mod 
calc_coeff.$(OBJ_EXT) : calc_coeff.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RXNS.mod \
            FUNITS.mod \
            COMPAR.mod \
            UR_FACS.mod \
            RUN.mod \
            DISCRETELEMENT.mod \
            DES_RXNS.mod 
calc_d.$(OBJ_EXT) : calc_d.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            SCALES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            CUTCELL.mod \
            QMOM_KINETIC_EQUATION.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_dif_g.$(OBJ_EXT) : calc_dif_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            CONSTANT.mod \
            SCALES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            RUN.mod \
            function.inc                                                
calc_dif_s.$(OBJ_EXT) : calc_dif_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            RUN.mod \
            function.inc                                                
calc_drag.$(OBJ_EXT) : calc_drag.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            DRAG.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod \
            QMOM_KINETIC_EQUATION.mod 
calc_e.$(OBJ_EXT) : calc_e.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            CONSTANT.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_gama.$(OBJ_EXT) : calc_gama.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            ENERGY.mod \
            RXNS.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            DISCRETELEMENT.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
calc_grbdry.$(OBJ_EXT) : calc_grbdry.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            TURB.mod \
            VISC_S.mod \
            GEOMETRY.mod \
            INDICES.mod \
            BC.mod \
            COMPAR.mod \
            TOLERANC.mod \
            MPI_UTILITY.mod \
            CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_h.$(OBJ_EXT) : calc_h.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            DES_RXNS.mod \
            DES_THERMO.mod \
            DISCRETELEMENT.mod 
calc_k_cp.$(OBJ_EXT) : calc_k_cp.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            PSCOR.mod \
            GEOMETRY.mod \
            CONSTANT.mod \
            RUN.mod \
            VISC_S.mod \
            TRACE.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
calc_k_g.$(OBJ_EXT) : calc_k_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            CONSTANT.mod \
            COMPAR.mod \
            RUN.mod \
            SENDRECV.mod \
            function.inc                                                
calc_k_s.$(OBJ_EXT) : calc_k_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            RUN.mod \
            function.inc                                                
calc_mflux.$(OBJ_EXT) : calc_mflux.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            MFLUX.mod \
            PHYSPROP.mod \
            RUN.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
calc_mu_g.$(OBJ_EXT) : calc_mu_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            VISC_G.mod \
            VISC_S.mod \
            INDICES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            COMPAR.mod \
            DRAG.mod \
            RUN.mod \
            TURB.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                
calc_mu_s.$(OBJ_EXT) : calc_mu_s.f \
            RUN.mod \
            VSHEAR.mod \
            VISC_S.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            COMPAR.mod \
            INDICES.mod \
            GEOMETRY.mod \
            QMOM_KINETIC_EQUATION.mod \
            PARAM.mod \
            PARAM1.mod \
            TRACE.mod \
            TOLERANC.mod \
            TURB.mod \
            DRAG.mod \
            KINTHEORY.mod \
            UR_FACS.mod \
            PARALLEL.mod \
            VISC_G.mod \
            IS.mod \
            SENDRECV.mod \
            CUTCELL.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            s_pr1.inc                                                    \
            s_pr2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
calc_mw.$(OBJ_EXT) : calc_mw.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod 
calc_outflow.$(OBJ_EXT) : calc_outflow.f \
            PARAM.mod \
            PARAM1.mod \
            BC.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
calc_p_star.$(OBJ_EXT) : calc_p_star.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            PGCOR.mod \
            PSCOR.mod \
            UR_FACS.mod \
            RESIDUAL.mod \
            COMPAR.mod \
            RUN.mod \
            VISC_S.mod \
            FLDVAR.mod \
            TOLERANC.mod \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_resid.$(OBJ_EXT) : calc_resid.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            RUN.mod \
            FLDVAR.mod \
            BC.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            RESIDUAL.mod \
            RXNS.mod \
            MFLUX.mod \
            function.inc                                                
calc_s_ddot_s.$(OBJ_EXT) : calc_s_ddot_s.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
calc_trd_g.$(OBJ_EXT) : calc_trd_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            BC.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            function.inc                                                
calc_trd_s.$(OBJ_EXT) : calc_trd_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            SENDRECV.mod \
            BC.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            function.inc                                                
calc_u_friction.$(OBJ_EXT) : calc_u_friction.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            TURB.mod \
            VISC_S.mod \
            GEOMETRY.mod \
            INDICES.mod \
            BC.mod \
            COMPAR.mod \
            TOLERANC.mod \
            MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_vol_fr.$(OBJ_EXT) : calc_vol_fr.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            VISC_S.mod \
            CONSTANT.mod \
            PGCOR.mod \
            PSCOR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_vol_fr.f 
calc_xsi.$(OBJ_EXT) : calc_xsi.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            VSHEAR.mod \
            CHISCHEME.mod \
            COMPAR.mod \
            SENDRECV.mod \
            xsi1.inc                                                     \
            function.inc                                                 \
            xsi2.inc                                                    
cal_d.$(OBJ_EXT) : cal_d.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            VSHEAR.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
check_ab_m.$(OBJ_EXT) : check_ab_m.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
check_convergence.$(OBJ_EXT) : check_convergence.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            RESIDUAL.mod \
            TOLERANC.mod \
            MPI_UTILITY.mod 
check_data_01.$(OBJ_EXT) : check_data_01.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            RUN.mod \
            PHYSPROP.mod \
            INDICES.mod \
            SCALARS.mod \
            FUNITS.mod 
check_data_02.$(OBJ_EXT) : check_data_02.f \
            PARAM.mod \
            PARAM1.mod \
            OUTPUT.mod \
            LEQSOL.mod \
            GEOMETRY.mod \
            RUN.mod \
            RXNS.mod 
check_data_03.$(OBJ_EXT) : check_data_03.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            BC.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
check_data_04.$(OBJ_EXT) : check_data_04.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            DISCRETELEMENT.mod \
            FUNITS.mod 
check_data_05.$(OBJ_EXT) : check_data_05.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FUNITS.mod \
            RUN.mod \
            INDICES.mod 
check_data_06.$(OBJ_EXT) : check_data_06.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            IC.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            INDICES.mod \
            FUNITS.mod \
            SCALARS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            RXNS.mod \
            function.inc                                                
check_data_07.$(OBJ_EXT) : check_data_07.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            BC.mod \
            INDICES.mod \
            FUNITS.mod \
            SCALARS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
check_data_08.$(OBJ_EXT) : check_data_08.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            IS.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
check_data_09.$(OBJ_EXT) : check_data_09.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            RXNS.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod 
check_data_20.$(OBJ_EXT) : check_data_20.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            VISC_G.mod \
            RXNS.mod \
            SCALARS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
check_data_30.$(OBJ_EXT) : check_data_30.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            RXNS.mod \
            VISC_S.mod \
            VISC_G.mod \
            GEOMETRY.mod \
            RUN.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            DISCRETELEMENT.mod \
            function.inc                                                
check_mass_balance.$(OBJ_EXT) : check_mass_balance.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            RXNS.mod \
            GEOMETRY.mod \
            RUN.mod \
            BC.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            OUTPUT.mod \
            CHECK.mod \
            MCHEM.mod \
            MFLUX.mod \
            XSI_ARRAY.mod \
            PARALLEL.mod \
            MATRIX.mod \
            function.inc                                                
check_one_axis.$(OBJ_EXT) : check_one_axis.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod 
check_plane.$(OBJ_EXT) : check_plane.f \
            FUNITS.mod \
            COMPAR.mod 
cn_extrapol.$(OBJ_EXT) : cn_extrapol.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            SCALARS.mod \
            TRACE.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
compare.$(OBJ_EXT) : compare.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
conv_dif_phi.$(OBJ_EXT) : conv_dif_phi.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            SENDRECV.mod \
            XSI_ARRAY.mod \
            MPI_UTILITY.mod \
            INDICES.mod \
            PARALLEL.mod \
            MATRIX.mod \
            TOLERANC.mod \
            CUTCELL.mod \
            SENDRECV3.mod \
            TMP_ARRAY.mod \
            VSHEAR.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            OUTPUT.mod \
            IS.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            function3.inc                                                \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
conv_dif_u_g.$(OBJ_EXT) : conv_dif_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            VISC_G.mod \
            COMPAR.mod \
            TOLERANC.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            MFLUX.mod \
            CUTCELL.mod \
            VSHEAR.mod \
            XSI_ARRAY.mod \
            TMP_ARRAY.mod \
            SENDRECV.mod \
            SENDRECV3.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_dif_u_s.$(OBJ_EXT) : conv_dif_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            PHYSPROP.mod \
            VISC_S.mod \
            COMPAR.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            MFLUX.mod \
            CUTCELL.mod \
            XSI_ARRAY.mod \
            TMP_ARRAY.mod \
            SENDRECV.mod \
            SENDRECV3.mod \
            VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_dif_v_g.$(OBJ_EXT) : conv_dif_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            VISC_G.mod \
            COMPAR.mod \
            TOLERANC.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            MFLUX.mod \
            CUTCELL.mod \
            XSI_ARRAY.mod \
            VSHEAR.mod \
            TMP_ARRAY.mod \
            SENDRECV.mod \
            SENDRECV3.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_dif_v_s.$(OBJ_EXT) : conv_dif_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            PHYSPROP.mod \
            VISC_S.mod \
            COMPAR.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            MFLUX.mod \
            CUTCELL.mod \
            XSI_ARRAY.mod \
            TMP_ARRAY.mod \
            SENDRECV.mod \
            SENDRECV3.mod \
            VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_dif_w_g.$(OBJ_EXT) : conv_dif_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            VISC_G.mod \
            COMPAR.mod \
            TOLERANC.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            MFLUX.mod \
            CUTCELL.mod \
            XSI_ARRAY.mod \
            TMP_ARRAY.mod \
            SENDRECV.mod \
            SENDRECV3.mod \
            VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_dif_w_s.$(OBJ_EXT) : conv_dif_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            PHYSPROP.mod \
            VISC_S.mod \
            COMPAR.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            MFLUX.mod \
            CUTCELL.mod \
            XSI_ARRAY.mod \
            TMP_ARRAY.mod \
            SENDRECV.mod \
            SENDRECV3.mod \
            VSHEAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_pp_g.$(OBJ_EXT) : conv_pp_g.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            PARALLEL.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            COMPAR.mod \
            MFLUX.mod \
            function.inc                                                
conv_rop.$(OBJ_EXT) : conv_rop.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            MFLUX.mod \
            PHYSPROP.mod \
            RUN.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            XSI_ARRAY.mod \
            function.inc                                                
conv_rop_g.$(OBJ_EXT) : conv_rop_g.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            COMPAR.mod \
            PARALLEL.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            XSI_ARRAY.mod \
            function.inc                                                
conv_rop_s.$(OBJ_EXT) : conv_rop_s.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            COMPAR.mod \
            PARALLEL.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            PSCOR.mod \
            XSI_ARRAY.mod \
            function.inc                                                
conv_source_epp.$(OBJ_EXT) : conv_source_epp.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            SENDRECV.mod \
            XSI_ARRAY.mod \
            PARALLEL.mod \
            MATRIX.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            RXNS.mod \
            INDICES.mod \
            PGCOR.mod \
            PSCOR.mod \
            VSHEAR.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
copy_a.$(OBJ_EXT) : copy_a.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            function.inc                                                
corner.$(OBJ_EXT) : corner.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            MATRIX.mod \
            CORNER.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
correct_0.$(OBJ_EXT) : correct_0.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PGCOR.mod \
            UR_FACS.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            function.inc                                                
correct_1.$(OBJ_EXT) : correct_1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            GEOMETRY.mod \
            PSCOR.mod \
            UR_FACS.mod \
            CONSTANT.mod \
            COMPAR.mod \
            SENDRECV.mod \
            CUTCELL.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
dgtsl.$(OBJ_EXT) : dgtsl.f 
dif_u_is.$(OBJ_EXT) : dif_u_is.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            OUTPUT.mod \
            INDICES.mod \
            IS.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
dif_v_is.$(OBJ_EXT) : dif_v_is.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            OUTPUT.mod \
            INDICES.mod \
            IS.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
dif_w_is.$(OBJ_EXT) : dif_w_is.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            OUTPUT.mod \
            INDICES.mod \
            IS.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
discretize.$(OBJ_EXT) : discretize.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod 
display_resid.$(OBJ_EXT) : display_resid.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RESIDUAL.mod \
            FLDVAR.mod \
            COMPAR.mod \
            GEOMETRY.mod 
drag_gs.$(OBJ_EXT) : drag_gs.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            CONSTANT.mod \
            COMPAR.mod \
            DRAG.mod \
            SENDRECV.mod \
            DISCRETELEMENT.mod \
            UR_FACS.mod \
            MFIX_PIC.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
drag_ss.$(OBJ_EXT) : drag_ss.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            SENDRECV.mod \
            DRAG.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
eosg.$(OBJ_EXT) : eosg.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            SCALES.mod \
            sc_p_g1.inc                                                  \
            sc_p_g2.inc                                                 
equal.$(OBJ_EXT) : equal.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            PHYSPROP.mod 
error_routine.$(OBJ_EXT) : error_routine.f \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
exchange.$(OBJ_EXT) : exchange.f \
            PARAM.mod \
            PARAM1.mod \
            COMPAR.mod 
exit.$(OBJ_EXT) : exit.f \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
flow_to_vel.$(OBJ_EXT) : flow_to_vel.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RUN.mod \
            BC.mod \
            SCALES.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod 
g_0.$(OBJ_EXT) : g_0.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            VISC_S.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
get_bc_area.$(OBJ_EXT) : get_bc_area.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            BC.mod \
            COMPAR.mod 
get_data.$(OBJ_EXT) : get_data.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            FUNITS.mod \
            COMPAR.mod \
            GRIDMAP.mod \
            DISCRETELEMENT.mod \
            LEQSOL.mod \
            PARALLEL.mod \
            QMOM_KINETIC_EQUATION.mod 
get_eq.$(OBJ_EXT) : get_eq.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod 
get_flow_bc.$(OBJ_EXT) : get_flow_bc.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            BC.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
get_hloss.$(OBJ_EXT) : get_hloss.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            BC.mod \
            INDICES.mod \
            ENERGY.mod 
get_is.$(OBJ_EXT) : get_is.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            IS.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod 
get_philoss.$(OBJ_EXT) : get_philoss.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            BC.mod \
            INDICES.mod \
            ENERGY.mod \
            COMPAR.mod \
            function.inc                                                
get_smass.$(OBJ_EXT) : get_smass.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
get_stats.$(OBJ_EXT) : get_stats.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            FUNITS.mod \
            RESIDUAL.mod \
            RUN.mod \
            COMPAR.mod \
            function.inc                                                
get_walls_bc.$(OBJ_EXT) : get_walls_bc.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            BC.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
in_bin_512.$(OBJ_EXT) : in_bin_512.f \
            MACHINE.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
in_bin_512i.$(OBJ_EXT) : in_bin_512i.f \
            MACHINE.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
init_ab_m.$(OBJ_EXT) : init_ab_m.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            PARALLEL.mod \
            COMPAR.mod 
init_fvars.$(OBJ_EXT) : init_fvars.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            SCALARS.mod \
            RXNS.mod \
            RUN.mod \
            COMPAR.mod 
init_namelist.$(OBJ_EXT) : init_namelist.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            IC.mod \
            BC.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            IS.mod \
            TOLERANC.mod \
            SCALES.mod \
            UR_FACS.mod \
            LEQSOL.mod \
            RESIDUAL.mod \
            RXNS.mod \
            SCALARS.mod \
            COMPAR.mod \
            PARALLEL.mod \
            CDIST.mod \
            namelist.inc                                                
init_resid.$(OBJ_EXT) : init_resid.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RESIDUAL.mod 
iterate.$(OBJ_EXT) : iterate.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            FUNITS.mod \
            TIME_CPU.mod \
            PSCOR.mod \
            COEFF.mod \
            LEQSOL.mod \
            VISC_G.mod \
            PGCOR.mod \
            CONT.mod \
            SCALARS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            DISCRETELEMENT.mod \
            RESIDUAL.mod \
            CUTCELL.mod \
            VTK.mod \
            DASHBOARD.mod \
            QMOM_KINETIC_EQUATION.mod \
            BC.mod \
            CONSTANT.mod 
k_epsilon_prop.$(OBJ_EXT) : k_epsilon_prop.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            DRAG.mod \
            RUN.mod \
            OUTPUT.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            VISC_G.mod \
            VISC_S.mod \
            TRACE.mod \
            INDICES.mod \
            CONSTANT.mod \
            VSHEAR.mod \
            TURB.mod \
            TOLERANC.mod \
            COMPAR.mod \
            TAU_G.mod \
            SENDRECV.mod \
            CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                
kintheory_drag_ss.$(OBJ_EXT) : kintheory_drag_ss.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            SENDRECV.mod \
            DRAG.mod \
            KINTHEORY.mod \
            function.inc                                                
kintheory_energy_dissipation_ss.$(OBJ_EXT) : kintheory_energy_dissipation_ss.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            KINTHEORY.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
kintheory_u_s.$(OBJ_EXT) : kintheory_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            KINTHEORY.mod \
            FLDVAR.mod \
            VISC_S.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
kintheory_v_s.$(OBJ_EXT) : kintheory_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            KINTHEORY.mod \
            FLDVAR.mod \
            VISC_S.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
kintheory_w_s.$(OBJ_EXT) : kintheory_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            KINTHEORY.mod \
            FLDVAR.mod \
            VISC_S.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
leq_bicgs.$(OBJ_EXT) : leq_bicgs.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            LEQSOL.mod \
            FUNITS.mod \
            PARALLEL.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            function.inc                                                
leq_bicgst.$(OBJ_EXT) : leq_bicgst.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            LEQSOL.mod \
            FUNITS.mod \
            PARALLEL.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            function.inc                                                
leq_cg.$(OBJ_EXT) : leq_cg.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            LEQSOL.mod \
            FUNITS.mod \
            PARALLEL.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            function.inc                                                
leq_gmres.$(OBJ_EXT) : leq_gmres.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            DEBUG.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PARALLEL.mod \
            FUNITS.mod \
            GRIDMAP.mod \
            function.inc                                                
leq_sor.$(OBJ_EXT) : leq_sor.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            LEQSOL.mod \
            function.inc                                                
line_too_big.$(OBJ_EXT) : line_too_big.f 
location_check.$(OBJ_EXT) : location_check.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod \
            GEOMETRY.mod 
location.$(OBJ_EXT) : location.f \
            PARAM.mod \
            PARAM1.mod 
machine.$(OBJ_EXT) : machine.f \
            MACHINE.mod \
            PARAM.mod \
            RUN.mod \
            FUNITS.mod 
make_upper_case.$(OBJ_EXT) : make_upper_case.f 
mark_phase_4_cor.$(OBJ_EXT) : mark_phase_4_cor.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            COMPAR.mod \
            VISC_S.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) mark_phase_4_cor.f 
mfix.$(OBJ_EXT) : mfix.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            TIME_CPU.mod \
            FUNITS.mod \
            OUTPUT.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PARALLEL_MPI.mod \
            DISCRETELEMENT.mod \
            MFIX_PIC.mod \
            CDIST.mod \
            FLDVAR.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            DASHBOARD.mod \
            PARALLEL.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            SENDRECV.mod \
            SENDRECV3.mod \
            INDICES.mod \
            LEQSOL.mod \
            function.inc                                                
mod_bc_i.$(OBJ_EXT) : mod_bc_i.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
mod_bc_j.$(OBJ_EXT) : mod_bc_j.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
mod_bc_k.$(OBJ_EXT) : mod_bc_k.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
open_file.$(OBJ_EXT) : open_file.f \
            CDIST.mod \
            COMPAR.mod 
open_files.$(OBJ_EXT) : open_files.f \
            MACHINE.mod \
            FUNITS.mod \
            COMPAR.mod \
            CDIST.mod \
            RUN.mod 
out_array_c.$(OBJ_EXT) : out_array_c.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
out_array.$(OBJ_EXT) : out_array.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
out_array_kc.$(OBJ_EXT) : out_array_kc.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            function.inc                                                
out_array_k.$(OBJ_EXT) : out_array_k.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            INDICES.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
out_bin_512.$(OBJ_EXT) : out_bin_512.f \
            MACHINE.mod 
out_bin_512i.$(OBJ_EXT) : out_bin_512i.f \
            MACHINE.mod 
out_bin_512r.$(OBJ_EXT) : out_bin_512r.f \
            MACHINE.mod 
out_bin_r.$(OBJ_EXT) : out_bin_r.f \
            PARAM.mod 
parse_line.$(OBJ_EXT) : parse_line.f \
            PARAM.mod \
            PARAM1.mod \
            PARSE.mod \
            COMPAR.mod 
parse_resid_string.$(OBJ_EXT) : parse_resid_string.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RESIDUAL.mod \
            FUNITS.mod \
            COMPAR.mod 
parse_rxn.$(OBJ_EXT) : parse_rxn.f \
            PARAM.mod \
            PARAM1.mod \
            PARSE.mod \
            RXNS.mod \
            COMPAR.mod 
partial_elim.$(OBJ_EXT) : partial_elim.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            DRAG.mod \
            FLDVAR.mod \
            RUN.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
physical_prop.$(OBJ_EXT) : physical_prop.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            TOLERANC.mod \
            CONSTANT.mod \
            SCALARS.mod \
            COMPAR.mod \
            FUNITS.mod \
            USR.mod \
            MPI_UTILITY.mod \
            DISCRETELEMENT.mod \
            CUTCELL.mod \
            usrnlst.inc                                                  \
            cp_fun1.inc                                                  \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
read_database.$(OBJ_EXT) : read_database.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            COMPAR.mod \
            RXNS.mod \
            FUNITS.mod \
            DISCRETELEMENT.mod \
            DES_RXNS.mod \
            mfix_directory_path.inc                                     
read_namelist.$(OBJ_EXT) : read_namelist.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            IC.mod \
            IS.mod \
            BC.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            TOLERANC.mod \
            FUNITS.mod \
            SCALES.mod \
            UR_FACS.mod \
            LEQSOL.mod \
            RESIDUAL.mod \
            RXNS.mod \
            SCALARS.mod \
            COMPAR.mod \
            PARALLEL.mod \
            DISCRETELEMENT.mod \
            MFIX_PIC.mod \
            USR.mod \
            DES_BC.mod \
            DES_IC.mod \
            DES_THERMO.mod \
            DES_RXNS.mod \
            CDIST.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            VTK.mod \
            POLYGON.mod \
            DASHBOARD.mod \
            STL.mod \
            QMOM_KINETIC_EQUATION.mod \
            usrnlst.inc                                                  \
            namelist.inc                                                 \
            des/desnamelist.inc                                          \
            cartesian_grid/cartesian_grid_namelist.inc                   \
            qmomk/qmomknamelist.inc                                     
read_res0.$(OBJ_EXT) : read_res0.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            IC.mod \
            BC.mod \
            IS.mod \
            CONSTANT.mod \
            FUNITS.mod \
            OUTPUT.mod \
            SCALES.mod \
            UR_FACS.mod \
            TOLERANC.mod \
            LEQSOL.mod \
            SCALARS.mod \
            RXNS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            FLDVAR.mod 
read_res1.$(OBJ_EXT) : read_res1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            RXNS.mod \
            SCALARS.mod \
            FUNITS.mod \
            ENERGY.mod \
            COMPAR.mod \
            CDIST.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            MFIX_NETCDF.mod 
remove_comment.$(OBJ_EXT) : remove_comment.f 
reset_new.$(OBJ_EXT) : reset_new.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            TRACE.mod \
            RUN.mod \
            SCALARS.mod 
rrates0.$(OBJ_EXT) : rrates0.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RXNS.mod \
            ENERGY.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            FUNITS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
rrates.$(OBJ_EXT) : rrates.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RXNS.mod \
            ENERGY.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            FUNITS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
rrates_init.$(OBJ_EXT) : rrates_init.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RXNS.mod \
            ENERGY.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
scalar_prop.$(OBJ_EXT) : scalar_prop.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            SCALARS.mod \
            TOLERANC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
seek_comment.$(OBJ_EXT) : seek_comment.f 
seek_end.$(OBJ_EXT) : seek_end.f 
set_bc0.$(OBJ_EXT) : set_bc0.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            BC.mod \
            FLDVAR.mod \
            INDICES.mod \
            RUN.mod \
            FUNITS.mod \
            SCALES.mod \
            SCALARS.mod \
            BOUNDFUNIJK.mod \
            TOLERANC.mod \
            SENDRECV.mod \
            sc_p_g1.inc                                                  \
            function.inc                                                 \
            sc_p_g2.inc                                                 
set_bc1.$(OBJ_EXT) : set_bc1.f \
            PARAM.mod \
            PARAM1.mod \
            BC.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
set_constants.$(OBJ_EXT) : set_constants.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            VISC_S.mod \
            ENERGY.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            RUN.mod \
            FUNITS.mod \
            DRAG.mod \
            COMPAR.mod 
set_constprop.$(OBJ_EXT) : set_constprop.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            VISC_S.mod \
            VISC_G.mod \
            ENERGY.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            RUN.mod \
            FUNITS.mod \
            DRAG.mod \
            COMPAR.mod \
            KINTHEORY.mod \
            function.inc                                                
set_flags.$(OBJ_EXT) : set_flags.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            BC.mod \
            IS.mod \
            INDICES.mod \
            PHYSPROP.mod \
            FUNITS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            SENDRECV3.mod \
            BOUNDFUNIJK.mod \
            MPI_UTILITY.mod \
            function.inc                                                 \
            function3.inc                                               
set_fluidbed_p.$(OBJ_EXT) : set_fluidbed_p.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            BC.mod \
            IC.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            FUNITS.mod \
            SCALES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            DISCRETELEMENT.mod \
            sc_p_g1.inc                                                  \
            b_force1.inc                                                 \
            function.inc                                                 \
            b_force2.inc                                                 \
            sc_p_g2.inc                                                 
set_geometry1.$(OBJ_EXT) : set_geometry1.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
set_geometry.$(OBJ_EXT) : set_geometry.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            GEOMETRY.mod \
            COMPAR.mod 
set_ic.$(OBJ_EXT) : set_ic.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            IC.mod \
            FLDVAR.mod \
            VISC_G.mod \
            INDICES.mod \
            SCALES.mod \
            ENERGY.mod \
            SCALARS.mod \
            COMPAR.mod \
            RUN.mod \
            SENDRECV.mod \
            sc_p_g1.inc                                                  \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            sc_p_g2.inc                                                 
set_increments3.$(OBJ_EXT) : set_increments3.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            FUNITS.mod \
            function.inc                                                 \
            function3.inc                                               
set_increments.$(OBJ_EXT) : set_increments.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            FUNITS.mod \
            function.inc                                                
set_index1a3.$(OBJ_EXT) : set_index1a3.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            BOUNDFUNIJK3.mod \
            function.inc                                                
set_index1a.$(OBJ_EXT) : set_index1a.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            BOUNDFUNIJK.mod \
            function.inc                                                
set_index1.$(OBJ_EXT) : set_index1.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            CONSTANT.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
set_l_scale.$(OBJ_EXT) : set_l_scale.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            VISC_G.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod 
set_max2.$(OBJ_EXT) : set_max2.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod 
set_mw_mix_g.$(OBJ_EXT) : set_mw_mix_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
set_outflow.$(OBJ_EXT) : set_outflow.f \
            PARAM.mod \
            PARAM1.mod \
            BC.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            SCALARS.mod \
            RUN.mod \
            COMPAR.mod \
            MFLUX.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
set_ro_g.$(OBJ_EXT) : set_ro_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
set_wall_bc.$(OBJ_EXT) : set_wall_bc.f \
            PARAM.mod \
            PARAM1.mod \
            BC.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            FUNITS.mod \
            COMPAR.mod \
            function.inc                                                
shift_dxyz.$(OBJ_EXT) : shift_dxyz.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod 
solve_continuity.$(OBJ_EXT) : solve_continuity.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            INDICES.mod \
            RESIDUAL.mod \
            CONT.mod \
            LEQSOL.mod \
            AMBM.mod 
solve_energy_eq.$(OBJ_EXT) : solve_energy_eq.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            BC.mod \
            ENERGY.mod \
            RXNS.mod \
            AMBM.mod \
            TMP_ARRAY.mod \
            TMP_ARRAY1.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod \
            DES_THERMO.mod \
            MFLUX.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
solve_epp.$(OBJ_EXT) : solve_epp.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PSCOR.mod \
            RESIDUAL.mod \
            LEQSOL.mod \
            PHYSPROP.mod \
            AMBM.mod \
            TMP_ARRAY1.mod 
solve_granular_energy.$(OBJ_EXT) : solve_granular_energy.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            VISC_S.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            BC.mod \
            ENERGY.mod \
            RXNS.mod \
            AMBM.mod \
            TMP_ARRAY.mod \
            COMPAR.mod \
            MFLUX.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
solve_k_epsilon_eq.$(OBJ_EXT) : solve_k_epsilon_eq.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            BC.mod \
            ENERGY.mod \
            RXNS.mod \
            TURB.mod \
            USR.mod \
            AMBM.mod \
            TMP_ARRAY.mod \
            COMPAR.mod \
            MFLUX.mod \
            CUTCELL.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
solve_lin_eq.$(OBJ_EXT) : solve_lin_eq.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            RESIDUAL.mod \
            TOLERANC.mod \
            LEQSOL.mod 
solve_pp_g.$(OBJ_EXT) : solve_pp_g.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            PGCOR.mod \
            RESIDUAL.mod \
            LEQSOL.mod \
            RUN.mod \
            AMBM.mod \
            TMP_ARRAY1.mod 
solve_scalar_eq.$(OBJ_EXT) : solve_scalar_eq.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            BC.mod \
            ENERGY.mod \
            RXNS.mod \
            SCALARS.mod \
            AMBM.mod \
            TMP_ARRAY.mod \
            COMPAR.mod \
            MFLUX.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
solve_species_eq.$(OBJ_EXT) : solve_species_eq.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            BC.mod \
            ENERGY.mod \
            RXNS.mod \
            AMBM.mod \
            MATRIX.mod \
            CHISCHEME.mod \
            TMP_ARRAY.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            MFLUX.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
solve_vel_star.$(OBJ_EXT) : solve_vel_star.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            GHDTHEORY.mod \
            OUTPUT.mod \
            INDICES.mod \
            DRAG.mod \
            RESIDUAL.mod \
            UR_FACS.mod \
            PGCOR.mod \
            PSCOR.mod \
            LEQSOL.mod \
            AMBM.mod \
            TMP_ARRAY1.mod \
            TMP_ARRAY.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod \
            QMOM_KINETIC_EQUATION.mod 
source_granular_energy.$(OBJ_EXT) : source_granular_energy.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            RUN.mod \
            DRAG.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            VISC_G.mod \
            VISC_S.mod \
            TRACE.mod \
            TURB.mod \
            INDICES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            COMPAR.mod \
            RESIDUAL.mod \
            KINTHEORY.mod \
            s_pr1.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            s_pr2.inc                                                   
source_phi.$(OBJ_EXT) : source_phi.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            COMPAR.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
source_pp_g.$(OBJ_EXT) : source_pp_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RXNS.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            BC.mod \
            VSHEAR.mod \
            XSI_ARRAY.mod \
            COMPAR.mod \
            UR_FACS.mod \
            CONSTANT.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            function.inc                                                
source_rop_g.$(OBJ_EXT) : source_rop_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            RXNS.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            COMPAR.mod \
            function.inc                                                
source_rop_s.$(OBJ_EXT) : source_rop_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            FLDVAR.mod \
            RXNS.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PGCOR.mod \
            PSCOR.mod \
            COMPAR.mod \
            function.inc                                                
source_u_g.$(OBJ_EXT) : source_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            TURB.mod \
            MPI_UTILITY.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_u_s.$(OBJ_EXT) : source_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            KINTHEORY.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_v_g.$(OBJ_EXT) : source_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            VSHEAR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_v_s.$(OBJ_EXT) : source_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            VSHEAR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            KINTHEORY.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_w_g.$(OBJ_EXT) : source_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_w_s.$(OBJ_EXT) : source_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            KINTHEORY.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
tau_u_g.$(OBJ_EXT) : tau_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            COMPAR.mod \
            SENDRECV.mod \
            BC.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_u_s.$(OBJ_EXT) : tau_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            VSHEAR.mod \
            SENDRECV.mod \
            COMPAR.mod \
            BC.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_v_g.$(OBJ_EXT) : tau_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            SENDRECV.mod \
            COMPAR.mod \
            BC.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_v_s.$(OBJ_EXT) : tau_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            SENDRECV.mod \
            COMPAR.mod \
            BC.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_w_g.$(OBJ_EXT) : tau_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            SENDRECV.mod \
            COMPAR.mod \
            BC.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_w_s.$(OBJ_EXT) : tau_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            SENDRECV.mod \
            COMPAR.mod \
            BC.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
test_lin_eq.$(OBJ_EXT) : test_lin_eq.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            function.inc                                                
time_march.$(OBJ_EXT) : time_march.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PGCOR.mod \
            PSCOR.mod \
            CONT.mod \
            COEFF.mod \
            TAU_G.mod \
            TAU_S.mod \
            VISC_G.mod \
            VISC_S.mod \
            FUNITS.mod \
            VSHEAR.mod \
            SCALARS.mod \
            TOLERANC.mod \
            DRAG.mod \
            RXNS.mod \
            COMPAR.mod \
            TIME_CPU.mod \
            DISCRETELEMENT.mod \
            MCHEM.mod \
            LEQSOL.mod \
            CDIST.mod \
            MFIX_NETCDF.mod \
            MPI_UTILITY.mod \
            CUTCELL.mod \
            VTK.mod \
            DASHBOARD.mod \
            QMOM_KINETIC_EQUATION.mod 
transfer.$(OBJ_EXT) : transfer.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod 
transport_prop.$(OBJ_EXT) : transport_prop.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            TOLERANC.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod 
undef_2_0.$(OBJ_EXT) : undef_2_0.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod 
under_relax.$(OBJ_EXT) : under_relax.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            function.inc                                                
update_old.$(OBJ_EXT) : update_old.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            TRACE.mod \
            VISC_S.mod \
            SCALARS.mod 
usr0.$(OBJ_EXT) : usr0.f \
            USR.mod 
usr1.$(OBJ_EXT) : usr1.f \
            USR.mod 
usr2.$(OBJ_EXT) : usr2.f \
            USR.mod 
usr3.$(OBJ_EXT) : usr3.f \
            USR.mod 
usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f \
            USR.mod 
usr_write_out0.$(OBJ_EXT) : usr_write_out0.f 
usr_write_out1.$(OBJ_EXT) : usr_write_out1.f 
utilities.$(OBJ_EXT) : utilities.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            BC.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            CONSTANT.mod \
            RUN.mod \
            COMPAR.mod \
            TOLERANC.mod \
            MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vavg_u_g.$(OBJ_EXT) : vavg_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            MFLUX.mod \
            function.inc                                                
vavg_u_s.$(OBJ_EXT) : vavg_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vavg_v_g.$(OBJ_EXT) : vavg_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            MFLUX.mod \
            function.inc                                                
vavg_v_s.$(OBJ_EXT) : vavg_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vavg_w_g.$(OBJ_EXT) : vavg_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            MFLUX.mod \
            function.inc                                                
vavg_w_s.$(OBJ_EXT) : vavg_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            BC.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vf_gs_x.$(OBJ_EXT) : vf_gs_x.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            DRAG.mod \
            DISCRETELEMENT.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
vf_gs_y.$(OBJ_EXT) : vf_gs_y.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            DRAG.mod \
            DISCRETELEMENT.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
vf_gs_z.$(OBJ_EXT) : vf_gs_z.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            DRAG.mod \
            DISCRETELEMENT.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
vtc_scalar.$(OBJ_EXT) : vtc_scalar.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            KINTHEORY.mod \
            function.inc                                                
write_ab_m.$(OBJ_EXT) : write_ab_m.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            INDICES.mod \
            function.inc                                                
write_ab_m_var.$(OBJ_EXT) : write_ab_m_var.f \
            PARAM.mod \
            PARAM1.mod \
            MATRIX.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            INDICES.mod \
            function.inc                                                
write_error.$(OBJ_EXT) : write_error.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod 
write_header.$(OBJ_EXT) : write_header.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            FUNITS.mod \
            COMPAR.mod 
write_out0.$(OBJ_EXT) : write_out0.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            IC.mod \
            BC.mod \
            IS.mod \
            FLDVAR.mod \
            CONSTANT.mod \
            INDICES.mod \
            FUNITS.mod \
            TOLERANC.mod \
            SCALES.mod \
            SCALARS.mod \
            UR_FACS.mod \
            LEQSOL.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            function.inc                                                
write_out1.$(OBJ_EXT) : write_out1.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            RXNS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
write_out3.$(OBJ_EXT) : write_out3.f \
            FUNITS.mod \
            COMPAR.mod 
write_res0.$(OBJ_EXT) : write_res0.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            IC.mod \
            IS.mod \
            BC.mod \
            CONSTANT.mod \
            FUNITS.mod \
            OUTPUT.mod \
            SCALES.mod \
            SCALARS.mod \
            RXNS.mod \
            UR_FACS.mod \
            LEQSOL.mod \
            TOLERANC.mod \
            CDIST.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod 
write_res1.$(OBJ_EXT) : write_res1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            SCALARS.mod \
            RXNS.mod \
            FUNITS.mod \
            OUTPUT.mod \
            ENERGY.mod \
            CDIST.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            MFIX_NETCDF.mod 
write_spx0.$(OBJ_EXT) : write_spx0.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            FUNITS.mod \
            CDIST.mod \
            COMPAR.mod \
            MPI_UTILITY.mod 
write_spx1.$(OBJ_EXT) : write_spx1.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            RUN.mod \
            FUNITS.mod \
            SCALARS.mod \
            OUTPUT.mod \
            RXNS.mod \
            CDIST.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            MFIX_NETCDF.mod 
write_table.$(OBJ_EXT) : write_table.f \
            PARAM.mod \
            PARAM1.mod \
            FUNITS.mod 
write_usr0.$(OBJ_EXT) : write_usr0.f 
write_usr1.$(OBJ_EXT) : write_usr1.f 
xerbla.$(OBJ_EXT) : xerbla.f \
            COMPAR.mod 
zero_array.$(OBJ_EXT) : zero_array.f \
            PARAM.mod \
            PARAM1.mod 
zero_norm_vel.$(OBJ_EXT) : zero_norm_vel.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            INDICES.mod \
            IS.mod \
            COMPAR.mod \
            function.inc                                                
allocate_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/allocate_cut_cell_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            CUTCELL.mod \
            STL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/allocate_cut_cell_arrays.f 
allocate_dummy_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/allocate_dummy_cut_cell_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/allocate_dummy_cut_cell_arrays.f 
calc_vort_out.$(OBJ_EXT) : ./cartesian_grid/calc_vort_out.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            FLDVAR.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/calc_vort_out.f 
cartesian_grid_init_namelist.$(OBJ_EXT) : ./cartesian_grid/cartesian_grid_init_namelist.f \
            PARAM1.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            POLYGON.mod \
            VTK.mod \
            PROGRESS_BAR.mod \
            DASHBOARD.mod \
            STL.mod \
            cartesian_grid/cartesian_grid_namelist.inc                  
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cartesian_grid_init_namelist.f 
CG_set_bc0.$(OBJ_EXT) : ./cartesian_grid/CG_set_bc0.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PHYSPROP.mod \
            BC.mod \
            FLDVAR.mod \
            INDICES.mod \
            RUN.mod \
            FUNITS.mod \
            SCALES.mod \
            SCALARS.mod \
            BOUNDFUNIJK.mod \
            TOLERANC.mod \
            SENDRECV.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            sc_p_g1.inc                                                  \
            function.inc                                                 \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_set_bc0.f 
CG_set_outflow.$(OBJ_EXT) : ./cartesian_grid/CG_set_outflow.f \
            PARAM.mod \
            PARAM1.mod \
            BC.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            SCALARS.mod \
            RUN.mod \
            COMPAR.mod \
            MFLUX.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_set_outflow.f 
CG_source_u_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_u_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_u_g.f 
CG_source_u_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_u_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            KINTHEORY.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_u_s.f 
CG_source_v_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_v_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            VSHEAR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_v_g.f 
CG_source_v_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_v_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            VSHEAR.mod \
            COMPAR.mod \
            SENDRECV.mod \
            KINTHEORY.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_v_s.f 
CG_source_w_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_w_g.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_w_g.f 
CG_source_w_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_w_s.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_S.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_S.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            KINTHEORY.mod \
            GHDTHEORY.mod \
            DRAG.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            OUTPUT.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_w_s.f 
check_data_cartesian.$(OBJ_EXT) : ./cartesian_grid/check_data_cartesian.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            RUN.mod \
            PHYSPROP.mod \
            INDICES.mod \
            SCALARS.mod \
            FUNITS.mod \
            LEQSOL.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            BC.mod \
            DISCRETELEMENT.mod \
            CUTCELL.mod \
            QUADRIC.mod \
            VTK.mod \
            POLYGON.mod \
            DASHBOARD.mod \
            STL.mod \
            FLDVAR.mod \
            SCALES.mod \
            PARALLEL.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            SENDRECV.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/check_data_cartesian.f 
cut_cell_preprocessing.$(OBJ_EXT) : ./cartesian_grid/cut_cell_preprocessing.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            VTK.mod \
            FLDVAR.mod \
            POLYGON.mod \
            STL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cut_cell_preprocessing.f 
deallocate_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/deallocate_cut_cell_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/deallocate_cut_cell_arrays.f 
define_quadrics.$(OBJ_EXT) : ./cartesian_grid/define_quadrics.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            VTK.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/define_quadrics.f 
dmp_cartesian.$(OBJ_EXT) : ./cartesian_grid/dmp_cartesian.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/dmp_cartesian.f 
eval_usr_fct.$(OBJ_EXT) : ./cartesian_grid/eval_usr_fct.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            FLDVAR.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/eval_usr_fct.f 
get_alpha.$(OBJ_EXT) : ./cartesian_grid/get_alpha.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            BC.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_alpha.f 
get_connectivity.$(OBJ_EXT) : ./cartesian_grid/get_connectivity.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            POLYGON.mod \
            STL.mod \
            FLDVAR.mod \
            VTK.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_connectivity.f 
get_cut_cell_flags.$(OBJ_EXT) : ./cartesian_grid/get_cut_cell_flags.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            VTK.mod \
            POLYGON.mod \
            STL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            SCALARS.mod \
            FUNITS.mod \
            RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_cut_cell_flags.f 
get_cut_cell_volume_area.$(OBJ_EXT) : ./cartesian_grid/get_cut_cell_volume_area.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            POLYGON.mod \
            STL.mod \
            BC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_cut_cell_volume_area.f 
get_delh.$(OBJ_EXT) : ./cartesian_grid/get_delh.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_delh.f 
get_master.$(OBJ_EXT) : ./cartesian_grid/get_master.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            BC.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_master.f 
get_poly_data.$(OBJ_EXT) : ./cartesian_grid/get_poly_data.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            RXNS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PROGRESS_BAR.mod \
            POLYGON.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_poly_data.f 
get_stl_data.$(OBJ_EXT) : ./cartesian_grid/get_stl_data.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            RXNS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PROGRESS_BAR.mod \
            STL.mod \
            VTK.mod \
            QUADRIC.mod \
            CONSTANT.mod \
            BC.mod \
            CUTCELL.mod \
            PARALLEL.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            SENDRECV.mod \
            STL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_stl_data.f 
set_Odxyz.$(OBJ_EXT) : ./cartesian_grid/set_Odxyz.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            VTK.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/set_Odxyz.f 
update_dashboard.$(OBJ_EXT) : ./cartesian_grid/update_dashboard.f \
            COMPAR.mod \
            PARALLEL.mod \
            SENDRECV.mod \
            RUN.mod \
            LEQSOL.mod \
            TIME_CPU.mod \
            RESIDUAL.mod \
            DASHBOARD.mod \
            VTK.mod \
            CONSTANT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/update_dashboard.f 
vtk_out.$(OBJ_EXT) : ./cartesian_grid/vtk_out.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            CONSTANT.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            QUADRIC.mod \
            CUTCELL.mod \
            FLDVAR.mod \
            VISC_S.mod \
            PHYSPROP.mod \
            PGCOR.mod \
            VTK.mod \
            RXNS.mod \
            OUTPUT.mod \
            SCALARS.mod \
            PSCOR.mod \
            DISCRETELEMENT.mod \
            MFIX_PIC.mod \
            MPI_UTILITY.mod \
            POLYGON.mod \
            STL.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/vtk_out.f 
write_progress_bar.$(OBJ_EXT) : ./cartesian_grid/write_progress_bar.f \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            SCALARS.mod \
            FUNITS.mod \
            RXNS.mod \
            COMPAR.mod \
            MPI_UTILITY.mod \
            PROGRESS_BAR.mod \
            PARALLEL.mod \
            SENDRECV.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/write_progress_bar.f 
calc_jacobian.$(OBJ_EXT) : ./chem/calc_jacobian.f \
            PARAM1.mod \
            MCHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/calc_jacobian.f 
check_data_chem.$(OBJ_EXT) : ./chem/check_data_chem.f \
            PARAM1.mod \
            RUN.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/check_data_chem.f 
dgpadm.$(OBJ_EXT) : ./chem/dgpadm.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/dgpadm.f 
exponential.$(OBJ_EXT) : ./chem/exponential.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/exponential.f 
fex.$(OBJ_EXT) : ./chem/fex.f \
            RUN.mod \
            PHYSPROP.mod \
            TOLERANC.mod \
            USR.mod \
            MCHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/fex.f 
g_derivs.$(OBJ_EXT) : ./chem/g_derivs.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/g_derivs.f 
jac.$(OBJ_EXT) : ./chem/jac.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/jac.f 
mchem_init.$(OBJ_EXT) : ./chem/mchem_init.f \
            PARAM1.mod \
            RUN.mod \
            PHYSPROP.mod \
            MCHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/mchem_init.f 
mchem_odepack_init.$(OBJ_EXT) : ./chem/mchem_odepack_init.f \
            MCHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/mchem_odepack_init.f 
mchem_time_march.$(OBJ_EXT) : ./chem/mchem_time_march.f \
            RUN.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RXNS.mod \
            MPI_UTILITY.mod \
            TOLERANC.mod \
            MCHEM.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/mchem_time_march.f 
misat_table_init.$(OBJ_EXT) : ./chem/misat_table_init.f \
            MCHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/misat_table_init.f 
react.$(OBJ_EXT) : ./chem/react.f \
            PARAM1.mod \
            TOLERANC.mod \
            FLDVAR.mod \
            PHYSPROP.mod \
            RXNS.mod \
            RUN.mod \
            MCHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/react.f 
usrfg.$(OBJ_EXT) : ./chem/usrfg.f \
            PARAM1.mod \
            RUN.mod \
            MCHEM.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/usrfg.f 
add_part_to_link_list.$(OBJ_EXT) : ./cohesion/add_part_to_link_list.f \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/add_part_to_link_list.f 
calc_app_coh_force.$(OBJ_EXT) : ./cohesion/calc_app_coh_force.f \
            DISCRETELEMENT.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_app_coh_force.f 
calc_cap_coh_force.$(OBJ_EXT) : ./cohesion/calc_cap_coh_force.f \
            DISCRETELEMENT.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_cap_coh_force.f 
calc_cohesive_forces.$(OBJ_EXT) : ./cohesion/calc_cohesive_forces.f \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_cohesive_forces.f 
calc_esc_coh_force.$(OBJ_EXT) : ./cohesion/calc_esc_coh_force.f \
            DISCRETELEMENT.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_esc_coh_force.f 
calc_square_well.$(OBJ_EXT) : ./cohesion/calc_square_well.f \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_square_well.f 
calc_van_der_waals.$(OBJ_EXT) : ./cohesion/calc_van_der_waals.f \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_van_der_waals.f 
check_link.$(OBJ_EXT) : ./cohesion/check_link.f \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_link.f 
check_sw_wall_interaction.$(OBJ_EXT) : ./cohesion/check_sw_wall_interaction.f \
            PARAM1.mod \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_sw_wall_interaction.f 
check_vdw_wall_interaction.$(OBJ_EXT) : ./cohesion/check_vdw_wall_interaction.f \
            PARAM1.mod \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_vdw_wall_interaction.f 
initialize_cohesion_parameters.$(OBJ_EXT) : ./cohesion/initialize_cohesion_parameters.f \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/initialize_cohesion_parameters.f 
initialize_coh_int_search.$(OBJ_EXT) : ./cohesion/initialize_coh_int_search.f \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/initialize_coh_int_search.f 
linked_interaction_eval.$(OBJ_EXT) : ./cohesion/linked_interaction_eval.f \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/linked_interaction_eval.f 
remove_part_from_link_list.$(OBJ_EXT) : ./cohesion/remove_part_from_link_list.f \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/remove_part_from_link_list.f 
unlinked_interaction_eval.$(OBJ_EXT) : ./cohesion/unlinked_interaction_eval.f \
            DISCRETELEMENT.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/unlinked_interaction_eval.f 
update_search_grids.$(OBJ_EXT) : ./cohesion/update_search_grids.f \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/update_search_grids.f 
calc_attrition_des.$(OBJ_EXT) : ./des/calc_attrition_des.f \
            RUN.mod \
            PARAM1.mod \
            DISCRETELEMENT.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            CONSTANT.mod
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_attrition_des.f
calc_force_des.$(OBJ_EXT) : ./des/calc_force_des.f \
            RUN.mod \
            PARAM1.mod \
            DISCRETELEMENT.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            CONSTANT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_force_des.f 
calc_rrate_des.$(OBJ_EXT) : ./des/calc_rrate_des.f \
            DISCRETELEMENT.mod \
            INTERPOLATION.mod \
            PARAM1.mod \
            RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_rrate_des.f 
calc_thermo_des.$(OBJ_EXT) : ./des/calc_thermo_des.f \
            DES_THERMO.mod \
            DISCRETELEMENT.mod \
            FLDVAR.mod \
            INTERPOLATION.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_thermo_des.f 
cell_near_wall.$(OBJ_EXT) : ./des/cell_near_wall.f \
            DISCRETELEMENT.mod \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            MATRIX.mod \
            INDICES.mod \
            PHYSPROP.mod \
            DRAG.mod \
            CONSTANT.mod \
            COMPAR.mod \
            SENDRECV.mod \
            TOLERANC.mod \
            IS.mod \
            BC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cell_near_wall.f 
cfassign.$(OBJ_EXT) : ./des/cfassign.f \
            PARAM1.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            CONSTANT.mod \
            COMPAR.mod \
            PARALLEL.mod \
            SENDRECV.mod \
            DISCRETELEMENT.mod \
            MFIX_PIC.mod \
            PARAM.mod \
            FLDVAR.mod \
            RUN.mod \
            INDICES.mod \
            BC.mod \
            CUTCELL.mod \
            b_force1.inc                                                 \
            b_force2.inc                                                 \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfassign.f 
cffctowall.$(OBJ_EXT) : ./des/cffctowall.f \
            PARAM1.mod \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffctowall.f 
cffctow.$(OBJ_EXT) : ./des/cffctow.f \
            PARAM1.mod \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffctow.f 
cfnewvalues.$(OBJ_EXT) : ./des/cfnewvalues.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            DISCRETELEMENT.mod \
            DES_BC.mod \
            MPI_UTILITY.mod \
            MFIX_PIC.mod \
            MPPIC_WALLBC.mod \
            RANDOMNO.mod \
            CUTCELL.mod \
            GEOMETRY.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfnewvalues.f 
cfrelvel.$(OBJ_EXT) : ./des/cfrelvel.f \
            DISCRETELEMENT.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfrelvel.f 
cfslide.$(OBJ_EXT) : ./des/cfslide.f \
            PARAM1.mod \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslide.f 
cfslidewall.$(OBJ_EXT) : ./des/cfslidewall.f \
            DISCRETELEMENT.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslidewall.f 
cfupdateold.$(OBJ_EXT) : ./des/cfupdateold.f \
            DISCRETELEMENT.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfupdateold.f 
cfwallcontact.$(OBJ_EXT) : ./des/cfwallcontact.f \
            PARAM1.mod \
            CONSTANT.mod \
            PARALLEL.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod \
            DES_BC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallcontact.f 
cfwallposvel.$(OBJ_EXT) : ./des/cfwallposvel.f \
            DISCRETELEMENT.mod \
            DES_BC.mod \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            MATRIX.mod \
            INDICES.mod \
            PHYSPROP.mod \
            DRAG.mod \
            CONSTANT.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallposvel.f 
check_des_bc.$(OBJ_EXT) : ./des/check_des_bc.f \
            COMPAR.mod \
            CONSTANT.mod \
            DES_BC.mod \
            DISCRETELEMENT.mod \
            FUNITS.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RUN.mod \
            MFIX_PIC.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_bc.f 
check_des_data.$(OBJ_EXT) : ./des/check_des_data.f \
            PARAM1.mod \
            GEOMETRY.mod \
            FUNITS.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod \
            RUN.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            DESGRID.mod \
            INDICES.mod \
            FLDVAR.mod \
            TOLERANC.mod \
            MPI_UTILITY.mod \
            OUTPUT.mod \
            MFIX_PIC.mod \
            CUTCELL.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_data.f 
check_des_ic.$(OBJ_EXT) : ./des/check_des_ic.f \
            DES_IC.mod \
            DISCRETELEMENT.mod \
            PARAM.mod \
            PARAM1.mod \
            DES_THERMO.mod \
            DES_RXNS.mod \
            COMPAR.mod \
            CONSTANT.mod \
            FUNITS.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_ic.f 
check_des_rxns.$(OBJ_EXT) : ./des/check_des_rxns.f \
            COMPAR.mod \
            DES_RXNS.mod \
            DISCRETELEMENT.mod \
            FUNITS.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_rxns.f 
check_des_thermo.$(OBJ_EXT) : ./des/check_des_thermo.f \
            COMPAR.mod \
            DES_THERMO.mod \
            DISCRETELEMENT.mod \
            FUNITS.mod \
            INTERPOLATION.mod \
            PHYSPROP.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_thermo.f 
des_allocate_arrays.$(OBJ_EXT) : ./des/des_allocate_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            DISCRETELEMENT.mod \
            INDICES.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            DES_BC.mod \
            FUNITS.mod \
            DESGRID.mod \
            DESMPI.mod \
            MFIX_PIC.mod \
            DES_THERMO.mod \
            DES_RXNS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_allocate_arrays.f 
des_check_particle.$(OBJ_EXT) : ./des/des_check_particle.f \
            COMPAR.mod \
            CONSTANT.mod \
            DES_BC.mod \
            DISCRETELEMENT.mod \
            FUNITS.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_check_particle.f 
des_functions.$(OBJ_EXT) : ./des/des_functions.f \
            PARAM.mod \
            PARAM1.mod \
            DISCRETELEMENT.mod \
            COMPAR.mod \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_functions.f 
des_granular_temperature.$(OBJ_EXT) : ./des/des_granular_temperature.f \
            DISCRETELEMENT.mod \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            PHYSPROP.mod \
            DES_BC.mod \
            FLDVAR.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_granular_temperature.f 
des_init_arrays.$(OBJ_EXT) : ./des/des_init_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            DISCRETELEMENT.mod \
            INDICES.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            DES_BC.mod \
            RUN.mod \
            DESGRID.mod \
            DESMPI.mod \
            DES_THERMO.mod \
            DES_RXNS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_arrays.f 
des_init_bc.$(OBJ_EXT) : ./des/des_init_bc.f \
            COMPAR.mod \
            CONSTANT.mod \
            DES_BC.mod \
            DISCRETELEMENT.mod \
            FUNITS.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_bc.f 
des_init_namelist.$(OBJ_EXT) : ./des/des_init_namelist.f \
            PARAM1.mod \
            DISCRETELEMENT.mod \
            MFIX_PIC.mod \
            DES_BC.mod \
            DES_IC.mod \
            DES_THERMO.mod \
            DES_RXNS.mod \
            des/desnamelist.inc                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_namelist.f 
des_mass_inlet.$(OBJ_EXT) : ./des/des_mass_inlet.f \
            COMPAR.mod \
            CONSTANT.mod \
            DES_BC.mod \
            DISCRETELEMENT.mod \
            FUNITS.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            DESGRID.mod \
            MPI_UTILITY.mod \
            DES_THERMO.mod \
            DES_RXNS.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_mass_inlet.f 
des_physical_prop.$(OBJ_EXT) : ./des/des_physical_prop.f \
            DES_RXNS.mod \
            DES_THERMO.mod \
            DISCRETELEMENT.mod \
            FUNITS.mod \
            PARAM.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_physical_prop.f 
des_reaction_model.$(OBJ_EXT) : ./des/des_reaction_model.f \
            CONSTANT.mod \
            DES_RXNS.mod \
            DES_THERMO.mod \
            DISCRETELEMENT.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_reaction_model.f 
des_rrates.$(OBJ_EXT) : ./des/des_rrates.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RXNS.mod \
            ENERGY.mod \
            GEOMETRY.mod \
            RUN.mod \
            INDICES.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            FUNITS.mod \
            TOLERANC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            USR.mod \
            DES_THERMO.mod \
            DES_RXNS.mod \
            DISCRETELEMENT.mod \
            INTERPOLATION.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_rrates.f 
des_set_ic.$(OBJ_EXT) : ./des/des_set_ic.f \
            COMPAR.mod \
            DES_THERMO.mod \
            DISCRETELEMENT.mod \
            DES_IC.mod \
            DES_RXNS.mod \
            FUNITS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_set_ic.f 
des_thermo_cond.$(OBJ_EXT) : ./des/des_thermo_cond.f \
            CONSTANT.mod \
            DES_THERMO.mod \
            DISCRETELEMENT.mod \
            FLDVAR.mod \
            FUNITS.mod \
            PARAM1.mod \
            PHYSPROP.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_cond.f 
des_thermo_conv.$(OBJ_EXT) : ./des/des_thermo_conv.f \
            CONSTANT.mod \
            DES_THERMO.mod \
            DISCRETELEMENT.mod \
            FLDVAR.mod \
            INTERPOLATION.mod \
            PARAM1.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_conv.f 
des_thermo_newvalues.$(OBJ_EXT) : ./des/des_thermo_newvalues.f \
            DES_THERMO.mod \
            DES_RXNS.mod \
            DISCRETELEMENT.mod \
            PARAM1.mod \
            PHYSPROP.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_newvalues.f 
des_thermo_rad.$(OBJ_EXT) : ./des/des_thermo_rad.f \
            CONSTANT.mod \
            DES_THERMO.mod \
            DISCRETELEMENT.mod \
            PARAM1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_rad.f 
des_time_march.$(OBJ_EXT) : ./des/des_time_march.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PGCOR.mod \
            PSCOR.mod \
            CONT.mod \
            COEFF.mod \
            TAU_G.mod \
            TAU_S.mod \
            VISC_G.mod \
            VISC_S.mod \
            FUNITS.mod \
            VSHEAR.mod \
            SCALARS.mod \
            DRAG.mod \
            RXNS.mod \
            COMPAR.mod \
            TIME_CPU.mod \
            DISCRETELEMENT.mod \
            CONSTANT.mod \
            SENDRECV.mod \
            DES_BC.mod \
            CUTCELL.mod \
            MPPIC_WALLBC.mod \
            MFIX_PIC.mod \
            DES_THERMO.mod \
            DES_RXNS.mod \
            INTERPOLATION.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_time_march.f 
des_wallbc_preprocessing.$(OBJ_EXT) : ./des/des_wallbc_preprocessing.f \
            PARAM1.mod \
            FUNITS.mod \
            RUN.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod \
            CUTCELL.mod \
            INDICES.mod \
            PHYSPROP.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            BC.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            MPI_UTILITY.mod \
            SENDRECV.mod \
            MFIX_PIC.mod \
            DISCRETELEMENT.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_wallbc_preprocessing.f 
drag_fgs.$(OBJ_EXT) : ./des/drag_fgs.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            DISCRETELEMENT.mod \
            CUTCELL.mod \
            CONSTANT.mod \
            DRAG.mod \
            INTERPOLATION.mod \
            DESMPI.mod \
            MFIX_PIC.mod \
            UR_FACS.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/drag_fgs.f 
gas_drag.$(OBJ_EXT) : ./des/gas_drag.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            DISCRETELEMENT.mod \
            DRAG.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/gas_drag.f 
generate_particle_config.$(OBJ_EXT) : ./des/generate_particle_config.f \
            PARAM1.mod \
            GEOMETRY.mod \
            FUNITS.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod \
            RUN.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            DESMPI.mod \
            CDIST.mod \
            MPI_UTILITY.mod \
            MFIX_PIC.mod \
            PARAM.mod \
            FLDVAR.mod \
            INDICES.mod \
            RANDOMNO.mod \
            SENDRECV.mod \
            INTERPOLATION.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/generate_particle_config.f 
grid_based_neighbor_search.$(OBJ_EXT) : ./des/grid_based_neighbor_search.f \
            PARAM1.mod \
            DISCRETELEMENT.mod \
            GEOMETRY.mod \
            DES_BC.mod \
            DES_THERMO.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/grid_based_neighbor_search.f 
make_arrays_des.$(OBJ_EXT) : ./des/make_arrays_des.f \
            PARAM1.mod \
            FUNITS.mod \
            RUN.mod \
            COMPAR.mod \
            DISCRETELEMENT.mod \
            CUTCELL.mod \
            DESMPI.mod \
            MPI_UTILITY.mod \
            GEOMETRY.mod \
            DES_IC.mod \
            DES_RXNS.mod \
            CDIST.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/make_arrays_des.f 
mppic_routines.$(OBJ_EXT) : ./des/mppic_routines.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            INDICES.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            DISCRETELEMENT.mod \
            CONSTANT.mod \
            CUTCELL.mod \
            INTERPOLATION.mod \
            MFIX_PIC.mod \
            DRAG.mod \
            DESMPI.mod \
            TOLERANC.mod \
            QUADRIC.mod \
            VTK.mod \
            MATRIX.mod \
            SCALES.mod \
            VISC_S.mod \
            RXNS.mod \
            IS.mod \
            TAU_S.mod \
            OUTPUT.mod \
            DES_BC.mod \
            MPI_UTILITY.mod \
            MPPIC_WALLBC.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mppic_routines.f 
neighbour.$(OBJ_EXT) : ./des/neighbour.f \
            PARAM1.mod \
            DISCRETELEMENT.mod \
            DESGRID.mod \
            DES_THERMO.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/neighbour.f 
nsquare.$(OBJ_EXT) : ./des/nsquare.f \
            PARAM1.mod \
            DISCRETELEMENT.mod \
            GEOMETRY.mod \
            DES_BC.mod \
            DES_THERMO.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/nsquare.f 
octree.$(OBJ_EXT) : ./des/octree.f \
            RUN.mod \
            PARAM1.mod \
            CONSTANT.mod \
            DISCRETELEMENT.mod \
            PARAM.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            MATRIX.mod \
            INDICES.mod \
            PHYSPROP.mod \
            DRAG.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/octree.f 
particles_in_cell.$(OBJ_EXT) : ./des/particles_in_cell.f \
            PARAM.mod \
            PARAM1.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            COMPAR.mod \
            PARALLEL.mod \
            SENDRECV.mod \
            DISCRETELEMENT.mod \
            DESGRID.mod \
            DESMPI.mod \
            CUTCELL.mod \
            MFIX_PIC.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/particles_in_cell.f 
quadtree.$(OBJ_EXT) : ./des/quadtree.f \
            RUN.mod \
            PARAM1.mod \
            CONSTANT.mod \
            DISCRETELEMENT.mod \
            PARAM.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            MATRIX.mod \
            INDICES.mod \
            PHYSPROP.mod \
            DRAG.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/quadtree.f 
read_des_restart.$(OBJ_EXT) : ./des/read_des_restart.f \
            PARAM1.mod \
            RUN.mod \
            DISCRETELEMENT.mod \
            DES_BC.mod \
            COMPAR.mod \
            DESMPI.mod \
            MACHINE.mod \
            CDIST.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/read_des_restart.f 
thermo_nbr.$(OBJ_EXT) : ./des/thermo_nbr.f \
            DES_THERMO.mod \
            DISCRETELEMENT.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/thermo_nbr.f 
walledgecontact.$(OBJ_EXT) : ./des/walledgecontact.f \
            DISCRETELEMENT.mod \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            MATRIX.mod \
            INDICES.mod \
            PHYSPROP.mod \
            DRAG.mod \
            CONSTANT.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/walledgecontact.f 
wallfacecontact.$(OBJ_EXT) : ./des/wallfacecontact.f \
            DISCRETELEMENT.mod \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            MATRIX.mod \
            INDICES.mod \
            PHYSPROP.mod \
            DRAG.mod \
            CONSTANT.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/wallfacecontact.f 
wallnodecontact.$(OBJ_EXT) : ./des/wallnodecontact.f \
            DISCRETELEMENT.mod \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            RUN.mod \
            GEOMETRY.mod \
            MATRIX.mod \
            INDICES.mod \
            PHYSPROP.mod \
            DRAG.mod \
            CONSTANT.mod \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/wallnodecontact.f 
write_des_data.$(OBJ_EXT) : ./des/write_des_data.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            DISCRETELEMENT.mod \
            RUN.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            SENDRECV.mod \
            DES_BC.mod \
            MPI_UTILITY.mod \
            COMPAR.mod \
            DESMPI.mod \
            CDIST.mod \
            DES_THERMO.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/write_des_data.f 
write_des_restart.$(OBJ_EXT) : ./des/write_des_restart.f \
            PARAM1.mod \
            DISCRETELEMENT.mod \
            RUN.mod \
            DES_BC.mod \
            COMPAR.mod \
            DESMPI.mod \
            MACHINE.mod \
            CDIST.mod \
            MPI_UTILITY.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/write_des_restart.f 
gaussj.$(OBJ_EXT) : ./dqmom/gaussj.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/gaussj.f 
odeint.$(OBJ_EXT) : ./dqmom/odeint.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/odeint.f 
rkck.$(OBJ_EXT) : ./dqmom/rkck.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/rkck.f 
rkqs.$(OBJ_EXT) : ./dqmom/rkqs.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/rkqs.f 
source_population_eq.$(OBJ_EXT) : ./dqmom/source_population_eq.f \
            PHYSPROP.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            SCALARS.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/source_population_eq.f 
usr_dqmom.$(OBJ_EXT) : ./dqmom/usr_dqmom.f \
            PARAM.mod \
            PARAM1.mod \
            RUN.mod \
            PHYSPROP.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            OUTPUT.mod \
            INDICES.mod \
            RXNS.mod \
            CONSTANT.mod \
            AMBM.mod \
            COMPAR.mod \
            SCALARS.mod \
            USR.mod \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/usr_dqmom.f 
adjust_eps_ghd.$(OBJ_EXT) : ./GhdTheory/adjust_eps_ghd.f \
            PARAM.mod \
            PARAM1.mod \
            TOLERANC.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            COMPAR.mod \
            SENDRECV.mod \
            GHDTHEORY.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/adjust_eps_ghd.f 
bulk_viscosity.$(OBJ_EXT) : ./GhdTheory/bulk_viscosity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/bulk_viscosity.f 
calc_d_ghd.$(OBJ_EXT) : ./GhdTheory/calc_d_ghd.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            INDICES.mod \
            PHYSPROP.mod \
            RUN.mod \
            SCALES.mod \
            COMPAR.mod \
            SENDRECV.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_d_ghd.f 
calc_external_forces.$(OBJ_EXT) : ./GhdTheory/calc_external_forces.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            GHDTHEORY.mod \
            PHYSPROP.mod \
            RUN.mod \
            CONSTANT.mod \
            DRAG.mod \
            BC.mod \
            SCALES.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            b_force1.inc                                                 \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_external_forces.f 
calc_nflux.$(OBJ_EXT) : ./GhdTheory/calc_nflux.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            CONSTANT.mod \
            INDICES.mod \
            MFLUX.mod \
            COMPAR.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_nflux.f 
chi_ij_GHD.$(OBJ_EXT) : ./GhdTheory/chi_ij_GHD.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/chi_ij_GHD.f 
cooling_rate.$(OBJ_EXT) : ./GhdTheory/cooling_rate.f \
            COMPAR.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/cooling_rate.f 
cooling_rate_tc.$(OBJ_EXT) : ./GhdTheory/cooling_rate_tc.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/cooling_rate_tc.f 
dufour_coeff.$(OBJ_EXT) : ./GhdTheory/dufour_coeff.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/dufour_coeff.f 
ghd.$(OBJ_EXT) : ./GhdTheory/ghd.f \
            DRAG.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ghd.f 
ghdmassflux.$(OBJ_EXT) : ./GhdTheory/ghdmassflux.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            VISC_S.mod \
            GHDTHEORY.mod \
            PHYSPROP.mod \
            RUN.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            DRAG.mod \
            IS.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ghdmassflux.f 
mass_mobility.$(OBJ_EXT) : ./GhdTheory/mass_mobility.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/mass_mobility.f 
ordinary_diff.$(OBJ_EXT) : ./GhdTheory/ordinary_diff.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ordinary_diff.f 
partial_elim_ghd.$(OBJ_EXT) : ./GhdTheory/partial_elim_ghd.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            GEOMETRY.mod \
            MATRIX.mod \
            PHYSPROP.mod \
            INDICES.mod \
            RUN.mod \
            COMPAR.mod \
            DRAG.mod \
            FLDVAR.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/partial_elim_ghd.f 
pressure.$(OBJ_EXT) : ./GhdTheory/pressure.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/pressure.f 
shear_viscosity.$(OBJ_EXT) : ./GhdTheory/shear_viscosity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/shear_viscosity.f 
source_ghd_granular_energy.$(OBJ_EXT) : ./GhdTheory/source_ghd_granular_energy.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            PHYSPROP.mod \
            RUN.mod \
            DRAG.mod \
            GEOMETRY.mod \
            FLDVAR.mod \
            VISC_S.mod \
            GHDTHEORY.mod \
            TRACE.mod \
            INDICES.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            COMPAR.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            b_force1.inc                                                 \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/source_ghd_granular_energy.f 
thermal_conductivity.$(OBJ_EXT) : ./GhdTheory/thermal_conductivity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/thermal_conductivity.f 
thermal_diffusivity.$(OBJ_EXT) : ./GhdTheory/thermal_diffusivity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/thermal_diffusivity.f 
thermal_mobility.$(OBJ_EXT) : ./GhdTheory/thermal_mobility.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/thermal_mobility.f 
transport_coeff_ghd.$(OBJ_EXT) : ./GhdTheory/transport_coeff_ghd.f \
            PARAM.mod \
            PARAM1.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            FLDVAR.mod \
            INDICES.mod \
            VISC_S.mod \
            GHDTHEORY.mod \
            PHYSPROP.mod \
            RUN.mod \
            CONSTANT.mod \
            TOLERANC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/transport_coeff_ghd.f 
qmomk_allocate_arrays.$(OBJ_EXT) : ./qmomk/qmomk_allocate_arrays.f \
            PARAM.mod \
            PARAM1.mod \
            INDICES.mod \
            GEOMETRY.mod \
            PHYSPROP.mod \
            QMOM_KINETIC_EQUATION.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_allocate_arrays.f 
qmomk_gas_drag.$(OBJ_EXT) : ./qmomk/qmomk_gas_drag.f \
            PARAM.mod \
            PARAM1.mod \
            PARALLEL.mod \
            MATRIX.mod \
            SCALES.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            VISC_G.mod \
            RXNS.mod \
            RUN.mod \
            TOLERANC.mod \
            GEOMETRY.mod \
            INDICES.mod \
            IS.mod \
            TAU_G.mod \
            BC.mod \
            COMPAR.mod \
            SENDRECV.mod \
            DISCRETELEMENT.mod \
            QMOM_KINETIC_EQUATION.mod \
            DRAG.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_gas_drag.f 
qmomk_init_bc.$(OBJ_EXT) : ./qmomk/qmomk_init_bc.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            BC.mod \
            QMOM_KINETIC_EQUATION.mod \
            QMOMK_QUADRATURE.mod \
            QMOMK_BC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_init_bc.f 
qmomk_initial_conditions.$(OBJ_EXT) : ./qmomk/qmomk_initial_conditions.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            QMOM_KINETIC_EQUATION.mod \
            QMOMK_QUADRATURE.mod \
            QMOMK_PARAMETERS.mod \
            IC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_initial_conditions.f 
qmomk_init_namelist.$(OBJ_EXT) : ./qmomk/qmomk_init_namelist.f \
            PARAM1.mod \
            QMOM_KINETIC_EQUATION.mod \
            qmomk/qmomknamelist.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_init_namelist.f 
qmomk_make_arrays.$(OBJ_EXT) : ./qmomk/qmomk_make_arrays.f \
            PARAM1.mod \
            GEOMETRY.mod \
            FUNITS.mod \
            COMPAR.mod \
            QMOM_KINETIC_EQUATION.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_make_arrays.f 
qmomk_read_restart.$(OBJ_EXT) : ./qmomk/qmomk_read_restart.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            FLDVAR.mod \
            CONT.mod \
            GEOMETRY.mod \
            INDICES.mod \
            RUN.mod \
            COMPAR.mod \
            PHYSPROP.mod \
            QMOM_KINETIC_EQUATION.mod \
            QMOMK_QUADRATURE.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_read_restart.f 
qmomk_set_bc.$(OBJ_EXT) : ./qmomk/qmomk_set_bc.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            COMPAR.mod \
            INDICES.mod \
            BC.mod \
            QMOM_KINETIC_EQUATION.mod \
            QMOMK_QUADRATURE.mod \
            QMOMK_BC.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_set_bc.f 
qmomk_time_march.$(OBJ_EXT) : ./qmomk/qmomk_time_march.f \
            PARAM.mod \
            PARAM1.mod \
            CONSTANT.mod \
            RUN.mod \
            OUTPUT.mod \
            PHYSPROP.mod \
            FLDVAR.mod \
            GEOMETRY.mod \
            CONT.mod \
            COEFF.mod \
            TAU_G.mod \
            TAU_S.mod \
            VISC_G.mod \
            VISC_S.mod \
            FUNITS.mod \
            VSHEAR.mod \
            SCALARS.mod \
            DRAG.mod \
            RXNS.mod \
            COMPAR.mod \
            TIME_CPU.mod \
            IS.mod \
            INDICES.mod \
            MATRIX.mod \
            SENDRECV.mod \
            QMOM_KINETIC_EQUATION.mod \
            QMOMK_FLUXES.mod \
            QMOMK_QUADRATURE.mod \
            QMOMK_COLLISION.mod \
            QMOMK_PARAMETERS.mod \
            UR_FACS.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_time_march.f 
qmomk_write_restart.$(OBJ_EXT) : ./qmomk/qmomk_write_restart.f \
            PARAM1.mod \
            QMOM_KINETIC_EQUATION.mod \
            RUN.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_write_restart.f 
get_values.$(OBJ_EXT) : ./thermochemical/get_values.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./thermochemical/get_values.f 
readTherm.$(OBJ_EXT) : ./thermochemical/readTherm.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./thermochemical/readTherm.f 
