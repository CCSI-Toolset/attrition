.$(FORTRAN_EXT).$(OBJ_EXT):
	$(FORTRAN_CMD) $(FORT_FLAGS) $<
  
mfix.exe : \
    ambm.mod \
    bc.mod \
    boundfunijk3.mod \
    boundfunijk.mod \
    cdist.mod \
    check.mod \
    chischeme.mod \
    coeff.mod \
    constant.mod \
    cont.mod \
    corner.mod \
    drag.mod \
    energy.mod \
    fldvar.mod \
    function.mod \
    funits.mod \
    geometry.mod \
    ic.mod \
    indices.mod \
    is.mod \
    kintheory.mod \
    leqsol.mod \
    machine.mod \
    matrix.mod \
    mfix_netcdf.mod \
    mflux.mod \
    output.mod \
    parallel.mod \
    param1.mod \
    param.mod \
    parse.mod \
    pgcor.mod \
    physprop.mod \
    pscor.mod \
    residual.mod \
    run.mod \
    rxns.mod \
    scalars.mod \
    scales.mod \
    tau_g.mod \
    tau_s.mod \
    time_cpu.mod \
    tmp_array1.mod \
    tmp_array.mod \
    toleranc.mod \
    trace.mod \
    turb.mod \
    ur_facs.mod \
    usr.mod \
    visc_g.mod \
    visc_s.mod \
    vshear.mod \
    xsi_array.mod \
    cutcell.mod \
    dashboard.mod \
    polygon.mod \
    progress_bar.mod \
    quadric.mod \
    stl.mod \
    vtk.mod \
    mchem.mod \
    des_bc.mod \
    desgrid.mod \
    des_ic.mod \
    desmpi.mod \
    desmpi_wrapper.mod \
    des_rxns.mod \
    des_thermo.mod \
    discretelement.mod \
    interpolation.mod \
    mfix_pic.mod \
    mppic_wallbc.mod \
    randomno.mod \
    sendrecvnode.mod \
    compar.mod \
    dbg_util.mod \
    debug.mod \
    gridmap.mod \
    mpi.mod \
    mpi_utility.mod \
    parallel_mpi.mod \
    sendrecv3.mod \
    sendrecv.mod \
    ghdtheory.mod \
    qmomk_bc.mod \
    qmomk_collision.mod \
    qmomk_fluxes.mod \
    qmom_kinetic_equation.mod \
    qmomk_parameters.mod \
    qmomk_quadrature.mod \
    qmomk_tools.mod \
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
ambm.mod : ambm_mod.f \
            param.mod \
            param1.mod \
            compar.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ambm_mod.f 
bc.mod : bc_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) bc_mod.f 
boundfunijk3.mod : boundfunijk3_mod.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            function3.inc                                               
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk3_mod.f 
boundfunijk.mod : boundfunijk_mod.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) boundfunijk_mod.f 
cdist.mod : cdist_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) cdist_mod.f 
check.mod : check_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) check_mod.f 
chischeme.mod : chischeme_mod.f \
            param.mod \
            param1.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) chischeme_mod.f 
coeff.mod : coeff_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) coeff_mod.f 
constant.mod : constant_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) constant_mod.f 
cont.mod : cont_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) cont_mod.f 
corner.mod : corner_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) corner_mod.f 
drag.mod : drag_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) drag_mod.f 
energy.mod : energy_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) energy_mod.f 
fldvar.mod : fldvar_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) fldvar_mod.f 
function.mod : function_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) function_mod.f 
funits.mod : funits_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) funits_mod.f 
geometry.mod : geometry_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) geometry_mod.f 
ic.mod : ic_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ic_mod.f 
indices.mod : indices_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) indices_mod.f 
is.mod : is_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) is_mod.f 
kintheory.mod : kintheory_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) kintheory_mod.f 
leqsol.mod : leqsol_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) leqsol_mod.f 
machine.mod : machine_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) machine_mod.f 
matrix.mod : matrix_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) matrix_mod.f 
mfix_netcdf.mod : mfix_netcdf_mod.f \
            MFIX_netcdf_constants.fi                                     \
            MFIX_netcdf_overloads.fi                                     \
            MFIX_netcdf_variables.fi                                     \
            MFIX_netcdf_misc.fi                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) mfix_netcdf_mod.f 
mflux.mod : mflux_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) mflux_mod.f 
output.mod : output_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) output_mod.f 
parallel.mod : parallel_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parallel_mod.f 
param1.mod : param1_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) param1_mod.f 
param.mod : param_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) param_mod.f 
parse.mod : parse_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) parse_mod.f 
pgcor.mod : pgcor_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pgcor_mod.f 
physprop.mod : physprop_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) physprop_mod.f 
pscor.mod : pscor_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) pscor_mod.f 
residual.mod : residual_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) residual_mod.f 
run.mod : run_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) run_mod.f 
rxns.mod : rxns_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) rxns_mod.f 
scalars.mod : scalars_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scalars_mod.f 
scales.mod : scales_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) scales_mod.f 
tau_g.mod : tau_g_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_g_mod.f 
tau_s.mod : tau_s_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tau_s_mod.f 
time_cpu.mod : time_cpu_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) time_cpu_mod.f 
tmp_array1.mod : tmp_array1_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array1_mod.f 
tmp_array.mod : tmp_array_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) tmp_array_mod.f 
toleranc.mod : toleranc_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) toleranc_mod.f 
trace.mod : trace_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) trace_mod.f 
turb.mod : turb_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) turb_mod.f 
ur_facs.mod : ur_facs_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ur_facs_mod.f 
usr.mod : usr_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) usr_mod.f 
visc_g.mod : visc_g_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_g_mod.f 
visc_s.mod : visc_s_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) visc_s_mod.f 
vshear.mod : vshear_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) vshear_mod.f 
xsi_array.mod : xsi_array_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) xsi_array_mod.f 
cutcell.mod : ./cartesian_grid/cutcell_mod.f \
            param.mod \
            param1.mod \
            progress_bar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cutcell_mod.f 
dashboard.mod : ./cartesian_grid/dashboard_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/dashboard_mod.f 
polygon.mod : ./cartesian_grid/polygon_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/polygon_mod.f 
progress_bar.mod : ./cartesian_grid/progress_bar_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/progress_bar_mod.f 
quadric.mod : ./cartesian_grid/quadric_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/quadric_mod.f 
stl.mod : ./cartesian_grid/stl_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/stl_mod.f 
vtk.mod : ./cartesian_grid/vtk_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/vtk_mod.f 
mchem.mod : ./chem/mchem_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/mchem_mod.f 
des_bc.mod : ./des/des_bc_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_bc_mod.f 
desgrid.mod : ./des/desgrid_mod.f \
            param1.mod \
            funits.mod \
            geometry.mod \
            compar.mod \
            discretelement.mod \
            constant.mod \
            desmpi_wrapper.mod \
            des_thermo.mod \
            des/desgrid_functions.inc                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desgrid_mod.f 
des_ic.mod : ./des/des_ic_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_ic_mod.f 
desmpi.mod : ./des/desmpi_mod.f \
            parallel_mpi.mod \
            mpi_utility.mod \
            discretelement.mod \
            desgrid.mod \
            compar.mod \
            physprop.mod \
            sendrecv.mod \
            des_bc.mod \
            desmpi_wrapper.mod \
            sendrecvnode.mod \
            mfix_pic.mod \
            des/desgrid_functions.inc                                    \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desmpi_mod.f 
desmpi_wrapper.mod : ./des/desmpi_wrapper_mod.f \
            parallel_mpi.mod \
            mpi_utility.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/desmpi_wrapper_mod.f 
des_rxns.mod : ./des/des_rxns_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_rxns_mod.f 
des_thermo.mod : ./des/des_thermo_mod.f \
            param.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_mod.f 
discretelement.mod : ./des/discretelement_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/discretelement_mod.f 
interpolation.mod : ./des/interpolation_mod.f \
            constant.mod \
            discretelement.mod \
            geometry.mod \
            param1.mod \
            compar.mod \
            indices.mod \
            param.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/interpolation_mod.f 
mfix_pic.mod : ./des/mfix_pic_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mfix_pic_mod.f 
mppic_wallbc.mod : ./des/mppic_wallbc_mod.f \
            param.mod \
            param1.mod \
            discretelement.mod \
            bc.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            funits.mod \
            mpi_utility.mod \
            constant.mod \
            physprop.mod \
            randomno.mod \
            cutcell.mod \
            fldvar.mod \
            mfix_pic.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mppic_wallbc_mod.f 
randomno.mod : ./des/randomno_mod.f \
            constant.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/randomno_mod.f 
sendrecvnode.mod : ./des/sendrecvnode_mod.f \
            parallel_mpi.mod \
            mpi_utility.mod \
            discretelement.mod \
            compar.mod \
            physprop.mod \
            sendrecv.mod \
            desmpi_wrapper.mod \
            desgrid.mod \
            function.inc                                                 \
            des/desgrid_functions.inc                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/sendrecvnode_mod.f 
compar.mod : ./dmp_modules/compar_mod.f \
            mpi.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/compar_mod.f 
dbg_util.mod : ./dmp_modules/dbg_util_mod.f \
            compar.mod \
            geometry.mod \
            parallel_mpi.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/dbg_util_mod.f 
debug.mod : ./dmp_modules/debug_mod.f \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/debug_mod.f 
gridmap.mod : ./dmp_modules/gridmap_mod.f \
            mpi_utility.mod \
            parallel_mpi.mod \
            geometry.mod \
            sendrecv.mod \
            compar.mod \
            run.mod \
            indices.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/gridmap_mod.f 
mpi.mod : ./dmp_modules/mpi_mod.f \
            mpif.h                                                      
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_mod.f 
mpi_utility.mod : ./dmp_modules/mpi_utility_mod.f \
            geometry.mod \
            compar.mod \
            parallel_mpi.mod \
            debug.mod \
            indices.mod \
            funits.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/mpi_utility_mod.f 
parallel_mpi.mod : ./dmp_modules/parallel_mpi_mod.f \
            geometry.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/parallel_mpi_mod.f 
sendrecv3.mod : ./dmp_modules/sendrecv3_mod.f \
            parallel_mpi.mod \
            debug.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            mpi.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv3_mod.f 
sendrecv.mod : ./dmp_modules/sendrecv_mod.f \
            parallel_mpi.mod \
            debug.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            mpi.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dmp_modules/sendrecv_mod.f 
ghdtheory.mod : ./GhdTheory/ghdtheory_mod.f \
            param.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ghdtheory_mod.f 
qmomk_bc.mod : ./qmomk/qmomk_bc_mod.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            bc.mod \
            qmom_kinetic_equation.mod \
            qmomk_quadrature.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_bc_mod.f 
qmomk_collision.mod : ./qmomk/qmomk_collision_mod.f \
            constant.mod \
            qmomk_parameters.mod \
            qmomk_quadrature.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_collision_mod.f 
qmomk_fluxes.mod : ./qmomk/qmomk_fluxes_mod.f \
            qmomk_parameters.mod \
            qmomk_collision.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_fluxes_mod.f 
qmom_kinetic_equation.mod : ./qmomk/qmom_kinetic_equation_mod.f \
            param.mod \
            param1.mod \
            qmomk_parameters.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmom_kinetic_equation_mod.f 
qmomk_parameters.mod : ./qmomk/qmomk_parameters_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_parameters_mod.f 
qmomk_quadrature.mod : ./qmomk/qmomk_quadrature_mod.f \
            qmomk_tools.mod \
            qmomk_parameters.mod \
            constant.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_quadrature_mod.f 
qmomk_tools.mod : ./qmomk/qmomk_tools_mod.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_tools_mod.f 
accum_resid.$(OBJ_EXT) : accum_resid.f \
            param.mod \
            param1.mod \
            matrix.mod \
            parallel.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            residual.mod \
            run.mod 
adjust_a_u_g.$(OBJ_EXT) : adjust_a_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            usr.mod \
            compar.mod \
            sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_u_s.$(OBJ_EXT) : adjust_a_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_v_g.$(OBJ_EXT) : adjust_a_v_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_v_s.$(OBJ_EXT) : adjust_a_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_w_g.$(OBJ_EXT) : adjust_a_w_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_a_w_s.$(OBJ_EXT) : adjust_a_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            sendrecv.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
adjust_dt.$(OBJ_EXT) : adjust_dt.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            compar.mod \
            mpi_utility.mod 
adjust_eps.$(OBJ_EXT) : adjust_eps.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            compar.mod \
            sendrecv.mod \
            function.inc                                                
adjust_leq.$(OBJ_EXT) : adjust_leq.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            leqsol.mod 
adjust_rop.$(OBJ_EXT) : adjust_rop.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
adjust_theta.$(OBJ_EXT) : adjust_theta.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            compar.mod \
            function.inc                                                
allocate_arrays.$(OBJ_EXT) : allocate_arrays.f \
            param.mod \
            param1.mod \
            ambm.mod \
            coeff.mod \
            cont.mod \
            drag.mod \
            energy.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            physprop.mod \
            pscor.mod \
            residual.mod \
            rxns.mod \
            run.mod \
            scalars.mod \
            turb.mod \
            tau_g.mod \
            tau_s.mod \
            tmp_array.mod \
            tmp_array1.mod \
            trace.mod \
            visc_g.mod \
            visc_s.mod \
            xsi_array.mod \
            vshear.mod \
            mflux.mod \
            mchem.mod \
            ghdtheory.mod \
            kintheory.mod \
            cdist.mod 
bc_phi.$(OBJ_EXT) : bc_phi.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            geometry.mod \
            output.mod \
            indices.mod \
            bc.mod \
            compar.mod \
            cutcell.mod \
            quadric.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
bc_theta.$(OBJ_EXT) : bc_theta.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            geometry.mod \
            output.mod \
            indices.mod \
            bc.mod \
            compar.mod \
            mpi_utility.mod \
            turb.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
b_m_p_star.$(OBJ_EXT) : b_m_p_star.f \
            param.mod \
            param1.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            rxns.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
bound_x.$(OBJ_EXT) : bound_x.f \
            param.mod \
            param1.mod 
calc_cell.$(OBJ_EXT) : calc_cell.f \
            param.mod \
            param1.mod 
calc_coeff.$(OBJ_EXT) : calc_coeff.f \
            param.mod \
            param1.mod \
            physprop.mod \
            rxns.mod \
            funits.mod \
            compar.mod \
            ur_facs.mod \
            run.mod \
            discretelement.mod \
            des_rxns.mod 
calc_d.$(OBJ_EXT) : calc_d.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            scales.mod \
            compar.mod \
            sendrecv.mod \
            cutcell.mod \
            qmom_kinetic_equation.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_dif_g.$(OBJ_EXT) : calc_dif_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            constant.mod \
            scales.mod \
            compar.mod \
            sendrecv.mod \
            run.mod \
            function.inc                                                
calc_dif_s.$(OBJ_EXT) : calc_dif_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            constant.mod \
            toleranc.mod \
            compar.mod \
            sendrecv.mod \
            run.mod \
            function.inc                                                
calc_drag.$(OBJ_EXT) : calc_drag.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            drag.mod \
            compar.mod \
            discretelement.mod \
            qmom_kinetic_equation.mod 
calc_e.$(OBJ_EXT) : calc_e.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            constant.mod \
            compar.mod \
            sendrecv.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_gama.$(OBJ_EXT) : calc_gama.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            energy.mod \
            rxns.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            discretelement.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
calc_grbdry.$(OBJ_EXT) : calc_grbdry.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            turb.mod \
            visc_s.mod \
            geometry.mod \
            indices.mod \
            bc.mod \
            compar.mod \
            toleranc.mod \
            mpi_utility.mod \
            cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_h.$(OBJ_EXT) : calc_h.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            constant.mod \
            des_rxns.mod \
            des_thermo.mod \
            discretelement.mod 
calc_k_cp.$(OBJ_EXT) : calc_k_cp.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            pscor.mod \
            geometry.mod \
            constant.mod \
            run.mod \
            visc_s.mod \
            trace.mod \
            compar.mod \
            sendrecv.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
calc_k_g.$(OBJ_EXT) : calc_k_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            constant.mod \
            compar.mod \
            run.mod \
            sendrecv.mod \
            function.inc                                                
calc_k_s.$(OBJ_EXT) : calc_k_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            constant.mod \
            toleranc.mod \
            compar.mod \
            sendrecv.mod \
            run.mod \
            function.inc                                                
calc_mflux.$(OBJ_EXT) : calc_mflux.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            mflux.mod \
            physprop.mod \
            run.mod \
            parallel.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
calc_mu_g.$(OBJ_EXT) : calc_mu_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            visc_g.mod \
            visc_s.mod \
            indices.mod \
            constant.mod \
            toleranc.mod \
            compar.mod \
            drag.mod \
            run.mod \
            turb.mod \
            sendrecv.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                
calc_mu_s.$(OBJ_EXT) : calc_mu_s.f \
            run.mod \
            vshear.mod \
            visc_s.mod \
            physprop.mod \
            constant.mod \
            fldvar.mod \
            compar.mod \
            indices.mod \
            geometry.mod \
            qmom_kinetic_equation.mod \
            param.mod \
            param1.mod \
            trace.mod \
            toleranc.mod \
            turb.mod \
            drag.mod \
            kintheory.mod \
            ur_facs.mod \
            parallel.mod \
            visc_g.mod \
            is.mod \
            sendrecv.mod \
            cutcell.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            s_pr1.inc                                                    \
            s_pr2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
calc_mw.$(OBJ_EXT) : calc_mw.f \
            param.mod \
            param1.mod \
            toleranc.mod 
calc_outflow.$(OBJ_EXT) : calc_outflow.f \
            param.mod \
            param1.mod \
            bc.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
calc_p_star.$(OBJ_EXT) : calc_p_star.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            pgcor.mod \
            pscor.mod \
            ur_facs.mod \
            residual.mod \
            compar.mod \
            run.mod \
            visc_s.mod \
            fldvar.mod \
            toleranc.mod \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_resid.$(OBJ_EXT) : calc_resid.f \
            param.mod \
            param1.mod \
            matrix.mod \
            parallel.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            run.mod \
            fldvar.mod \
            bc.mod \
            constant.mod \
            physprop.mod \
            residual.mod \
            rxns.mod \
            mflux.mod \
            function.inc                                                
calc_s_ddot_s.$(OBJ_EXT) : calc_s_ddot_s.f \
            param.mod \
            param1.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
calc_trd_g.$(OBJ_EXT) : calc_trd_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            bc.mod \
            cutcell.mod \
            quadric.mod \
            function.inc                                                
calc_trd_s.$(OBJ_EXT) : calc_trd_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            sendrecv.mod \
            bc.mod \
            cutcell.mod \
            quadric.mod \
            function.inc                                                
calc_u_friction.$(OBJ_EXT) : calc_u_friction.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            turb.mod \
            visc_s.mod \
            geometry.mod \
            indices.mod \
            bc.mod \
            compar.mod \
            toleranc.mod \
            mpi_utility.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
calc_vol_fr.$(OBJ_EXT) : calc_vol_fr.f \
            param.mod \
            param1.mod \
            parallel.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            visc_s.mod \
            constant.mod \
            pgcor.mod \
            pscor.mod \
            compar.mod \
            sendrecv.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) calc_vol_fr.f 
calc_xsi.$(OBJ_EXT) : calc_xsi.f \
            param.mod \
            param1.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            vshear.mod \
            chischeme.mod \
            compar.mod \
            sendrecv.mod \
            xsi1.inc                                                     \
            function.inc                                                 \
            xsi2.inc                                                    
cal_d.$(OBJ_EXT) : cal_d.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            bc.mod \
            vshear.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
check_ab_m.$(OBJ_EXT) : check_ab_m.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
check_convergence.$(OBJ_EXT) : check_convergence.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            residual.mod \
            toleranc.mod \
            mpi_utility.mod 
check_data_01.$(OBJ_EXT) : check_data_01.f \
            param.mod \
            param1.mod \
            constant.mod \
            run.mod \
            physprop.mod \
            indices.mod \
            scalars.mod \
            funits.mod 
check_data_02.$(OBJ_EXT) : check_data_02.f \
            param.mod \
            param1.mod \
            output.mod \
            leqsol.mod \
            geometry.mod \
            run.mod \
            rxns.mod 
check_data_03.$(OBJ_EXT) : check_data_03.f \
            param.mod \
            param1.mod \
            geometry.mod \
            bc.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod 
check_data_04.$(OBJ_EXT) : check_data_04.f \
            param.mod \
            param1.mod \
            run.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            discretelement.mod \
            funits.mod 
check_data_05.$(OBJ_EXT) : check_data_05.f \
            param.mod \
            param1.mod \
            physprop.mod \
            funits.mod \
            run.mod \
            indices.mod 
check_data_06.$(OBJ_EXT) : check_data_06.f \
            param.mod \
            param1.mod \
            geometry.mod \
            ic.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            indices.mod \
            funits.mod \
            scalars.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            rxns.mod \
            function.inc                                                
check_data_07.$(OBJ_EXT) : check_data_07.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            bc.mod \
            indices.mod \
            funits.mod \
            scalars.mod \
            compar.mod \
            sendrecv.mod \
            function.inc                                                
check_data_08.$(OBJ_EXT) : check_data_08.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            is.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
check_data_09.$(OBJ_EXT) : check_data_09.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            rxns.mod \
            indices.mod \
            funits.mod \
            compar.mod 
check_data_20.$(OBJ_EXT) : check_data_20.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            constant.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            visc_g.mod \
            rxns.mod \
            scalars.mod \
            compar.mod \
            sendrecv.mod \
            function.inc                                                
check_data_30.$(OBJ_EXT) : check_data_30.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            fldvar.mod \
            rxns.mod \
            visc_s.mod \
            visc_g.mod \
            geometry.mod \
            run.mod \
            constant.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod \
            discretelement.mod \
            function.inc                                                
check_mass_balance.$(OBJ_EXT) : check_mass_balance.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            fldvar.mod \
            rxns.mod \
            geometry.mod \
            run.mod \
            bc.mod \
            constant.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod \
            output.mod \
            check.mod \
            mchem.mod \
            mflux.mod \
            xsi_array.mod \
            parallel.mod \
            matrix.mod \
            function.inc                                                
check_one_axis.$(OBJ_EXT) : check_one_axis.f \
            param.mod \
            param1.mod \
            funits.mod 
check_plane.$(OBJ_EXT) : check_plane.f \
            funits.mod \
            compar.mod 
cn_extrapol.$(OBJ_EXT) : cn_extrapol.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            scalars.mod \
            trace.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            function.inc                                                
compare.$(OBJ_EXT) : compare.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
conv_dif_phi.$(OBJ_EXT) : conv_dif_phi.f \
            param.mod \
            param1.mod \
            run.mod \
            geometry.mod \
            compar.mod \
            sendrecv.mod \
            xsi_array.mod \
            mpi_utility.mod \
            indices.mod \
            parallel.mod \
            matrix.mod \
            toleranc.mod \
            cutcell.mod \
            sendrecv3.mod \
            tmp_array.mod \
            vshear.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            output.mod \
            is.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            function3.inc                                                \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
conv_dif_u_g.$(OBJ_EXT) : conv_dif_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            visc_g.mod \
            compar.mod \
            toleranc.mod \
            physprop.mod \
            fldvar.mod \
            output.mod \
            mflux.mod \
            cutcell.mod \
            vshear.mod \
            xsi_array.mod \
            tmp_array.mod \
            sendrecv.mod \
            sendrecv3.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_dif_u_s.$(OBJ_EXT) : conv_dif_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            physprop.mod \
            visc_s.mod \
            compar.mod \
            toleranc.mod \
            fldvar.mod \
            output.mod \
            mflux.mod \
            cutcell.mod \
            xsi_array.mod \
            tmp_array.mod \
            sendrecv.mod \
            sendrecv3.mod \
            vshear.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_dif_v_g.$(OBJ_EXT) : conv_dif_v_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            visc_g.mod \
            compar.mod \
            toleranc.mod \
            physprop.mod \
            fldvar.mod \
            output.mod \
            mflux.mod \
            cutcell.mod \
            xsi_array.mod \
            vshear.mod \
            tmp_array.mod \
            sendrecv.mod \
            sendrecv3.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_dif_v_s.$(OBJ_EXT) : conv_dif_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            physprop.mod \
            visc_s.mod \
            compar.mod \
            toleranc.mod \
            fldvar.mod \
            output.mod \
            mflux.mod \
            cutcell.mod \
            xsi_array.mod \
            tmp_array.mod \
            sendrecv.mod \
            sendrecv3.mod \
            vshear.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_dif_w_g.$(OBJ_EXT) : conv_dif_w_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            visc_g.mod \
            compar.mod \
            toleranc.mod \
            physprop.mod \
            fldvar.mod \
            output.mod \
            mflux.mod \
            cutcell.mod \
            xsi_array.mod \
            tmp_array.mod \
            sendrecv.mod \
            sendrecv3.mod \
            vshear.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_dif_w_s.$(OBJ_EXT) : conv_dif_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            physprop.mod \
            visc_s.mod \
            compar.mod \
            toleranc.mod \
            fldvar.mod \
            output.mod \
            mflux.mod \
            cutcell.mod \
            xsi_array.mod \
            tmp_array.mod \
            sendrecv.mod \
            sendrecv3.mod \
            vshear.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            function3.inc                                               
conv_pp_g.$(OBJ_EXT) : conv_pp_g.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            parallel.mod \
            matrix.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            compar.mod \
            mflux.mod \
            function.inc                                                
conv_rop.$(OBJ_EXT) : conv_rop.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            mflux.mod \
            physprop.mod \
            run.mod \
            parallel.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            xsi_array.mod \
            function.inc                                                
conv_rop_g.$(OBJ_EXT) : conv_rop_g.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            compar.mod \
            parallel.mod \
            matrix.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            xsi_array.mod \
            function.inc                                                
conv_rop_s.$(OBJ_EXT) : conv_rop_s.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            compar.mod \
            parallel.mod \
            matrix.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            pscor.mod \
            xsi_array.mod \
            function.inc                                                
conv_source_epp.$(OBJ_EXT) : conv_source_epp.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            compar.mod \
            sendrecv.mod \
            xsi_array.mod \
            parallel.mod \
            matrix.mod \
            constant.mod \
            physprop.mod \
            rxns.mod \
            indices.mod \
            pgcor.mod \
            pscor.mod \
            vshear.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
copy_a.$(OBJ_EXT) : copy_a.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            physprop.mod \
            function.inc                                                
corner.$(OBJ_EXT) : corner.f \
            param.mod \
            param1.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            matrix.mod \
            corner.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
correct_0.$(OBJ_EXT) : correct_0.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            pgcor.mod \
            ur_facs.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            cutcell.mod \
            quadric.mod \
            function.inc                                                
correct_1.$(OBJ_EXT) : correct_1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            geometry.mod \
            pscor.mod \
            ur_facs.mod \
            constant.mod \
            compar.mod \
            sendrecv.mod \
            cutcell.mod \
            ep_s1.inc                                                    \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            ep_s2.inc                                                   
dgtsl.$(OBJ_EXT) : dgtsl.f 
dif_u_is.$(OBJ_EXT) : dif_u_is.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            output.mod \
            indices.mod \
            is.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
dif_v_is.$(OBJ_EXT) : dif_v_is.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            output.mod \
            indices.mod \
            is.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
dif_w_is.$(OBJ_EXT) : dif_w_is.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            output.mod \
            indices.mod \
            is.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
discretize.$(OBJ_EXT) : discretize.f \
            param.mod \
            param1.mod \
            run.mod 
display_resid.$(OBJ_EXT) : display_resid.f \
            param.mod \
            param1.mod \
            physprop.mod \
            residual.mod \
            fldvar.mod \
            compar.mod \
            geometry.mod 
drag_gs.$(OBJ_EXT) : drag_gs.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            constant.mod \
            compar.mod \
            drag.mod \
            sendrecv.mod \
            discretelement.mod \
            ur_facs.mod \
            mfix_pic.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
drag_ss.$(OBJ_EXT) : drag_ss.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            sendrecv.mod \
            drag.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
eosg.$(OBJ_EXT) : eosg.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            scales.mod \
            sc_p_g1.inc                                                  \
            sc_p_g2.inc                                                 
equal.$(OBJ_EXT) : equal.f \
            param.mod \
            param1.mod \
            indices.mod \
            physprop.mod 
error_routine.$(OBJ_EXT) : error_routine.f \
            funits.mod \
            compar.mod \
            mpi_utility.mod 
exchange.$(OBJ_EXT) : exchange.f \
            param.mod \
            param1.mod \
            compar.mod 
exit.$(OBJ_EXT) : exit.f \
            funits.mod \
            compar.mod \
            mpi_utility.mod 
flow_to_vel.$(OBJ_EXT) : flow_to_vel.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            run.mod \
            bc.mod \
            scales.mod \
            indices.mod \
            funits.mod \
            compar.mod 
g_0.$(OBJ_EXT) : g_0.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            visc_s.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
get_bc_area.$(OBJ_EXT) : get_bc_area.f \
            param.mod \
            param1.mod \
            geometry.mod \
            bc.mod \
            compar.mod 
get_data.$(OBJ_EXT) : get_data.f \
            param.mod \
            param1.mod \
            run.mod \
            funits.mod \
            compar.mod \
            gridmap.mod \
            discretelement.mod \
            leqsol.mod \
            parallel.mod \
            qmom_kinetic_equation.mod 
get_eq.$(OBJ_EXT) : get_eq.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod 
get_flow_bc.$(OBJ_EXT) : get_flow_bc.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            bc.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            sendrecv.mod \
            function.inc                                                
get_hloss.$(OBJ_EXT) : get_hloss.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            bc.mod \
            indices.mod \
            energy.mod 
get_is.$(OBJ_EXT) : get_is.f \
            param.mod \
            param1.mod \
            geometry.mod \
            is.mod \
            indices.mod \
            funits.mod \
            compar.mod 
get_philoss.$(OBJ_EXT) : get_philoss.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            bc.mod \
            indices.mod \
            energy.mod \
            compar.mod \
            function.inc                                                
get_smass.$(OBJ_EXT) : get_smass.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            function.inc                                                
get_stats.$(OBJ_EXT) : get_stats.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            funits.mod \
            residual.mod \
            run.mod \
            compar.mod \
            function.inc                                                
get_walls_bc.$(OBJ_EXT) : get_walls_bc.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            bc.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            sendrecv.mod \
            function.inc                                                
in_bin_512.$(OBJ_EXT) : in_bin_512.f \
            machine.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
in_bin_512i.$(OBJ_EXT) : in_bin_512i.f \
            machine.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
init_ab_m.$(OBJ_EXT) : init_ab_m.f \
            param.mod \
            param1.mod \
            matrix.mod \
            parallel.mod \
            compar.mod 
init_fvars.$(OBJ_EXT) : init_fvars.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            scalars.mod \
            rxns.mod \
            run.mod \
            compar.mod 
init_namelist.$(OBJ_EXT) : init_namelist.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            geometry.mod \
            ic.mod \
            bc.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            is.mod \
            toleranc.mod \
            scales.mod \
            ur_facs.mod \
            leqsol.mod \
            residual.mod \
            rxns.mod \
            scalars.mod \
            compar.mod \
            parallel.mod \
            cdist.mod \
            namelist.inc                                                
init_resid.$(OBJ_EXT) : init_resid.f \
            param.mod \
            param1.mod \
            physprop.mod \
            residual.mod 
iterate.$(OBJ_EXT) : iterate.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            output.mod \
            indices.mod \
            funits.mod \
            time_cpu.mod \
            pscor.mod \
            coeff.mod \
            leqsol.mod \
            visc_g.mod \
            pgcor.mod \
            cont.mod \
            scalars.mod \
            compar.mod \
            mpi_utility.mod \
            discretelement.mod \
            residual.mod \
            cutcell.mod \
            vtk.mod \
            dashboard.mod \
            qmom_kinetic_equation.mod \
            bc.mod \
            constant.mod 
k_epsilon_prop.$(OBJ_EXT) : k_epsilon_prop.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            drag.mod \
            run.mod \
            output.mod \
            geometry.mod \
            fldvar.mod \
            visc_g.mod \
            visc_s.mod \
            trace.mod \
            indices.mod \
            constant.mod \
            vshear.mod \
            turb.mod \
            toleranc.mod \
            compar.mod \
            tau_g.mod \
            sendrecv.mod \
            cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                
kintheory_drag_ss.$(OBJ_EXT) : kintheory_drag_ss.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            sendrecv.mod \
            drag.mod \
            kintheory.mod \
            function.inc                                                
kintheory_energy_dissipation_ss.$(OBJ_EXT) : kintheory_energy_dissipation_ss.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            constant.mod \
            toleranc.mod \
            kintheory.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
kintheory_u_s.$(OBJ_EXT) : kintheory_u_s.f \
            param.mod \
            param1.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            kintheory.mod \
            fldvar.mod \
            visc_s.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
kintheory_v_s.$(OBJ_EXT) : kintheory_v_s.f \
            param.mod \
            param1.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            kintheory.mod \
            fldvar.mod \
            visc_s.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
kintheory_w_s.$(OBJ_EXT) : kintheory_w_s.f \
            param.mod \
            param1.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            kintheory.mod \
            fldvar.mod \
            visc_s.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
leq_bicgs.$(OBJ_EXT) : leq_bicgs.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            leqsol.mod \
            funits.mod \
            parallel.mod \
            mpi_utility.mod \
            sendrecv.mod \
            function.inc                                                
leq_bicgst.$(OBJ_EXT) : leq_bicgst.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            leqsol.mod \
            funits.mod \
            parallel.mod \
            mpi_utility.mod \
            sendrecv.mod \
            function.inc                                                
leq_cg.$(OBJ_EXT) : leq_cg.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            leqsol.mod \
            funits.mod \
            parallel.mod \
            mpi_utility.mod \
            sendrecv.mod \
            function.inc                                                
leq_gmres.$(OBJ_EXT) : leq_gmres.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            debug.mod \
            compar.mod \
            mpi_utility.mod \
            parallel.mod \
            funits.mod \
            gridmap.mod \
            function.inc                                                
leq_sor.$(OBJ_EXT) : leq_sor.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            leqsol.mod \
            function.inc                                                
line_too_big.$(OBJ_EXT) : line_too_big.f 
location_check.$(OBJ_EXT) : location_check.f \
            param.mod \
            param1.mod \
            funits.mod \
            geometry.mod 
location.$(OBJ_EXT) : location.f \
            param.mod \
            param1.mod 
machine.$(OBJ_EXT) : machine.f \
            machine.mod \
            param.mod \
            run.mod \
            funits.mod 
make_upper_case.$(OBJ_EXT) : make_upper_case.f 
mark_phase_4_cor.$(OBJ_EXT) : mark_phase_4_cor.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            fldvar.mod \
            physprop.mod \
            constant.mod \
            compar.mod \
            visc_s.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) mark_phase_4_cor.f 
mfix.$(OBJ_EXT) : mfix.f \
            param.mod \
            param1.mod \
            run.mod \
            time_cpu.mod \
            funits.mod \
            output.mod \
            compar.mod \
            mpi_utility.mod \
            parallel_mpi.mod \
            discretelement.mod \
            mfix_pic.mod \
            cdist.mod \
            fldvar.mod \
            cutcell.mod \
            quadric.mod \
            dashboard.mod \
            parallel.mod \
            matrix.mod \
            geometry.mod \
            sendrecv.mod \
            sendrecv3.mod \
            indices.mod \
            leqsol.mod \
            function.inc                                                
mod_bc_i.$(OBJ_EXT) : mod_bc_i.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod \
            function.inc                                                
mod_bc_j.$(OBJ_EXT) : mod_bc_j.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod \
            function.inc                                                
mod_bc_k.$(OBJ_EXT) : mod_bc_k.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod \
            function.inc                                                
open_file.$(OBJ_EXT) : open_file.f \
            cdist.mod \
            compar.mod 
open_files.$(OBJ_EXT) : open_files.f \
            machine.mod \
            funits.mod \
            compar.mod \
            cdist.mod \
            run.mod 
out_array_c.$(OBJ_EXT) : out_array_c.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
out_array.$(OBJ_EXT) : out_array.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
out_array_kc.$(OBJ_EXT) : out_array_kc.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            mpi_utility.mod \
            function.inc                                                
out_array_k.$(OBJ_EXT) : out_array_k.f \
            param.mod \
            param1.mod \
            geometry.mod \
            fldvar.mod \
            physprop.mod \
            indices.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
out_bin_512.$(OBJ_EXT) : out_bin_512.f \
            machine.mod 
out_bin_512i.$(OBJ_EXT) : out_bin_512i.f \
            machine.mod 
out_bin_512r.$(OBJ_EXT) : out_bin_512r.f \
            machine.mod 
out_bin_r.$(OBJ_EXT) : out_bin_r.f \
            param.mod 
parse_line.$(OBJ_EXT) : parse_line.f \
            param.mod \
            param1.mod \
            parse.mod \
            compar.mod 
parse_resid_string.$(OBJ_EXT) : parse_resid_string.f \
            param.mod \
            param1.mod \
            physprop.mod \
            residual.mod \
            funits.mod \
            compar.mod 
parse_rxn.$(OBJ_EXT) : parse_rxn.f \
            param.mod \
            param1.mod \
            parse.mod \
            rxns.mod \
            compar.mod 
partial_elim.$(OBJ_EXT) : partial_elim.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            matrix.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            drag.mod \
            fldvar.mod \
            run.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
physical_prop.$(OBJ_EXT) : physical_prop.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            toleranc.mod \
            constant.mod \
            scalars.mod \
            compar.mod \
            funits.mod \
            usr.mod \
            mpi_utility.mod \
            discretelement.mod \
            cutcell.mod \
            usrnlst.inc                                                  \
            cp_fun1.inc                                                  \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
read_database.$(OBJ_EXT) : read_database.f \
            param.mod \
            param1.mod \
            physprop.mod \
            constant.mod \
            compar.mod \
            rxns.mod \
            funits.mod \
            discretelement.mod \
            des_rxns.mod \
            mfix_directory_path.inc                                     
read_namelist.$(OBJ_EXT) : read_namelist.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            geometry.mod \
            ic.mod \
            is.mod \
            bc.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            toleranc.mod \
            funits.mod \
            scales.mod \
            ur_facs.mod \
            leqsol.mod \
            residual.mod \
            rxns.mod \
            scalars.mod \
            compar.mod \
            parallel.mod \
            discretelement.mod \
            mfix_pic.mod \
            usr.mod \
            des_bc.mod \
            des_ic.mod \
            des_thermo.mod \
            des_rxns.mod \
            cdist.mod \
            quadric.mod \
            cutcell.mod \
            vtk.mod \
            polygon.mod \
            dashboard.mod \
            stl.mod \
            qmom_kinetic_equation.mod \
            usrnlst.inc                                                  \
            namelist.inc                                                 \
            des/desnamelist.inc                                          \
            cartesian_grid/cartesian_grid_namelist.inc                   \
            qmomk/qmomknamelist.inc                                     
read_res0.$(OBJ_EXT) : read_res0.f \
            param.mod \
            param1.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            ic.mod \
            bc.mod \
            is.mod \
            constant.mod \
            funits.mod \
            output.mod \
            scales.mod \
            ur_facs.mod \
            toleranc.mod \
            leqsol.mod \
            scalars.mod \
            rxns.mod \
            compar.mod \
            mpi_utility.mod \
            fldvar.mod 
read_res1.$(OBJ_EXT) : read_res1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            rxns.mod \
            scalars.mod \
            funits.mod \
            energy.mod \
            compar.mod \
            cdist.mod \
            mpi_utility.mod \
            sendrecv.mod \
            mfix_netcdf.mod 
remove_comment.$(OBJ_EXT) : remove_comment.f 
reset_new.$(OBJ_EXT) : reset_new.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            trace.mod \
            run.mod \
            scalars.mod 
rrates0.$(OBJ_EXT) : rrates0.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            rxns.mod \
            energy.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            funits.mod \
            compar.mod \
            sendrecv.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
rrates.$(OBJ_EXT) : rrates.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            rxns.mod \
            energy.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            funits.mod \
            compar.mod \
            sendrecv.mod \
            function.inc                                                
rrates_init.$(OBJ_EXT) : rrates_init.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            rxns.mod \
            energy.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
scalar_prop.$(OBJ_EXT) : scalar_prop.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            scalars.mod \
            toleranc.mod \
            compar.mod \
            sendrecv.mod \
            function.inc                                                
seek_comment.$(OBJ_EXT) : seek_comment.f 
seek_end.$(OBJ_EXT) : seek_end.f 
set_bc0.$(OBJ_EXT) : set_bc0.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod \
            mpi_utility.mod \
            physprop.mod \
            constant.mod \
            bc.mod \
            fldvar.mod \
            indices.mod \
            run.mod \
            funits.mod \
            scales.mod \
            scalars.mod \
            boundfunijk.mod \
            toleranc.mod \
            sendrecv.mod \
            sc_p_g1.inc                                                  \
            function.inc                                                 \
            sc_p_g2.inc                                                 
set_bc1.$(OBJ_EXT) : set_bc1.f \
            param.mod \
            param1.mod \
            bc.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
set_constants.$(OBJ_EXT) : set_constants.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            visc_s.mod \
            energy.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            run.mod \
            funits.mod \
            drag.mod \
            compar.mod 
set_constprop.$(OBJ_EXT) : set_constprop.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            visc_s.mod \
            visc_g.mod \
            energy.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            run.mod \
            funits.mod \
            drag.mod \
            compar.mod \
            kintheory.mod \
            function.inc                                                
set_flags.$(OBJ_EXT) : set_flags.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            bc.mod \
            is.mod \
            indices.mod \
            physprop.mod \
            funits.mod \
            compar.mod \
            sendrecv.mod \
            sendrecv3.mod \
            boundfunijk.mod \
            mpi_utility.mod \
            function.inc                                                 \
            function3.inc                                               
set_fluidbed_p.$(OBJ_EXT) : set_fluidbed_p.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            bc.mod \
            ic.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            funits.mod \
            scales.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            discretelement.mod \
            sc_p_g1.inc                                                  \
            b_force1.inc                                                 \
            function.inc                                                 \
            b_force2.inc                                                 \
            sc_p_g2.inc                                                 
set_geometry1.$(OBJ_EXT) : set_geometry1.f \
            param.mod \
            param1.mod \
            parallel.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
set_geometry.$(OBJ_EXT) : set_geometry.f \
            param.mod \
            param1.mod \
            run.mod \
            geometry.mod \
            compar.mod 
set_ic.$(OBJ_EXT) : set_ic.f \
            param.mod \
            param1.mod \
            geometry.mod \
            constant.mod \
            physprop.mod \
            ic.mod \
            fldvar.mod \
            visc_g.mod \
            indices.mod \
            scales.mod \
            energy.mod \
            scalars.mod \
            compar.mod \
            run.mod \
            sendrecv.mod \
            sc_p_g1.inc                                                  \
            s_pr1.inc                                                    \
            function.inc                                                 \
            s_pr2.inc                                                    \
            sc_p_g2.inc                                                 
set_increments3.$(OBJ_EXT) : set_increments3.f \
            param.mod \
            param1.mod \
            indices.mod \
            geometry.mod \
            compar.mod \
            physprop.mod \
            fldvar.mod \
            funits.mod \
            function.inc                                                 \
            function3.inc                                               
set_increments.$(OBJ_EXT) : set_increments.f \
            param.mod \
            param1.mod \
            indices.mod \
            geometry.mod \
            compar.mod \
            physprop.mod \
            fldvar.mod \
            funits.mod \
            function.inc                                                
set_index1a3.$(OBJ_EXT) : set_index1a3.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            boundfunijk3.mod \
            function.inc                                                
set_index1a.$(OBJ_EXT) : set_index1a.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            boundfunijk.mod \
            function.inc                                                
set_index1.$(OBJ_EXT) : set_index1.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            constant.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
set_l_scale.$(OBJ_EXT) : set_l_scale.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            visc_g.mod \
            geometry.mod \
            indices.mod \
            compar.mod 
set_max2.$(OBJ_EXT) : set_max2.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod 
set_mw_mix_g.$(OBJ_EXT) : set_mw_mix_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
set_outflow.$(OBJ_EXT) : set_outflow.f \
            param.mod \
            param1.mod \
            bc.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            scalars.mod \
            run.mod \
            compar.mod \
            mflux.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
set_ro_g.$(OBJ_EXT) : set_ro_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
set_wall_bc.$(OBJ_EXT) : set_wall_bc.f \
            param.mod \
            param1.mod \
            bc.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            compar.mod \
            function.inc                                                
shift_dxyz.$(OBJ_EXT) : shift_dxyz.f \
            param.mod \
            param1.mod \
            geometry.mod 
solve_continuity.$(OBJ_EXT) : solve_continuity.f \
            param.mod \
            param1.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            indices.mod \
            residual.mod \
            cont.mod \
            leqsol.mod \
            ambm.mod 
solve_energy_eq.$(OBJ_EXT) : solve_energy_eq.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            output.mod \
            indices.mod \
            drag.mod \
            residual.mod \
            ur_facs.mod \
            pgcor.mod \
            pscor.mod \
            leqsol.mod \
            bc.mod \
            energy.mod \
            rxns.mod \
            ambm.mod \
            tmp_array.mod \
            tmp_array1.mod \
            compar.mod \
            discretelement.mod \
            des_thermo.mod \
            mflux.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
solve_epp.$(OBJ_EXT) : solve_epp.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            pscor.mod \
            residual.mod \
            leqsol.mod \
            physprop.mod \
            ambm.mod \
            tmp_array1.mod 
solve_granular_energy.$(OBJ_EXT) : solve_granular_energy.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            visc_s.mod \
            geometry.mod \
            fldvar.mod \
            constant.mod \
            output.mod \
            indices.mod \
            drag.mod \
            residual.mod \
            ur_facs.mod \
            pgcor.mod \
            pscor.mod \
            leqsol.mod \
            bc.mod \
            energy.mod \
            rxns.mod \
            ambm.mod \
            tmp_array.mod \
            compar.mod \
            mflux.mod \
            radtn1.inc                                                   \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            radtn2.inc                                                  
solve_k_epsilon_eq.$(OBJ_EXT) : solve_k_epsilon_eq.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            constant.mod \
            output.mod \
            indices.mod \
            drag.mod \
            residual.mod \
            ur_facs.mod \
            pgcor.mod \
            pscor.mod \
            leqsol.mod \
            bc.mod \
            energy.mod \
            rxns.mod \
            turb.mod \
            usr.mod \
            ambm.mod \
            tmp_array.mod \
            compar.mod \
            mflux.mod \
            cutcell.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
solve_lin_eq.$(OBJ_EXT) : solve_lin_eq.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod \
            residual.mod \
            toleranc.mod \
            leqsol.mod 
solve_pp_g.$(OBJ_EXT) : solve_pp_g.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            pgcor.mod \
            residual.mod \
            leqsol.mod \
            run.mod \
            ambm.mod \
            tmp_array1.mod 
solve_scalar_eq.$(OBJ_EXT) : solve_scalar_eq.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            output.mod \
            indices.mod \
            drag.mod \
            residual.mod \
            ur_facs.mod \
            pgcor.mod \
            pscor.mod \
            leqsol.mod \
            bc.mod \
            energy.mod \
            rxns.mod \
            scalars.mod \
            ambm.mod \
            tmp_array.mod \
            compar.mod \
            mflux.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
solve_species_eq.$(OBJ_EXT) : solve_species_eq.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            output.mod \
            indices.mod \
            drag.mod \
            residual.mod \
            ur_facs.mod \
            pgcor.mod \
            pscor.mod \
            leqsol.mod \
            bc.mod \
            energy.mod \
            rxns.mod \
            ambm.mod \
            matrix.mod \
            chischeme.mod \
            tmp_array.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            mflux.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
solve_vel_star.$(OBJ_EXT) : solve_vel_star.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            ghdtheory.mod \
            output.mod \
            indices.mod \
            drag.mod \
            residual.mod \
            ur_facs.mod \
            pgcor.mod \
            pscor.mod \
            leqsol.mod \
            ambm.mod \
            tmp_array1.mod \
            tmp_array.mod \
            compar.mod \
            discretelement.mod \
            qmom_kinetic_equation.mod 
source_granular_energy.$(OBJ_EXT) : source_granular_energy.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            run.mod \
            drag.mod \
            geometry.mod \
            fldvar.mod \
            visc_g.mod \
            visc_s.mod \
            trace.mod \
            turb.mod \
            indices.mod \
            constant.mod \
            toleranc.mod \
            compar.mod \
            residual.mod \
            kintheory.mod \
            s_pr1.inc                                                    \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            ep_s2.inc                                                    \
            fun_avg2.inc                                                 \
            s_pr2.inc                                                   
source_phi.$(OBJ_EXT) : source_phi.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            compar.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
source_pp_g.$(OBJ_EXT) : source_pp_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            physprop.mod \
            fldvar.mod \
            rxns.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            bc.mod \
            vshear.mod \
            xsi_array.mod \
            compar.mod \
            ur_facs.mod \
            constant.mod \
            cutcell.mod \
            quadric.mod \
            function.inc                                                
source_rop_g.$(OBJ_EXT) : source_rop_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            rxns.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            compar.mod \
            function.inc                                                
source_rop_s.$(OBJ_EXT) : source_rop_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            fldvar.mod \
            rxns.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            pgcor.mod \
            pscor.mod \
            compar.mod \
            function.inc                                                
source_u_g.$(OBJ_EXT) : source_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            turb.mod \
            mpi_utility.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_u_s.$(OBJ_EXT) : source_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            kintheory.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_v_g.$(OBJ_EXT) : source_v_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            vshear.mod \
            compar.mod \
            sendrecv.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_v_s.$(OBJ_EXT) : source_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            bc.mod \
            vshear.mod \
            compar.mod \
            sendrecv.mod \
            kintheory.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_w_g.$(OBJ_EXT) : source_w_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
source_w_s.$(OBJ_EXT) : source_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            kintheory.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
tau_u_g.$(OBJ_EXT) : tau_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            compar.mod \
            sendrecv.mod \
            bc.mod \
            quadric.mod \
            cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_u_s.$(OBJ_EXT) : tau_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            vshear.mod \
            sendrecv.mod \
            compar.mod \
            bc.mod \
            quadric.mod \
            cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_v_g.$(OBJ_EXT) : tau_v_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            sendrecv.mod \
            compar.mod \
            bc.mod \
            quadric.mod \
            cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_v_s.$(OBJ_EXT) : tau_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            sendrecv.mod \
            compar.mod \
            bc.mod \
            quadric.mod \
            cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_w_g.$(OBJ_EXT) : tau_w_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            sendrecv.mod \
            compar.mod \
            bc.mod \
            quadric.mod \
            cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
tau_w_s.$(OBJ_EXT) : tau_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            sendrecv.mod \
            compar.mod \
            bc.mod \
            quadric.mod \
            cutcell.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
test_lin_eq.$(OBJ_EXT) : test_lin_eq.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            function.inc                                                
time_march.$(OBJ_EXT) : time_march.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            pgcor.mod \
            pscor.mod \
            cont.mod \
            coeff.mod \
            tau_g.mod \
            tau_s.mod \
            visc_g.mod \
            visc_s.mod \
            funits.mod \
            vshear.mod \
            scalars.mod \
            toleranc.mod \
            drag.mod \
            rxns.mod \
            compar.mod \
            time_cpu.mod \
            discretelement.mod \
            mchem.mod \
            leqsol.mod \
            cdist.mod \
            mfix_netcdf.mod \
            mpi_utility.mod \
            cutcell.mod \
            vtk.mod \
            dashboard.mod \
            qmom_kinetic_equation.mod 
transfer.$(OBJ_EXT) : transfer.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod 
transport_prop.$(OBJ_EXT) : transport_prop.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            physprop.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            toleranc.mod \
            compar.mod \
            discretelement.mod 
undef_2_0.$(OBJ_EXT) : undef_2_0.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod 
under_relax.$(OBJ_EXT) : under_relax.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            function.inc                                                
update_old.$(OBJ_EXT) : update_old.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            trace.mod \
            visc_s.mod \
            scalars.mod 
usr0.$(OBJ_EXT) : usr0.f \
            usr.mod 
usr1.$(OBJ_EXT) : usr1.f \
            usr.mod 
usr2.$(OBJ_EXT) : usr2.f \
            usr.mod 
usr3.$(OBJ_EXT) : usr3.f \
            usr.mod 
usr_init_namelist.$(OBJ_EXT) : usr_init_namelist.f \
            usr.mod 
usr_write_out0.$(OBJ_EXT) : usr_write_out0.f 
usr_write_out1.$(OBJ_EXT) : usr_write_out1.f 
utilities.$(OBJ_EXT) : utilities.f \
            param.mod \
            param1.mod \
            parallel.mod \
            bc.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            constant.mod \
            run.mod \
            compar.mod \
            toleranc.mod \
            mpi_utility.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vavg_u_g.$(OBJ_EXT) : vavg_u_g.f \
            param.mod \
            param1.mod \
            run.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            mflux.mod \
            function.inc                                                
vavg_u_s.$(OBJ_EXT) : vavg_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vavg_v_g.$(OBJ_EXT) : vavg_v_g.f \
            param.mod \
            param1.mod \
            run.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            mflux.mod \
            function.inc                                                
vavg_v_s.$(OBJ_EXT) : vavg_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vavg_w_g.$(OBJ_EXT) : vavg_w_g.f \
            param.mod \
            param1.mod \
            run.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            mflux.mod \
            function.inc                                                
vavg_w_s.$(OBJ_EXT) : vavg_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            bc.mod \
            geometry.mod \
            physprop.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
vf_gs_x.$(OBJ_EXT) : vf_gs_x.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            drag.mod \
            discretelement.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
vf_gs_y.$(OBJ_EXT) : vf_gs_y.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            drag.mod \
            discretelement.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
vf_gs_z.$(OBJ_EXT) : vf_gs_z.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            drag.mod \
            discretelement.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
vtc_scalar.$(OBJ_EXT) : vtc_scalar.f \
            param.mod \
            param1.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            kintheory.mod \
            function.inc                                                
write_ab_m.$(OBJ_EXT) : write_ab_m.f \
            param.mod \
            param1.mod \
            matrix.mod \
            compar.mod \
            mpi_utility.mod \
            indices.mod \
            function.inc                                                
write_ab_m_var.$(OBJ_EXT) : write_ab_m_var.f \
            param.mod \
            param1.mod \
            matrix.mod \
            geometry.mod \
            compar.mod \
            mpi_utility.mod \
            indices.mod \
            function.inc                                                
write_error.$(OBJ_EXT) : write_error.f \
            param.mod \
            param1.mod \
            funits.mod 
write_header.$(OBJ_EXT) : write_header.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            funits.mod \
            compar.mod 
write_out0.$(OBJ_EXT) : write_out0.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            geometry.mod \
            ic.mod \
            bc.mod \
            is.mod \
            fldvar.mod \
            constant.mod \
            indices.mod \
            funits.mod \
            toleranc.mod \
            scales.mod \
            scalars.mod \
            ur_facs.mod \
            leqsol.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            function.inc                                                
write_out1.$(OBJ_EXT) : write_out1.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            rxns.mod \
            compar.mod \
            mpi_utility.mod 
write_out3.$(OBJ_EXT) : write_out3.f \
            funits.mod \
            compar.mod 
write_res0.$(OBJ_EXT) : write_res0.f \
            param.mod \
            param1.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            ic.mod \
            is.mod \
            bc.mod \
            constant.mod \
            funits.mod \
            output.mod \
            scales.mod \
            scalars.mod \
            rxns.mod \
            ur_facs.mod \
            leqsol.mod \
            toleranc.mod \
            cdist.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod 
write_res1.$(OBJ_EXT) : write_res1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            scalars.mod \
            rxns.mod \
            funits.mod \
            output.mod \
            energy.mod \
            cdist.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            mfix_netcdf.mod 
write_spx0.$(OBJ_EXT) : write_spx0.f \
            param.mod \
            param1.mod \
            run.mod \
            funits.mod \
            cdist.mod \
            compar.mod \
            mpi_utility.mod 
write_spx1.$(OBJ_EXT) : write_spx1.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            run.mod \
            funits.mod \
            scalars.mod \
            output.mod \
            rxns.mod \
            cdist.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            mfix_netcdf.mod 
write_table.$(OBJ_EXT) : write_table.f \
            param.mod \
            param1.mod \
            funits.mod 
write_usr0.$(OBJ_EXT) : write_usr0.f 
write_usr1.$(OBJ_EXT) : write_usr1.f 
xerbla.$(OBJ_EXT) : xerbla.f \
            compar.mod 
zero_array.$(OBJ_EXT) : zero_array.f \
            param.mod \
            param1.mod 
zero_norm_vel.$(OBJ_EXT) : zero_norm_vel.f \
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            physprop.mod \
            fldvar.mod \
            indices.mod \
            is.mod \
            compar.mod \
            function.inc                                                
allocate_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/allocate_cut_cell_arrays.f \
            param.mod \
            param1.mod \
            indices.mod \
            cutcell.mod \
            stl.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/allocate_cut_cell_arrays.f 
allocate_dummy_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/allocate_dummy_cut_cell_arrays.f \
            param.mod \
            param1.mod \
            indices.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/allocate_dummy_cut_cell_arrays.f 
calc_vort_out.$(OBJ_EXT) : ./cartesian_grid/calc_vort_out.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            fldvar.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/calc_vort_out.f 
cartesian_grid_init_namelist.$(OBJ_EXT) : ./cartesian_grid/cartesian_grid_init_namelist.f \
            param1.mod \
            quadric.mod \
            cutcell.mod \
            polygon.mod \
            vtk.mod \
            progress_bar.mod \
            dashboard.mod \
            stl.mod \
            cartesian_grid/cartesian_grid_namelist.inc                  
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cartesian_grid_init_namelist.f 
CG_set_bc0.$(OBJ_EXT) : ./cartesian_grid/CG_set_bc0.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod \
            mpi_utility.mod \
            physprop.mod \
            bc.mod \
            fldvar.mod \
            indices.mod \
            run.mod \
            funits.mod \
            scales.mod \
            scalars.mod \
            boundfunijk.mod \
            toleranc.mod \
            sendrecv.mod \
            cutcell.mod \
            quadric.mod \
            sc_p_g1.inc                                                  \
            function.inc                                                 \
            sc_p_g2.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_set_bc0.f 
CG_set_outflow.$(OBJ_EXT) : ./cartesian_grid/CG_set_outflow.f \
            param.mod \
            param1.mod \
            bc.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            scalars.mod \
            run.mod \
            compar.mod \
            mflux.mod \
            cutcell.mod \
            quadric.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_set_outflow.f 
CG_source_u_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_u_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_u_g.f 
CG_source_u_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_u_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            kintheory.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_u_s.f 
CG_source_v_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_v_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            vshear.mod \
            compar.mod \
            sendrecv.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_v_g.f 
CG_source_v_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_v_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            bc.mod \
            vshear.mod \
            compar.mod \
            sendrecv.mod \
            kintheory.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_v_s.f 
CG_source_w_g.$(OBJ_EXT) : ./cartesian_grid/CG_source_w_g.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_w_g.f 
CG_source_w_s.$(OBJ_EXT) : ./cartesian_grid/CG_source_w_s.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_s.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_s.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            kintheory.mod \
            ghdtheory.mod \
            drag.mod \
            cutcell.mod \
            quadric.mod \
            output.mod \
            b_force1.inc                                                 \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                    \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/CG_source_w_s.f 
check_data_cartesian.$(OBJ_EXT) : ./cartesian_grid/check_data_cartesian.f \
            param.mod \
            param1.mod \
            constant.mod \
            run.mod \
            physprop.mod \
            indices.mod \
            scalars.mod \
            funits.mod \
            leqsol.mod \
            compar.mod \
            mpi_utility.mod \
            bc.mod \
            discretelement.mod \
            cutcell.mod \
            quadric.mod \
            vtk.mod \
            polygon.mod \
            dashboard.mod \
            stl.mod \
            fldvar.mod \
            scales.mod \
            parallel.mod \
            toleranc.mod \
            geometry.mod \
            sendrecv.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/check_data_cartesian.f 
cut_cell_preprocessing.$(OBJ_EXT) : ./cartesian_grid/cut_cell_preprocessing.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            vtk.mod \
            fldvar.mod \
            polygon.mod \
            stl.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/cut_cell_preprocessing.f 
deallocate_cut_cell_arrays.$(OBJ_EXT) : ./cartesian_grid/deallocate_cut_cell_arrays.f \
            param.mod \
            param1.mod \
            indices.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/deallocate_cut_cell_arrays.f 
define_quadrics.$(OBJ_EXT) : ./cartesian_grid/define_quadrics.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            vtk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/define_quadrics.f 
dmp_cartesian.$(OBJ_EXT) : ./cartesian_grid/dmp_cartesian.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/dmp_cartesian.f 
eval_usr_fct.$(OBJ_EXT) : ./cartesian_grid/eval_usr_fct.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            fldvar.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/eval_usr_fct.f 
get_alpha.$(OBJ_EXT) : ./cartesian_grid/get_alpha.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            bc.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_alpha.f 
get_connectivity.$(OBJ_EXT) : ./cartesian_grid/get_connectivity.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            polygon.mod \
            stl.mod \
            fldvar.mod \
            vtk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_connectivity.f 
get_cut_cell_flags.$(OBJ_EXT) : ./cartesian_grid/get_cut_cell_flags.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            vtk.mod \
            polygon.mod \
            stl.mod \
            physprop.mod \
            fldvar.mod \
            scalars.mod \
            funits.mod \
            rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_cut_cell_flags.f 
get_cut_cell_volume_area.$(OBJ_EXT) : ./cartesian_grid/get_cut_cell_volume_area.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            polygon.mod \
            stl.mod \
            bc.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_cut_cell_volume_area.f 
get_delh.$(OBJ_EXT) : ./cartesian_grid/get_delh.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_delh.f 
get_master.$(OBJ_EXT) : ./cartesian_grid/get_master.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            bc.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_master.f 
get_poly_data.$(OBJ_EXT) : ./cartesian_grid/get_poly_data.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            rxns.mod \
            compar.mod \
            mpi_utility.mod \
            progress_bar.mod \
            polygon.mod \
            parallel.mod \
            constant.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_poly_data.f 
get_stl_data.$(OBJ_EXT) : ./cartesian_grid/get_stl_data.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            rxns.mod \
            compar.mod \
            mpi_utility.mod \
            progress_bar.mod \
            stl.mod \
            vtk.mod \
            quadric.mod \
            constant.mod \
            bc.mod \
            cutcell.mod \
            parallel.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            sendrecv.mod \
            stl.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/get_stl_data.f 
set_Odxyz.$(OBJ_EXT) : ./cartesian_grid/set_Odxyz.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            vtk.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/set_Odxyz.f 
update_dashboard.$(OBJ_EXT) : ./cartesian_grid/update_dashboard.f \
            compar.mod \
            parallel.mod \
            sendrecv.mod \
            run.mod \
            leqsol.mod \
            time_cpu.mod \
            residual.mod \
            dashboard.mod \
            vtk.mod \
            constant.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/update_dashboard.f 
vtk_out.$(OBJ_EXT) : ./cartesian_grid/vtk_out.f \
            param.mod \
            param1.mod \
            parallel.mod \
            constant.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            quadric.mod \
            cutcell.mod \
            fldvar.mod \
            visc_s.mod \
            physprop.mod \
            pgcor.mod \
            vtk.mod \
            rxns.mod \
            output.mod \
            scalars.mod \
            pscor.mod \
            discretelement.mod \
            mfix_pic.mod \
            mpi_utility.mod \
            polygon.mod \
            stl.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/vtk_out.f 
write_progress_bar.$(OBJ_EXT) : ./cartesian_grid/write_progress_bar.f \
            param.mod \
            param1.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            scalars.mod \
            funits.mod \
            rxns.mod \
            compar.mod \
            mpi_utility.mod \
            progress_bar.mod \
            parallel.mod \
            sendrecv.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cartesian_grid/write_progress_bar.f 
calc_jacobian.$(OBJ_EXT) : ./chem/calc_jacobian.f \
            param1.mod \
            mchem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/calc_jacobian.f 
check_data_chem.$(OBJ_EXT) : ./chem/check_data_chem.f \
            param1.mod \
            run.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/check_data_chem.f 
dgpadm.$(OBJ_EXT) : ./chem/dgpadm.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/dgpadm.f 
exponential.$(OBJ_EXT) : ./chem/exponential.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/exponential.f 
fex.$(OBJ_EXT) : ./chem/fex.f \
            run.mod \
            physprop.mod \
            toleranc.mod \
            usr.mod \
            mchem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/fex.f 
g_derivs.$(OBJ_EXT) : ./chem/g_derivs.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/g_derivs.f 
jac.$(OBJ_EXT) : ./chem/jac.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/jac.f 
mchem_init.$(OBJ_EXT) : ./chem/mchem_init.f \
            param1.mod \
            run.mod \
            physprop.mod \
            mchem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/mchem_init.f 
mchem_odepack_init.$(OBJ_EXT) : ./chem/mchem_odepack_init.f \
            mchem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/mchem_odepack_init.f 
mchem_time_march.$(OBJ_EXT) : ./chem/mchem_time_march.f \
            run.mod \
            physprop.mod \
            fldvar.mod \
            rxns.mod \
            mpi_utility.mod \
            toleranc.mod \
            mchem.mod \
            ep_s1.inc                                                    \
            function.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/mchem_time_march.f 
misat_table_init.$(OBJ_EXT) : ./chem/misat_table_init.f \
            mchem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/misat_table_init.f 
react.$(OBJ_EXT) : ./chem/react.f \
            param1.mod \
            toleranc.mod \
            fldvar.mod \
            physprop.mod \
            rxns.mod \
            run.mod \
            mchem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/react.f 
usrfg.$(OBJ_EXT) : ./chem/usrfg.f \
            param1.mod \
            run.mod \
            mchem.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./chem/usrfg.f 
add_part_to_link_list.$(OBJ_EXT) : ./cohesion/add_part_to_link_list.f \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/add_part_to_link_list.f 
calc_app_coh_force.$(OBJ_EXT) : ./cohesion/calc_app_coh_force.f \
            discretelement.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_app_coh_force.f 
calc_cap_coh_force.$(OBJ_EXT) : ./cohesion/calc_cap_coh_force.f \
            discretelement.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_cap_coh_force.f 
calc_cohesive_forces.$(OBJ_EXT) : ./cohesion/calc_cohesive_forces.f \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_cohesive_forces.f 
calc_esc_coh_force.$(OBJ_EXT) : ./cohesion/calc_esc_coh_force.f \
            discretelement.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_esc_coh_force.f 
calc_square_well.$(OBJ_EXT) : ./cohesion/calc_square_well.f \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_square_well.f 
calc_van_der_waals.$(OBJ_EXT) : ./cohesion/calc_van_der_waals.f \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/calc_van_der_waals.f 
check_link.$(OBJ_EXT) : ./cohesion/check_link.f \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_link.f 
check_sw_wall_interaction.$(OBJ_EXT) : ./cohesion/check_sw_wall_interaction.f \
            param1.mod \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_sw_wall_interaction.f 
check_vdw_wall_interaction.$(OBJ_EXT) : ./cohesion/check_vdw_wall_interaction.f \
            param1.mod \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/check_vdw_wall_interaction.f 
initialize_cohesion_parameters.$(OBJ_EXT) : ./cohesion/initialize_cohesion_parameters.f \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/initialize_cohesion_parameters.f 
initialize_coh_int_search.$(OBJ_EXT) : ./cohesion/initialize_coh_int_search.f \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/initialize_coh_int_search.f 
linked_interaction_eval.$(OBJ_EXT) : ./cohesion/linked_interaction_eval.f \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/linked_interaction_eval.f 
remove_part_from_link_list.$(OBJ_EXT) : ./cohesion/remove_part_from_link_list.f \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/remove_part_from_link_list.f 
unlinked_interaction_eval.$(OBJ_EXT) : ./cohesion/unlinked_interaction_eval.f \
            discretelement.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/unlinked_interaction_eval.f 
update_search_grids.$(OBJ_EXT) : ./cohesion/update_search_grids.f \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./cohesion/update_search_grids.f 
calc_attrition_des.$(OBJ_EXT) : ./des/calc_attrition_des.f \
            run.mod \
            param1.mod \
            discretelement.mod \
            geometry.mod \
            compar.mod \
            constant.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_attrition_des.f
calc_force_des.$(OBJ_EXT) : ./des/calc_force_des.f \
            run.mod \
            param1.mod \
            discretelement.mod \
            geometry.mod \
            compar.mod \
            constant.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_force_des.f
calc_rrate_des.$(OBJ_EXT) : ./des/calc_rrate_des.f \
            discretelement.mod \
            interpolation.mod \
            param1.mod \
            rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_rrate_des.f 
calc_thermo_des.$(OBJ_EXT) : ./des/calc_thermo_des.f \
            des_thermo.mod \
            discretelement.mod \
            fldvar.mod \
            interpolation.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/calc_thermo_des.f 
cell_near_wall.$(OBJ_EXT) : ./des/cell_near_wall.f \
            discretelement.mod \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            matrix.mod \
            indices.mod \
            physprop.mod \
            drag.mod \
            constant.mod \
            compar.mod \
            sendrecv.mod \
            toleranc.mod \
            is.mod \
            bc.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cell_near_wall.f 
cfassign.$(OBJ_EXT) : ./des/cfassign.f \
            param1.mod \
            physprop.mod \
            geometry.mod \
            constant.mod \
            compar.mod \
            parallel.mod \
            sendrecv.mod \
            discretelement.mod \
            mfix_pic.mod \
            param.mod \
            fldvar.mod \
            run.mod \
            indices.mod \
            bc.mod \
            cutcell.mod \
            b_force1.inc                                                 \
            b_force2.inc                                                 \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfassign.f 
cffctowall.$(OBJ_EXT) : ./des/cffctowall.f \
            param1.mod \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffctowall.f 
cffctow.$(OBJ_EXT) : ./des/cffctow.f \
            param1.mod \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cffctow.f 
cfnewvalues.$(OBJ_EXT) : ./des/cfnewvalues.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            constant.mod \
            fldvar.mod \
            discretelement.mod \
            des_bc.mod \
            mpi_utility.mod \
            mfix_pic.mod \
            mppic_wallbc.mod \
            randomno.mod \
            cutcell.mod \
            geometry.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfnewvalues.f 
cfrelvel.$(OBJ_EXT) : ./des/cfrelvel.f \
            discretelement.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfrelvel.f 
cfslide.$(OBJ_EXT) : ./des/cfslide.f \
            param1.mod \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslide.f 
cfslidewall.$(OBJ_EXT) : ./des/cfslidewall.f \
            discretelement.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfslidewall.f 
cfupdateold.$(OBJ_EXT) : ./des/cfupdateold.f \
            discretelement.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfupdateold.f 
cfwallcontact.$(OBJ_EXT) : ./des/cfwallcontact.f \
            param1.mod \
            constant.mod \
            parallel.mod \
            compar.mod \
            discretelement.mod \
            des_bc.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallcontact.f 
cfwallposvel.$(OBJ_EXT) : ./des/cfwallposvel.f \
            discretelement.mod \
            des_bc.mod \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            matrix.mod \
            indices.mod \
            physprop.mod \
            drag.mod \
            constant.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/cfwallposvel.f 
check_des_bc.$(OBJ_EXT) : ./des/check_des_bc.f \
            compar.mod \
            constant.mod \
            des_bc.mod \
            discretelement.mod \
            funits.mod \
            geometry.mod \
            indices.mod \
            param.mod \
            param1.mod \
            physprop.mod \
            run.mod \
            mfix_pic.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_bc.f 
check_des_data.$(OBJ_EXT) : ./des/check_des_data.f \
            param1.mod \
            geometry.mod \
            funits.mod \
            compar.mod \
            discretelement.mod \
            run.mod \
            constant.mod \
            physprop.mod \
            desgrid.mod \
            indices.mod \
            fldvar.mod \
            toleranc.mod \
            mpi_utility.mod \
            output.mod \
            mfix_pic.mod \
            cutcell.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_data.f 
check_des_ic.$(OBJ_EXT) : ./des/check_des_ic.f \
            des_ic.mod \
            discretelement.mod \
            param.mod \
            param1.mod \
            des_thermo.mod \
            des_rxns.mod \
            compar.mod \
            constant.mod \
            funits.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_ic.f 
check_des_rxns.$(OBJ_EXT) : ./des/check_des_rxns.f \
            compar.mod \
            des_rxns.mod \
            discretelement.mod \
            funits.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_rxns.f 
check_des_thermo.$(OBJ_EXT) : ./des/check_des_thermo.f \
            compar.mod \
            des_thermo.mod \
            discretelement.mod \
            funits.mod \
            interpolation.mod \
            physprop.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/check_des_thermo.f 
des_allocate_arrays.$(OBJ_EXT) : ./des/des_allocate_arrays.f \
            param.mod \
            param1.mod \
            constant.mod \
            discretelement.mod \
            indices.mod \
            geometry.mod \
            compar.mod \
            physprop.mod \
            des_bc.mod \
            funits.mod \
            desgrid.mod \
            desmpi.mod \
            mfix_pic.mod \
            des_thermo.mod \
            des_rxns.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_allocate_arrays.f 
des_check_particle.$(OBJ_EXT) : ./des/des_check_particle.f \
            compar.mod \
            constant.mod \
            des_bc.mod \
            discretelement.mod \
            funits.mod \
            geometry.mod \
            indices.mod \
            param1.mod \
            physprop.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_check_particle.f 
des_functions.$(OBJ_EXT) : ./des/des_functions.f \
            param.mod \
            param1.mod \
            discretelement.mod \
            compar.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_functions.f 
des_granular_temperature.$(OBJ_EXT) : ./des/des_granular_temperature.f \
            discretelement.mod \
            param.mod \
            param1.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            compar.mod \
            sendrecv.mod \
            physprop.mod \
            des_bc.mod \
            fldvar.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_granular_temperature.f 
des_init_arrays.$(OBJ_EXT) : ./des/des_init_arrays.f \
            param.mod \
            param1.mod \
            discretelement.mod \
            indices.mod \
            geometry.mod \
            compar.mod \
            physprop.mod \
            des_bc.mod \
            run.mod \
            desgrid.mod \
            desmpi.mod \
            des_thermo.mod \
            des_rxns.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_arrays.f 
des_init_bc.$(OBJ_EXT) : ./des/des_init_bc.f \
            compar.mod \
            constant.mod \
            des_bc.mod \
            discretelement.mod \
            funits.mod \
            geometry.mod \
            indices.mod \
            param.mod \
            param1.mod \
            physprop.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_bc.f 
des_init_namelist.$(OBJ_EXT) : ./des/des_init_namelist.f \
            param1.mod \
            discretelement.mod \
            mfix_pic.mod \
            des_bc.mod \
            des_ic.mod \
            des_thermo.mod \
            des_rxns.mod \
            des/desnamelist.inc                                         
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_init_namelist.f 
des_mass_inlet.$(OBJ_EXT) : ./des/des_mass_inlet.f \
            compar.mod \
            constant.mod \
            des_bc.mod \
            discretelement.mod \
            funits.mod \
            geometry.mod \
            indices.mod \
            param1.mod \
            physprop.mod \
            desgrid.mod \
            mpi_utility.mod \
            des_thermo.mod \
            des_rxns.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_mass_inlet.f 
des_physical_prop.$(OBJ_EXT) : ./des/des_physical_prop.f \
            des_rxns.mod \
            des_thermo.mod \
            discretelement.mod \
            funits.mod \
            param.mod \
            param1.mod \
            physprop.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_physical_prop.f 
des_reaction_model.$(OBJ_EXT) : ./des/des_reaction_model.f \
            constant.mod \
            des_rxns.mod \
            des_thermo.mod \
            discretelement.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_reaction_model.f 
des_rrates.$(OBJ_EXT) : ./des/des_rrates.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            rxns.mod \
            energy.mod \
            geometry.mod \
            run.mod \
            indices.mod \
            physprop.mod \
            constant.mod \
            funits.mod \
            toleranc.mod \
            compar.mod \
            sendrecv.mod \
            usr.mod \
            des_thermo.mod \
            des_rxns.mod \
            discretelement.mod \
            interpolation.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            usrnlst.inc                                                 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_rrates.f 
des_set_ic.$(OBJ_EXT) : ./des/des_set_ic.f \
            compar.mod \
            des_thermo.mod \
            discretelement.mod \
            des_ic.mod \
            des_rxns.mod \
            funits.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_set_ic.f 
des_thermo_cond.$(OBJ_EXT) : ./des/des_thermo_cond.f \
            constant.mod \
            des_thermo.mod \
            discretelement.mod \
            fldvar.mod \
            funits.mod \
            param1.mod \
            physprop.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_cond.f 
des_thermo_conv.$(OBJ_EXT) : ./des/des_thermo_conv.f \
            constant.mod \
            des_thermo.mod \
            discretelement.mod \
            fldvar.mod \
            interpolation.mod \
            param1.mod \
            physprop.mod \
            compar.mod \
            geometry.mod \
            indices.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_conv.f 
des_thermo_newvalues.$(OBJ_EXT) : ./des/des_thermo_newvalues.f \
            des_thermo.mod \
            des_rxns.mod \
            discretelement.mod \
            param1.mod \
            physprop.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_newvalues.f 
des_thermo_rad.$(OBJ_EXT) : ./des/des_thermo_rad.f \
            constant.mod \
            des_thermo.mod \
            discretelement.mod \
            param1.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_thermo_rad.f 
des_time_march.$(OBJ_EXT) : ./des/des_time_march.f \
            param.mod \
            param1.mod \
            run.mod \
            output.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            pgcor.mod \
            pscor.mod \
            cont.mod \
            coeff.mod \
            tau_g.mod \
            tau_s.mod \
            visc_g.mod \
            visc_s.mod \
            funits.mod \
            vshear.mod \
            scalars.mod \
            drag.mod \
            rxns.mod \
            compar.mod \
            time_cpu.mod \
            discretelement.mod \
            constant.mod \
            sendrecv.mod \
            des_bc.mod \
            cutcell.mod \
            mppic_wallbc.mod \
            mfix_pic.mod \
            des_thermo.mod \
            des_rxns.mod \
            interpolation.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_time_march.f 
des_wallbc_preprocessing.$(OBJ_EXT) : ./des/des_wallbc_preprocessing.f \
            param1.mod \
            funits.mod \
            run.mod \
            compar.mod \
            discretelement.mod \
            cutcell.mod \
            indices.mod \
            physprop.mod \
            parallel.mod \
            geometry.mod \
            bc.mod \
            constant.mod \
            fldvar.mod \
            mpi_utility.mod \
            sendrecv.mod \
            mfix_pic.mod \
            discretelement.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/des_wallbc_preprocessing.f 
drag_fgs.$(OBJ_EXT) : ./des/drag_fgs.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            discretelement.mod \
            cutcell.mod \
            constant.mod \
            drag.mod \
            interpolation.mod \
            desmpi.mod \
            mfix_pic.mod \
            ur_facs.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/drag_fgs.f 
gas_drag.$(OBJ_EXT) : ./des/gas_drag.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            discretelement.mod \
            drag.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/gas_drag.f 
generate_particle_config.$(OBJ_EXT) : ./des/generate_particle_config.f \
            param1.mod \
            geometry.mod \
            funits.mod \
            compar.mod \
            discretelement.mod \
            run.mod \
            constant.mod \
            physprop.mod \
            desmpi.mod \
            cdist.mod \
            mpi_utility.mod \
            mfix_pic.mod \
            param.mod \
            fldvar.mod \
            indices.mod \
            randomno.mod \
            sendrecv.mod \
            interpolation.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/generate_particle_config.f 
grid_based_neighbor_search.$(OBJ_EXT) : ./des/grid_based_neighbor_search.f \
            param1.mod \
            discretelement.mod \
            geometry.mod \
            des_bc.mod \
            des_thermo.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/grid_based_neighbor_search.f 
make_arrays_des.$(OBJ_EXT) : ./des/make_arrays_des.f \
            param1.mod \
            funits.mod \
            run.mod \
            compar.mod \
            discretelement.mod \
            cutcell.mod \
            desmpi.mod \
            mpi_utility.mod \
            geometry.mod \
            des_ic.mod \
            des_rxns.mod \
            cdist.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/make_arrays_des.f 
mppic_routines.$(OBJ_EXT) : ./des/mppic_routines.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            indices.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            discretelement.mod \
            constant.mod \
            cutcell.mod \
            interpolation.mod \
            mfix_pic.mod \
            drag.mod \
            desmpi.mod \
            toleranc.mod \
            quadric.mod \
            vtk.mod \
            matrix.mod \
            scales.mod \
            visc_s.mod \
            rxns.mod \
            is.mod \
            tau_s.mod \
            output.mod \
            des_bc.mod \
            mpi_utility.mod \
            mppic_wallbc.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/mppic_routines.f 
neighbour.$(OBJ_EXT) : ./des/neighbour.f \
            param1.mod \
            discretelement.mod \
            desgrid.mod \
            des_thermo.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/neighbour.f 
nsquare.$(OBJ_EXT) : ./des/nsquare.f \
            param1.mod \
            discretelement.mod \
            geometry.mod \
            des_bc.mod \
            des_thermo.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/nsquare.f 
octree.$(OBJ_EXT) : ./des/octree.f \
            run.mod \
            param1.mod \
            constant.mod \
            discretelement.mod \
            param.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            matrix.mod \
            indices.mod \
            physprop.mod \
            drag.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/octree.f 
particles_in_cell.$(OBJ_EXT) : ./des/particles_in_cell.f \
            param.mod \
            param1.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            compar.mod \
            parallel.mod \
            sendrecv.mod \
            discretelement.mod \
            desgrid.mod \
            desmpi.mod \
            cutcell.mod \
            mfix_pic.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/particles_in_cell.f 
quadtree.$(OBJ_EXT) : ./des/quadtree.f \
            run.mod \
            param1.mod \
            constant.mod \
            discretelement.mod \
            param.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            matrix.mod \
            indices.mod \
            physprop.mod \
            drag.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/quadtree.f 
read_des_restart.$(OBJ_EXT) : ./des/read_des_restart.f \
            param1.mod \
            run.mod \
            discretelement.mod \
            des_bc.mod \
            compar.mod \
            desmpi.mod \
            machine.mod \
            cdist.mod \
            mpi_utility.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/read_des_restart.f 
thermo_nbr.$(OBJ_EXT) : ./des/thermo_nbr.f \
            des_thermo.mod \
            discretelement.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/thermo_nbr.f 
walledgecontact.$(OBJ_EXT) : ./des/walledgecontact.f \
            discretelement.mod \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            matrix.mod \
            indices.mod \
            physprop.mod \
            drag.mod \
            constant.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/walledgecontact.f 
wallfacecontact.$(OBJ_EXT) : ./des/wallfacecontact.f \
            discretelement.mod \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            matrix.mod \
            indices.mod \
            physprop.mod \
            drag.mod \
            constant.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/wallfacecontact.f 
wallnodecontact.$(OBJ_EXT) : ./des/wallnodecontact.f \
            discretelement.mod \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            run.mod \
            geometry.mod \
            matrix.mod \
            indices.mod \
            physprop.mod \
            drag.mod \
            constant.mod \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/wallnodecontact.f 
write_des_data.$(OBJ_EXT) : ./des/write_des_data.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            discretelement.mod \
            run.mod \
            geometry.mod \
            physprop.mod \
            sendrecv.mod \
            des_bc.mod \
            mpi_utility.mod \
            compar.mod \
            desmpi.mod \
            cdist.mod \
            des_thermo.mod \
            function.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./des/write_des_data.f 
write_des_restart.$(OBJ_EXT) : ./des/write_des_restart.f \
            param1.mod \
            discretelement.mod \
            run.mod \
            des_bc.mod \
            compar.mod \
            desmpi.mod \
            machine.mod \
            cdist.mod \
            mpi_utility.mod 
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
            physprop.mod \
            constant.mod \
            fldvar.mod \
            scalars.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/source_population_eq.f 
usr_dqmom.$(OBJ_EXT) : ./dqmom/usr_dqmom.f \
            param.mod \
            param1.mod \
            run.mod \
            physprop.mod \
            geometry.mod \
            fldvar.mod \
            output.mod \
            indices.mod \
            rxns.mod \
            constant.mod \
            ambm.mod \
            compar.mod \
            scalars.mod \
            usr.mod \
            ep_s1.inc                                                    \
            ep_s2.inc                                                    \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./dqmom/usr_dqmom.f 
adjust_eps_ghd.$(OBJ_EXT) : ./GhdTheory/adjust_eps_ghd.f \
            param.mod \
            param1.mod \
            toleranc.mod \
            constant.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            compar.mod \
            sendrecv.mod \
            ghdtheory.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s1.inc                                                    \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/adjust_eps_ghd.f 
bulk_viscosity.$(OBJ_EXT) : ./GhdTheory/bulk_viscosity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/bulk_viscosity.f 
calc_d_ghd.$(OBJ_EXT) : ./GhdTheory/calc_d_ghd.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            indices.mod \
            physprop.mod \
            run.mod \
            scales.mod \
            compar.mod \
            sendrecv.mod \
            ep_s1.inc                                                    \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                 \
            ep_s2.inc                                                   
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_d_ghd.f 
calc_external_forces.$(OBJ_EXT) : ./GhdTheory/calc_external_forces.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            ghdtheory.mod \
            physprop.mod \
            run.mod \
            constant.mod \
            drag.mod \
            bc.mod \
            scales.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                 \
            b_force1.inc                                                 \
            b_force2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_external_forces.f 
calc_nflux.$(OBJ_EXT) : ./GhdTheory/calc_nflux.f \
            param.mod \
            param1.mod \
            parallel.mod \
            fldvar.mod \
            geometry.mod \
            physprop.mod \
            constant.mod \
            indices.mod \
            mflux.mod \
            compar.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/calc_nflux.f 
chi_ij_GHD.$(OBJ_EXT) : ./GhdTheory/chi_ij_GHD.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/chi_ij_GHD.f 
cooling_rate.$(OBJ_EXT) : ./GhdTheory/cooling_rate.f \
            compar.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/cooling_rate.f 
cooling_rate_tc.$(OBJ_EXT) : ./GhdTheory/cooling_rate_tc.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/cooling_rate_tc.f 
dufour_coeff.$(OBJ_EXT) : ./GhdTheory/dufour_coeff.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/dufour_coeff.f 
ghd.$(OBJ_EXT) : ./GhdTheory/ghd.f \
            drag.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/ghd.f 
ghdmassflux.$(OBJ_EXT) : ./GhdTheory/ghdmassflux.f \
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            visc_s.mod \
            ghdtheory.mod \
            physprop.mod \
            run.mod \
            constant.mod \
            toleranc.mod \
            drag.mod \
            is.mod \
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
            param.mod \
            param1.mod \
            parallel.mod \
            geometry.mod \
            matrix.mod \
            physprop.mod \
            indices.mod \
            run.mod \
            compar.mod \
            drag.mod \
            fldvar.mod \
            fun_avg1.inc                                                 \
            function.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/partial_elim_ghd.f 
pressure.$(OBJ_EXT) : ./GhdTheory/pressure.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/pressure.f 
shear_viscosity.$(OBJ_EXT) : ./GhdTheory/shear_viscosity.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/shear_viscosity.f 
source_ghd_granular_energy.$(OBJ_EXT) : ./GhdTheory/source_ghd_granular_energy.f \
            param.mod \
            param1.mod \
            parallel.mod \
            physprop.mod \
            run.mod \
            drag.mod \
            geometry.mod \
            fldvar.mod \
            visc_s.mod \
            ghdtheory.mod \
            trace.mod \
            indices.mod \
            constant.mod \
            toleranc.mod \
            compar.mod \
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
            param.mod \
            param1.mod \
            geometry.mod \
            compar.mod \
            fldvar.mod \
            indices.mod \
            visc_s.mod \
            ghdtheory.mod \
            physprop.mod \
            run.mod \
            constant.mod \
            toleranc.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./GhdTheory/transport_coeff_ghd.f 
qmomk_allocate_arrays.$(OBJ_EXT) : ./qmomk/qmomk_allocate_arrays.f \
            param.mod \
            param1.mod \
            indices.mod \
            geometry.mod \
            physprop.mod \
            qmom_kinetic_equation.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_allocate_arrays.f 
qmomk_gas_drag.$(OBJ_EXT) : ./qmomk/qmomk_gas_drag.f \
            param.mod \
            param1.mod \
            parallel.mod \
            matrix.mod \
            scales.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            visc_g.mod \
            rxns.mod \
            run.mod \
            toleranc.mod \
            geometry.mod \
            indices.mod \
            is.mod \
            tau_g.mod \
            bc.mod \
            compar.mod \
            sendrecv.mod \
            discretelement.mod \
            qmom_kinetic_equation.mod \
            drag.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_gas_drag.f 
qmomk_init_bc.$(OBJ_EXT) : ./qmomk/qmomk_init_bc.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            bc.mod \
            qmom_kinetic_equation.mod \
            qmomk_quadrature.mod \
            qmomk_bc.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_init_bc.f 
qmomk_initial_conditions.$(OBJ_EXT) : ./qmomk/qmomk_initial_conditions.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            qmom_kinetic_equation.mod \
            qmomk_quadrature.mod \
            qmomk_parameters.mod \
            ic.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_initial_conditions.f 
qmomk_init_namelist.$(OBJ_EXT) : ./qmomk/qmomk_init_namelist.f \
            param1.mod \
            qmom_kinetic_equation.mod \
            qmomk/qmomknamelist.inc                                     
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_init_namelist.f 
qmomk_make_arrays.$(OBJ_EXT) : ./qmomk/qmomk_make_arrays.f \
            param1.mod \
            geometry.mod \
            funits.mod \
            compar.mod \
            qmom_kinetic_equation.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_make_arrays.f 
qmomk_read_restart.$(OBJ_EXT) : ./qmomk/qmomk_read_restart.f \
            param.mod \
            param1.mod \
            constant.mod \
            fldvar.mod \
            cont.mod \
            geometry.mod \
            indices.mod \
            run.mod \
            compar.mod \
            physprop.mod \
            qmom_kinetic_equation.mod \
            qmomk_quadrature.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_read_restart.f 
qmomk_set_bc.$(OBJ_EXT) : ./qmomk/qmomk_set_bc.f \
            param.mod \
            param1.mod \
            constant.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            compar.mod \
            indices.mod \
            bc.mod \
            qmom_kinetic_equation.mod \
            qmomk_quadrature.mod \
            qmomk_bc.mod \
            function.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_set_bc.f 
qmomk_time_march.$(OBJ_EXT) : ./qmomk/qmomk_time_march.f \
            param.mod \
            param1.mod \
            constant.mod \
            run.mod \
            output.mod \
            physprop.mod \
            fldvar.mod \
            geometry.mod \
            cont.mod \
            coeff.mod \
            tau_g.mod \
            tau_s.mod \
            visc_g.mod \
            visc_s.mod \
            funits.mod \
            vshear.mod \
            scalars.mod \
            drag.mod \
            rxns.mod \
            compar.mod \
            time_cpu.mod \
            is.mod \
            indices.mod \
            matrix.mod \
            sendrecv.mod \
            qmom_kinetic_equation.mod \
            qmomk_fluxes.mod \
            qmomk_quadrature.mod \
            qmomk_collision.mod \
            qmomk_parameters.mod \
            ur_facs.mod \
            function.inc                                                 \
            fun_avg1.inc                                                 \
            fun_avg2.inc                                                
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_time_march.f 
qmomk_write_restart.$(OBJ_EXT) : ./qmomk/qmomk_write_restart.f \
            param1.mod \
            qmom_kinetic_equation.mod \
            run.mod 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./qmomk/qmomk_write_restart.f 
get_values.$(OBJ_EXT) : ./thermochemical/get_values.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./thermochemical/get_values.f 
readTherm.$(OBJ_EXT) : ./thermochemical/readTherm.f 
	$(FORTRAN_CMD) $(FORT_FLAGS) ./thermochemical/readTherm.f 
