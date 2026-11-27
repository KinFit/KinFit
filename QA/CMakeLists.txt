# KinFit QA test suite (CTest).
#
#   cmake -S . -B build            # tests are on by default (KINFIT_BUILD_TESTS)
#   cmake --build build
#   ctest --test-dir build --output-on-failure
#
# Every test case is its own CTest entry: <executable>.<case>.
# Cases listed under DISABLED encode the *correct* behaviour of a known bug;
# CTest reports them as "Not Run (Disabled)". Remove a case from DISABLED once
# the bug is fixed (see TODO.md).

# kinfit_add_test(<exe> CASES <case>... [DISABLED <case>...] [LABELS <l>...])
function(kinfit_add_test exe)
  cmake_parse_arguments(T "" "" "CASES;DISABLED;LABELS" ${ARGN})
  add_executable(${exe} ${exe}.cxx)
  target_link_libraries(${exe} PRIVATE KinFit)
  target_include_directories(${exe} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  foreach(case ${T_CASES} ${T_DISABLED})
    add_test(NAME ${exe}.${case} COMMAND ${exe} ${case}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${exe}.${case} PROPERTIES LABELS "${T_LABELS}")
  endforeach()
  foreach(case ${T_DISABLED})
    set_tests_properties(${exe}.${case} PROPERTIES DISABLED TRUE)
  endforeach()
endfunction()

kinfit_add_test(test_particle
  LABELS core io
  CASES construct_from_R_Z construct_from_point_gives_PCA sign_of_R
        covariance_5x5 covariance_wrong_size_ignored io_roundtrip_tclonesarray)

kinfit_add_test(test_vertexfinder
  LABELS core vertex
  CASES lines_exact_intersection lines_parallel_tracks_fail_cleanly
        helix_exact_vertex helix_low_pt_loopers helix_negative_field
        helix_zero_field_equals_lines helix_neutral_tracks_are_lines
        helix_refit_after_setters helix_parallel_tracks_fail_cleanly)

kinfit_add_test(test_kinfitter
  LABELS core fitter
  CASES exact_4C exact_Mass exact_Vertex exact_MassVtx exact_MissingParticle
        exact_3C
        pulls_4C pulls_Mass pulls_Vertex pulls_MissingParticle pulls_3C
        jacobian_4C jacobian_Mass jacobian_Vertex jacobian_MissingParticle
        jacobian_3C
        pull_of_unconstrained_parameter mass_exclusive_candidates
        missing_particle_single_track_fails_cleanly
        single_iteration_returns_true
  DISABLED
        # MissingMass: f uses P_init + sum p_i, the Jacobian P_init - sum p_i
        exact_MissingMass pulls_MissingMass jacobian_MissingMass
        # MassVtx: vertex row of f has the opposite sign of its Jacobian
        pulls_MassVtx jacobian_MassVtx
        # relative convergence criterion divides by measured values (R=0)
        converges_with_zero_parameters)

kinfit_add_test(test_decaycandfinder
  LABELS core vertex
  CASES momentum_from_daughter_energies direction_from_vertices
        covariance_accessors_agree covariance_symmetric_positive
        usable_as_3C_mother
  DISABLED
        # flight line anchored at (decay - primary) instead of decay vertex
        R_Z_of_flight_line)

kinfit_add_test(test_coordinateconversion
  LABELS core
  CASES sph_to_cart cart_to_sph_forward
  DISABLED
        # theta = atan(pT/pz) and sign of R wrong for pz < 0
        cart_to_sph_backward
        # cross terms added as if the inputs were fully correlated
        errors_sph_to_cart)

kinfit_add_test(test_analyzer
  LABELS analysis io
  CASES builder_selects_true_combination builder_no_double_use
        builder_unknown_task_no_abort
        analyzer_mass_fit_output analyzer_event_limit
        analyzer_processes_first_event)

# The self-contained examples double as smoke tests.
if(KINFIT_BUILD_EXAMPLES)
  add_test(NAME example.toy_lambda_fit COMMAND toy_lambda_fit 200)
  add_test(NAME example.vertex_constB COMMAND vertex_constB 3)
  set_tests_properties(example.toy_lambda_fit example.vertex_constB
                       PROPERTIES LABELS example)
endif()
