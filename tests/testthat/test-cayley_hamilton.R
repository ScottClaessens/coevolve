# Helper: get rendered functions block and append test wrapper functions
build_stan_code <- function(extra_functions) {
  fn_block <- coevolve:::write_functions_block(
    lon_lat = NULL, dist_k = NA, dist_cov = "exp_quad"
  )
  # Insert extra functions before the closing brace
  sub("\\}\\s*$", paste0("\n", extra_functions, "\n}\n"), fn_block)
}

test_that("Cayley-Hamilton matrix_exp matches Stan matrix_exp for 2x2", {
  skip_on_cran()
  stan_code <- build_stan_code("
  matrix ch_exp_2(matrix A, real t) {
    real tr_At = (A[1,1] + A[2,2]) * t;
    real det_At = (A[1,1]*A[2,2] - A[1,2]*A[2,1]) * t * t;
    vector[2] ch = cayley_hamilton_coeffs(tr_At, det_At);
    matrix[2,2] result = ch[2] * t * A;
    result[1,1] += ch[1];
    result[2,2] += ch[1];
    return result;
  }

  matrix ref_exp_2(matrix A, real t) {
    return matrix_exp(A * t);
  }
")
  mod <- cmdstanr::cmdstan_model(
    stan_file = cmdstanr::write_stan_file(stan_code),
    compile = TRUE,
    force_recompile = TRUE
  )
  mod$expose_functions(global = TRUE)

  test_cases <- list(
    # real distinct eigenvalues
    list(A = matrix(c(-0.5, 0.3, 0.2, -0.8), 2, 2), t = 0.5),
    list(A = matrix(c(-1.0, 0.5, 0.5, -1.0), 2, 2), t = 2.0),
    # near-repeated eigenvalues
    list(A = matrix(c(-0.5, 0.0, 0.0, -0.5), 2, 2), t = 1.0),
    # complex eigenvalues
    list(A = matrix(c(-0.1, 2.0, -2.0, -0.1), 2, 2), t = 0.3),
    list(A = matrix(c(-0.5, 1.5, -1.0, -0.5), 2, 2), t = 1.0),
    # very small time
    list(A = matrix(c(-0.5, 0.3, 0.2, -0.8), 2, 2), t = 0.001),
    # larger values
    list(A = matrix(c(-2.0, 1.0, 0.5, -3.0), 2, 2), t = 0.1)
  )

  for (i in seq_along(test_cases)) {
    tc <- test_cases[[i]]
    ch_result <- ch_exp_2(tc$A, tc$t)
    ref_result <- ref_exp_2(tc$A, tc$t)
    expect_equal(
      ch_result, ref_result,
      tolerance = 1e-10,
      info = sprintf("2x2 test case %d", i)
    )
  }
})

test_that("Cayley-Hamilton matrix_exp matches Stan matrix_exp for 3x3", {
  skip_on_cran()
  stan_code <- build_stan_code("
  matrix ch_exp_3(matrix A, real t) {
    matrix[3,3] M = A * t;
    vector[3] ch = cayley_hamilton_coeffs_3(M);
    matrix[3,3] A2 = A * A;
    real c1t = ch[2] * t;
    real c2t2 = ch[3] * t * t;
    matrix[3,3] result = c1t * A + c2t2 * A2;
    result[1,1] += ch[1];
    result[2,2] += ch[1];
    result[3,3] += ch[1];
    return result;
  }

  matrix ref_exp_3(matrix A, real t) {
    return matrix_exp(A * t);
  }
")
  mod <- cmdstanr::cmdstan_model(
    stan_file = cmdstanr::write_stan_file(stan_code),
    compile = TRUE,
    force_recompile = TRUE
  )
  mod$expose_functions(global = TRUE)

  test_cases <- list(
    # three distinct real eigenvalues
    list(
      A = matrix(c(-1, 0.3, 0.1, 0.2, -0.8, 0.4, 0.1, 0.2, -0.6), 3, 3),
      t = 0.5
    ),
    # complex eigenvalues likely
    list(
      A = matrix(c(-0.5, 1.5, 0.1, -1.0, -0.5, 0.2, 0.3, 0.1, -0.3), 3, 3),
      t = 1.0
    ),
    # near-degenerate (diagonal)
    list(
      A = matrix(c(-0.5, 0, 0, 0, -0.5, 0, 0, 0, -0.5), 3, 3),
      t = 1.0
    ),
    # small time
    list(
      A = matrix(c(-1, 0.5, 0.2, 0.3, -1.5, 0.1, 0.4, 0.6, -2), 3, 3),
      t = 0.01
    ),
    # larger time
    list(
      A = matrix(c(-0.3, 0.1, 0.1, 0.2, -0.4, 0.1, 0.1, 0.2, -0.5), 3, 3),
      t = 3.0
    )
  )

  for (i in seq_along(test_cases)) {
    tc <- test_cases[[i]]
    ch_result <- ch_exp_3(tc$A, tc$t)
    ref_result <- ref_exp_3(tc$A, tc$t)
    expect_equal(
      ch_result, ref_result,
      tolerance = 1e-8,
      info = sprintf("3x3 test case %d", i)
    )
  }
})
