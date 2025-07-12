# reload fit from cmdstan csv in fixtures folder
reload_fit <- function(coevfit, filename) {
  coevfit$fit <-
    cmdstanr::as_cmdstan_fit(
      testthat::test_path("fixtures", filename)
    )
  coevfit
}

# manually fix parameters in stan code
manually_fix_parameters <- function(scode) {
  scode %>%
    stringr::str_remove(
      stringr::fixed(
        paste0(
          "  vector<upper=0>[J] A_diag; // autoregressive terms of A\n",
          "  vector[num_effects - J] A_offdiag; // cross-lagged terms of A\n",
          "  vector<lower=0>[J] Q_sigma; // std deviation parameters of the ",
          "Q mat\n",
          "  vector[J] b; // SDE intercepts\n",
          "  array[N_tree] vector[J] eta_anc; // ancestral states\n"
        )
      )
    ) %>%
    stringr::str_replace(
      pattern = stringr::fixed(
        paste0(
          "  matrix[J,J] A = diag_matrix(A_diag); // selection matrix\n",
          "  matrix[J,J] Q = diag_matrix(Q_sigma^2); // drift matrix\n"
        )
      ),
      replacement = paste0(
        "  array[N_tree] vector[J] eta_anc;\n",
        "  vector[J] b = rep_vector(0.0, J);\n",
        "  matrix[J,J] A = diag_matrix(rep_vector(-0.5, J));\n",
        "  matrix[J,J] Q = diag_matrix(rep_vector(1.5, J));\n"
      )
    ) %>%
    stringr::str_replace(
      pattern = stringr::fixed(
        paste0(
          "  // fill off diagonal of A matrix\n",
          "  {\n",
          "    int ticker = 1;\n",
          "    for (i in 1:J) {\n",
          "      for (j in 1:J) {\n",
          "        if (i != j) {\n",
          "          if (effects_mat[i,j] == 1) {\n",
          "            A[i,j] = A_offdiag[ticker];\n",
          "            ticker += 1;\n",
          "          } else if (effects_mat[i,j] == 0) {\n",
          "            A[i,j] = 0;\n",
          "          }\n",
          "        }\n",
          "      }\n",
          "    }\n",
          "  }\n"
        )
      ),
      replacement = paste0(
        "  A[2,1] = 3;\n",
        "  A[1,2] = 0;\n",
        "  for (t in 1:N_tree) eta_anc[t] = rep_vector(0.0, J);\n"
      )
    ) %>%
    stringr::str_remove(stringr::fixed("  b ~ std_normal();\n")) %>%
    stringr::str_remove(stringr::fixed("    eta_anc[t] ~ std_normal();\n")) %>%
    stringr::str_remove(stringr::fixed("  A_offdiag ~ std_normal();\n")) %>%
    stringr::str_remove(stringr::fixed("  A_diag ~ std_normal();\n")) %>%
    stringr::str_remove(stringr::fixed("  Q_sigma ~ std_normal();\n"))
}
