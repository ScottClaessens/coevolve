functions {

  // Charles Driver's solver for the asymptotic Q matrix
  matrix ksolve (matrix A, matrix Q) {
    int d = rows(A);
    int d2 = (d * d - d) %/% 2;
    matrix [d + d2, d + d2] O;
    vector [d + d2] triQ;
    matrix[d,d] AQ;
    int z = 0;         // z is row of output
    for (j in 1:d) {   // for column reference of solution vector
      for (i in 1:j) { // and row reference...
        if (j >= i) {  // if i and j denote a covariance parameter
          int y = 0;   // start new output row
          z += 1;      // shift current output row down
          for (ci in 1:d) {   // for columns and
            for (ri in 1:d) { // rows of solution
              if (ci >= ri) { // when in upper tri (inc diag)
                y += 1;       // move to next column of output
                if (i == j) { // if output row is a diag element
                  if (ri == i) O[z, y] = 2 * A[ri, ci];
                  if (ci == i) O[z, y] = 2 * A[ci, ri];
                }
                if (i != j) { // if output row is not a diag element
                  //if column matches row, sum both A diags
                  if (y == z) O[z, y] = A[ri, ri] + A[ci, ci];
                  if (y != z) { // otherwise...
                    // if solution element is related to output row...
                    if (ci == ri) { // if solution element is variance
                      // if variance of solution corresponds to row
                      if (ci == i) O[z, y] = A[j, ci];
                      // if variance of solution corresponds to col
                      if (ci == j) O[z, y] = A[i, ci];
                    }
                    //if solution element is a related covariance
                    if (ci != ri && (ri == i || ri == j || ci == i || ci == j )) {
                      // for row 1,2 / 2,1 of output,
                      // if solution row ri 1 (match)
                      // and column ci 3, we need A[2,3]
                      if (ri == i) O[z, y] = A[j, ci];
                      if (ri == j) O[z, y] = A[i, ci];
                      if (ci == i) O[z, y] = A[j, ri];
                      if (ci == j) O[z, y] = A[i, ri];
                    }
                  }
                }
                if (is_nan(O[z, y])) O[z, y] = 0;
              }
            }
          }
        }
      }
    }
    z = 0; // get upper tri of Q
    for (j in 1:d) {
      for (i in 1:j) {
        z += 1;
        triQ[z] = Q[i, j];
      }
    }
    triQ = -O \ triQ; // get upper tri of asymQ
    z = 0; // put upper tri of asymQ into matrix
    for (j in 1:d) {
      for (i in 1:j) {
        z += 1;
        AQ[i, j] = triQ[z];
        if (i != j) AQ[j, i] = triQ[z];
      }
    }
    return AQ;
  }

  // Cayley-Hamilton matrix exponential for 2x2 matrices
  // exp(M) = c0*I + c1*M where c0,c1 derived from eigenvalues
  // Returns {c0, c1} coefficients
  vector cayley_hamilton_coeffs(real tr_M, real det_M) {
    real disc = tr_M * tr_M - 4 * det_M;
    real c0;
    real c1;
    if (disc >= 0) {
      real sq = sqrt(disc);
      real lam1 = (tr_M + sq) * 0.5;
      real lam2 = (tr_M - sq) * 0.5;
      real e1 = exp(lam1);
      real e2 = exp(lam2);
      if (abs(lam1 - lam2) > 1e-10) {
        c1 = (e1 - e2) / (lam1 - lam2);
        c0 = e1 - lam1 * c1;
      } else {
        c1 = e1;
        c0 = (1 - lam1) * e1;
      }
    } else {
      real alpha = tr_M * 0.5;
      real beta = sqrt(-disc) * 0.5;
      real ea = exp(alpha);
      real sb = sin(beta);
      real cb = cos(beta);
      c1 = ea * sb / beta;
      c0 = ea * (cb - alpha * sb / beta);
    }
    return [c0, c1]';
  }

  // Cayley-Hamilton matrix exponential for 3x3 matrices
  // exp(M) = c0*I + c1*M + c2*M^2
  // Returns {c0, c1, c2} coefficients
  vector cayley_hamilton_coeffs_3(matrix M) {
    int d = 3;
    // characteristic polynomial: lambda^3 - p*lambda^2 + q*lambda - r = 0
    real p = trace(M);
    real p2 = p * p;
    real tr_M2 = trace(M * M);
    real q = (p2 - tr_M2) * 0.5;
    real r = determinant(M);
    // depressed cubic: mu^3 + a*mu + b = 0  (substituting lambda = mu + p/3)
    real a = q - p2 / 3.0;
    real b = (p * q / 3.0) - (2.0 * p2 * p / 27.0) - r;
    // discriminant: delta = -(4a^3 + 27b^2)
    real delta = -(4.0 * a * a * a + 27.0 * b * b);
    real c0;
    real c1;
    real c2;
    real p3 = p / 3.0;
    if (delta > 1e-12 && a < 0) {
      // three distinct real eigenvalues — trigonometric method
      real m = sqrt(-a / 3.0);
      real theta = acos(
        fmin(fmax(-b / (2.0 * m * m * m), -1.0), 1.0)
      ) / 3.0;
      real lam1 = 2.0 * m * cos(theta) + p3;
      real lam2 = 2.0 * m * cos(theta - 2.0 * pi() / 3.0) + p3;
      real lam3 = 2.0 * m * cos(theta - 4.0 * pi() / 3.0) + p3;
      real e1 = exp(lam1);
      real e2 = exp(lam2);
      real e3 = exp(lam3);
      // solve 3x3 Vandermonde system for c0, c1, c2
      real d12 = lam1 - lam2;
      real d13 = lam1 - lam3;
      real d23 = lam2 - lam3;
      // Lagrange interpolation coefficients
      real L1 = e1 / (d12 * d13);
      real L2 = -e2 / (d12 * d23);
      real L3 = e3 / (d13 * d23);
      c2 = L1 + L2 + L3;
      c1 = -(L1 * (lam2 + lam3) + L2 * (lam1 + lam3) + L3 * (lam1 + lam2));
      c0 = L1 * lam2 * lam3 + L2 * lam1 * lam3 + L3 * lam1 * lam2;
    } else if (delta < -1e-12) {
      // one real + two complex conjugate eigenvalues
      real h = sqrt(b * b / 4.0 + a * a * a / 27.0);
      real A_cr = -b / 2.0 + h;
      real B_cr = -b / 2.0 - h;
      // cube roots (preserving sign)
      real S = (A_cr >= 0) ? pow(A_cr, 1.0/3.0) : -pow(-A_cr, 1.0/3.0);
      real T = (B_cr >= 0) ? pow(B_cr, 1.0/3.0) : -pow(-B_cr, 1.0/3.0);
      real lam1 = S + T + p3;  // real eigenvalue
      real alpha = -(S + T) / 2.0 + p3;  // real part of complex pair
      real beta = (S - T) * sqrt(3.0) / 2.0;  // imaginary part
      real e1 = exp(lam1);
      real ea = exp(alpha);
      real cb = cos(beta);
      real sb = sin(beta);
      real e_re = ea * cb;  // Re(exp(alpha + i*beta))
      real e_im = ea * sb;  // Im(exp(alpha + i*beta))
      // solve Vandermonde system with one real + complex conjugate pair
      // Using real arithmetic throughout
      real d_r = lam1 - alpha;  // lam1 - Re(complex)
      real denom = d_r * d_r + beta * beta;
      // from Lagrange: c2 = e1/((lam1-z)(lam1-z*)) + ...
      real L1_r = e1 / denom;
      // complex Lagrange terms combine to real:
      // L2 + L3 = 2*Re(exp(z) / ((z-lam1)(z-z*)))
      // where z = alpha + i*beta, z-lam1 = -d_r + i*beta, z-z* = 2i*beta
      // (z-lam1)(z-z*) = (-d_r + i*beta)(2i*beta) = -2*beta^2 - 2*d_r*beta*i
      // 1/(...) = (-2*beta^2 + 2*d_r*beta*i) / (4*beta^4 + 4*d_r^2*beta^2)
      //         = (-1 + d_r*i/beta) / (2*(beta^2 + d_r^2))
      //         = (-1 + d_r*i/beta) / (2*denom)
      // exp(z) * 1/(...) = (e_re + i*e_im)*(-1 + i*d_r/beta) / (2*denom)
      //                   = (-e_re - e_im*d_r/beta + i*(-e_im + e_re*d_r/beta)) / (2*denom)
      // 2*Re(...) = 2*(-e_re - e_im*d_r/beta) / (2*denom)
      //           = (-e_re - e_im*d_r/beta) / denom
      real L23_r = (-e_re - e_im * d_r / beta) / denom;
      c2 = L1_r + L23_r;
      // c1 from Lagrange: -sum of Li*(sum of other eigenvalues)
      // L1 contributes: -L1_r * 2*alpha  (sum of complex pair = 2*alpha)
      // L2+L3 contribute: -2*Re(L2*(lam1 + z*))
      // lam1 + z* = lam1 + alpha - i*beta
      // L2*(lam1+z*): need Re(L2_complex * (lam1+alpha-i*beta))
      // L2_complex = (-e_re - e_im*d_r/beta + i*(-e_im + e_re*d_r/beta)) / (2*denom)
      // Let Lr = (-e_re - e_im*d_r/beta)/(2*denom), Li = (-e_im + e_re*d_r/beta)/(2*denom)
      // L2*(lam1+alpha-i*beta) = (Lr+iLi)*(lam1+alpha-i*beta)
      // Re = Lr*(lam1+alpha) + Li*beta
      real Lr = (-e_re - e_im * d_r / beta) / (2.0 * denom);
      real Li = (-e_im + e_re * d_r / beta) / (2.0 * denom);
      // 2*Re(L2*(lam1+z*)) = 2*(Lr*(lam1+alpha) + Li*beta)
      c1 = -(L1_r * 2.0 * alpha + 2.0 * (Lr * (lam1 + alpha) + Li * beta));
      // c0 from Lagrange: sum of Li*(product of other eigenvalues)
      // L1 contributes: L1_r * |z|^2 = L1_r * (alpha^2 + beta^2)
      // L2+L3 contribute: 2*Re(L2 * lam1 * z*)
      // L2*lam1*z* = (Lr+iLi)*lam1*(alpha-i*beta)
      // Re = lam1*(Lr*alpha + Li*beta)
      real z_mod2 = alpha * alpha + beta * beta;
      c0 = L1_r * z_mod2
           + 2.0 * lam1 * (Lr * alpha + Li * beta);
    } else {
      // near-degenerate: all eigenvalues close together
      // Taylor: exp(M) = exp(p/3)*(I + D + D²/2 + D³/6) where D = M - p/3*I
      // By Cayley-Hamilton: D³ = -a*M + (r - p3³)*I (no M² term)
      real e_avg = exp(p3);
      real p3_3 = p3 * p3 * p3;
      c0 = e_avg * (1.0 - p3 + p3 * p3 / 2.0 + (r - p3_3) / 6.0);
      c1 = e_avg * (1.0 - p3 - a / 6.0);
      c2 = e_avg * 0.5;
    }
    return [c0, c1, c2]';
  }

  // return number of matches of y in vector x
  int num_matches(vector x, real y) {
    int n = 0;
    for (i in 1:rows(x))
      if (x[i] == y)
        n += 1;
    return n;
  }

  // return indices in vector x where x == y
  array[] int which_equal(vector x, real y) {
    array [num_matches(x, y)] int match_positions;
    int pos = 1;
    for (i in 1:rows(x)) {
      if (x[i] == y) {
        match_positions[pos] = i;
        pos += 1;
      }
    }
    return match_positions;
  }

  {{#approximate_gps}}
  {{#exp_quad}}
  // spectral density function of Gaussian process with exp_quad kernel
  vector spd_gp_exp_quad(data array[] vector x, real sdgp, real lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    real constant = square(sdgp) * sqrt(2 * pi())^D;
    vector[NB] out;
    real neg_half_lscale2 = -0.5 * square(lscale);
    constant = constant * lscale^D;
    for (m in 1:NB) {
      out[m] = constant * exp(neg_half_lscale2 * dot_self(x[m]));
    }
    return out;
  }
  {{/exp_quad}}
  {{#exponential}}
  // spectral density function of Gaussian process with exponential kernel
  vector spd_gp_exponential(data array[] vector x, real sdgp, real lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    real constant = square(sdgp) *
      (2^D * pi()^(D / 2.0) * tgamma((D + 1.0) / 2)) / sqrt(pi());
    real expo = -(D + 1.0) / 2;
    vector[NB] out;
    real lscale2 = square(lscale);
    constant = constant * lscale^D;
    for (m in 1:NB) {
      out[m] = constant * (1 + lscale2 * dot_self(x[m]))^expo;
    }
    return out;
  }
  {{/exponential}}
  {{#matern32}}
  // spectral density function of Gaussian process with matern32 kernel
  vector spd_gp_matern32(data array[] vector x, real sdgp, real lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    real constant = square(sdgp) *
      (2^D * pi()^(D / 2.0) * tgamma((D + 3.0) / 2) * 3^(3.0 / 2)) /
      (0.5 * sqrt(pi()));
    real expo = -(D + 3.0) / 2;
    vector[NB] out;
    real lscale2 = square(lscale);
    constant = constant * lscale^D;
    for (m in 1:NB) {
      out[m] = constant * (3 + lscale2 * dot_self(x[m]))^expo;
    }
    return out;
  }
  {{/matern32}}
  {{/approximate_gps}}

}
