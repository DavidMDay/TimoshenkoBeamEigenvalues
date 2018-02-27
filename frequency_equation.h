enum boundary_condition { freefree, clampedclamped, clampedfree };
double frequency_equation(double a,
                          double b,
                          double gamma2,
                          bool is_sub_critical,
                          boundary_condition bc,
                          double* F_a,
                          double* F_b);


double euler_wave_number( unsigned id, boundary_condition bc );

