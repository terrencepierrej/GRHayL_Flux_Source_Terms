
// Free parameters related to physical system:

// Set free-parameter values.

const REAL domain_size    = 10.;
const REAL sinh_width     = 0.4;
const REAL sinhv2_const_dr= 0.05;
const REAL SymTP_bScale   = 1.0;

griddata.params.xmin = -domain_size, griddata.params.xmax = domain_size;
griddata.params.ymin = -domain_size, griddata.params.ymax = domain_size;
griddata.params.zmin = -domain_size, griddata.params.zmax = domain_size;

