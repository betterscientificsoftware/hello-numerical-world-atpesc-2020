#include "heat.H"

// Double class' statics
int         Double::nadds  = 0;
int         Double::nmults = 0;
int         Double::ndivs  = 0;
std::size_t Double::nbytes = 0;

// Command-line argument variables
int noout        = 0;
int savi         = 0;
int outi         = 100;
int save         = 0;
char const *runame = "heat_results";
char const *alg  = "ftcs";
char const *prec = "double";
char const *ic   = "const(1)";
Double lenx      = 1.0;
Double alpha     = 0.2;
Double dt        = 0.004;
Double dx        = 0.1;
Double bc0       = 0;
Double bc1       = 1;
Double maxt      = 2.0;
Double min_change = 1e-8*1e-8;

// Various arrays of numerical data
Double *curr           = 0; // current solution
Double *last           = 0; // last solution
Double *exact          = 0; // exact solution (when available)
Double *change_history = 0; // solution l2norm change history
Double *error_history  = 0; // solution error history (when available)
Double *cn_Amat        = 0; // A matrix for Crank-Nicholson

// Number of points in space, x, and time, t.
int Nx = (int) (lenx/dx);
int Nt = (int) (maxt/dt);

// Utilities
Double
l2_norm(int n, Double const *a, Double const *b)
{
    int i;
    Double sum = 0;
    for (i = 0; i < n; i++)
    {
        Double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sum / n;
}

void
copy(int n, Double *dst, Double const *src)
{
    int i;
    for (i = 0; i < n; i++)
        dst[i] = src[i];
}

void
write_array(int t, int n, Double dx, Double const *a)
{
    int i;
    char fname[128];
    char vname[64];
    FILE *outf;

    if (noout) return;

    if (t == TSTART)
    {
        snprintf(fname, sizeof(fname), "%s/%s_soln_00000.curve", runame, runame);
        snprintf(vname, sizeof(vname), "Temperature");
    }
    else if (t == TFINAL)
    {
        snprintf(fname, sizeof(fname), "%s/%s_soln_final.curve", runame, runame);
        snprintf(vname, sizeof(vname), "Temperature");
    }
    else if (t == RESIDUAL)
    {
        snprintf(fname, sizeof(fname), "%s/%s_change.curve", runame, runame);
        snprintf(vname, sizeof(vname), "%s/%s_l2_change", runame, runame);
    }
    else if (t == ERROR)
    {
        snprintf(fname, sizeof(fname), "%s/%s_error.curve", runame, runame);
        snprintf(vname, sizeof(vname), "%s/%s_l2", runame, runame);
    }
    else
    {
        if (a == exact)
        {
            snprintf(fname, sizeof(fname), "%s/%s_exact_%05d.curve", runame, runame, t);
            snprintf(vname, sizeof(vname), "exact_temperature");
        } 
        else
        {
            snprintf(fname, sizeof(fname), "%s/%s_soln_%05d.curve", runame, runame, t);
            snprintf(vname, sizeof(vname), "Temperature");
        }
    }


    outf = fopen(fname,"w");
    fprintf(outf, "# %s\n", vname);
    for (i = 0; i < n; i++)
        fprintf(outf, "%8.4g %8.4g\n", i*((double)dx), (double) a[i]);
    fclose(outf);
}

void
set_initial_condition(int n, Double *a, Double dx, char const *ic)
{
    int i;
    Double x;

    if (!strncmp(ic, "const(", 6)) /* const(val) */
    {
        double cval = strtod(ic+6, 0);
        for (i = 0; i < n; i++)
            a[i] = cval;
    }
    else if (!strncmp(ic, "step(", 5)) /* step(left,xmid,right) */
    {
        char *p;
        double left = strtod(ic+5, &p);
        double xmid = strtod(p+1, &p);
        double right = strtod(p+1, 0);
        for (i = 0, x = 0; i < n; i++, x+=dx)
        {
            if (x < xmid) a[i] = left;
            else          a[i] = right;
        }
    }
    else if (!strncmp(ic, "ramp(", 5)) /* ramp(left,right) */
    {
        char *p;
        double left = strtod(ic+5, &p);
        double right = strtod(p+1, 0);
        double dv = (right-left)/(n-1);
        for (i = 0, x = left; i < n; i++, x+=dv)
            a[i] = x;
    }
    else if (!strncmp(ic, "rand(", 5)) /* rand(seed,base,amp) */
    {
        char *p, *ep;
        int seed = (int) strtol(ic+5,&p,10);
        double base = strtod(p+1, &p);
        double amp = strtod(p+1, 0);
        const double maxr = ((long long)1<<31)-1;
        srandom(seed);
        for (i = 0; i < n; i++)
            a[i] = base + amp * (2*random()/maxr - 1);
    }
    else if (!strncmp(ic, "sin(Pi*x)", 9)) /* sin(PI*x) */
    {
        for (i = 0, x = 0; i < n; i++, x+=dx)
            a[i] = sin(M_PI*x);
    }
    else if (!strncmp(ic, "spikes(", 7)) /* spikes(Const,Amp,Loc,Amp,Loc,...) */
    {
        char *next;
        double cval = strtod(ic+7, &next);
        char const *p = next;
        for (i = 0, x = 0; i < n; i++)
            a[i] = cval;
        while (*p != ')')
        {
            char *ep_amp, *ep_idx;
            double amp = strtod(p+1, &ep_amp);
            int idx = (int) strtod(ep_amp+1, &ep_idx);
            assert(idx<n);
            a[idx] = amp;
            p = ep_idx;
        }

    }

    write_array(TSTART, Nx, dx, a);
}

static void
r83_np_fa(int n, Double *a)
{
    int i;

    for ( i = 1; i <= n-1; i++ )
    {
        assert ( a[1+(i-1)*3] != 0.0 );

        // Store the multiplier in L.
        a[2+(i-1)*3] = a[2+(i-1)*3] / a[1+(i-1)*3];

        // Modify the diagonal entry in the next column.
        a[1+i*3] = a[1+i*3] - a[2+(i-1)*3] * a[0+i*3];
    }

    assert( a[1+(n-1)*3] != 0.0 );
}

void
initialize_crankn(int n,
    Double alpha, Double dx, Double dt,
    Double **_cn_Amat)
{
    int i;
    Double const w = alpha * dt / dx / dx;

    // Build a tri-diagonal matrix
    Double *cn_Amat = new Double[3*n]();

    cn_Amat[0+0*3] = 0.0;
    cn_Amat[1+0*3] = 1.0;
    cn_Amat[0+1*3] = 0.0;

    for ( i = 1; i < n - 1; i++ )
    {
        cn_Amat[2+(i-1)*3] =           - w;
        cn_Amat[1+ i   *3] = 1.0 + 2.0 * w;
        cn_Amat[0+(i+1)*3] =           - w;
    }
        
    cn_Amat[2+(n-2)*3] = 0.0;
    cn_Amat[1+(n-1)*3] = 1.0;
    cn_Amat[2+(n-1)*3] = 0.0;

    // Factor the matrix.
    r83_np_fa(n, cn_Amat);

    // Return the generated matrix
    *_cn_Amat = cn_Amat;
}

// Licensing: This code is distributed under the GNU LGPL license. 
// Modified: 30 May 2009 Author: John Burkardt
// Modified by Mark C. Miller, miller86@llnl.gov, July 23, 2017
static void 
r83_np_sl ( int n, Double const *a_lu, Double const *b, Double *x)
{
    int i;

    for ( i = 0; i < n; i++ )
        x[i] = b[i];

    // Solve L * Y = B.
    for ( i = 1; i < n; i++ )
        x[i] = x[i] - a_lu[2+(i-1)*3] * x[i-1];

    // Solve U * X = Y.
    for ( i = n; 1 <= i; i-- )
    {
        x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
        if ( 1 < i )
            x[i-2] = x[i-2] - a_lu[0+(i-1)*3] * x[i-1];
    }
}

bool
update_solution_crankn(int n,
    Double *curr, Double const *last,
    Double const *cn_Amat,
    Double bc_0, Double bc_1)
{
    // Do the solve
    r83_np_sl (n, cn_Amat, last, curr);
    curr[0] = bc_0;
    curr[n-1] = bc_1;

    return true;
}


static char clargs[2048];

#define HANDLE_ARG(VAR, TYPE, STYLE, HELP) \
{ \
    char const *style = #STYLE; \
    char const *q = style[1]=='s'?"\"":""; \
    void *valp = (void*) &VAR; \
    int const len = strlen(#VAR)+1; \
    std::stringstream strmvar; \
    for (i = 1; i < argc; i++) \
    {\
        int valid_style = style[1]=='d'||style[1]=='g'||style[1]=='s'; \
        if (strncmp(argv[i], #VAR"=", len)) \
            continue; \
        assert(valid_style); \
	if (strlen(argv[i]+len)) \
        {\
            if      (style[1] == 'd') /* int */ \
                *((int*) valp) = (int) strtol(argv[i]+len,0,10); \
            else if (style[1] == 'g') /* double */ \
                *((Double*) valp) = (double) strtod(argv[i]+len,0); \
            else if (style[1] == 's') /* char* */ \
                *((char**) valp) = (char*) strdup(argv[i]+len); \
        }\
    }\
    strmvar << VAR; \
    if (help) \
    {\
        char tmp[256]; \
        int len = snprintf(tmp, sizeof(tmp), "    %s=%s%s%s", \
            #VAR, q, strmvar.str().c_str(), q);\
        snprintf(tmp, sizeof(tmp), "%s (%s)", #HELP, #TYPE); \
        fprintf(stderr, "    %s=%s%s%s%*s\n", \
            #VAR, q, strmvar.str().c_str(), q, 80-len, tmp);\
    }\
    else \
    { \
        char tmp[64]; \
        fprintf(stderr, "    %s=%s%s%s\n", \
            #VAR, q, strmvar.str().c_str(), q);\
        snprintf(tmp, sizeof(tmp), "    %s=%s%s%s\n", \
            #VAR, q, strmvar.str().c_str(), q);\
        strcat(clargs, tmp); \
    } \
}

void
process_args(int argc, char **argv)
{
    int i;
    int help = 0;

    // quick pass for 'help' anywhere on command line
    clargs[0] ='\0';
    for (i = 0; i < argc && !help; i++)
        help = 0!=strcasestr(argv[i], "help");
    
    if (help)
        fprintf(stderr, "Usage: ./heat <arg>=<value> <arg>=<value>...\n");

    HANDLE_ARG(runame, char*, %s, name to give run and results dir);
    HANDLE_ARG(prec, char*, %s, precision half|float|double|quad);
    HANDLE_ARG(alpha, double, %g, material thermal diffusivity (sq-meters/second));
    HANDLE_ARG(lenx, double, %g, material length (meters));
    HANDLE_ARG(dx, double, %g, x-incriment. Best if lenx/dx==int. (meters));
    HANDLE_ARG(dt, double, %g, t-incriment (seconds));
    HANDLE_ARG(maxt, double, %g, >0:max sim time (seconds) | <0:min l2 change in soln);
    HANDLE_ARG(bc0, double, %g, boundary condition @ x=0: u(0,t) (Kelvin));
    HANDLE_ARG(bc1, double, %g, boundary condition @ x=lenx: u(lenx,t) (Kelvin));
    HANDLE_ARG(ic, char*, %s, initial condition @ t=0: u(x,0) (Kelvin));
    HANDLE_ARG(alg, char*, %s, algorithm ftcs|upwind15|crankn);
    HANDLE_ARG(savi, int, %d, save every i-th solution step);
    HANDLE_ARG(save, int, %d, save error in every saved solution);
    HANDLE_ARG(outi, int, %d, output progress every i-th solution step);
    HANDLE_ARG(noout, int, %d, disable all file outputs);

    if (help)
    {
        fprintf(stderr, "Examples...\n");
        fprintf(stderr, "    ./heat dx=0.01 dt=0.0002 alg=ftcs\n");
        fprintf(stderr, "    ./heat dx=0.1 bc0=273 bc1=273 ic=\"spikes(273,5,373)\"\n");
        exit(1);
    }

    // Handle possible termination by change threshold criterion
    if (maxt < 0)
    {
        min_change = -maxt * -maxt;
        maxt = INT_MAX;
    }

    // Handle output results dir creation and save of command-line
    if (access(runame, F_OK) == 0)
    {
        fprintf(stderr, "An entry \"%s\" already exists\n", runame);
        exit(1);
    } 

    // Make the output dir and save clargs there too
    mkdir(runame, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
    char fname[128];
    sprintf(fname, "%s/clargs.out", runame);
    FILE *outf = fopen(fname, "w");
    fprintf(outf, "%s", clargs);
    fclose(outf);
}

void 
compute_exact_solution(int n, Double *a, Double dx, char const *ic,
    Double alpha, Double t, Double bc0, Double bc1)
{
    int i;
    Double x;
    
    // For any time t for Sin(Pi*x) initial condition
    // and zero boundary condition
    if (bc0 == 0 && bc1 == 0 && !strncmp(ic, "sin(Pi*x)", 9))
    {
        for (i = 0, x = 0; i < n; i++, x+=dx)
            a[i] = sin(M_PI*x)*exp(-alpha*M_PI*M_PI*t);
    }
    // For any time t for constant initial condition
    // and zero boundary condition
    else if (bc0 == 0 && bc1 == 0 && !strncmp(ic, "const(", 6))
    {
        Double cval = strtod(ic+6, 0);
        for (i = 0, x = 0; i < n; i++, x+=dx)
        {
            int n;
            Double fsum = 0;

            // Brute force sum first 1000 terms of Fourier series
            for (n = 1; n < 1000; n++)
            {
                Double coeff = 2*cval*(1-pow(-1.0,(double)n))/(n*M_PI);
                Double func = sin(n*M_PI*x)*exp(((double)-alpha)*n*n*M_PI*M_PI*((double)t));
                fsum += coeff * func;
            }
            a[i] = fsum;
        }
    }
    // For t>>0 (steady state) for any initial and boundary conditions
    else
    {
        for (i = 0, x = 0; i < n; i++, x+=dx)
            a[i] = bc0 + (bc1-bc0)*x;
    }
}

bool                          // false if unstable, true otherwise
update_solution_ftcs(
    int n,                    // number of samples
    Double *uk1,              // new array of u(x,k+1) to compute/return
    Double const *uk0,        // old/last array u(x,k) of samples computed
    Double alpha,             // thermal diffusivity
    Double dx, Double dt,     // spacing in space, x, and time, t.
    Double bc0, Double bc1)   // boundary conditions @ x=0 & x=Lx
{
    Double r = alpha * dt / (dx * dx);

    // sanity check for stability
    if (r > 0.5) return false; 

    // FTCS update algorithm
    for (int i = 1; i < n-1; i++)
        uk1[i] = r*uk0[i+1] + (1-2*r)*uk0[i] + r*uk0[i-1];

    // enforce boundary conditions
    uk1[0  ] = bc0;
    uk1[n-1] = bc1;

    return true;
}

bool
update_solution_upwind15(int n, Double *curr, Double const *last,
    Double alpha, Double dx, Double dt,
    Double bc_0, Double bc_1)
{
    Double const f2 = 1.0/24;
    Double const f1 = 1.0/6;
    Double const f0 = 1.0/4;
    Double const k = alpha * alpha * dt / (dx * dx);
    Double const k2 = k*k;

    int i;
    curr[0  ] = bc_0;
    curr[1  ] = last[1  ] + k * (last[0  ] - 2 * last[1  ] + last[2  ]);
    curr[n-2] = last[n-2] + k * (last[n-3] - 2 * last[n-2] + last[n-1]);
    curr[n-1] = bc_1;
    for (i = 2; i < n-2; i++)
        curr[i] =  f2*(12*k2  -2*k    )*last[i-2]
                  +f2*(12*k2  -2*k    )*last[i+2]
                  -f1*(12*k2  -8*k    )*last[i-1]
                  -f1*(12*k2  -8*k    )*last[i+1]
                  +f0*(12*k2 -10*k  +4)*last[i  ];

    return true;
}


static void
initialize(void)
{
    Nx = (int) (lenx/dx)+1;
    Nt = (int) (maxt/dt);
    dx = lenx/(Nx-1);

    curr = new Double[Nx]();
    last = new Double[Nx]();
    if (save)
    {
        exact = new Double[Nx]();
        change_history = new Double[Nx]();
        error_history = new Double[Nx]();
    }

    assert(strncmp(alg, "ftcs", 4)==0 ||
           strncmp(alg, "upwind15", 8)==0 ||
           strncmp(alg, "crankn", 6)==0);

#ifdef HAVE_FEENABLEEXCEPT
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
#endif

    if (!strncmp(alg, "crankn", 6))
        initialize_crankn(Nx, alpha, dx, dt, &cn_Amat);

    /* Initial condition */
    set_initial_condition(Nx, last, dx, ic);
}

int finalize(int ti, Double maxt, Double change)
{
    int retval = 0;

    write_array(TFINAL, Nx, dx, curr);
    if (save)
    {
        write_array(RESIDUAL, ti, dt, change_history);
        write_array(ERROR, ti, dt, error_history);
    }

    if (outi)
    {
        printf("Iteration %04d: last change l2=%g\n", ti, (double) change);
        printf("Counts: %s\n", Double::counts_string());
    }

    delete [] curr;
    delete [] last;
    if (exact) delete [] exact;
    if (change_history) delete [] change_history;
    if (error_history) delete [] error_history;
    if (cn_Amat) delete [] cn_Amat;
    if (strncmp(alg, "ftcs", 4)) free((void*)alg);
    if (strncmp(prec, "double", 6)) free((void*)prec);
    if (strncmp(ic, "const(1)", 8)) free((void*)ic);

    return retval;
}

static bool
update_solution()
{
    if (!strcmp(alg, "ftcs"))
        return update_solution_ftcs(Nx, curr, last, alpha, dx, dt, bc0, bc1);
    else if (!strcmp(alg, "upwind15"))
        return update_solution_upwind15(Nx, curr, last, alpha, dx, dt, bc0, bc1);
    else if (!strcmp(alg, "crankn"))
        return update_solution_crankn(Nx, curr, last, cn_Amat, bc0, bc1);
    return false;
}

static Double
update_output_files(int ti)
{
    Double change;

    if (ti>0 && save)
    {
        compute_exact_solution(Nx, exact, dx, ic, alpha, ti*dt, bc0, bc1);
        if (savi && ti%savi==0)
            write_array(ti, Nx, dx, exact);
    }

    if (ti>0 && savi && ti%savi==0)
        write_array(ti, Nx, dx, curr);

    change = l2_norm(Nx, curr, last);
    if (save)
    {
        change_history[ti] = change;
        error_history[ti] = l2_norm(Nx, curr, exact);
    }

    return change;
}

int main(int argc, char **argv)
{
    int ti;
    Double change;

    // Read command-line args and set values
    process_args(argc, argv);

    // Allocate arrays and set initial conditions
    initialize();

    // Iterate to max iterations or solution change is below threshold
    for (ti = 0; ti*dt < maxt; ti++)
    {
        // compute the next solution step
        if (!update_solution())
        {
            fprintf(stderr, "Solution criteria violated. Make better choices\n");
            exit(1);
        }

        // compute amount of change in solution
        change = update_output_files(ti);

        // Handle possible termination by change threshold
        if (maxt == INT_MAX && change < min_change)
        {
            printf("Stopped after %06d iterations for threshold %g\n",
                ti, (double) change);
            break;
        }

        // Output progress
        if (outi && ti%outi==0)
            printf("Iteration %04d: last change l2=%g\n", ti, (double) change);

        // Copy current solution to last
        copy(Nx, last, curr);
    }

    // Delete storage and output final results
    return finalize(ti, maxt, change);
}
