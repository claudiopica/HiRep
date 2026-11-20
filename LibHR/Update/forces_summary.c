#include "libhr_core.h"
#include "update.h"
#include "inverters.h"
#include "memory.h"

static scalar_field *la = NULL; /* local action field for Metropolis test */

static void force_to_zero(mon_type type, monomial_data const *data) {
    switch (type) {
    case PureGauge: {
        mon_pg_par *par_pg = (mon_pg_par *)(data->par);
        zero(*par_pg->force_par.momenta);
        break;
    }

    case LuscherWeisz: {
        mon_lw_par *par_lw = (mon_lw_par *)(data->par);
        zero(*par_lw->force_par.momenta);
        break;
    }

    case FourFermion: {
        /* Do nothing */
        break;
    }

    case HMC: {
        mon_hmc_par *par_hmc = (mon_hmc_par *)(data->par);
        zero(*par_hmc->fpar.momenta);
        break;
    }

    case RHMC: {
        mon_rhmc_par *par_rhmc = (mon_rhmc_par *)(data->par);
        zero(*par_rhmc->fpar.momenta);
        break;
    }

    case TM: {
        mon_tm_par *par_tm = (mon_tm_par *)(data->par);
        zero(*par_tm->fpar.momenta);
        break;
    }

    case TM_alt: {
        mon_tm_par *par_tm_alt = (mon_tm_par *)(data->par);
        zero(*par_tm_alt->fpar.momenta);
        break;
    }

    case Hasenbusch: {
        mon_hasenbusch_par *par_hb = (mon_hasenbusch_par *)(data->par);
        zero(*par_hb->fpar.momenta);
        break;
    }

    case Hasenbusch_tm: {
        mon_hasenbusch_tm_par *par_hb_tm = (mon_hasenbusch_tm_par *)(data->par);
        zero(*par_hb_tm->fpar.momenta);
        break;
    }

    case Hasenbusch_tm_alt: {
        mon_hasenbusch_tm_par *par_hb_tm_alt = (mon_hasenbusch_tm_par *)(data->par);
        zero(*par_hb_tm_alt->fpar.momenta);
        break;
    }

    case HMC_ff: {
        mon_hmc_par *par_hmc_ff = (mon_hmc_par *)(data->par);
        zero(*par_hmc_ff->fpar.momenta);
        break;
    }

    case Hasenbusch_ff: {
        mon_hasenbusch_par *par_hb_ff = (mon_hasenbusch_par *)(data->par);
        zero(*par_hb_ff->fpar.momenta);
        break;
    }

    case Scalar: {
        mon_scalar_par *par_scalar = (mon_scalar_par *)(data->par);
        zero(*par_scalar->force_par.momenta);
        break;
    }

    default: {
        error(1, 1, __func__, "Unsupported monomial\n");
        break;
    }
    }
}

static void square_norm_force(mon_type type, monomial_data const *data) {
    switch (type) {
    case PureGauge: {
        mon_pg_par *par_pg = (mon_pg_par *)(data->par);
        lprintf("GAUGE FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_pg->force_par.momenta),
                sqnorm(*par_pg->force_par.momenta));
        break;
    }

    case LuscherWeisz: {
        mon_lw_par *par_lw = (mon_lw_par *)(data->par);
        lprintf("LW FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_lw->force_par.momenta),
                sqnorm(*par_lw->force_par.momenta));
        break;
    }

    case FourFermion: {
        lprintf("FOUR FERMION FORCE", 0, "Skipping force measurement\n");
        break;
    }

    case HMC: {
        mon_hmc_par *par_hmc = (mon_hmc_par *)(data->par);
        lprintf("HMC FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_hmc->fpar.momenta),
                sqnorm(*par_hmc->fpar.momenta));
        break;
    }

    case RHMC: {
        mon_rhmc_par *par_rhmc = (mon_rhmc_par *)(data->par);
        lprintf("RHMC FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_rhmc->fpar.momenta),
                sqnorm(*par_rhmc->fpar.momenta));
        break;
    }

    case TM: {
        mon_tm_par *par_tm = (mon_tm_par *)(data->par);
        lprintf("TM FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_tm->fpar.momenta),
                sqnorm(*par_tm->fpar.momenta));
        break;
    }

    case TM_alt: {
        mon_tm_par *par_tm_alt = (mon_tm_par *)(data->par);
        lprintf("TM ALT FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_tm_alt->fpar.momenta),
                sqnorm(*par_tm_alt->fpar.momenta));
        break;
    }

    case Hasenbusch: {
        mon_hasenbusch_par *par_hb = (mon_hasenbusch_par *)(data->par);
        lprintf("HASENBUSCH FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_hb->fpar.momenta),
                sqnorm(*par_hb->fpar.momenta));
        break;
    }

    case Hasenbusch_tm: {
        mon_hasenbusch_tm_par *par_hb_tm = (mon_hasenbusch_tm_par *)(data->par);
        lprintf("HASENBUSCH TM FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_hb_tm->fpar.momenta),
                sqnorm(*par_hb_tm->fpar.momenta));
        break;
    }

    case Hasenbusch_tm_alt: {
        mon_hasenbusch_tm_par *par_hb_tm_alt = (mon_hasenbusch_tm_par *)(data->par);
        lprintf("HASENBUSCH TM ALT FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_hb_tm_alt->fpar.momenta),
                sqnorm(*par_hb_tm_alt->fpar.momenta));
        break;
    }

    case HMC_ff: {
        mon_hmc_par *par_hmc_ff = (mon_hmc_par *)(data->par);
        lprintf("HMC FF FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_hmc_ff->fpar.momenta),
                sqnorm(*par_hmc_ff->fpar.momenta));
        break;
    }

    case Hasenbusch_ff: {
        mon_hasenbusch_par *par_hb_ff = (mon_hasenbusch_par *)(data->par);
        lprintf("HASENBUSCH FF FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_hb_ff->fpar.momenta),
                sqnorm(*par_hb_ff->fpar.momenta));
        break;
    }

    case Scalar: {
        mon_scalar_par *par_scalar = (mon_scalar_par *)(data->par);
        lprintf("SCALAR FORCE", 0, "Max norm: %0.15e, square norm: %0.15e\n", max(*par_scalar->force_par.momenta),
                sqnorm(*par_scalar->force_par.momenta));
        break;
    }

    default: {
        error(1, 1, __func__, "Unsupported monomial\n");
        break;
    }
    }
}

void print_force_summary() {
    lprintf("INFO", 0, "Evaluating force summary ... \n");

    if (la == NULL) { la = alloc_scalar_field(1, &glattice); }

    /* generate new pseudofermions */
    for (int i = 0; i < num_mon(); ++i) {
        monomial const *m = mon_n(i);
        m->gaussian_pf(m);
    }

    /* compute starting action */
    local_hmc_action(NEW, la, suN_momenta, scalar_momenta);

    /* correct pseudofermion distribution */
    for (int i = 0; i < num_mon(); ++i) {
        monomial const *m = mon_n(i);
        m->correct_pf(m);
    }

    for (int n = 0; n < num_mon(); n++) {
        monomial const *m = mon_n(n);
        mon_type type = m->data.type;
        force_to_zero(type, &m->data);
        m->update_force(0.1, m->force_par);
        square_norm_force(type, &m->data);
    }
}