#include "SCGS.h"

typedef struct CFD_t CFD_t;

#include <stdlib.h>
#include <math.h>
#include "../../CFD.h"
#include "../schemes/schemes.h"
#include "../../utils/cALGEBRA/cMAT.h"
#include "../../utils/cALGEBRA/cVEC.h"
#include "../../utils/cLOG/cLOG.h"

void CFD_SCGS(CFD_t *cfd)
{
    SCGS_t *scgs;

    scgs = CFD_SCGS_Allocate(cfd);

    for (int iteration = 0; iteration < cfd->engine->method->maxIter; iteration++)
    {
        CFD_SCGS_Reset(scgs);

        for (cfd->engine->method->index->j = cfd->engine->mesh->n_ghosts;
             cfd->engine->method->index->j < cfd->engine->mesh->nodes->Ny + cfd->engine->mesh->n_ghosts;
             cfd->engine->method->index->j++)
        {
            for (cfd->engine->method->index->i = cfd->engine->mesh->n_ghosts;
                 cfd->engine->method->index->i < cfd->engine->mesh->nodes->Nx + cfd->engine->mesh->n_ghosts;
                 cfd->engine->method->index->i++)
            {

                CFD_SCGS_System_Compose(cfd, scgs);

                CFD_SCGS_BC_NoSlip_Normal(cfd, scgs);

                CFD_SCGS_System_Solve(scgs);

                CFD_SCGS_Apply_Correction(cfd, scgs);

                CFD_SCGS_Update_Residual(scgs);
            }
        }

        CFD_SCGS_BC_NoSlip_Tangetial(cfd, scgs);

        if (fmax(scgs->residual->u, fmax(scgs->residual->v, scgs->residual->p)) < cfd->engine->method->tolerance &&
            iteration > 1)
        {
            log_info("Algorithm converged in %d iterations", iteration);
            break;
        }
    }
}