#include "SCGS_BC.h"

#include "SCGS.h"
#include "../../CFD.h"
#include "../../utils/cALGEBRA/cMAT.h"

void CFD_SCGS_BC_NoSlip_Normal(CFD_t *cfd, SCGS_t *scgs)
{
    if (cfd->engine->method->index->i == cfd->engine->mesh->n_ghosts)
    {
        scgs->vanka->A->data[0][0] = 1.0;
        scgs->vanka->A->data[0][4] = 0.0;
        scgs->vanka->R->data[0] = 0.0;

        scgs->vanka->A->data[4][0] = 0.0;
    }

    if (cfd->engine->method->index->i == cfd->engine->mesh->nodes->Nx + cfd->engine->mesh->n_ghosts - 1)
    {
        scgs->vanka->A->data[1][1] = 1.0;
        scgs->vanka->A->data[1][4] = 0.0;
        scgs->vanka->R->data[1] = 0.0;

        scgs->vanka->A->data[4][1] = 0.0;
    }

    if (cfd->engine->method->index->j == cfd->engine->mesh->n_ghosts)
    {
        scgs->vanka->A->data[2][2] = 1.0;
        scgs->vanka->A->data[2][4] = 0.0;
        scgs->vanka->R->data[2] = 0.0;

        scgs->vanka->A->data[4][2] = 0.0;
    }

    if (cfd->engine->method->index->j == cfd->engine->mesh->nodes->Ny + cfd->engine->mesh->n_ghosts - 1)
    {
        scgs->vanka->A->data[3][3] = 1.0;
        scgs->vanka->A->data[3][4] = 0.0;
        scgs->vanka->R->data[3] = 0.0;

        scgs->vanka->A->data[4][3] = 0.0;
    }
}

void CFD_SCGS_BC_NoSlip_Tangetial(CFD_t *cfd)
{
    // Horizontal walls
    for (uint16_t i = 0;
         i < cfd->engine->mesh->nodes->Nx + 2 * cfd->engine->mesh->n_ghosts;
         i++)
    {
        cfd->engine->method->state->u->data[cfd->engine->mesh->n_ghosts - 1][i] = 2.0 * 0.0 - cfd->engine->method->state->u->data[cfd->engine->mesh->n_ghosts][i];
        cfd->engine->method->state->u->data[cfd->engine->mesh->n_ghosts + cfd->engine->mesh->nodes->Ny][i] = 2.0 * cfd->in->uLid - cfd->engine->method->state->u->data[cfd->engine->mesh->n_ghosts + cfd->engine->mesh->nodes->Ny - 1][i];
    }

    // Vertical walls
    for (uint16_t j = 0;
         j < cfd->engine->mesh->nodes->Ny + 2 * cfd->engine->mesh->n_ghosts;
         j++)
    {
        cfd->engine->method->state->v->data[j][cfd->engine->mesh->n_ghosts - 1] = 2.0 * 0.0 - cfd->engine->method->state->v->data[j][cfd->engine->mesh->n_ghosts];
        cfd->engine->method->state->v->data[j][cfd->engine->mesh->n_ghosts + cfd->engine->mesh->nodes->Nx] = 2.0 * 0.0 - cfd->engine->method->state->v->data[j][cfd->engine->mesh->n_ghosts + cfd->engine->mesh->nodes->Nx - 1];
    }
}
