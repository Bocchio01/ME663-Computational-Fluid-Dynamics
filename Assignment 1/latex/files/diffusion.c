#include "diffusion.h"

#include "../../CFD.h"
#include "../../utils/cLOG/cLOG.h"
#include "../../utils/cALGEBRA/cVEC.h"
#include "schemes.h"

void CFD_Setup_Diffusion(CFD_t *cfd)
{
    switch (cfd->engine->schemes->diffusion->type)
    {
    case SECOND:
        cfd->engine->schemes->diffusion->callable = CFD_Scheme_Diffusion_SECOND;
        break;
    case FOURTH:
        cfd->engine->schemes->diffusion->callable = CFD_Scheme_Diffusion_FOURTH;
        break;
    }
}

void CFD_Scheme_Diffusion_SECOND(CFD_t *cfd)
{
    cVEC_t *Ap = cfd->engine->schemes->diffusion->coefficients;

    double dx = cfd->engine->mesh->element->size->dx;
    double dy = cfd->engine->mesh->element->size->dy;
    double nu = cfd->in->fluid->nu;

    Ap->data[WWSS] = 0.0;
    Ap->data[WWS] = 0.0;
    Ap->data[WWP] = 0.0;
    Ap->data[WWN] = 0.0;
    Ap->data[WWNN] = 0.0;
    Ap->data[WSS] = 0.0;
    Ap->data[WS] = 0.0;
    Ap->data[WP] = (dy * nu) / (dx);
    Ap->data[WN] = 0.0;
    Ap->data[WNN] = 0.0;
    Ap->data[PSS] = 0.0;
    Ap->data[PS] = (dx * nu) / (dy);
    Ap->data[PP] = -1.0 * dx * dy * nu * (-2.0 / (dx * dx) - 2.0 / (dy * dy));
    Ap->data[PN] = (dx * nu) / (dy);
    Ap->data[PNN] = 0.0;
    Ap->data[ESS] = 0.0;
    Ap->data[ES] = 0.0;
    Ap->data[EP] = (dy * nu) / (dx);
    Ap->data[EN] = 0.0;
    Ap->data[ENN] = 0.0;
    Ap->data[EESS] = 0.0;
    Ap->data[EES] = 0.0;
    Ap->data[EEP] = 0.0;
    Ap->data[EEN] = 0.0;
    Ap->data[EENN] = 0.0;
}

void CFD_Scheme_Diffusion_FOURTH(CFD_t *cfd)
{
    cVEC_t *Ap = cfd->engine->schemes->diffusion->coefficients;

    double dx = cfd->engine->mesh->element->size->dx;
    double dy = cfd->engine->mesh->element->size->dy;
    double nu = cfd->in->fluid->nu;

    Ap->data[WWSS] = 0.0;
    Ap->data[WWS] = 0.0;
    Ap->data[WWP] = -(dy * nu) / (12.0 * dx);
    Ap->data[WWN] = 0.0;
    Ap->data[WWNN] = 0.0;
    Ap->data[WSS] = 0.0;
    Ap->data[WS] = 0.0;
    Ap->data[WP] = (4.0 * dy * nu) / (3.0 * dx);
    Ap->data[WN] = 0.0;
    Ap->data[WNN] = 0.0;
    Ap->data[PSS] = -(dx * nu) / (12.0 * dy);
    Ap->data[PS] = (4.0 * dx * nu) / (3.0 * dy);
    Ap->data[PP] = -dx * dy * nu * (-5.0 / (2.0 * dx * dx) - 5.0 / (2.0 * dy * dy));
    Ap->data[PN] = (4.0 * dx * nu) / (3.0 * dy);
    Ap->data[PNN] = -(dx * nu) / (12.0 * dy);
    Ap->data[ESS] = 0.0;
    Ap->data[ES] = 0.0;
    Ap->data[EP] = (4.0 * dy * nu) / (3.0 * dx);
    Ap->data[EN] = 0.0;
    Ap->data[ENN] = 0.0;
    Ap->data[EESS] = 0.0;
    Ap->data[EES] = 0.0;
    Ap->data[EEP] = -(dy * nu) / (12.0 * dx);
    Ap->data[EEN] = 0.0;
    Ap->data[EENN] = 0.0;
}