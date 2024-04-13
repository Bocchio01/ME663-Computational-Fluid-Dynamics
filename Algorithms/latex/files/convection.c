#include "convection.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../CFD.h"
#include "../../utils/cLOG/cLOG.h"
#include "../../utils/cALGEBRA/cVEC.h"
#include "schemes.h"

void CFD_Setup_Convection(CFD_t *cfd)
{
    switch (cfd->engine->schemes->convection->type)
    {
    case UDS:
        cfd->engine->schemes->convection->callable = CFD_Scheme_Convection_UDS;
        break;
    case CDS:
        cfd->engine->schemes->convection->callable = CFD_Scheme_Convection_CDS;
        break;
    case QUICK:
        cfd->engine->schemes->convection->callable = CFD_Scheme_Convection_QUICK;
        break;
    case HYBRID:
        cfd->engine->schemes->convection->callable = CFD_Scheme_Convection_HYBRID;
        break;
    }
}

void CFD_Scheme_Get_Flux(CFD_t *cfd, int i, int j, phi_t phi)
{

    double uPP = CFD_Get_State(cfd, u, i + 0, j + 0);
    double vPP = CFD_Get_State(cfd, v, i + 0, j + 0);

    double uPN = CFD_Get_State(cfd, u, i + 0, j + 1);
    double vPN = CFD_Get_State(cfd, v, i + 0, j + 1);

    // double uEN = CFD_Get_State(cfd, u, i + 1, j + 1);
    // double vEN = CFD_Get_State(cfd, v, i + 1, j + 1);

    double uEP = CFD_Get_State(cfd, u, i + 1, j + 0);
    double vEP = CFD_Get_State(cfd, v, i + 1, j + 0);

    // double uES = CFD_Get_State(cfd, u, i + 1, j - 1);
    double vES = CFD_Get_State(cfd, v, i + 1, j - 1);

    // double uPS = CFD_Get_State(cfd, u, i + 0, j - 1);
    double vPS = CFD_Get_State(cfd, v, i + 0, j - 1);

    // double uWS = CFD_Get_State(cfd, u, i - 1, j - 1);
    // double vWS = CFD_Get_State(cfd, v, i - 1, j - 1);

    double uWP = CFD_Get_State(cfd, u, i - 1, j + 0);
    // double vWP = CFD_Get_State(cfd, v, i - 1, j + 0);

    double uWN = CFD_Get_State(cfd, u, i - 1, j + 1);
    // double vWN = CFD_Get_State(cfd, v, i - 1, j + 1);

    double dx = cfd->engine->mesh->element->size->dx;
    double dy = cfd->engine->mesh->element->size->dy;

    switch (phi)
    {
    case u:

        cfd->engine->schemes->convection->F->w = 1.0 / 2.0 * (uWP + uPP) * dy;
        cfd->engine->schemes->convection->F->e = 1.0 / 2.0 * (uPP + uEP) * dy;
        cfd->engine->schemes->convection->F->s = 1.0 / 2.0 * (vPS + vES) * dx;
        cfd->engine->schemes->convection->F->n = 1.0 / 2.0 * (vPP + vEP) * dx;

        break;

    case v:

        cfd->engine->schemes->convection->F->w = 1.0 / 2.0 * (uWP + uWN) * dy;
        cfd->engine->schemes->convection->F->e = 1.0 / 2.0 * (uPP + uPN) * dy;
        cfd->engine->schemes->convection->F->s = 1.0 / 2.0 * (vPS + vPP) * dx;
        cfd->engine->schemes->convection->F->n = 1.0 / 2.0 * (vPP + vPN) * dx;

        break;

    default:
        log_fatal("Error: CFD_Scheme_Get_Flux phi not found");
        exit(EXIT_FAILURE);
        break;
    }
}

void CFD_Scheme_Convection_UDS(CFD_t *cfd, uint16_t i, uint16_t j, phi_t phi)
{

    CFD_Scheme_Get_Flux(cfd, i, j, phi);

    F_coefficients_t *F = cfd->engine->schemes->convection->F;
    cVEC_t *Ap = cfd->engine->schemes->convection->coefficients;

    Ap->data[WWSS] = 0.0;
    Ap->data[WWS] = 0.0;
    Ap->data[WWP] = 0.0;
    Ap->data[WWN] = 0.0;
    Ap->data[WWNN] = 0.0;
    Ap->data[WSS] = 0.0;
    Ap->data[WS] = 0.0;
    Ap->data[WP] = fmax(0.0, F->w);
    Ap->data[WN] = 0.0;
    Ap->data[WNN] = 0.0;
    Ap->data[PSS] = 0.0;
    Ap->data[PS] = fmax(0.0, F->s);
    Ap->data[PP] = fmax(0.0, F->e) + fmax(0.0, F->n) - fmin(0.0, F->s) - fmin(0.0, F->w);
    Ap->data[PN] = -fmin(0.0, F->n);
    Ap->data[PNN] = 0.0;
    Ap->data[ESS] = 0.0;
    Ap->data[ES] = 0.0;
    Ap->data[EP] = -fmin(0.0, F->e);
    Ap->data[EN] = 0.0;
    Ap->data[ENN] = 0.0;
    Ap->data[EESS] = 0.0;
    Ap->data[EES] = 0.0;
    Ap->data[EEP] = 0.0;
    Ap->data[EEN] = 0.0;
    Ap->data[EENN] = 0.0;
}

void CFD_Scheme_Convection_CDS(CFD_t *cfd, uint16_t i, uint16_t j, phi_t phi)
{
    CFD_Scheme_Get_Flux(cfd, i, j, phi);

    F_coefficients_t *F = cfd->engine->schemes->convection->F;
    cVEC_t *Ap = cfd->engine->schemes->convection->coefficients;

    Ap->data[WWSS] = 0.0;
    Ap->data[WWS] = 0.0;
    Ap->data[WWP] = 0.0;
    Ap->data[WWN] = 0.0;
    Ap->data[WWNN] = 0.0;
    Ap->data[WSS] = 0.0;
    Ap->data[WS] = 0.0;
    Ap->data[WP] = F->w / 2.0;
    Ap->data[WN] = 0.0;
    Ap->data[WNN] = 0.0;
    Ap->data[PSS] = 0.0;
    Ap->data[PS] = F->s / 2.0;
    Ap->data[PP] = F->e / 2.0 + F->n / 2.0 - F->s / 2.0 - F->w / 2.0;
    Ap->data[PN] = -F->n / 2.0;
    Ap->data[PNN] = 0.0;
    Ap->data[ESS] = 0.0;
    Ap->data[ES] = 0.0;
    Ap->data[EP] = -F->e / 2.0;
    Ap->data[EN] = 0.0;
    Ap->data[ENN] = 0.0;
    Ap->data[EESS] = 0.0;
    Ap->data[EES] = 0.0;
    Ap->data[EEP] = 0.0;
    Ap->data[EEN] = 0.0;
    Ap->data[EENN] = 0.0;
}

void CFD_Scheme_Convection_QUICK(CFD_t *cfd, uint16_t i, uint16_t j, phi_t phi)
{
    CFD_Scheme_Get_Flux(cfd, i, j, phi);

    F_coefficients_t *F = cfd->engine->schemes->convection->F;
    cVEC_t *Ap = cfd->engine->schemes->convection->coefficients;

    Ap->data[WWSS] = 0.0;
    Ap->data[WWS] = 0.0;
    Ap->data[WWP] = -1.0 / 8.0 * fmax(0.0, F->w);
    Ap->data[WWN] = 0.0;
    Ap->data[WWNN] = 0.0;
    Ap->data[WSS] = 0.0;
    Ap->data[WS] = 0.0;
    Ap->data[WP] = 1.0 / 8.0 * fmax(0.0, F->e) + 3.0 / 4.0 * fmax(0.0, F->w) + 3.0 / 8.0 * fmin(0.0, F->w);
    Ap->data[WN] = 0.0;
    Ap->data[WNN] = 0.0;
    Ap->data[PSS] = -1.0 / 8.0 * fmax(0.0, F->s);
    Ap->data[PS] = 1.0 / 8.0 * fmax(0.0, F->n) + 3.0 / 4.0 * fmax(0.0, F->s) + 3.0 / 8.0 * fmin(0.0, F->s);
    Ap->data[PP] = 3.0 / 4.0 * fmax(0.0, F->e) + 3.0 / 8.0 * fmin(0.0, F->e) + 3.0 / 4.0 * fmax(0.0, F->n) + 3.0 / 8.0 * fmin(0.0, F->n) - 3.0 / 8.0 * fmax(0.0, F->s) - 3.0 / 4.0 * fmin(0.0, F->s) - 3.0 / 8.0 * fmax(0.0, F->w) - 3.0 / 4.0 * fmin(0.0, F->w);
    Ap->data[PN] = 1.0 / 8.0 * (-3.0) * fmax(0.0, F->n) - 3.0 / 4.0 * fmin(0.0, F->n) - 1.0 / 8.0 * fmin(0.0, F->s);
    Ap->data[PNN] = 1.0 / 8.0 * fmin(0.0, F->n);
    Ap->data[ESS] = 0.0;
    Ap->data[ES] = 0.0;
    Ap->data[EP] = 1.0 / 8.0 * (-3.0) * fmax(0.0, F->e) - 3.0 / 4.0 * fmin(0.0, F->e) - 1.0 / 8.0 * fmin(0.0, F->w);
    Ap->data[EN] = 0.0;
    Ap->data[ENN] = 0.0;
    Ap->data[EESS] = 0.0;
    Ap->data[EES] = 0.0;
    Ap->data[EEP] = 1.0 / 8.0 * fmin(0.0, F->e);
    Ap->data[EEN] = 0.0;
    Ap->data[EENN] = 0.0;
}

void CFD_Scheme_Convection_HYBRID(CFD_t *cfd, uint16_t i, uint16_t j, phi_t phi)
{
    CFD_Scheme_Get_Flux(cfd, i, j, phi);

    // F_coefficients_t *F = cfd->engine->schemes->convection->F;
    // cVEC_t *Ap = cfd->engine->schemes->convection->coefficients;

    log_fatal("HYBRID not implemented yet");
    exit(EXIT_FAILURE);
}