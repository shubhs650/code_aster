#ifndef PHYSICSANDMODELISATIONS_H_
#define PHYSICSANDMODELISATIONS_H_

/**
 * @file PhysicsAndModelings.h
 * @brief Fichier definissant les physiques et les modelisations disponibles
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @enum Physics
 * @brief Physiques existantes dans Code_Aster
 * @author Nicolas Sellenet
 */
enum Physics { Mechanics, Thermal, Acoustic };
const int nbPhysics = 3;
/**
 * @var PhysicNames
 * @brief Nom Aster des differentes physiques disponibles
 */
extern const char *const PhysicNames[nbPhysics];

/**
 * @enum Modelings
 * @brief Modelisations existantes dans Code_Aster
 * @author Nicolas Sellenet
 */
enum Modelings {
    PlanarBar,
    BARRE_2D_NL,
    DIS_T_2D,
    DIS_TR_2D,
    FLUIDE_2D,
    FLUI_ABSO_2D,
    FLUI_PESA_2D,
    FLUI_STRU_2D,
    Tridimensional,
    TridimensionalAbsorbingBoundary,
    DIAG_3D,
    DIL_3D,
    FAISCEAU_3D,
    FLUIDE_3D,
    FLUI_ABSO_3D,
    GRAD_HHO_3D,
    GRAD_INCO_3D,
    GRAD_VARI_3D,
    HH2D_3D,
    HH2MD_3D,
    HH2MS_3D,
    HH2MS_DIL_3D,
    HH2M_SI_3D,
    HH2S_3D,
    HH2SUDA_3D,
    HHD_3D,
    HHM_3D,
    HHMD_3D,
    HHMS_3D,
    HHO_3D,
    HHS_3D,
    HM_3D,
    HMD_3D,
    HMS_3D,
    HMS_DIL_3D,
    HM_SI_3D,
    HM_SI_DIL_3D,
    HS_3D,
    INCO_UP_3D,
    INCO_UPG_3D,
    INCO_UPO_3D,
    INTERFACE_3D,
    INTERFACE_S_3D,
    INTERF_POU_3D,
    JOINT_3D,
    JOINT_HYME_3D,
    SECH_3D,
    SECH_3D_DIAG,
    MIX_STA_3D,
    SI_3D,
    THH2D_3D,
    THH2MD_3D,
    THH2MS_3D,
    THH2S_3D,
    THHD_3D,
    THHM_3D,
    THHMD_3D,
    THHMS_3D,
    THHS_3D,
    THM_3D,
    THMD_3D,
    THMS_3D,
    THMS_DIL_3D,
    THVD_3D,
    THVS_3D,
    Axisymmetrical,
    AXIS_DIAG,
    AXIS_FLUIDE,
    AXIS_FLUI_ABSO,
    AXIS_FLUI_STRU,
    AXIS_FOURIER,
    AXIS_GRAD_INCO,
    AXIS_GRAD_HHO,
    AXIS_GRAD_VARI,
    AXIS_HH2D,
    AXIS_HH2MD,
    AXIS_HH2MS,
    AXIS_HH2S,
    AXIS_HHD,
    AXIS_HHM,
    AXIS_HHMD,
    AXIS_HHMS,
    AXIS_HHO,
    AXIS_HHS,
    AXIS_HM,
    AXIS_HMD,
    AXIS_HMS,
    AXIS_INCO_UP,
    AXIS_INCO_UPG,
    AXIS_INCO_UPO,
    AXIS_INTERFACE,
    AXIS_INTERFACE_S,
    AXIS_JHMS,
    AXIS_JOINT,
    AXIS_SECH,
    AXIS_SECH_DIAG,
    AXIS_SI,
    AXIS_THH2D,
    AXIS_THH2MD,
    AXIS_THH2MS,
    AXIS_THH2S,
    AXIS_THHD,
    AXIS_THHMD,
    AXIS_THHMS,
    AXIS_THHS,
    AXIS_THM,
    AXIS_THMD,
    AXIS_THMS,
    AXIS_THVD,
    AXIS_THVS,
    BARRE,
    BARRE_3D_NL,
    CABLE,
    CABLE_GAINE,
    CABLE_POULIE,
    COQUE,
    COQUE_3D,
    COQUE_AXIS,
    COQUE_PLAN,
    COQUE_SOLIDE,
    PlaneStress,
    C_PLAN_SI,
    DIS_T,
    DIS_TR,
    DKT,
    DKTG,
    DST,
    PlaneStrain,
    D_PLAN_2DG,
    D_PLAN_ABSO,
    D_PLAN_DIL,
    D_PLAN_GRAD_HHO,
    D_PLAN_GRAD_INCO,
    D_PLAN_GRAD_SIGM,
    D_PLAN_GRAD_VARI,
    D_PLAN_HH2D,
    D_PLAN_HH2MD,
    D_PLAN_HH2MS,
    D_PLAN_HH2MS_DIL,
    D_PLAN_HH2M_SI,
    D_PLAN_HH2S,
    D_PLAN_HH2SUDA,
    D_PLAN_HHD,
    D_PLAN_HHM,
    D_PLAN_HHMD,
    D_PLAN_HHMS,
    D_PLAN_HHO,
    D_PLAN_HHS,
    D_PLAN_HM,
    D_PLAN_HMD,
    D_PLAN_HMS,
    D_PLAN_HMS_DIL,
    D_PLAN_HM_SI,
    D_PLAN_HM_SI_DIL,
    D_PLAN_HS,
    D_PLAN_INCO_UP,
    D_PLAN_INCO_UPG,
    D_PLAN_INCO_UPO,
    D_PLAN_MIX_STA,
    D_PLAN_SI,
    D_PLAN_THH2D,
    D_PLAN_THH2MD,
    D_PLAN_THH2MS,
    D_PLAN_THH2S,
    D_PLAN_THHD,
    D_PLAN_THHMD,
    D_PLAN_THHMS,
    D_PLAN_THHS,
    D_PLAN_THM,
    D_PLAN_THMD,
    D_PLAN_THMS,
    D_PLAN_THMS_DIL,
    D_PLAN_THVD,
    D_PLAN_THVS,
    PLAQ_MITC,
    FLUI_STRU,
    GRILLE_EXCENTRE,
    GRILLE_MEMBRANE,
    MEMBRANE,
    Planar,
    PLAN_ABSO,
    PLAN_DIAG,
    PLAN_HHO,
    PLAN_INTERFACE,
    PLAN_INTERFACE_S,
    PLAN_JHMS,
    PLAN_JOINT,
    PLAN_JOINT_HYME,
    POU_D_E,
    POU_D_EM,
    POU_D_SQUE,
    POU_D_T,
    POU_D_TG,
    POU_D_TGM,
    POU_D_T_GD,
    POU_FLUI_STRU,
    Q4G,
    Q4GG,
    TUYAU_3M,
    TUYAU_6M,
};
const int nbModelings = 200;
/**
 * @var ModelingNames
 * @brief Nom Aster des differentes modelisations disponibles
 */
extern const char *const ModelingNames[nbModelings];

const int nbModelingsMechanics = 188;
extern const Modelings MechanicsModelings[nbModelingsMechanics];

const int nbModelingsThermal = 17;
extern const Modelings ThermalModelings[nbModelingsThermal];

const int nbModelingsAcoustic = 4;
extern const Modelings AcousticModelings[nbModelingsAcoustic];

/**
 * @enum Formulation
 * @author Nicolas Sellenet
 */
enum Formulation { NoFormulation, Linear, Quadratic, UPPhi, UP, UPsi, Dil, DilInco, Sta, StaInco };
const int nbFormulation = 10;
/**
 * @var FormulationNames
 * @brief Nom Aster des differentes physiques disponibles
 */
extern const char *const FormulationNames[nbFormulation];

#endif /* PHYSICSANDMODELISATIONS_H_ */
