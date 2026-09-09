/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2026 - EDF - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#include "astercxx.h"

#include "tria_mitc.h"

#include <cmath>
#include <iomanip> // Pour setprecision
#include <iostream>
#include <vector>

#include <math.h>
#include <stdalign.h>
#include <stdlib.h>
#include <string.h>

VectorReal B_p1_tr6( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 21, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)
    // Quadrature rules
    static const double weights_39d[6] = { 0.054975871827661,  0.054975871827661,
                                           0.054975871827661,  0.1116907948390055,
                                           0.1116907948390055, 0.1116907948390055 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE1_C0_D10_Q39d[1][1][1][3] = { { { { -1.0, 1.0, 0.0 } } } };
    static const double FE1_C1_D01_Q39d[1][1][1][3] = { { { { -1.0, 0.0, 1.0 } } } };
    static const double FE5_C1_D01_Q39d[1][1][6][6] = {
        { { { 0.6336951459609189, 0.0, -0.6336951459609159, 3.267390291921835, 0.0,
              -3.267390291921836 },
            { 0.6336951459609192, 0.0, 2.267390291921836, 0.3663048540390839, -2.901085437882756,
              -0.3663048540390842 },
            { -2.267390291921832, 0.0, -0.6336951459609156, 0.366304854039084, 2.901085437882748,
              -0.3663048540390843 },
            { -0.7837939636638601, 0.0, 0.78379396366386, 0.43241207267228, 0.0,
              -0.4324120726722802 },
            { -0.7837939636638601, 0.0, -0.5675879273277196, 1.78379396366386, 1.35138189099158,
              -1.78379396366386 },
            { 0.5675879273277193, 0.0, 0.7837939636638599, 1.78379396366386, -1.35138189099158,
              -1.78379396366386 } } }
    };
    static const double FE5_C1_D10_Q39d[1][1][6][6] = {
        { { { 0.6336951459609198, 2.267390291921836, 0.0, 0.3663048540390836, -0.3663048540390846,
              -2.901085437882756 },
            { 0.6336951459609188, -0.633695145960915, 0.0, 3.267390291921836, -3.267390291921836,
              0.0 },
            { -2.267390291921832, -0.6336951459609163, 0.0, 0.3663048540390844, -0.3663048540390835,
              2.901085437882748 },
            { -0.7837939636638603, -0.5675879273277198, 0.0, 1.78379396366386, -1.78379396366386,
              1.35138189099158 },
            { -0.78379396366386, 0.78379396366386, 0.0, 0.43241207267228, -0.4324120726722801,
              0.0 },
            { 0.5675879273277195, 0.7837939636638604, 0.0, 1.78379396366386, -1.78379396366386,
              -1.35138189099158 } } }
    };
    static const double FE5_C2_Q39d[1][1][6][3] = {
        { { { 0.09157621350977008, 0.816847572980459, 0.09157621350977098 },
            { 0.09157621350976999, 0.09157621350977097, 0.816847572980459 },
            { 0.8168475729804581, 0.09157621350977091, 0.09157621350977101 },
            { 0.445948490915965, 0.10810301816807, 0.445948490915965 },
            { 0.4459484909159651, 0.445948490915965, 0.10810301816807 },
            { 0.10810301816807, 0.445948490915965, 0.445948490915965 } } }
    };
    static const double FE5_C3_Q39d[1][1][6][3] = {
        { { { -0.09157621350977084, 0.09157621350977084, 0.9084237864902291 },
            { -0.8168475729804592, 0.8168475729804592, 0.1831524270195409 },
            { -0.09157621350977109, 0.09157621350977109, 0.9084237864902289 },
            { -0.4459484909159651, 0.4459484909159651, 0.5540515090840348 },
            { -0.10810301816807, 0.10810301816807, 0.89189698183193 },
            { -0.445948490915965, 0.445948490915965, 0.554051509084035 } } }
    };
    static const double FE5_C4_Q39d[1][1][6][3] = {
        { { { 0.8168475729804592, 0.1831524270195409, 0.8168475729804592 },
            { 0.09157621350977091, 0.9084237864902291, 0.09157621350977091 },
            { 0.09157621350977108, 0.9084237864902289, 0.09157621350977108 },
            { 0.10810301816807, 0.89189698183193, 0.10810301816807 },
            { 0.4459484909159651, 0.5540515090840349, 0.4459484909159651 },
            { 0.445948490915965, 0.554051509084035, 0.445948490915965 } } }
    };
    // ------------------------
    // Section: Jacobian
    // Inputs: coordinate_dofs, FE1_C0_D10_Q39d, FE1_C1_D01_Q39d
    // Outputs: J_c2, J_c3, J_c0, J_c1
    double J_c0 = 0.0;
    double J_c3 = 0.0;
    double J_c1 = 0.0;
    double J_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            J_c0 += coordinate_dofs[(ic)*3] * FE1_C0_D10_Q39d[0][0][0][ic];
            J_c3 += coordinate_dofs[(ic)*3 + 1] * FE1_C1_D01_Q39d[0][0][0][ic];
            J_c1 += coordinate_dofs[(ic)*3] * FE1_C1_D01_Q39d[0][0][0][ic];
            J_c2 += coordinate_dofs[(ic)*3 + 1] * FE1_C0_D10_Q39d[0][0][0][ic];
        }
    }
    // ------------------------
    double sp_39d_0 = pow( c[3], 3 );
    double sp_39d_1 = c[4] * sp_39d_0;
    double sp_39d_2 = -sp_39d_1;
    double sp_39d_3 = J_c0 * J_c3;
    double sp_39d_4 = J_c1 * J_c2;
    double sp_39d_5 = -sp_39d_4;
    double sp_39d_6 = sp_39d_3 + sp_39d_5;
    double sp_39d_7 = J_c0 / sp_39d_6;
    double sp_39d_8 = -J_c1;
    double sp_39d_9 = sp_39d_8 / sp_39d_6;
    double sp_39d_10 = sp_39d_7 + sp_39d_7;
    double sp_39d_11 = sp_39d_9 + sp_39d_9;
    double sp_39d_12 = sp_39d_10 / 2;
    double sp_39d_13 = sp_39d_11 / 2;
    double sp_39d_14 = J_c3 / sp_39d_6;
    double sp_39d_15 = -J_c2;
    double sp_39d_16 = sp_39d_15 / sp_39d_6;
    double sp_39d_17 = sp_39d_16 / 2;
    double sp_39d_18 = sp_39d_14 / 2;
    double sp_39d_19 = sp_39d_7 / 2;
    double sp_39d_20 = sp_39d_9 / 2;
    double sp_39d_21 = sp_39d_16 + sp_39d_16;
    double sp_39d_22 = sp_39d_14 + sp_39d_14;
    double sp_39d_23 = sp_39d_21 / 2;
    double sp_39d_24 = sp_39d_22 / 2;
    double sp_39d_25 = -c[1];
    double sp_39d_26 = 1.0 + sp_39d_25;
    double sp_39d_27 = 2 * sp_39d_12;
    double sp_39d_28 = 2 * sp_39d_13;
    double sp_39d_29 = 2 * sp_39d_23;
    double sp_39d_30 = 2 * sp_39d_24;
    double sp_39d_31 = c[0] * sp_39d_0;
    double sp_39d_32 = pow( c[1], 2 );
    double sp_39d_33 = -sp_39d_32;
    double sp_39d_34 = 1.0 + sp_39d_33;
    double sp_39d_35 = 24.0 * sp_39d_34;
    double sp_39d_36 = sp_39d_31 / sp_39d_35;
    double sp_39d_37 = c[0] * c[2];
    double sp_39d_38 = c[3] * sp_39d_37;
    double sp_39d_39 = 1.0 + c[1];
    double sp_39d_40 = 4.0 * sp_39d_39;
    double sp_39d_41 = sp_39d_38 / sp_39d_40;
    double sp_39d_42 = fabs( sp_39d_6 );
    double sp_39d_43 = sp_39d_2 * sp_39d_42;
    for ( int iq = 0; iq < 6; ++iq ) {
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C1_D01_Q39d
        // Outputs: w0_d01_c1
        double w0_d01_c1 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_d01_c1 += w[(ic)*2 + 1] * FE5_C1_D01_Q39d[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C1_D10_Q39d
        // Outputs: w0_d10_c1
        double w0_d10_c1 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_d10_c1 += w[(ic)*2 + 1] * FE5_C1_D10_Q39d[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C1_D01_Q39d
        // Outputs: w0_d01_c0
        double w0_d01_c0 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_d01_c0 += w[(ic)*2] * FE5_C1_D01_Q39d[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C1_D10_Q39d
        // Outputs: w0_d10_c0
        double w0_d10_c0 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_d10_c0 += w[(ic)*2] * FE5_C1_D10_Q39d[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C4_Q39d
        // Outputs: w0_c4
        double w0_c4 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c4 += w[( ic ) + 15] * FE5_C4_Q39d[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C3_Q39d
        // Outputs: w0_c3
        double w0_c3 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c3 += w[( ic ) + 15] * FE5_C3_Q39d[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Intermediates
        // Inputs: w0_d01_c1, w0_d10_c1, w0_d01_c0, w0_d10_c0, w0_c4, w0_c3
        // Outputs: fw0, fw1, fw2, fw3, fw4, fw5, fw6
        double fw0 = 0;
        double fw1 = 0;
        double fw2 = 0;
        double fw3 = 0;
        double fw4 = 0;
        double fw5 = 0;
        double fw6 = 0;
        {
            double sv_39d_0 = w0_d01_c1 * sp_39d_7;
            double sv_39d_1 = w0_d10_c1 * sp_39d_9;
            double sv_39d_2 = sv_39d_0 + sv_39d_1;
            double sv_39d_3 = sv_39d_2 + sv_39d_2;
            double sv_39d_4 = sv_39d_3 / 2;
            double sv_39d_5 = sv_39d_4 * sp_39d_12;
            double sv_39d_6 = sv_39d_4 * sp_39d_13;
            double sv_39d_7 = sv_39d_5 + sv_39d_5;
            double sv_39d_8 = sv_39d_6 + sv_39d_6;
            double sv_39d_9 = w0_d01_c0 * sp_39d_7;
            double sv_39d_10 = w0_d10_c0 * sp_39d_9;
            double sv_39d_11 = sv_39d_9 + sv_39d_10;
            double sv_39d_12 = w0_d10_c1 * sp_39d_14;
            double sv_39d_13 = w0_d01_c1 * sp_39d_16;
            double sv_39d_14 = sv_39d_12 + sv_39d_13;
            double sv_39d_15 = sv_39d_11 + sv_39d_14;
            double sv_39d_16 = sv_39d_15 / 2;
            double sv_39d_17 = sv_39d_16 * sp_39d_17;
            double sv_39d_18 = sv_39d_16 * sp_39d_18;
            double sv_39d_19 = sv_39d_16 * sp_39d_19;
            double sv_39d_20 = sv_39d_16 * sp_39d_20;
            double sv_39d_21 = sv_39d_17 + sv_39d_17;
            double sv_39d_22 = sv_39d_18 + sv_39d_18;
            double sv_39d_23 = sv_39d_19 + sv_39d_19;
            double sv_39d_24 = sv_39d_20 + sv_39d_20;
            double sv_39d_25 = sv_39d_7 + sv_39d_21;
            double sv_39d_26 = sv_39d_8 + sv_39d_22;
            double sv_39d_27 = w0_d10_c0 * sp_39d_14;
            double sv_39d_28 = w0_d01_c0 * sp_39d_16;
            double sv_39d_29 = sv_39d_27 + sv_39d_28;
            double sv_39d_30 = sv_39d_29 + sv_39d_29;
            double sv_39d_31 = sv_39d_30 / 2;
            double sv_39d_32 = sv_39d_31 * sp_39d_23;
            double sv_39d_33 = sv_39d_31 * sp_39d_24;
            double sv_39d_34 = sv_39d_32 + sv_39d_32;
            double sv_39d_35 = sv_39d_33 + sv_39d_33;
            double sv_39d_36 = sv_39d_34 + sv_39d_23;
            double sv_39d_37 = sv_39d_35 + sv_39d_24;
            double sv_39d_38 = sv_39d_25 + sv_39d_21;
            double sv_39d_39 = sv_39d_26 + sv_39d_22;
            double sv_39d_40 = sv_39d_36 + sv_39d_23;
            double sv_39d_41 = sv_39d_37 + sv_39d_24;
            double sv_39d_42 = sv_39d_38 * sp_39d_26;
            double sv_39d_43 = sv_39d_39 * sp_39d_26;
            double sv_39d_44 = sv_39d_40 * sp_39d_26;
            double sv_39d_45 = sv_39d_41 * sp_39d_26;
            double sv_39d_46 = sv_39d_4 + sv_39d_31;
            double sv_39d_47 = sv_39d_46 * sp_39d_27;
            double sv_39d_48 = sv_39d_46 * sp_39d_28;
            double sv_39d_49 = sv_39d_46 * sp_39d_29;
            double sv_39d_50 = sv_39d_46 * sp_39d_30;
            double sv_39d_51 = c[1] * sv_39d_47;
            double sv_39d_52 = c[1] * sv_39d_48;
            double sv_39d_53 = c[1] * sv_39d_49;
            double sv_39d_54 = c[1] * sv_39d_50;
            double sv_39d_55 = sv_39d_42 + sv_39d_51;
            double sv_39d_56 = sv_39d_43 + sv_39d_52;
            double sv_39d_57 = sv_39d_44 + sv_39d_53;
            double sv_39d_58 = sv_39d_45 + sv_39d_54;
            double sv_39d_59 = sv_39d_55 * sp_39d_36;
            double sv_39d_60 = sv_39d_56 * sp_39d_36;
            double sv_39d_61 = sv_39d_57 * sp_39d_36;
            double sv_39d_62 = sv_39d_58 * sp_39d_36;
            double sv_39d_63 = w0_c4 * sp_39d_7;
            double sv_39d_64 = w0_c3 * sp_39d_9;
            double sv_39d_65 = sv_39d_63 + sv_39d_64;
            double sv_39d_66 = sv_39d_65 * sp_39d_7;
            double sv_39d_67 = sv_39d_65 * sp_39d_9;
            double sv_39d_68 = sv_39d_66 + sv_39d_66;
            double sv_39d_69 = sv_39d_67 + sv_39d_67;
            double sv_39d_70 = w0_c3 * sp_39d_14;
            double sv_39d_71 = w0_c4 * sp_39d_16;
            double sv_39d_72 = sv_39d_70 + sv_39d_71;
            double sv_39d_73 = sv_39d_72 * sp_39d_16;
            double sv_39d_74 = sv_39d_72 * sp_39d_14;
            double sv_39d_75 = sv_39d_73 + sv_39d_73;
            double sv_39d_76 = sv_39d_74 + sv_39d_74;
            double sv_39d_77 = sv_39d_68 + sv_39d_75;
            double sv_39d_78 = sv_39d_76 + sv_39d_69;
            double sv_39d_79 = sv_39d_77 * sp_39d_41;
            double sv_39d_80 = sv_39d_78 * sp_39d_41;
            double sv_39d_81 = sv_39d_59 * sp_39d_42;
            double sv_39d_82 = sv_39d_60 * sp_39d_42;
            double sv_39d_83 = sv_39d_61 * sp_39d_42;
            double sv_39d_84 = sv_39d_62 * sp_39d_42;
            double sv_39d_85 = sv_39d_79 * sp_39d_42;
            double sv_39d_86 = sv_39d_80 * sp_39d_42;
            fw0 = sv_39d_84 * weights_39d[iq];
            fw1 = sv_39d_83 * weights_39d[iq];
            fw2 = sv_39d_82 * weights_39d[iq];
            fw3 = sv_39d_81 * weights_39d[iq];
            fw4 = sp_39d_43 * weights_39d[iq];
            fw5 = sv_39d_86 * weights_39d[iq];
            fw6 = sv_39d_85 * weights_39d[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw3, fw2, FE5_C1_D10_Q39d, fw0, fw1, FE5_C1_D01_Q39d
        // Outputs: A
        {
            for ( int i = 0; i < 6; ++i ) {
                A[2 * ( i )] += fw0 * FE5_C1_D10_Q39d[0][0][iq][i];
                A[2 * ( i )] += fw1 * FE5_C1_D01_Q39d[0][0][iq][i];
                A[( 2 * ( i ) + 1 )] += fw2 * FE5_C1_D10_Q39d[0][0][iq][i];
                A[( 2 * ( i ) + 1 )] += fw3 * FE5_C1_D01_Q39d[0][0][iq][i];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: FE5_C4_Q39d, FE5_C2_Q39d, fw4, fw5, fw6, FE5_C3_Q39d
        // Outputs: A
        {
            for ( int i = 0; i < 3; ++i ) {
                A[( ( i ) + 12 )] += fw4 * FE5_C2_Q39d[0][0][iq][i];
                A[( ( i ) + 15 )] += fw5 * FE5_C3_Q39d[0][0][iq][i];
                A[( ( i ) + 15 )] += fw6 * FE5_C4_Q39d[0][0][iq][i];
            }
        }
        // ------------------------
    }
    return A;
}

VectorReal B_p2_tr6( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 21, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)
    // Quadrature rules
    static const double weights_8e7[3] = { 0.2777777777777777, 0.4444444444444444,
                                           0.2777777777777777 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE1_C1_D01_F_Q8e7[1][1][1][3] = { { { { -1.0, 0.0, 1.0 } } } };
    static const double FE5_C0_F_Q8e7[1][3][3][6] = {
        { { { 0.0, 0.6872983346207417, -0.08729833462074169, 0.4, 0.0, 0.0 },
            { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
            { 0.0, -0.08729833462074159, 0.6872983346207417, 0.4, 0.0, 0.0 } },
          { { 0.6872983346207417, 0.0, -0.08729833462074171, 0.0, 0.3999999999999999, 0.0 },
            { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
            { -0.08729833462074173, 0.0, 0.6872983346207417, 0.0, 0.4, 0.0 } },
          { { 0.6872983346207417, -0.0872983346207416, 0.0, 0.0, 0.0, 0.4 },
            { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
            { -0.08729833462074155, 0.6872983346207416, 0.0, 0.0, 0.0, 0.4 } } }
    };
    static const double FE5_C2_D10_F_Q8e7[1][1][1][3] = { { { { -1.0, 1.0, 0.0 } } } };
    static const double FE5_C3_F_Q8e7[1][3][3][3] = {
        { { { -0.1127016653792581, 0.1127016653792581, 0.8872983346207418 },
            { -0.5, 0.5, 0.5 },
            { -0.8872983346207419, 0.8872983346207419, 0.1127016653792581 } },
          { { -0.1127016653792584, 0.1127016653792584, 0.8872983346207416 },
            { -0.5000000000000002, 0.5000000000000002, 0.4999999999999998 },
            { -0.8872983346207419, 0.8872983346207419, 0.1127016653792581 } },
          { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } } }
    };
    static const double FE5_C4_F_Q8e7[1][3][3][3] = {
        { { { 0.8872983346207418, 0.1127016653792582, 0.8872983346207418 },
            { 0.5, 0.5, 0.5 },
            { 0.1127016653792582, 0.8872983346207418, 0.1127016653792582 } },
          { { 0.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 } },
          { { 0.1127016653792584, 0.8872983346207417, 0.1127016653792584 },
            { 0.5000000000000001, 0.4999999999999999, 0.5000000000000001 },
            { 0.8872983346207419, 0.1127016653792581, 0.8872983346207419 } } }
    };
    static const double triangle_reference_facet_jacobian[3][2][1] = { { { -1.0 }, { 1.0 } },
                                                                       { { 0.0 }, { 1.0 } },
                                                                       { { 1.0 }, { 0.0 } } };
    static const double triangle_reference_facet_normals[3][2] = {
        { 0.7071067811865475, 0.7071067811865475 }, { -1.0, -0.0 }, { 0.0, -1.0 }
    };
    // ------------------------
    // Section: Function
    // Inputs: w, FE5_C2_D10_F_Q8e7
    // Outputs: w0_d10_c2
    double w0_d10_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            w0_d10_c2 += w[( ic ) + 12] * FE5_C2_D10_F_Q8e7[0][0][0][ic];
        }
    }
    // ------------------------
    // ------------------------
    // Section: Jacobian
    // Inputs: FE5_C2_D10_F_Q8e7, coordinate_dofs, FE1_C1_D01_F_Q8e7
    // Outputs: J_c2, J_c3, J_c0, J_c1
    double J_c3 = 0.0;
    double J_c0 = 0.0;
    double J_c1 = 0.0;
    double J_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            J_c3 += coordinate_dofs[(ic)*3 + 1] * FE1_C1_D01_F_Q8e7[0][0][0][ic];
            J_c0 += coordinate_dofs[(ic)*3] * FE5_C2_D10_F_Q8e7[0][0][0][ic];
            J_c1 += coordinate_dofs[(ic)*3] * FE1_C1_D01_F_Q8e7[0][0][0][ic];
            J_c2 += coordinate_dofs[(ic)*3 + 1] * FE5_C2_D10_F_Q8e7[0][0][0][ic];
        }
    }
    // ------------------------
    // ------------------------
    // Section: Function
    // Inputs: w, FE1_C1_D01_F_Q8e7
    // Outputs: w0_d01_c2
    double w0_d01_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            w0_d01_c2 += w[( ic ) + 12] * FE1_C1_D01_F_Q8e7[0][0][0][ic];
        }
    }
    // ------------------------
    double sp_8e7_0 = J_c0 * J_c3;
    double sp_8e7_1 = J_c1 * J_c2;
    double sp_8e7_2 = -sp_8e7_1;
    double sp_8e7_3 = sp_8e7_0 + sp_8e7_2;
    double sp_8e7_4 = J_c3 / sp_8e7_3;
    double sp_8e7_5 = w0_d10_c2 * sp_8e7_4;
    double sp_8e7_6 = -J_c2;
    double sp_8e7_7 = sp_8e7_6 / sp_8e7_3;
    double sp_8e7_8 = w0_d01_c2 * sp_8e7_7;
    double sp_8e7_9 = sp_8e7_5 + sp_8e7_8;
    double sp_8e7_10 = J_c0 / sp_8e7_3;
    double sp_8e7_11 = triangle_reference_facet_normals[entity_local_index[0]][1] * sp_8e7_10;
    double sp_8e7_12 = -J_c1;
    double sp_8e7_13 = sp_8e7_12 / sp_8e7_3;
    double sp_8e7_14 = triangle_reference_facet_normals[entity_local_index[0]][0] * sp_8e7_13;
    double sp_8e7_15 = sp_8e7_11 + sp_8e7_14;
    double sp_8e7_16 = sp_8e7_15 * sp_8e7_15;
    double sp_8e7_17 = triangle_reference_facet_normals[entity_local_index[0]][0] * sp_8e7_4;
    double sp_8e7_18 = triangle_reference_facet_normals[entity_local_index[0]][1] * sp_8e7_7;
    double sp_8e7_19 = sp_8e7_17 + sp_8e7_18;
    double sp_8e7_20 = sp_8e7_19 * sp_8e7_19;
    double sp_8e7_21 = sp_8e7_16 + sp_8e7_20;
    double sp_8e7_22 = sqrt( sp_8e7_21 );
    double sp_8e7_23 = sp_8e7_15 / sp_8e7_22;
    double sp_8e7_24 = -sp_8e7_23;
    double sp_8e7_25 = w0_d01_c2 * sp_8e7_10;
    double sp_8e7_26 = w0_d10_c2 * sp_8e7_13;
    double sp_8e7_27 = sp_8e7_25 + sp_8e7_26;
    double sp_8e7_28 = sp_8e7_19 / sp_8e7_22;
    double sp_8e7_29 = sp_8e7_24 * sp_8e7_4;
    double sp_8e7_30 = sp_8e7_24 * sp_8e7_7;
    double sp_8e7_31 = sp_8e7_13 * sp_8e7_28;
    double sp_8e7_32 = sp_8e7_10 * sp_8e7_28;
    double sp_8e7_33 = sp_8e7_29 + sp_8e7_31;
    double sp_8e7_34 = sp_8e7_30 + sp_8e7_32;
    double sp_8e7_35 = -sp_8e7_4;
    double sp_8e7_36 = -sp_8e7_7;
    double sp_8e7_37 = -sp_8e7_24;
    double sp_8e7_38 = sp_8e7_35 * sp_8e7_24;
    double sp_8e7_39 = sp_8e7_36 * sp_8e7_24;
    double sp_8e7_40 = -sp_8e7_13;
    double sp_8e7_41 = -sp_8e7_10;
    double sp_8e7_42 = sp_8e7_40 * sp_8e7_28;
    double sp_8e7_43 = sp_8e7_41 * sp_8e7_28;
    double sp_8e7_44 = -sp_8e7_28;
    double sp_8e7_45 = sp_8e7_38 + sp_8e7_42;
    double sp_8e7_46 = sp_8e7_39 + sp_8e7_43;
    double sp_8e7_47 = J_c0 * triangle_reference_facet_jacobian[entity_local_index[0]][0][0];
    double sp_8e7_48 = J_c1 * triangle_reference_facet_jacobian[entity_local_index[0]][1][0];
    double sp_8e7_49 = sp_8e7_47 + sp_8e7_48;
    double sp_8e7_50 = sp_8e7_49 * sp_8e7_49;
    double sp_8e7_51 = triangle_reference_facet_jacobian[entity_local_index[0]][0][0] * J_c2;
    double sp_8e7_52 = triangle_reference_facet_jacobian[entity_local_index[0]][1][0] * J_c3;
    double sp_8e7_53 = sp_8e7_51 + sp_8e7_52;
    double sp_8e7_54 = sp_8e7_53 * sp_8e7_53;
    double sp_8e7_55 = sp_8e7_50 + sp_8e7_54;
    double sp_8e7_56 = sqrt( sp_8e7_55 );
    for ( int iq = 0; iq < 3; ++iq ) {
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C0_F_Q8e7
        // Outputs: w0_c0
        double w0_c0 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_c0 += w[(ic)*2] * FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C3_F_Q8e7
        // Outputs: w0_c3
        double w0_c3 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c3 += w[( ic ) + 15] * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C4_F_Q8e7
        // Outputs: w0_c4
        double w0_c4 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c4 += w[( ic ) + 15] * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C0_F_Q8e7
        // Outputs: w0_c1
        double w0_c1 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_c1 += w[(ic)*2 + 1] * FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C3_F_Q8e7
        // Outputs: w0_c5
        double w0_c5 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c5 += w[( ic ) + 18] * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C4_F_Q8e7
        // Outputs: w0_c6
        double w0_c6 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c6 += w[( ic ) + 18] * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Intermediates
        // Inputs: w0_c0, w0_c3, w0_c4, w0_c1, w0_c5, w0_c6
        // Outputs: fw0, fw1, fw2, fw3, fw4, fw5, fw6, fw7
        double fw0 = 0;
        double fw1 = 0;
        double fw2 = 0;
        double fw3 = 0;
        double fw4 = 0;
        double fw5 = 0;
        double fw6 = 0;
        double fw7 = 0;
        {
            double sv_8e7_0 = -w0_c0;
            double sv_8e7_1 = sp_8e7_9 + sv_8e7_0;
            double sv_8e7_2 = w0_c3 * sp_8e7_4;
            double sv_8e7_3 = w0_c4 * sp_8e7_7;
            double sv_8e7_4 = sv_8e7_2 + sv_8e7_3;
            double sv_8e7_5 = -sv_8e7_4;
            double sv_8e7_6 = sv_8e7_1 + sv_8e7_5;
            double sv_8e7_7 = sv_8e7_6 * sp_8e7_24;
            double sv_8e7_8 = -w0_c1;
            double sv_8e7_9 = sp_8e7_27 + sv_8e7_8;
            double sv_8e7_10 = w0_c4 * sp_8e7_10;
            double sv_8e7_11 = w0_c3 * sp_8e7_13;
            double sv_8e7_12 = sv_8e7_10 + sv_8e7_11;
            double sv_8e7_13 = -sv_8e7_12;
            double sv_8e7_14 = sv_8e7_9 + sv_8e7_13;
            double sv_8e7_15 = sv_8e7_14 * sp_8e7_28;
            double sv_8e7_16 = sv_8e7_7 + sv_8e7_15;
            double sv_8e7_17 = sv_8e7_16 * sp_8e7_33;
            double sv_8e7_18 = sv_8e7_16 * sp_8e7_34;
            double sv_8e7_19 = w0_c5 * sp_8e7_4;
            double sv_8e7_20 = w0_c6 * sp_8e7_7;
            double sv_8e7_21 = sv_8e7_19 + sv_8e7_20;
            double sv_8e7_22 = sv_8e7_21 * sp_8e7_24;
            double sv_8e7_23 = w0_c6 * sp_8e7_10;
            double sv_8e7_24 = w0_c5 * sp_8e7_13;
            double sv_8e7_25 = sv_8e7_23 + sv_8e7_24;
            double sv_8e7_26 = sv_8e7_25 * sp_8e7_28;
            double sv_8e7_27 = sv_8e7_22 + sv_8e7_26;
            double sv_8e7_28 = sv_8e7_27 * sp_8e7_33;
            double sv_8e7_29 = sv_8e7_27 * sp_8e7_34;
            double sv_8e7_30 = sv_8e7_27 * sp_8e7_37;
            double sv_8e7_31 = sv_8e7_27 * sp_8e7_45;
            double sv_8e7_32 = sv_8e7_27 * sp_8e7_46;
            double sv_8e7_33 = sv_8e7_27 * sp_8e7_44;
            double sv_8e7_34 = sv_8e7_17 * sp_8e7_56;
            double sv_8e7_35 = sv_8e7_18 * sp_8e7_56;
            double sv_8e7_36 = sv_8e7_28 * sp_8e7_56;
            double sv_8e7_37 = sv_8e7_29 * sp_8e7_56;
            double sv_8e7_38 = sv_8e7_30 * sp_8e7_56;
            double sv_8e7_39 = sv_8e7_31 * sp_8e7_56;
            double sv_8e7_40 = sv_8e7_32 * sp_8e7_56;
            double sv_8e7_41 = sv_8e7_33 * sp_8e7_56;
            fw0 = sv_8e7_38 * weights_8e7[iq];
            fw1 = sv_8e7_41 * weights_8e7[iq];
            fw2 = sv_8e7_36 * weights_8e7[iq];
            fw3 = sv_8e7_37 * weights_8e7[iq];
            fw4 = sv_8e7_39 * weights_8e7[iq];
            fw5 = sv_8e7_40 * weights_8e7[iq];
            fw6 = sv_8e7_34 * weights_8e7[iq];
            fw7 = sv_8e7_35 * weights_8e7[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw1, fw0, FE5_C0_F_Q8e7
        // Outputs: A
        {
            for ( int i = 0; i < 6; ++i ) {
                A[2 * ( i )] += fw0 * FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][i];
                A[( 2 * ( i ) + 1 )] += fw1 * FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][i];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: FE5_C4_F_Q8e7, fw7, fw3, FE1_C1_D01_F_Q8e7, fw2, fw6, FE5_C3_F_Q8e7, fw5, fw4,
        // FE5_C2_D10_F_Q8e7 Outputs: A
        {
            for ( int i = 0; i < 3; ++i ) {
                A[( ( i ) + 12 )] += fw2 * FE5_C2_D10_F_Q8e7[0][0][0][i];
                A[( ( i ) + 12 )] += fw3 * FE1_C1_D01_F_Q8e7[0][0][0][i];
                A[( ( i ) + 15 )] += fw4 * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 15 )] += fw5 * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 18 )] += fw6 * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 18 )] += fw7 * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][i];
            }
        }
        // ------------------------
    }
    return A;
}

VectorReal B_p4_tr6( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 21 * 21, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)
    // Quadrature rules
    static const double weights_39d[6] = { 0.054975871827661,  0.054975871827661,
                                           0.054975871827661,  0.1116907948390055,
                                           0.1116907948390055, 0.1116907948390055 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE1_C0_D10_Q39d[1][1][1][3] = { { { { -1.0, 1.0, 0.0 } } } };
    static const double FE1_C1_D01_Q39d[1][1][1][3] = { { { { -1.0, 0.0, 1.0 } } } };
    static const double FE5_C1_D01_Q39d[1][1][6][6] = {
        { { { 0.6336951459609189, 0.0, -0.6336951459609159, 3.267390291921835, 0.0,
              -3.267390291921836 },
            { 0.6336951459609192, 0.0, 2.267390291921836, 0.3663048540390839, -2.901085437882756,
              -0.3663048540390842 },
            { -2.267390291921832, 0.0, -0.6336951459609156, 0.366304854039084, 2.901085437882748,
              -0.3663048540390843 },
            { -0.7837939636638601, 0.0, 0.78379396366386, 0.43241207267228, 0.0,
              -0.4324120726722802 },
            { -0.7837939636638601, 0.0, -0.5675879273277196, 1.78379396366386, 1.35138189099158,
              -1.78379396366386 },
            { 0.5675879273277193, 0.0, 0.7837939636638599, 1.78379396366386, -1.35138189099158,
              -1.78379396366386 } } }
    };
    static const double FE5_C1_D10_Q39d[1][1][6][6] = {
        { { { 0.6336951459609198, 2.267390291921836, 0.0, 0.3663048540390836, -0.3663048540390846,
              -2.901085437882756 },
            { 0.6336951459609188, -0.633695145960915, 0.0, 3.267390291921836, -3.267390291921836,
              0.0 },
            { -2.267390291921832, -0.6336951459609163, 0.0, 0.3663048540390844, -0.3663048540390835,
              2.901085437882748 },
            { -0.7837939636638603, -0.5675879273277198, 0.0, 1.78379396366386, -1.78379396366386,
              1.35138189099158 },
            { -0.78379396366386, 0.78379396366386, 0.0, 0.43241207267228, -0.4324120726722801,
              0.0 },
            { 0.5675879273277195, 0.7837939636638604, 0.0, 1.78379396366386, -1.78379396366386,
              -1.35138189099158 } } }
    };
    static const double FE5_C3_Q39d[1][1][6][3] = {
        { { { -0.09157621350977084, 0.09157621350977084, 0.9084237864902291 },
            { -0.8168475729804592, 0.8168475729804592, 0.1831524270195409 },
            { -0.09157621350977109, 0.09157621350977109, 0.9084237864902289 },
            { -0.4459484909159651, 0.4459484909159651, 0.5540515090840348 },
            { -0.10810301816807, 0.10810301816807, 0.89189698183193 },
            { -0.445948490915965, 0.445948490915965, 0.554051509084035 } } }
    };
    static const double FE5_C4_Q39d[1][1][6][3] = {
        { { { 0.8168475729804592, 0.1831524270195409, 0.8168475729804592 },
            { 0.09157621350977091, 0.9084237864902291, 0.09157621350977091 },
            { 0.09157621350977108, 0.9084237864902289, 0.09157621350977108 },
            { 0.10810301816807, 0.89189698183193, 0.10810301816807 },
            { 0.4459484909159651, 0.5540515090840349, 0.4459484909159651 },
            { 0.445948490915965, 0.554051509084035, 0.445948490915965 } } }
    };
    // ------------------------
    // Section: Jacobian
    // Inputs: FE1_C0_D10_Q39d, FE1_C1_D01_Q39d, coordinate_dofs
    // Outputs: J_c2, J_c1, J_c0, J_c3
    double J_c0 = 0.0;
    double J_c3 = 0.0;
    double J_c1 = 0.0;
    double J_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            J_c0 += coordinate_dofs[(ic)*3] * FE1_C0_D10_Q39d[0][0][0][ic];
            J_c3 += coordinate_dofs[(ic)*3 + 1] * FE1_C1_D01_Q39d[0][0][0][ic];
            J_c1 += coordinate_dofs[(ic)*3] * FE1_C1_D01_Q39d[0][0][0][ic];
            J_c2 += coordinate_dofs[(ic)*3 + 1] * FE1_C0_D10_Q39d[0][0][0][ic];
        }
    }
    // ------------------------
    double sp_39d_0 = J_c0 * J_c3;
    double sp_39d_1 = J_c1 * J_c2;
    double sp_39d_2 = -sp_39d_1;
    double sp_39d_3 = sp_39d_0 + sp_39d_2;
    double sp_39d_4 = J_c0 / sp_39d_3;
    double sp_39d_5 = -J_c1;
    double sp_39d_6 = sp_39d_5 / sp_39d_3;
    double sp_39d_7 = sp_39d_4 * sp_39d_4;
    double sp_39d_8 = sp_39d_4 * sp_39d_6;
    double sp_39d_9 = sp_39d_6 * sp_39d_6;
    double sp_39d_10 = sp_39d_7 + sp_39d_7;
    double sp_39d_11 = sp_39d_8 + sp_39d_8;
    double sp_39d_12 = sp_39d_9 + sp_39d_9;
    double sp_39d_13 = J_c3 / sp_39d_3;
    double sp_39d_14 = -J_c2;
    double sp_39d_15 = sp_39d_14 / sp_39d_3;
    double sp_39d_16 = sp_39d_15 * sp_39d_15;
    double sp_39d_17 = sp_39d_13 * sp_39d_15;
    double sp_39d_18 = sp_39d_13 * sp_39d_13;
    double sp_39d_19 = sp_39d_16 + sp_39d_16;
    double sp_39d_20 = sp_39d_17 + sp_39d_17;
    double sp_39d_21 = sp_39d_18 + sp_39d_18;
    double sp_39d_22 = sp_39d_10 + sp_39d_19;
    double sp_39d_23 = sp_39d_11 + sp_39d_20;
    double sp_39d_24 = sp_39d_21 + sp_39d_12;
    double sp_39d_25 = c[0] * c[2];
    double sp_39d_26 = c[3] * sp_39d_25;
    double sp_39d_27 = 1.0 + c[1];
    double sp_39d_28 = 4.0 * sp_39d_27;
    double sp_39d_29 = sp_39d_26 / sp_39d_28;
    double sp_39d_30 = sp_39d_22 * sp_39d_29;
    double sp_39d_31 = sp_39d_23 * sp_39d_29;
    double sp_39d_32 = sp_39d_24 * sp_39d_29;
    double sp_39d_33 = sp_39d_4 + sp_39d_4;
    double sp_39d_34 = sp_39d_6 + sp_39d_6;
    double sp_39d_35 = sp_39d_33 / 2;
    double sp_39d_36 = sp_39d_34 / 2;
    double sp_39d_37 = sp_39d_35 * sp_39d_35;
    double sp_39d_38 = sp_39d_35 * sp_39d_36;
    double sp_39d_39 = sp_39d_36 * sp_39d_36;
    double sp_39d_40 = sp_39d_37 + sp_39d_37;
    double sp_39d_41 = sp_39d_38 + sp_39d_38;
    double sp_39d_42 = sp_39d_39 + sp_39d_39;
    double sp_39d_43 = sp_39d_15 / 2;
    double sp_39d_44 = sp_39d_13 / 2;
    double sp_39d_45 = sp_39d_4 / 2;
    double sp_39d_46 = sp_39d_6 / 2;
    double sp_39d_47 = sp_39d_43 * sp_39d_43;
    double sp_39d_48 = sp_39d_44 * sp_39d_43;
    double sp_39d_49 = sp_39d_45 * sp_39d_43;
    double sp_39d_50 = sp_39d_46 * sp_39d_43;
    double sp_39d_51 = sp_39d_44 * sp_39d_44;
    double sp_39d_52 = sp_39d_45 * sp_39d_44;
    double sp_39d_53 = sp_39d_44 * sp_39d_46;
    double sp_39d_54 = sp_39d_45 * sp_39d_45;
    double sp_39d_55 = sp_39d_45 * sp_39d_46;
    double sp_39d_56 = sp_39d_46 * sp_39d_46;
    double sp_39d_57 = sp_39d_47 + sp_39d_47;
    double sp_39d_58 = sp_39d_48 + sp_39d_48;
    double sp_39d_59 = sp_39d_49 + sp_39d_49;
    double sp_39d_60 = sp_39d_50 + sp_39d_50;
    double sp_39d_61 = sp_39d_51 + sp_39d_51;
    double sp_39d_62 = sp_39d_52 + sp_39d_52;
    double sp_39d_63 = sp_39d_53 + sp_39d_53;
    double sp_39d_64 = sp_39d_54 + sp_39d_54;
    double sp_39d_65 = sp_39d_55 + sp_39d_55;
    double sp_39d_66 = sp_39d_56 + sp_39d_56;
    double sp_39d_67 = sp_39d_40 + sp_39d_57;
    double sp_39d_68 = sp_39d_41 + sp_39d_58;
    double sp_39d_69 = sp_39d_42 + sp_39d_61;
    double sp_39d_70 = sp_39d_15 + sp_39d_15;
    double sp_39d_71 = sp_39d_13 + sp_39d_13;
    double sp_39d_72 = sp_39d_70 / 2;
    double sp_39d_73 = sp_39d_71 / 2;
    double sp_39d_74 = sp_39d_72 * sp_39d_72;
    double sp_39d_75 = sp_39d_73 * sp_39d_72;
    double sp_39d_76 = sp_39d_73 * sp_39d_73;
    double sp_39d_77 = sp_39d_74 + sp_39d_74;
    double sp_39d_78 = sp_39d_75 + sp_39d_75;
    double sp_39d_79 = sp_39d_76 + sp_39d_76;
    double sp_39d_80 = sp_39d_77 + sp_39d_64;
    double sp_39d_81 = sp_39d_78 + sp_39d_65;
    double sp_39d_82 = sp_39d_79 + sp_39d_66;
    double sp_39d_83 = sp_39d_67 + sp_39d_57;
    double sp_39d_84 = sp_39d_68 + sp_39d_58;
    double sp_39d_85 = sp_39d_59 + sp_39d_59;
    double sp_39d_86 = sp_39d_60 + sp_39d_60;
    double sp_39d_87 = sp_39d_69 + sp_39d_61;
    double sp_39d_88 = sp_39d_62 + sp_39d_62;
    double sp_39d_89 = sp_39d_63 + sp_39d_63;
    double sp_39d_90 = sp_39d_80 + sp_39d_64;
    double sp_39d_91 = sp_39d_81 + sp_39d_65;
    double sp_39d_92 = sp_39d_82 + sp_39d_66;
    double sp_39d_93 = -c[1];
    double sp_39d_94 = 1.0 + sp_39d_93;
    double sp_39d_95 = sp_39d_83 * sp_39d_94;
    double sp_39d_96 = sp_39d_84 * sp_39d_94;
    double sp_39d_97 = sp_39d_85 * sp_39d_94;
    double sp_39d_98 = sp_39d_86 * sp_39d_94;
    double sp_39d_99 = sp_39d_87 * sp_39d_94;
    double sp_39d_100 = sp_39d_88 * sp_39d_94;
    double sp_39d_101 = sp_39d_89 * sp_39d_94;
    double sp_39d_102 = sp_39d_90 * sp_39d_94;
    double sp_39d_103 = sp_39d_91 * sp_39d_94;
    double sp_39d_104 = sp_39d_92 * sp_39d_94;
    double sp_39d_105 = 2 * sp_39d_35;
    double sp_39d_106 = 2 * sp_39d_36;
    double sp_39d_107 = 2 * sp_39d_72;
    double sp_39d_108 = 2 * sp_39d_73;
    double sp_39d_109 = sp_39d_105 * sp_39d_35;
    double sp_39d_110 = sp_39d_106 * sp_39d_35;
    double sp_39d_111 = sp_39d_107 * sp_39d_35;
    double sp_39d_112 = sp_39d_108 * sp_39d_35;
    double sp_39d_113 = sp_39d_105 * sp_39d_36;
    double sp_39d_114 = sp_39d_106 * sp_39d_36;
    double sp_39d_115 = sp_39d_107 * sp_39d_36;
    double sp_39d_116 = sp_39d_108 * sp_39d_36;
    double sp_39d_117 = sp_39d_105 * sp_39d_72;
    double sp_39d_118 = sp_39d_106 * sp_39d_72;
    double sp_39d_119 = sp_39d_107 * sp_39d_72;
    double sp_39d_120 = sp_39d_108 * sp_39d_72;
    double sp_39d_121 = sp_39d_105 * sp_39d_73;
    double sp_39d_122 = sp_39d_106 * sp_39d_73;
    double sp_39d_123 = sp_39d_107 * sp_39d_73;
    double sp_39d_124 = sp_39d_108 * sp_39d_73;
    double sp_39d_125 = c[1] * sp_39d_109;
    double sp_39d_126 = c[1] * sp_39d_113;
    double sp_39d_127 = c[1] * sp_39d_117;
    double sp_39d_128 = c[1] * sp_39d_121;
    double sp_39d_129 = c[1] * sp_39d_110;
    double sp_39d_130 = c[1] * sp_39d_114;
    double sp_39d_131 = c[1] * sp_39d_118;
    double sp_39d_132 = c[1] * sp_39d_122;
    double sp_39d_133 = c[1] * sp_39d_111;
    double sp_39d_134 = c[1] * sp_39d_112;
    double sp_39d_135 = c[1] * sp_39d_115;
    double sp_39d_136 = c[1] * sp_39d_116;
    double sp_39d_137 = c[1] * sp_39d_119;
    double sp_39d_138 = c[1] * sp_39d_123;
    double sp_39d_139 = c[1] * sp_39d_120;
    double sp_39d_140 = c[1] * sp_39d_124;
    double sp_39d_141 = sp_39d_95 + sp_39d_125;
    double sp_39d_142 = sp_39d_96 + sp_39d_126;
    double sp_39d_143 = sp_39d_97 + sp_39d_127;
    double sp_39d_144 = sp_39d_98 + sp_39d_128;
    double sp_39d_145 = sp_39d_96 + sp_39d_129;
    double sp_39d_146 = sp_39d_99 + sp_39d_130;
    double sp_39d_147 = sp_39d_100 + sp_39d_131;
    double sp_39d_148 = sp_39d_101 + sp_39d_132;
    double sp_39d_149 = sp_39d_97 + sp_39d_133;
    double sp_39d_150 = sp_39d_98 + sp_39d_134;
    double sp_39d_151 = sp_39d_100 + sp_39d_135;
    double sp_39d_152 = sp_39d_101 + sp_39d_136;
    double sp_39d_153 = sp_39d_102 + sp_39d_137;
    double sp_39d_154 = sp_39d_103 + sp_39d_138;
    double sp_39d_155 = sp_39d_103 + sp_39d_139;
    double sp_39d_156 = sp_39d_104 + sp_39d_140;
    double sp_39d_157 = pow( c[3], 3 );
    double sp_39d_158 = c[0] * sp_39d_157;
    double sp_39d_159 = pow( c[1], 2 );
    double sp_39d_160 = -sp_39d_159;
    double sp_39d_161 = 1.0 + sp_39d_160;
    double sp_39d_162 = 24.0 * sp_39d_161;
    double sp_39d_163 = sp_39d_158 / sp_39d_162;
    double sp_39d_164 = sp_39d_141 * sp_39d_163;
    double sp_39d_165 = sp_39d_142 * sp_39d_163;
    double sp_39d_166 = sp_39d_143 * sp_39d_163;
    double sp_39d_167 = sp_39d_144 * sp_39d_163;
    double sp_39d_168 = sp_39d_145 * sp_39d_163;
    double sp_39d_169 = sp_39d_146 * sp_39d_163;
    double sp_39d_170 = sp_39d_147 * sp_39d_163;
    double sp_39d_171 = sp_39d_148 * sp_39d_163;
    double sp_39d_172 = sp_39d_149 * sp_39d_163;
    double sp_39d_173 = sp_39d_150 * sp_39d_163;
    double sp_39d_174 = sp_39d_151 * sp_39d_163;
    double sp_39d_175 = sp_39d_152 * sp_39d_163;
    double sp_39d_176 = sp_39d_153 * sp_39d_163;
    double sp_39d_177 = sp_39d_154 * sp_39d_163;
    double sp_39d_178 = sp_39d_155 * sp_39d_163;
    double sp_39d_179 = sp_39d_156 * sp_39d_163;
    double sp_39d_180 = fabs( sp_39d_3 );
    double sp_39d_181 = sp_39d_30 * sp_39d_180;
    double sp_39d_182 = sp_39d_31 * sp_39d_180;
    double sp_39d_183 = sp_39d_32 * sp_39d_180;
    double sp_39d_184 = sp_39d_164 * sp_39d_180;
    double sp_39d_185 = sp_39d_165 * sp_39d_180;
    double sp_39d_186 = sp_39d_166 * sp_39d_180;
    double sp_39d_187 = sp_39d_167 * sp_39d_180;
    double sp_39d_188 = sp_39d_168 * sp_39d_180;
    double sp_39d_189 = sp_39d_169 * sp_39d_180;
    double sp_39d_190 = sp_39d_170 * sp_39d_180;
    double sp_39d_191 = sp_39d_171 * sp_39d_180;
    double sp_39d_192 = sp_39d_172 * sp_39d_180;
    double sp_39d_193 = sp_39d_173 * sp_39d_180;
    double sp_39d_194 = sp_39d_174 * sp_39d_180;
    double sp_39d_195 = sp_39d_175 * sp_39d_180;
    double sp_39d_196 = sp_39d_176 * sp_39d_180;
    double sp_39d_197 = sp_39d_177 * sp_39d_180;
    double sp_39d_198 = sp_39d_178 * sp_39d_180;
    double sp_39d_199 = sp_39d_179 * sp_39d_180;
    for ( int iq = 0; iq < 6; ++iq ) {
        // ------------------------
        // Section: Intermediates
        // Inputs:
        // Outputs: fw0, fw1, fw2, fw3, fw4, fw5, fw6, fw7, fw8, fw9, fw10, fw11, fw12, fw13, fw14,
        // fw15, fw16, fw17, fw18
        double fw0 = 0;
        double fw1 = 0;
        double fw2 = 0;
        double fw3 = 0;
        double fw4 = 0;
        double fw5 = 0;
        double fw6 = 0;
        double fw7 = 0;
        double fw8 = 0;
        double fw9 = 0;
        double fw10 = 0;
        double fw11 = 0;
        double fw12 = 0;
        double fw13 = 0;
        double fw14 = 0;
        double fw15 = 0;
        double fw16 = 0;
        double fw17 = 0;
        double fw18 = 0;
        {
            fw0 = sp_39d_199 * weights_39d[iq];
            fw1 = sp_39d_198 * weights_39d[iq];
            fw2 = sp_39d_197 * weights_39d[iq];
            fw3 = sp_39d_196 * weights_39d[iq];
            fw4 = sp_39d_195 * weights_39d[iq];
            fw5 = sp_39d_193 * weights_39d[iq];
            fw6 = sp_39d_194 * weights_39d[iq];
            fw7 = sp_39d_192 * weights_39d[iq];
            fw8 = sp_39d_191 * weights_39d[iq];
            fw9 = sp_39d_190 * weights_39d[iq];
            fw10 = sp_39d_187 * weights_39d[iq];
            fw11 = sp_39d_186 * weights_39d[iq];
            fw12 = sp_39d_189 * weights_39d[iq];
            fw13 = sp_39d_188 * weights_39d[iq];
            fw14 = sp_39d_185 * weights_39d[iq];
            fw15 = sp_39d_184 * weights_39d[iq];
            fw16 = sp_39d_183 * weights_39d[iq];
            fw17 = sp_39d_182 * weights_39d[iq];
            fw18 = sp_39d_181 * weights_39d[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw12, fw0, fw13, fw6, FE5_C1_D01_Q39d, fw14, fw2, fw3, fw5, fw1, fw9, fw10,
        // FE5_C1_D10_Q39d, fw7, fw8, fw15, fw4, fw11 Outputs: A
        {
            double temp_0[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_0[j] = fw0 * FE5_C1_D10_Q39d[0][0][iq][j];
            }
            double temp_1[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_1[j] = fw1 * FE5_C1_D01_Q39d[0][0][iq][j];
            }
            double temp_2[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_2[j] = fw2 * FE5_C1_D10_Q39d[0][0][iq][j];
            }
            double temp_3[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_3[j] = fw3 * FE5_C1_D01_Q39d[0][0][iq][j];
            }
            double temp_4[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_4[j] = fw4 * FE5_C1_D10_Q39d[0][0][iq][j];
            }
            double temp_5[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_5[j] = fw5 * FE5_C1_D01_Q39d[0][0][iq][j];
            }
            double temp_6[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_6[j] = fw6 * FE5_C1_D10_Q39d[0][0][iq][j];
            }
            double temp_7[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_7[j] = fw7 * FE5_C1_D01_Q39d[0][0][iq][j];
            }
            double temp_8[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_8[j] = fw8 * FE5_C1_D10_Q39d[0][0][iq][j];
            }
            double temp_9[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_9[j] = fw9 * FE5_C1_D01_Q39d[0][0][iq][j];
            }
            double temp_10[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_10[j] = fw10 * FE5_C1_D10_Q39d[0][0][iq][j];
            }
            double temp_11[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_11[j] = fw11 * FE5_C1_D01_Q39d[0][0][iq][j];
            }
            double temp_12[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_12[j] = fw12 * FE5_C1_D10_Q39d[0][0][iq][j];
            }
            double temp_13[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_13[j] = fw13 * FE5_C1_D01_Q39d[0][0][iq][j];
            }
            double temp_14[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_14[j] = fw14 * FE5_C1_D10_Q39d[0][0][iq][j];
            }
            double temp_15[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_15[j] = fw15 * FE5_C1_D01_Q39d[0][0][iq][j];
            }
            for ( int j = 0; j < 6; ++j ) {
                for ( int i = 0; i < 6; ++i ) {
                    A[21 * ( 2 * ( i ) ) + 2 * ( j )] += FE5_C1_D10_Q39d[0][0][iq][i] * temp_0[j];
                    A[21 * ( 2 * ( i ) ) + 2 * ( j )] += FE5_C1_D10_Q39d[0][0][iq][i] * temp_1[j];
                    A[21 * ( 2 * ( i ) ) + 2 * ( j )] += FE5_C1_D01_Q39d[0][0][iq][i] * temp_2[j];
                    A[21 * ( 2 * ( i ) ) + 2 * ( j )] += FE5_C1_D01_Q39d[0][0][iq][i] * temp_3[j];
                    A[21 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D10_Q39d[0][0][iq][i] * temp_4[j];
                    A[21 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D10_Q39d[0][0][iq][i] * temp_5[j];
                    A[21 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D01_Q39d[0][0][iq][i] * temp_6[j];
                    A[21 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D01_Q39d[0][0][iq][i] * temp_7[j];
                    A[21 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE5_C1_D10_Q39d[0][0][iq][i] * temp_8[j];
                    A[21 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE5_C1_D10_Q39d[0][0][iq][i] * temp_9[j];
                    A[21 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE5_C1_D01_Q39d[0][0][iq][i] * temp_10[j];
                    A[21 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE5_C1_D01_Q39d[0][0][iq][i] * temp_11[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D10_Q39d[0][0][iq][i] * temp_12[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D10_Q39d[0][0][iq][i] * temp_13[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D01_Q39d[0][0][iq][i] * temp_14[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D01_Q39d[0][0][iq][i] * temp_15[j];
                }
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw17, fw18, FE5_C4_Q39d, fw16, FE5_C3_Q39d
        // Outputs: A
        {
            double temp_0[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_0[j] = fw16 * FE5_C3_Q39d[0][0][iq][j];
            }
            double temp_1[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_1[j] = fw17 * FE5_C4_Q39d[0][0][iq][j];
            }
            double temp_2[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_2[j] = fw17 * FE5_C3_Q39d[0][0][iq][j];
            }
            double temp_3[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_3[j] = fw18 * FE5_C4_Q39d[0][0][iq][j];
            }
            for ( int j = 0; j < 3; ++j ) {
                for ( int i = 0; i < 3; ++i ) {
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 15 )] += FE5_C3_Q39d[0][0][iq][i] * temp_0[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 15 )] += FE5_C3_Q39d[0][0][iq][i] * temp_1[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 15 )] += FE5_C4_Q39d[0][0][iq][i] * temp_2[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 15 )] += FE5_C4_Q39d[0][0][iq][i] * temp_3[j];
                }
            }
        }
        // ------------------------
    }

    return A;
}

VectorReal B_p5_tr6( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 21 * 21, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)
    // Quadrature rules
    static const double weights_8e7[3] = { 0.2777777777777777, 0.4444444444444444,
                                           0.2777777777777777 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE1_C1_D01_F_Q8e7[1][1][1][3] = { { { { -1.0, 0.0, 1.0 } } } };
    static const double FE5_C0_F_Q8e7[1][3][3][6] = {
        { { { 0.0, 0.6872983346207417, -0.08729833462074169, 0.4, 0.0, 0.0 },
            { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
            { 0.0, -0.08729833462074159, 0.6872983346207417, 0.4, 0.0, 0.0 } },
          { { 0.6872983346207417, 0.0, -0.08729833462074171, 0.0, 0.3999999999999999, 0.0 },
            { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
            { -0.08729833462074173, 0.0, 0.6872983346207417, 0.0, 0.4, 0.0 } },
          { { 0.6872983346207417, -0.0872983346207416, 0.0, 0.0, 0.0, 0.4 },
            { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
            { -0.08729833462074155, 0.6872983346207416, 0.0, 0.0, 0.0, 0.4 } } }
    };
    static const double FE5_C2_D10_F_Q8e7[1][1][1][3] = { { { { -1.0, 1.0, 0.0 } } } };
    static const double FE5_C3_F_Q8e7[1][3][3][3] = {
        { { { -0.1127016653792581, 0.1127016653792581, 0.8872983346207418 },
            { -0.5, 0.5, 0.5 },
            { -0.8872983346207419, 0.8872983346207419, 0.1127016653792581 } },
          { { -0.1127016653792584, 0.1127016653792584, 0.8872983346207416 },
            { -0.5000000000000002, 0.5000000000000002, 0.4999999999999998 },
            { -0.8872983346207419, 0.8872983346207419, 0.1127016653792581 } },
          { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } } }
    };
    static const double FE5_C4_F_Q8e7[1][3][3][3] = {
        { { { 0.8872983346207418, 0.1127016653792582, 0.8872983346207418 },
            { 0.5, 0.5, 0.5 },
            { 0.1127016653792582, 0.8872983346207418, 0.1127016653792582 } },
          { { 0.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 } },
          { { 0.1127016653792584, 0.8872983346207417, 0.1127016653792584 },
            { 0.5000000000000001, 0.4999999999999999, 0.5000000000000001 },
            { 0.8872983346207419, 0.1127016653792581, 0.8872983346207419 } } }
    };
    static const double triangle_reference_facet_jacobian[3][2][1] = { { { -1.0 }, { 1.0 } },
                                                                       { { 0.0 }, { 1.0 } },
                                                                       { { 1.0 }, { 0.0 } } };
    static const double triangle_reference_facet_normals[3][2] = {
        { 0.7071067811865475, 0.7071067811865475 }, { -1.0, -0.0 }, { 0.0, -1.0 }
    };
    // ------------------------
    // Section: Jacobian
    // Inputs: FE1_C1_D01_F_Q8e7, coordinate_dofs, FE5_C2_D10_F_Q8e7
    // Outputs: J_c2, J_c1, J_c0, J_c3
    double J_c3 = 0.0;
    double J_c0 = 0.0;
    double J_c1 = 0.0;
    double J_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            J_c3 += coordinate_dofs[(ic)*3 + 1] * FE1_C1_D01_F_Q8e7[0][0][0][ic];
            J_c0 += coordinate_dofs[(ic)*3] * FE5_C2_D10_F_Q8e7[0][0][0][ic];
            J_c1 += coordinate_dofs[(ic)*3] * FE1_C1_D01_F_Q8e7[0][0][0][ic];
            J_c2 += coordinate_dofs[(ic)*3 + 1] * FE5_C2_D10_F_Q8e7[0][0][0][ic];
        }
    }
    // ------------------------
    double sp_8e7_0 = J_c0 * J_c3;
    double sp_8e7_1 = J_c1 * J_c2;
    double sp_8e7_2 = -sp_8e7_1;
    double sp_8e7_3 = sp_8e7_0 + sp_8e7_2;
    double sp_8e7_4 = J_c3 / sp_8e7_3;
    double sp_8e7_5 = -J_c2;
    double sp_8e7_6 = sp_8e7_5 / sp_8e7_3;
    double sp_8e7_7 = -sp_8e7_4;
    double sp_8e7_8 = -sp_8e7_6;
    double sp_8e7_9 = J_c0 / sp_8e7_3;
    double sp_8e7_10 = triangle_reference_facet_normals[entity_local_index[0]][1] * sp_8e7_9;
    double sp_8e7_11 = -J_c1;
    double sp_8e7_12 = sp_8e7_11 / sp_8e7_3;
    double sp_8e7_13 = triangle_reference_facet_normals[entity_local_index[0]][0] * sp_8e7_12;
    double sp_8e7_14 = sp_8e7_10 + sp_8e7_13;
    double sp_8e7_15 = sp_8e7_14 * sp_8e7_14;
    double sp_8e7_16 = triangle_reference_facet_normals[entity_local_index[0]][0] * sp_8e7_4;
    double sp_8e7_17 = triangle_reference_facet_normals[entity_local_index[0]][1] * sp_8e7_6;
    double sp_8e7_18 = sp_8e7_16 + sp_8e7_17;
    double sp_8e7_19 = sp_8e7_18 * sp_8e7_18;
    double sp_8e7_20 = sp_8e7_15 + sp_8e7_19;
    double sp_8e7_21 = sqrt( sp_8e7_20 );
    double sp_8e7_22 = sp_8e7_14 / sp_8e7_21;
    double sp_8e7_23 = -sp_8e7_22;
    double sp_8e7_24 = sp_8e7_23 * sp_8e7_4;
    double sp_8e7_25 = sp_8e7_23 * sp_8e7_6;
    double sp_8e7_26 = -sp_8e7_23;
    double sp_8e7_27 = sp_8e7_7 * sp_8e7_23;
    double sp_8e7_28 = sp_8e7_8 * sp_8e7_23;
    double sp_8e7_29 = -sp_8e7_12;
    double sp_8e7_30 = -sp_8e7_9;
    double sp_8e7_31 = sp_8e7_18 / sp_8e7_21;
    double sp_8e7_32 = sp_8e7_12 * sp_8e7_31;
    double sp_8e7_33 = sp_8e7_9 * sp_8e7_31;
    double sp_8e7_34 = sp_8e7_29 * sp_8e7_31;
    double sp_8e7_35 = sp_8e7_30 * sp_8e7_31;
    double sp_8e7_36 = -sp_8e7_31;
    double sp_8e7_37 = sp_8e7_24 + sp_8e7_32;
    double sp_8e7_38 = sp_8e7_25 + sp_8e7_33;
    double sp_8e7_39 = sp_8e7_27 + sp_8e7_34;
    double sp_8e7_40 = sp_8e7_28 + sp_8e7_35;
    double sp_8e7_41 = sp_8e7_37 * sp_8e7_37;
    double sp_8e7_42 = sp_8e7_38 * sp_8e7_37;
    double sp_8e7_43 = sp_8e7_38 * sp_8e7_38;
    double sp_8e7_44 = sp_8e7_37 * sp_8e7_26;
    double sp_8e7_45 = sp_8e7_38 * sp_8e7_26;
    double sp_8e7_46 = sp_8e7_39 * sp_8e7_37;
    double sp_8e7_47 = sp_8e7_39 * sp_8e7_38;
    double sp_8e7_48 = sp_8e7_40 * sp_8e7_37;
    double sp_8e7_49 = sp_8e7_40 * sp_8e7_38;
    double sp_8e7_50 = sp_8e7_37 * sp_8e7_36;
    double sp_8e7_51 = sp_8e7_38 * sp_8e7_36;
    double sp_8e7_52 = J_c0 * triangle_reference_facet_jacobian[entity_local_index[0]][0][0];
    double sp_8e7_53 = J_c1 * triangle_reference_facet_jacobian[entity_local_index[0]][1][0];
    double sp_8e7_54 = sp_8e7_52 + sp_8e7_53;
    double sp_8e7_55 = sp_8e7_54 * sp_8e7_54;
    double sp_8e7_56 = triangle_reference_facet_jacobian[entity_local_index[0]][0][0] * J_c2;
    double sp_8e7_57 = triangle_reference_facet_jacobian[entity_local_index[0]][1][0] * J_c3;
    double sp_8e7_58 = sp_8e7_56 + sp_8e7_57;
    double sp_8e7_59 = sp_8e7_58 * sp_8e7_58;
    double sp_8e7_60 = sp_8e7_55 + sp_8e7_59;
    double sp_8e7_61 = sqrt( sp_8e7_60 );
    double sp_8e7_62 = sp_8e7_41 * sp_8e7_61;
    double sp_8e7_63 = sp_8e7_42 * sp_8e7_61;
    double sp_8e7_64 = sp_8e7_43 * sp_8e7_61;
    double sp_8e7_65 = sp_8e7_44 * sp_8e7_61;
    double sp_8e7_66 = sp_8e7_45 * sp_8e7_61;
    double sp_8e7_67 = sp_8e7_46 * sp_8e7_61;
    double sp_8e7_68 = sp_8e7_47 * sp_8e7_61;
    double sp_8e7_69 = sp_8e7_48 * sp_8e7_61;
    double sp_8e7_70 = sp_8e7_49 * sp_8e7_61;
    double sp_8e7_71 = sp_8e7_50 * sp_8e7_61;
    double sp_8e7_72 = sp_8e7_51 * sp_8e7_61;
    for ( int iq = 0; iq < 3; ++iq ) {
        // ------------------------
        // Section: Intermediates
        // Inputs:
        // Outputs: fw0, fw1, fw2, fw3, fw4, fw5, fw6, fw7, fw8, fw9, fw10
        double fw0 = 0;
        double fw1 = 0;
        double fw2 = 0;
        double fw3 = 0;
        double fw4 = 0;
        double fw5 = 0;
        double fw6 = 0;
        double fw7 = 0;
        double fw8 = 0;
        double fw9 = 0;
        double fw10 = 0;
        {
            fw0 = sp_8e7_65 * weights_8e7[iq];
            fw1 = sp_8e7_66 * weights_8e7[iq];
            fw2 = sp_8e7_71 * weights_8e7[iq];
            fw3 = sp_8e7_72 * weights_8e7[iq];
            fw4 = sp_8e7_62 * weights_8e7[iq];
            fw5 = sp_8e7_63 * weights_8e7[iq];
            fw6 = sp_8e7_64 * weights_8e7[iq];
            fw7 = sp_8e7_67 * weights_8e7[iq];
            fw8 = sp_8e7_68 * weights_8e7[iq];
            fw9 = sp_8e7_69 * weights_8e7[iq];
            fw10 = sp_8e7_70 * weights_8e7[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw0, fw2, fw3, FE5_C3_F_Q8e7, fw1, FE5_C0_F_Q8e7, FE5_C4_F_Q8e7
        // Outputs: A
        {
            double temp_0[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_0[j] = fw0 * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_1[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_1[j] = fw1 * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_2[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_2[j] = fw2 * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_3[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_3[j] = fw3 * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            for ( int j = 0; j < 3; ++j ) {
                for ( int i = 0; i < 6; ++i ) {
                    A[21 * ( 2 * ( i ) ) + ( ( j ) + 18 )] +=
                        FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_0[j];
                    A[21 * ( 2 * ( i ) ) + ( ( j ) + 18 )] +=
                        FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_1[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( ( j ) + 18 )] +=
                        FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_2[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( ( j ) + 18 )] +=
                        FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_3[j];
                }
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw6, FE5_C3_F_Q8e7, FE5_C2_D10_F_Q8e7, fw5, fw9, fw10, fw7, fw8, fw4,
        // FE1_C1_D01_F_Q8e7, FE5_C4_F_Q8e7 Outputs: A
        {
            double temp_0[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_0[j] = fw4 * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_1[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_1[j] = fw5 * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_2[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_2[j] = fw5 * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_3[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_3[j] = fw6 * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_4[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_4[j] = fw7 * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_5[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_5[j] = fw8 * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_6[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_6[j] = fw9 * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_7[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_7[j] = fw10 * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_8[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_8[j] = fw4 * FE5_C2_D10_F_Q8e7[0][0][0][j];
            }
            double temp_9[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_9[j] = fw5 * FE1_C1_D01_F_Q8e7[0][0][0][j];
            }
            double temp_10[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_10[j] = fw5 * FE5_C2_D10_F_Q8e7[0][0][0][j];
            }
            double temp_11[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_11[j] = fw6 * FE1_C1_D01_F_Q8e7[0][0][0][j];
            }
            double temp_12[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_12[j] = fw7 * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_13[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_13[j] = fw9 * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_14[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_14[j] = fw8 * FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_15[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_15[j] = fw10 * FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            for ( int j = 0; j < 3; ++j ) {
                for ( int i = 0; i < 3; ++i ) {
                    A[21 * ( ( i ) + 12 ) + ( ( j ) + 18 )] +=
                        FE5_C2_D10_F_Q8e7[0][0][0][i] * temp_0[j];
                    A[21 * ( ( i ) + 12 ) + ( ( j ) + 18 )] +=
                        FE5_C2_D10_F_Q8e7[0][0][0][i] * temp_1[j];
                    A[21 * ( ( i ) + 12 ) + ( ( j ) + 18 )] +=
                        FE1_C1_D01_F_Q8e7[0][0][0][i] * temp_2[j];
                    A[21 * ( ( i ) + 12 ) + ( ( j ) + 18 )] +=
                        FE1_C1_D01_F_Q8e7[0][0][0][i] * temp_3[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 18 )] +=
                        FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_4[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 18 )] +=
                        FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_5[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 18 )] +=
                        FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_6[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 18 )] +=
                        FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_7[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 12 )] +=
                        FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_8[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 12 )] +=
                        FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_9[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 12 )] +=
                        FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_10[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 12 )] +=
                        FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_11[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 15 )] +=
                        FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_12[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 15 )] +=
                        FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_13[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 15 )] +=
                        FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_14[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 15 )] +=
                        FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_15[j];
                }
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw0, fw2, fw3, FE5_C3_F_Q8e7, fw1, FE5_C0_F_Q8e7, FE5_C4_F_Q8e7
        // Outputs: A
        {
            double temp_0[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_0[j] = fw0 * FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_1[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_1[j] = fw1 * FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_2[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_2[j] = fw2 * FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            double temp_3[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_3[j] = fw3 * FE5_C0_F_Q8e7[0][entity_local_index[0]][iq][j];
            }
            for ( int j = 0; j < 6; ++j ) {
                for ( int i = 0; i < 3; ++i ) {
                    A[21 * ( ( i ) + 18 ) + 2 * ( j )] +=
                        FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_0[j];
                    A[21 * ( ( i ) + 18 ) + 2 * ( j )] +=
                        FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_1[j];
                    A[21 * ( ( i ) + 18 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C3_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_2[j];
                    A[21 * ( ( i ) + 18 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C4_F_Q8e7[0][entity_local_index[0]][iq][i] * temp_3[j];
                }
            }
        }
        // ------------------------
    }
    return A;
}

extern "C" {
void BP1_tr6_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A ) {
    if ( nw <= 0 || ncd <= 0 || ncst <= 0 || A == nullptr ) {
        std::cerr << "Erreur : tailles invalides dans K_elem_fortran." << std::endl;
        return;
    }

    // Création des vecteurs à partir des pointeurs C
    VectorReal w_vec( w, w + nw );
    VectorInt entities_vec( entities, entities + ne );
    VectorReal coordinate_dofs_vec( coordinate_dofs, coordinate_dofs + ncd );
    VectorReal cst_vec( cst, cst + ncst );

    // Appel de la fonction C++ de calcul
    VectorReal A_vec = B_p1_tr6( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
void BP2_tr6_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A ) {
    if ( nw <= 0 || ncd <= 0 || ncst <= 0 || A == nullptr ) {
        std::cerr << "Erreur : tailles invalides dans K_elem_fortran." << std::endl;
        return;
    }

    // Création des vecteurs à partir des pointeurs C
    VectorReal w_vec( w, w + nw );
    VectorInt entities_vec( entities, entities + ne );
    VectorReal coordinate_dofs_vec( coordinate_dofs, coordinate_dofs + ncd );
    VectorReal cst_vec( cst, cst + ncst );

    // Appel de la fonction C++ de calcul
    VectorReal A_vec = B_p2_tr6( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
void BP4_tr6_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A ) {
    if ( nw <= 0 || ncd <= 0 || ncst <= 0 || A == nullptr ) {
        std::cerr << "Erreur : tailles invalides dans K_elem_fortran." << std::endl;
        return;
    }

    // Création des vecteurs à partir des pointeurs C
    VectorReal w_vec( w, w + nw );
    VectorInt entities_vec( entities, entities + ne );
    VectorReal coordinate_dofs_vec( coordinate_dofs, coordinate_dofs + ncd );
    VectorReal cst_vec( cst, cst + ncst );

    // Appel de la fonction C++ de calcul
    VectorReal A_vec = B_p4_tr6( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
void BP5_tr6_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A ) {
    if ( nw <= 0 || ncd <= 0 || ncst <= 0 || A == nullptr ) {
        std::cerr << "Erreur : tailles invalides dans K_elem_fortran." << std::endl;
        return;
    }

    // Création des vecteurs à partir des pointeurs C
    VectorReal w_vec( w, w + nw );
    VectorInt entities_vec( entities, entities + ne );
    VectorReal coordinate_dofs_vec( coordinate_dofs, coordinate_dofs + ncd );
    VectorReal cst_vec( cst, cst + ncst );

    // Appel de la fonction C++ de calcul
    VectorReal A_vec = B_p5_tr6( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
}
