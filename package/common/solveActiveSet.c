/**
 * Copyright (C) Till Blaha 2022-2023
 * MAVLab -- Faculty of Aerospace Engineering -- Delft University of Techology
 */

/**
 * @file solveActiveSet.c
 * 
 * @brief Implementations of solveActiveSet.h
*/

#include "solveActiveSet.h"

activeSetAlgo solveActiveSet(activeSetAlgoChoice choice) {
    switch (choice) {
        case AS_QR_NAIVE:
            return &solveActiveSet_qr_naive;
        case AS_CHOL:
            return &solveActiveSet_chol;
        case AS_CG:
            return &solveActiveSet_cg;
        default:
            return &solveActiveSet_qr;
    }
};
