#pragma once
/* NSUBGRID is passed by the Makefile (-DNSUBGRID=...) */
#ifndef NSUBGRID
#error "NSUBGRID must be defined (-DNSUBGRID=...)"
#endif
#ifdef MHD
#if NSUBGRID != 2
#error "Metal MHD requires NSUBGRID=2"
#endif
#define MR_UCT_PRODUCT_FLUX_FIELDS 5
#define MR_UCT_PRODUCT_RECORD_FIELDS 5
#define MR_UCT_PRODUCT_FIELDS (MR_UCT_PRODUCT_FLUX_FIELDS + MR_UCT_PRODUCT_RECORD_FIELDS)
#define MR_UCT_PRODUCT_FACES (3 * (2 * NSUBGRID) * (2 * NSUBGRID) * (2 * NSUBGRID))
#define MR_SOLVER_LLF 1
#define MR_SOLVER_HLLD 4
#define MR_SOLVER_UCT_HLLD 7
#define MR_SOLVER2D_LLF 1
#define MR_SOLVER2D_HLLD 4
#else
#if NSUBGRID != 1
#error "Metal hydro requires NSUBGRID=1"
#endif
#endif
