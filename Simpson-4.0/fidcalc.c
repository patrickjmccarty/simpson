/*
 * fidcalc.c
 *
 *  Created on: Jun 24, 2010
 *      Author: zdenek
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "fidcalc.h"
#include "defs.h"
#include "cm.h"
#include "blockdiag.h"
#include "spinach.h"
#include "ham.h"
#include "pulse.h"
#include "fftw3.h"

	/* for acurate timings on windows */
#define TIMING
#include "timing.h"

// Result = +1 if there is symmetry, -1 otherwise
int is_rhosymmetry(mat_complx *fstart, mat_complx *fdetect)
{
	mat_complx *tmp;
	double sum;

	tmp = cm_adjoint(fdetect);
	cm_addto(tmp,fdetect);
	cm_muld(tmp,-0.5);
	cm_addto(tmp,fstart);
	sum = cm_sumnorm1(tmp);
	free_complx_matrix(tmp);
	return (sum < 0.0001 ? 1 : -1);
}

/* exit program if command not allowed inside acq_block */
void acqblock_disable_command(Sim_wsp *wsp, char *cmd)
{
	if (wsp->evalmode == EM_ACQBLOCK) {
		fprintf(stderr,"Error: acq_block - command '%s' not allowed\n",cmd);
		exit(1);
	}
}

void acqblock_sto_incr(Sim_wsp *wsp)
{
	if ( ++(wsp->acqblock_sto) >= ACQBLOCK_STO_END) {
		fprintf(stderr,"Error: acq_block - overflow of stored propagators\n");
		exit(1);
	}
	DEBUGPRINT(" - prop stored at time %g\n",wsp->t);
	return;
}

void modnumbers(double taur,double dt, int *k, int *l)
{
   double a,b;
   int m,n;

   if (fabs(taur-dt) < 1e-6) {
	   *k = *l = 1;
	   return;
   }

   a = taur;
   b = dt;
   m = 1;
   n = 1;
   if (taur>dt) {
      do {
         b += dt;
	 n++;
	 if (b-a > 1e-6) {
	    a += taur;
	    m++;
	 }
      } while ( fabs(a-b) > 1e-6 );
   } else {
      do {
         a += taur;
	 m++;
	 if (a-b > 1e-6) {
	    b += dt;
	    n++;
	 }
      } while ( fabs(a-b) > 1e-6 );
   }
   *k = m;
   *l = n;
}

void modnumbers_mas(double a, double b, double c, int *k, int *l, int *m)
{
	double v[3] = {a, b, c}, incr[3], d;
	int i, pos[3] = {0, 1, 2}, mul[3] = {1, 1, 1};

	/* sort the numbers first */
	if (v[0] > v[1]) {
		d = v[1]; v[1] = v[0]; v[0] = d;
		i = pos[1]; pos[1] = pos[0]; pos[0] = i;
	}
	if (v[1] > v[2]) {
		d = v[2]; v[2] = v[1]; v[1] = d;
		i = pos[2]; pos[2] = pos[1]; pos[1] = i;
	}
	if (v[0] > v[1]) {
		d = v[1]; v[1] = v[0]; v[0] = d;
		i = pos[1]; pos[1] = pos[0]; pos[0] = i;
	}
	for (i=0; i<3; i++) incr[i] = v[i];

	//printf("modnum: v = [ %f, %f, %f], pos = [ %d, %d, %d]\n",v[0],v[1],v[2],pos[0],pos[1],pos[2]);

	while ( fabs(v[0]-v[1])>TINY || fabs(v[1]-v[2])>TINY) {
		v[0] += incr[0];
		mul[0]++;
		if ( v[0]-v[1] > TINY) {
			v[1] += incr[1];
			mul[1]++;
			if ( v[1]-v[2] > TINY) {
				v[2] += incr[2];
				mul[2]++;
			}
		}
	}
	//printf("modnum: v = [ %f, %f, %f], pos = [ %d, %d, %d]\n",v[0],v[1],v[2],pos[0],pos[1],pos[2]);
	//printf("modnum: mul = [ %d, %d, %d]\n",mul[0],mul[1],mul[2]);
	//*k = mul[pos[0]];
	//*l = mul[pos[1]];
	//*m = mul[pos[2]];
	for (i=0; i<3; i++) {
		switch (pos[i]) {
		case 0: *k = mul[i]; break;
		case 1: *l = mul[i]; break;
		case 2: *m = mul[i]; break;
		}
	}
	//printf("modnum: k = %d, l= %d, m = %d\n",*k, *l, *m);
	//exit(1);
}


void scan_contrib_direct(mat_complx *rho, mat_complx **dets, int Ng, int *irow, int *icol)
{
	int i, j, k, contrib, nnz = 0;
	int matdim = rho->row;
	complx z1, z2;

	irow[0] = 1;
	for (i=1; i<=matdim; i++) {
		irow[i] = irow[i-1];
		for (j=1; j<=matdim; j++) {
			contrib = 0;
			z1 = cm_getelem(rho,j,i);
			if (fabs(z1.re) < TINY && fabs(z1.im) < TINY) continue;
			for (k=0; k<Ng; k++) {
				z2 = cm_getelem(dets[k],i,j);
				//printf("[%d,%d] = (%g,%g) * (%g,%g)\n",i,j,z1.re,z1.im,z2.re,z2.im);
				if ( fabs(z2.re) > TINY || fabs(z2.im) > TINY ) {
					contrib = 1;
					break;
				}
			}
			if (contrib != 0) {
				irow[i]++;
				icol[nnz] = j;
				nnz++;
				//printf("contrib [%d, %d]\n",i,j);
			}
		}
	}
}

void direct_acqblock(Tcl_Interp *interp,Tcl_Obj *obj,Sim_info *sim,Sim_wsp *wsp)
{
	int i, k=0, l=0, m=0, phase;
	double cosph, sinph, t0;
	complx *fidptr, z;
	mat_complx *rho, *det;
	blk_mat_complx *Ud, *T;

	/* acquire first point now */
	_evolve_with_prop(sim,wsp);
	_reset_prop(sim,wsp);
	fidptr = &(wsp->fid[++(wsp->curr_nsig)]);
	if ( fabs(wsp->acqphase) < TINY ) {
		phase = 0;
		cosph = 1;
		sinph = 0;
	} else {
		phase = 1;
		cosph = cos(wsp->acqphase*DEG2RAD);
		sinph = sin(wsp->acqphase*DEG2RAD);
	}
	if (sim->acq_adjoint == 0) {
		z = cm_trace(wsp->fdetect,wsp->sigma);
	} else {
		z = cm_trace_adjoint(wsp->fdetect,wsp->sigma);
	}
	if (phase) {
		fidptr->re += cosph*z.re+sinph*z.im;
		fidptr->im += -sinph*z.re+cosph*z.im;
	} else {
		fidptr->re += z.re;
		fidptr->im += z.im;
	}

	/* initialize propagator counters */
	wsp->acqblock_sto = ACQBLOCK_STO_INI;  // points to free position
	wsp->acqblock_t0 = t0 = wsp->t;
	wsp->evalmode = EM_ACQBLOCK;

	/* evaluate acq_block to get its propagator(s) */
	if (Tcl_EvalObjEx(interp,obj,0) != TCL_OK) {
		fprintf(stderr,"acq_block error: (1) can not execute block:\n'%s'\n\n",Tcl_GetString(obj));
		exit(1);
	}
	t0 = wsp->t - t0;
	if (t0 > wsp->dw*(wsp->Nacq-1)) {
		fprintf(stderr,"Error: acq_block is too long: duration %g > acquisition time %g\n",t0,wsp->dw*(wsp->Nacq-1));
		exit(1);
	}
	l = wsp->acqblock_sto - ACQBLOCK_STO_INI;
	DEBUGPRINT("acq_block: FIRST run, elapsed time %g us, %d dwelltime propagators created\n",t0,l);

	/* repeat acq_block to fill multiple of dwelltimes */
	if (sim->taur < 0) {
		/* static case */
		modnumbers(t0,wsp->dw,&k,&l);
		if (verbose & VERBOSE_ACQBLOCK) {
			printf("acq_block synchronization info: %d * acq_block(%g us) = %d * dwelltime(%g us)\n",k,t0,l,wsp->dw);
		}
	} else {
		/* we have MAS */
		modnumbers_mas(t0,wsp->dw,sim->taur,&k,&l,&m);
		if (verbose & VERBOSE_ACQBLOCK) {
			printf("acq_block synchronization info: %d * acq_block(%g us) = %d * dwelltime(%g us) = %d * tauR(%g us)\n",k,t0,l,wsp->dw,m,sim->taur);
		}
	}
	if (k*t0 > wsp->dw*(wsp->Nacq-1)) {
		fprintf(stderr,"\n\nWARNING: acq_block is NOT properly synchronized: Hamiltonian period %g > acquisition time %g\n\n",k*t0,wsp->dw*(wsp->Nacq-1));
		exit(1);
	}
	for (i=1; i<k; i++) {
		if (Tcl_EvalObjEx(interp,obj,0) != TCL_OK) {
			fprintf(stderr,"acq_block error: (%d) can not execute block:\n'%s'\n\n",i+1,Tcl_GetString(obj));
			exit(1);
		}
	}
	/* we should be done now with propagators */
	k = ACQBLOCK_STO_INI;       // first propagator
	m = wsp->acqblock_sto - 1;  // final propagator
	DEBUGPRINT("direct_acqblock: propagators from %d to %d\n",k,m);
	if ( wsp->acqblock_sto - ACQBLOCK_STO_INI != l) {
		fprintf(stderr,"Error: direct_acqblock - propagator count mismatch\n");
		exit(1);
	}

	/* make all propagators to be in the same basis as the final one */
	Ud = wsp->STO[m];
	for (i=k; i<m; i++) {
		if (wsp->STO[i]->basis == Ud->basis) continue;
		T = create_blk_mat_complx_copy(Ud);
		blk_cm_change_basis(T,wsp->STO[i],sim);
		free_blk_mat_complx(wsp->STO[i]);
		wsp->STO[i] = T;
	}
	//blk_cm_print(T,"U total i basis 0");
	//exit(1);

	/* prepare transformed matrices */
	//DEBUGPRINT("\t pointers wsp->fdetect = %p\n",wsp->fdetect);
	det = wsp->matrix[k] = cm_change_basis_2(wsp->fdetect,Ud->basis,sim);
	rho = cm_change_basis_2(wsp->sigma,Ud->basis,sim);
	if (blk_cm_isdiag(Ud) || (sim->sparse && Ud->dim > MAXDIMDIAGONALIZE) ) {
		DEBUGPRINT("\t\t prop is already diagonal, or DO NOT diagonalize big (>%d) sparse matrices\n",MAXDIMDIAGONALIZE);
		for (i=k; i<m; i++) {
			wsp->matrix[i+1] = cm_dup(det);
			blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
			free_blk_mat_complx(wsp->STO[i]);
			wsp->STO[i] = NULL;
		}
	} else {
		/* need to be diagonalized */
		DEBUGPRINT("\t\t need to diagonalize prop\n");
		//blk_cm_print(Ud,"will diagonalize this");
		T = blk_cm_diag(Ud);
		//blk_cm_print(Ud,"result"); blk_cm_print(T,"transformation");
		DEBUGPRINT("\t 333\n");
		for (i=k; i<m; i++) {
			DEBUGPRINT("\t transforming %d\n",i);
			wsp->matrix[i+1] = cm_dup(det);
			blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
			blk_simtrans_adj(wsp->matrix[i+1],T,sim);
			free_blk_mat_complx(wsp->STO[i]);
			wsp->STO[i] = NULL;
		}
		blk_simtrans_adj(rho,T,sim);
		blk_simtrans_adj(det,T,sim);
		free_blk_mat_complx(T);
	}

	// control print out
	//cm_print(rho,"rho");
	//for (i=k; i<=m; i++) cm_print(wsp->matrix[i],"dets");


	/* run the acquisition */
	DEBUGPRINT("\n----> acquisition\n");
	l = k+1;
	for (i=1; i<wsp->Nacq; i++) {
		if ( l > m ) {
			DEBUGPRINT("\t\t   period finished\n");
			l = k;
			blk_simtrans(rho,Ud,sim);
		}
		if (sim->acq_adjoint) {
			z = cm_trace_adjoint(wsp->matrix[l],rho);
		} else {
			z = cm_trace(wsp->matrix[l],rho);
		}
		fidptr++;
		if (phase) {
			fidptr->re += cosph*z.re+sinph*z.im;
			fidptr->im += -sinph*z.re+cosph*z.im;
		} else {
			fidptr->re += z.re;
			fidptr->im += z.im;
		}
		l++;
		DEBUGPRINT("\t\t %i is done\n",i);
	}
	wsp->curr_nsig += wsp->Nacq - 1;
	free_complx_matrix(rho);
	free_blk_mat_complx(wsp->STO[m]); wsp->STO[m] = NULL;
	for (i=k; i<=m; i++) {
		free_complx_matrix(wsp->matrix[i]);
		wsp->matrix[i] = NULL;
	}

	//cv_print(wsp->fid,"acq_block FID");
}

/****
 * this is for simulation in time domain with FWT interpolation
 ****/
void direct_acqblock_time_FWT(Tcl_Interp *interp,Tcl_Obj *obj,Sim_info *sim,Sim_wsp *wsp)
{
	int i, j, matdim, r, c, N, k=0, l=0, m=0;
	double t0;
	complx  z, zz;
	mat_complx *rho, *det;
	blk_mat_complx *Ud, *T, *TT;

	assert(sim->interpolation == 2); // run this ONLY for FWT interpol

	_evolve_with_prop(sim,wsp);
	_reset_prop(sim,wsp);

	// initialize propagator counters
	wsp->acqblock_sto = ACQBLOCK_STO_INI;  // points to free position
	wsp->acqblock_t0 = t0 = wsp->t;
	wsp->evalmode = EM_ACQBLOCK;

	// evaluate acq_block to get its propagator(s)
	if (Tcl_EvalObjEx(interp,obj,0) != TCL_OK) {
		fprintf(stderr,"acq_block error: (1) can not execute block:\n'%s'\n\n",Tcl_GetString(obj));
		exit(1);
	}
	t0 = wsp->t - t0;
	if (t0 > wsp->dw*(wsp->Nacq-1)) {
		fprintf(stderr,"Error: acq_block is too long: duration %g > acquisition time %g\n",t0,wsp->dw*(wsp->Nacq-1));
		exit(1);
	}
	l = wsp->acqblock_sto - ACQBLOCK_STO_INI;
	DEBUGPRINT("acq_block: FIRST run, elapsed time %g us, %d dwelltime propagators created\n",t0,l);

	// repeat acq_block to fill multiple of dwelltimes
	if (sim->taur < 0) {
		// static case
		modnumbers(t0,wsp->dw,&k,&l);
		if (verbose & VERBOSE_ACQBLOCK) {
			printf("acq_block synchronization info: %d * acq_block(%g us) = %d * dwelltime(%g us)\n",k,t0,l,wsp->dw);
		}
	} else {
		// we have MAS
		modnumbers_mas(t0,wsp->dw,sim->taur,&k,&l,&m);
		if (verbose & VERBOSE_ACQBLOCK) {
			printf("acq_block synchronization info: %d * acq_block(%g us) = %d * dwelltime(%g us) = %d * tauR(%g us)\n",k,t0,l,wsp->dw,m,sim->taur);
		}
	}
	t0 *= k;
	if (t0 > wsp->dw*(wsp->Nacq-1)) {
		fprintf(stderr,"\n\nWARNING: acq_block is NOT properly synchronized: Hamiltonian period %g > acquisition time %g\n\n",k*t0,wsp->dw*(wsp->Nacq-1));
		exit(1);
	}
	for (i=1; i<k; i++) {
		if (Tcl_EvalObjEx(interp,obj,0) != TCL_OK) {
			fprintf(stderr,"acq_block error: (%d) can not execute block:\n'%s'\n\n",i+1,Tcl_GetString(obj));
			exit(1);
		}
	}
	// we should be done now with propagators
	k = ACQBLOCK_STO_INI;       // first propagator
	m = wsp->acqblock_sto - 1;  // final propagator
	DEBUGPRINT("direct_acqblock: propagators from %d to %d\n",k,m);
	if ( wsp->acqblock_sto - ACQBLOCK_STO_INI != l) {
		fprintf(stderr,"Error: direct_acqblock - propagator count mismatch\n");
		exit(1);
	}

	// make all propagators to be in the same basis as the final one
	Ud = wsp->STO[m];
	for (i=k; i<m; i++) {
		if (wsp->STO[i]->basis == Ud->basis) continue;
		T = create_blk_mat_complx_copy(Ud);
		blk_cm_change_basis(T,wsp->STO[i],sim);
		free_blk_mat_complx(wsp->STO[i]);
		wsp->STO[i] = T;
	}

	// prepare transformed matrices
	//DEBUGPRINT("\t pointers wsp->fdetect = %p\n",wsp->fdetect);
	// prepare detection operator
	if (sim->acq_adjoint) {
		if (wsp->fdetect != sim->fdetect)
			cm_adjointi(wsp->fdetect);
		else
			wsp->fdetect = cm_adjoint(sim->fdetect);
	}
	det = wsp->matrix[k] = cm_change_basis_2(wsp->fdetect,Ud->basis,sim);
	rho = cm_change_basis_2(wsp->sigma,Ud->basis,sim);
	if (blk_cm_isdiag(Ud)) {
		T = NULL;
		for (i=k; i<m; i++) {
			wsp->matrix[i+1] = cm_dup(det);
			blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
			free_blk_mat_complx(wsp->STO[i]);
			wsp->STO[i] = NULL;
		}
	} else {
		DEBUGPRINT("\t\t need to diagonalize prop\n");
		//blk_cm_print(Ud,"will diagonalize this");
		T = blk_cm_diag(Ud);
		for (i=k; i<m; i++) {
			TT = blk_cm_mul(wsp->STO[i],T);
			free_blk_mat_complx(wsp->STO[i]);
			wsp->STO[i] = NULL;
			wsp->matrix[i+1] = cm_dup(det);
			blk_simtrans_adj(wsp->matrix[i+1],TT,sim);
			free_blk_mat_complx(TT);
		}
		blk_simtrans_adj(rho,T,sim);
		blk_simtrans_adj(det,T,sim);
		free_blk_mat_complx(T);
	}
	matdim = Ud->dim;

	// control print out
	//cm_print(rho,"rho");
	//for (i=k; i<=m; i++) cm_print(wsp->matrix[i],"dets");

    int *irow, *icol, *ic;
    irow = (int*)malloc((matdim+1)*sizeof(int));
    icol = ic = (int*)malloc(matdim*matdim*sizeof(int)); // just large enough
    N = m-k+1;
    scan_contrib_direct(rho,&wsp->matrix[k],N,irow,icol);

   	// now I am sure the file exists
   	if (ftell(wsp->interpol_file) == 0) {
   		// write header
   		r = irow[matdim]-1;
   		fwrite(&r,sizeof(int),1,wsp->interpol_file);
   		fwrite(&N,sizeof(int),1,wsp->interpol_file);
   		fwrite(&t0,sizeof(double),1,wsp->interpol_file);
   		fwrite(irow,sizeof(int),matdim+1,wsp->interpol_file);
   		fwrite(icol,sizeof(int),r,wsp->interpol_file);
   	}
   	// write crystallite index
   	if (wsp->ig == 0) fwrite(&wsp->cryst_idx,sizeof(int),1,wsp->interpol_file);
   	// write diagonalized period propagator
   	for (i=0; i<Ud->Nblocks; i++) {
   		fwrite(Ud->m[i].data,sizeof(complx),Ud->blk_dims[i],wsp->interpol_file);
   	}

	ic = icol;
	for (r=1; r<=matdim; r++) {
		int nc = irow[r] - irow[r-1];
		for (i=0; i<nc; i++) {
			c = *ic;
			zz = cm_getelem(rho,r,c);
			for (j=k; j<=m; j++) {
				z = Cmul(zz,cm_getelem(wsp->matrix[j],c,r));
				fwrite(&z,sizeof(complx),1,wsp->interpol_file);
			}
			ic++;
		}
	}

	free(irow);
	free(icol);
	free_complx_matrix(rho);
	free_blk_mat_complx(wsp->STO[m]); wsp->STO[m] = NULL;
	for (i=k; i<=m; i++) {
		free_complx_matrix(wsp->matrix[i]);
		wsp->matrix[i] = NULL;
	}
}


/****
 * this is for simulation in frequency domain
 ****/
void direct_acqblock_freq(Tcl_Interp *interp,Tcl_Obj *obj,Sim_info *sim,Sim_wsp *wsp)
{
	int i, j, k=0, m=0, phase, matdim, bin;
	double cosph, sinph, t0, period, *freq, freqT, *dptr, diff;
	complx *z1, z, zz;
	mat_complx *rho, *det;
	blk_mat_complx *Ud, *T, *TT;

	// prepare phase factor
	if ( fabs(wsp->acqphase) < TINY ) {
		phase = 0;
		cosph = 1;
		sinph = 0;
	} else {
		phase = 1;
		cosph = cos(wsp->acqphase*DEG2RAD);
		sinph = sin(wsp->acqphase*DEG2RAD);
	}

	// prepare initial density matrix
	_evolve_with_prop(sim,wsp);
	_reset_prop(sim,wsp);
	// prepare detection operator
	if (sim->acq_adjoint) {
		if (wsp->fdetect != sim->fdetect)
			cm_adjointi(wsp->fdetect);
		else
			wsp->fdetect = cm_adjoint(sim->fdetect);
	}

	// measure hamiltonian period - should be all included inside acq_block {}
	wsp->evalmode = EM_MEASURE;
	t0 = wsp->t;
	if (Tcl_EvalObjEx(interp,obj,0) != TCL_OK) {
		fprintf(stderr,"acq_block error: Measure time - can not execute block:\n'%s'\n\n",Tcl_GetString(obj));
		fprintf(stderr,"Error: %s\n",Tcl_GetStringResult(interp));
		exit(1);
	}
	period = wsp->t - t0;
	freqT = 2*M_PI*1e6/period;
	wsp->t = t0; // reset time counter
	// adjust dwelltime
	wsp->dw = period/sim->points_per_cycle;
	printf("period = %g, freqT = %g, N = %d, dw = %g\n",period,freqT,sim->points_per_cycle,wsp->dw);

	/* initialize propagator counters */
	wsp->acqblock_sto = ACQBLOCK_STO_INI;  // points to free position
	wsp->acqblock_t0 = t0 = wsp->t;
	wsp->evalmode = EM_ACQBLOCK;

	/* evaluate acq_block to get its propagator(s) */
	if (Tcl_EvalObjEx(interp,obj,0) != TCL_OK) {
		fprintf(stderr,"acq_block error: (1) can not execute block:\n'%s'\n\n",Tcl_GetString(obj));
		exit(1);
	}
	t0 = wsp->t - t0;
	if (fabs(t0-period) > TINY) {
		fprintf(stderr,"Error: acq_block : mismatch in durations, measured %g <> executed %g\n",period,t0);
		exit(1);
	}
	/* we should be done now with propagators */
	k = ACQBLOCK_STO_INI;       // first propagator
	m = wsp->acqblock_sto - 1;  // final propagator
	DEBUGPRINT("direct_acqblock: propagators from %d to %d\n",k,m);
	if ( wsp->acqblock_sto - ACQBLOCK_STO_INI != sim->points_per_cycle) {
		fprintf(stderr,"Error: direct_acqblock - freq - propagator count mismatch\n");
		exit(1);
	}

	/* make all propagators to be in the same basis as the final one */
	Ud = wsp->STO[m];
	for (i=k; i<m; i++) {
		if (wsp->STO[i]->basis == Ud->basis) continue;
		T = create_blk_mat_complx_copy(Ud);
		blk_cm_change_basis(T,wsp->STO[i],sim);
		free_blk_mat_complx(wsp->STO[i]);
		wsp->STO[i] = T;
	}
	//blk_cm_print(T,"U total i basis 0");
	//exit(1);

	/* prepare transformed matrices */
	//DEBUGPRINT("\t pointers wsp->fdetect = %p\n",wsp->fdetect);
	det = wsp->matrix[k] = cm_change_basis_2(wsp->fdetect,Ud->basis,sim);
	rho = cm_change_basis_2(wsp->sigma,Ud->basis,sim);
	if (blk_cm_isdiag(Ud)) {
		T = NULL;
	} else {
		DEBUGPRINT("\t\t need to diagonalize prop\n");
		//blk_cm_print(Ud,"will diagonalize this");
		T = blk_cm_diag(Ud);
	}
	matdim = Ud->dim;
	freq = dptr = (double*)malloc(matdim*sizeof(double));
	for (i=0; i<Ud->Nblocks; i++) {
		z1 = Ud->m[i].data;
		assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
		for (j=0; j<Ud->blk_dims[i]; j++) {
			*dptr = -Carg((*z1))/period*1.0e6;
			printf("freq[%d] = %g\n",j,*dptr);
			//*z1 = CRpow(*z1,1.0/(double)sim->points_per_cycle);
			//z1->re = cos(*dptr*period*1.0e-6/sim->points_per_cycle);
			//z1->im = -sin(*dptr*period*1.0e-6/sim->points_per_cycle);
			z1++; dptr++;
		}
	}

	if (T == NULL) {
		// Ud was diagonal
		for (i=k; i<m; i++) {
			//TT = blk_cm_mul(wsp->STO[i],Ud);
			//free_blk_mat_complx(wsp->STO[i]);
			//wsp->STO[i] = NULL;
			//wsp->matrix[i+1] = cm_dup(det);
			//blk_simtrans_adj(wsp->matrix[i+1],TT,sim);
			//free_blk_mat_complx(TT);
			wsp->matrix[i+1] = cm_dup(det);
			blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
			free_blk_mat_complx(wsp->STO[i]);
			wsp->STO[i] = NULL;
		}
	} else {
		// Ud was diagonalized
		for (i=k; i<m; i++) {
			//TT = blk_cm_mul(wsp->STO[i],T);
			//free_blk_mat_complx(wsp->STO[i]);
			//wsp->STO[i] = blk_cm_mul(TT,Ud);
			//free_blk_mat_complx(TT);
			//wsp->matrix[i+1] = cm_dup(det);
			//blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
			//free_blk_mat_complx(wsp->STO[i]);
			//wsp->STO[i] = NULL;
			TT = blk_cm_mul(wsp->STO[i],T);
			free_blk_mat_complx(wsp->STO[i]);
			wsp->STO[i] = NULL;
			wsp->matrix[i+1] = cm_dup(det);
			blk_simtrans_adj(wsp->matrix[i+1],TT,sim);
			free_blk_mat_complx(TT);
		}
		blk_simtrans_adj(rho,T,sim);
		blk_simtrans_adj(det,T,sim);
		free_blk_mat_complx(T);
	}

	// control print out
	//cm_print(rho,"rho");
	//for (i=k; i<=m; i++) cm_print(wsp->matrix[i],"dets");

	// construct the spectrum
	int r, c, N = sim->points_per_cycle;
	double binsize = sim->sw*2*M_PI/sim->np;
	printf("binsize = %g\n",binsize);
    fftw_complex *fftin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *fftout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p = fftw_plan_dft_1d(N, fftin, fftout, FFTW_FORWARD, FFTW_ESTIMATE);

    int *irow, *icol, *ic;
    irow = (int*)malloc((matdim+1)*sizeof(int));
    icol = ic = (int*)malloc(matdim*matdim*sizeof(int)); // just large enough
    scan_contrib_direct(rho,&wsp->matrix[k],N,irow,icol);

    if (sim->interpolation == 1) { // FWT interpolation
    	// now I am sure the file exists
    	if (ftell(wsp->interpol_file) == 0) {
    		// write header
    		r = irow[matdim]-1;
    		fwrite(&r,sizeof(int),1,wsp->interpol_file);
    		fwrite(&N,sizeof(int),1,wsp->interpol_file);
    		fwrite(&period,sizeof(double),1,wsp->interpol_file);
    		fwrite(irow,sizeof(int),matdim+1,wsp->interpol_file);
    		fwrite(icol,sizeof(int),r,wsp->interpol_file);
    	}
    	// write crystallite index
    	if (wsp->ig == 0) fwrite(&wsp->cryst_idx,sizeof(int),1,wsp->interpol_file);
    	// write diagonalized period propagator
    	for (i=0; i<Ud->Nblocks; i++) {
    		fwrite(Ud->m[i].data,sizeof(complx),Ud->blk_dims[i],wsp->interpol_file);
    	}
    } else if (sim->interpolation == 2) { // ASG interpolation
    	if (ftell(wsp->interpol_file) == 0) {
    		// write header
    		r = irow[matdim]-1;
    		fwrite(&r,sizeof(int),1,wsp->interpol_file);
    		fwrite(&N,sizeof(int),1,wsp->interpol_file);
    		fwrite(&period,sizeof(double),1,wsp->interpol_file);
    	}
    	// write crystallite index
    	if (wsp->ig == 0) fwrite(&wsp->cryst_idx,sizeof(int),1,wsp->interpol_file);
    }

	ic = icol;
	for (r=1; r<=matdim; r++) {
		int nc = irow[r] - irow[r-1];
		for (i=0; i<nc; i++) {
			c = *ic;
			zz = cm_getelem(rho,r,c);
			if (sim->interpolation == 1) { // FWT interpolation storage
				for (j=0; j<N; j++) {
					z = Cmul(zz,cm_getelem(wsp->matrix[k+j],c,r));
					fwrite(&z,sizeof(complx),1,wsp->interpol_file);
				}
				ic++;
				continue;
			}
			diff = freq[c-1] - freq[r-1];
			if (sim->interpolation == 2) { // ASG interpol
				fwrite(&diff,sizeof(double),1,wsp->interpol_file);
			}
			complx ph = Cexpi(-diff*period*1.0e-6/N);
			complx phmul = Complx(1.0,0.0);
			for (j=0; j<N; j++) {
				z = Cmul(phmul,cm_getelem(wsp->matrix[k+j],c,r));
				fftin[j][0] = zz.re*z.re - zz.im*z.im;
				fftin[j][1] = zz.im*z.re + zz.re*z.im;
				phmul = Cmul(phmul,ph);
			}
			fftw_execute(p);
			if (sim->interpolation ==  2) { // ASG interpol
				fwrite(fftout,sizeof(fftw_complex),N,wsp->interpol_file);
			}
			for (j=0; j<N; j++) {
				bin = (int)(1.5-(diff + freqT*(j-N/2+1))/binsize);
				printf("index %d -> freq = %g, bin %d -> ",j,diff+freqT*(j-N/2+1),bin);
				if (bin < 1) bin += sim->np;
				if (bin > sim->np) bin -= sim->np;
				printf("%d\n",bin);
				assert(bin >= 1 && bin <= sim->np);
				wsp->fid[bin].re += fftout[j][0]*cosph - fftout[j][1]*sinph;
				wsp->fid[bin].im += fftout[j][1]*cosph + fftout[j][0]*sinph;
			}
			ic++;
		}
	}

	free(irow);
	free(icol);
    fftw_destroy_plan(p);
    fftw_free(fftin); fftw_free(fftout);


	wsp->curr_nsig += wsp->Nacq - 1;
	free_complx_matrix(rho);
	for (i=k; i<=m; i++) {
		free_complx_matrix(wsp->matrix[i]);
		wsp->matrix[i] = NULL;
	}


}





void gcompute_acqblock(Tcl_Interp *interp,Tcl_Obj *obj,Sim_info *sim,Sim_wsp *wsp)
{
	double t0;
	int q = wsp->Uisunit;

	_evolve_with_prop(sim,wsp);
	_reset_prop(sim,wsp);

	t0 = wsp->t;
	/* evaluate acq_block to get its propagator */
	if (Tcl_EvalObjEx(interp,obj,0) != TCL_OK) {
		fprintf(stderr,"acq_block error: (igcompute %d) can not execute block:\n'%s'\n\n",wsp->acqblock_sto-ACQBLOCK_STO_INI+1,Tcl_GetString(obj));
		exit(1);
	}
	if (fabs(wsp->t - t0 - sim->taur/(double)sim->ngamma) > TINY) {
		fprintf(stderr,"Error: igcompute - acq_block event not synchronized with gamma-averaging\n");
		fprintf(stderr,"                   acq_block duration = %g us, need to be %g\n",wsp->t-t0,sim->taur/(double)sim->ngamma);
		exit(1);
	}
	/* store cumulative propagators, all should have the same basis */
	if (wsp->acqblock_sto == ACQBLOCK_STO_INI) {
		wsp->STO[wsp->acqblock_sto] = blk_cm_dup(wsp->U);
	} else {
		wsp->STO[wsp->acqblock_sto] = blk_cm_mul(wsp->U,wsp->STO[wsp->acqblock_sto - 1]);
	}
	//blk_cm_print(wsp->U,"gcompute_acqblock step U");
	//blk_cm_print(wsp->STO[wsp->acqblock_sto],"gcompute_acqblock U(i)");
	// store hopefully unchanged sigma in the basis of acq_block propagators
	//      but only if it changed...
	if (!q) wsp->matrix[wsp->acqblock_sto] = cm_change_basis_2(wsp->sigma,wsp->U->basis,sim);
	acqblock_sto_incr(wsp);
}

/****************
 * subroutine of new_gcompute - algorithm without diagonalization
 *   DO NOT call directly, it assumes some tasks already performed
 *   in new_gcompute
 */
void new_gcompute_sparse(Sim_info *sim,Sim_wsp *wsp,int k,int Ng,int m,double cosph,double sinph)
{
	blk_mat_complx *Uf, *Upow = NULL;
	int i, l, ig, F_idx, F_last, phase=0;
	complx *fidptr, z, zval;
	int j, jj, m1, m2;

	DEBUGPRINT("\n-->inside new_gcompute_sparse <---\n\n");
	if (sim->interpolation > 0) {
		fprintf(stderr,"Error: entering algorithm without diagonalization where interpolation is not allowed.\n");
		exit(1);
	}

	if (fabs(cosph) + fabs(sinph) > TINY) phase = 1;

	/* test for unitarity */
	//Uf = cm_adjoint(wsp->STO[k+Ng-1]);
	//Upow = cm_mul(Uf,wsp->STO[k+Ng-1]);
	//cm_print(Upow,"Utot+ * Utot");
	//exit(1);
	/* end of the test */

	l = k+Ng-1;
	/* continue with transforming the sigmas and dets */
	Uf = wsp->STO[l];
	mat_complx *det = wsp->matrix[k+Ng];
	for (i=k; i<l; i++) {
		blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
		wsp->matrix[i+Ng+1] = cm_dup(det);
		blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
		free_blk_mat_complx(wsp->STO[i]);
		wsp->STO[i] = NULL;
	}

	/* acquisition */
	F_idx = k+Ng;
	F_last = k+Ng+Ng-1;
	m1 = m / Ng; m2 = m % Ng;
	if (m1) Upow = blk_cm_power(Uf,m1);
	for (i=1; i<wsp->Nacq; i++) {
		if (m1 != 0) {
			for (j=F_idx; j<=F_last; j++) blk_simtrans_adj(wsp->matrix[j],Upow,sim);
		}
		jj = (i*m) % Ng + Ng;
		for (j=m2; j>0; j--) {
			blk_simtrans_adj(wsp->matrix[F_idx + (jj-j)%Ng],Uf,sim);
		}
		Czero(zval);
		for (ig=0; ig<Ng; ig++) {
			z = cm_trace(wsp->matrix[F_idx+(i*m+ig)%Ng],wsp->matrix[k+ig]);
			zval.re += z.re;
			zval.im += z.im;
		}
		fidptr = &(wsp->fid[++(wsp->curr_nsig)]);
		//printf("\nsparse gcompute: add (%f,%f) to fid (%f,%f)\n",zval.re,zval.im,fidptr->re,fidptr->im);
		if (phase) {
			fidptr->re += cosph*zval.re+sinph*zval.im;
			fidptr->im += -sinph*zval.re+cosph*zval.im;
		} else {
			fidptr->re += zval.re;
			fidptr->im += zval.im;
		}
	}
	cv_muld(wsp->fid, 1.0/(double)Ng);

	/* cleaning */
	if (Upow != NULL) free_blk_mat_complx(Upow);
	free_blk_mat_complx(wsp->STO[l]);
	wsp->STO[l] = NULL;
	for (i=k; i<=l; i++) {
		free_complx_matrix(wsp->matrix[i]);
		wsp->matrix[i] = NULL;
		free_complx_matrix(wsp->matrix[Ng+i]);
		wsp->matrix[Ng+i] = NULL;
	}
}

/**********
 * subroutine for new_gcompute
 * average over gamma angles of transformed matrices
 * res += rho^T .* det
 ***/
void gamma_ave(mat_complx *res, mat_complx **A, int ll, int Ng)
{
	// res is always dense
	int ig, i;
	const int dim = res->row;
	complx *z1, *z2, *z3;

	for (ig=0; ig<Ng; ig++) {
		mat_complx *rho = A[ig];
		mat_complx *det = A[Ng+((ig+ll)%Ng)];
		switch (100*rho->type+det->type) {
			case 100*MAT_DENSE+MAT_DENSE: {
				const int len = dim*dim;
				z1 = res->data;
				z3 = det->data;
				for (i=0; i<len; i++) {
					if (i%dim == 0) z2 = rho->data+i/dim;
					z1->re += z2->re*z3->re - z2->im*z3->im;
					z1->im += z2->re*z3->im + z2->im*z3->re;
					z1++;
					z3++;
					z2 += dim;
				}
				break;}
			case 100*MAT_DENSE+MAT_DENSE_DIAG: {
				z1 = res->data;
				z2 = rho->data;
				z3 = det->data;
				for (i=0; i<dim; i++) {
					z1->re += z2->re*z3->re - z2->im*z3->im;
					z1->im += z2->re*z3->im + z2->im*z3->re;
					z1 += dim+1;
					z2 += dim+1;
					z3++;
				}
				break; }
			case 100*MAT_DENSE+MAT_SPARSE_DIAG: {
				const int len = det->irow[dim]-1;
				z3 = det->data;
				for (i=0; i<len; i++) {
					int k = (det->icol[i]-1)*(dim+1);
					z1 = res->data + k;
					z2 = rho->data + k;
					z1->re += z2->re*z3->re - z2->im*z3->im;
					z1->im += z2->re*z3->im + z2->im*z3->re;
					z3++;
				}
				break; }
			case 100*MAT_DENSE_DIAG+MAT_DENSE: {
				z1 = res->data;
				z2 = rho->data;
				z3 = det->data;
				for (i=0; i<dim; i++) {
					z1->re += z2->re*z3->re - z2->im*z3->im;
					z1->im += z2->re*z3->im + z2->im*z3->re;
					z1 += dim+1;
					z2++;
					z3 += dim+1;
				}
				break; }
			case 100*MAT_DENSE_DIAG+MAT_DENSE_DIAG: {
				z1 = res->data;
				z2 = rho->data;
				z3 = det->data;
				for (i=0; i<dim; i++) {
					z1->re += z2->re*z3->re - z2->im*z3->im;
					z1->im += z2->re*z3->im + z2->im*z3->re;
					z1 += dim+1;
					z2++;
					z3++;
				}
				break; }
			case 100*MAT_DENSE_DIAG+MAT_SPARSE_DIAG: {
				const int len = det->irow[dim]-1;
				int *k = det->icol;
				z3 = det->data;
				for (i=0; i<len; i++) {
					z1 = res->data + (*k -1)*(dim+1);
					z2 = rho->data + *k-1;
					z1->re += z2->re*z3->re - z2->im*z3->im;
					z1->im += z2->re*z3->im + z2->im*z3->re;
					z3++;
					k++;
				}
				break; }
			case 100*MAT_SPARSE_DIAG+MAT_DENSE: {
				const int len = rho->irow[dim]-1;
				z2 = rho->data;
				for (i=0; i<len; i++) {
					int k = (rho->icol[i]-1)*(dim+1);
					z1 = res->data + k;
					z3 = det->data + k;
					z1->re += z2->re*z3->re - z2->im*z3->im;
					z1->im += z2->re*z3->im + z2->im*z3->re;
					z3++;
				}
				break; }
			case 100*MAT_SPARSE_DIAG+MAT_DENSE_DIAG: {
				const int len = rho->irow[dim]-1;
				z2 = rho->data;
				int *k = rho->icol;
				for (i=0; i<len; i++) {
					z1 = res->data + (*k -1)*(dim+1);
					z3 = det->data + *k -1;
					z1->re += z2->re*z3->re - z2->im*z3->im;
					z1->im += z2->re*z3->im + z2->im*z3->re;
					z2++;
					k++;
				}
				break; }
			case 100*MAT_SPARSE_DIAG+MAT_SPARSE_DIAG: {
				const int len = rho->irow[dim]-1;
				z2 = rho->data;
				int *k = rho->icol;
				for (i=0; i<len; i++) {
					z1 = res->data + (*k -1)*(dim+1);
					complx zz = cm_getelem(det,*k,*k);
					z1->re += z2->re*zz.re - z2->im*zz.im;
					z1->im += z2->re*zz.im + z2->im*zz.re;
					z2++;
					k++;
				}
				break; }
			default:
				fprintf(stderr,"gamma_ave error: invalid matrix types: rho=%s, det=%s\n",matrix_type(rho->type),matrix_type(det->type));
				exit(1);
		}
	}
	return;
}

void new_gcompute(Sim_info *sim, Sim_wsp *wsp)
{
	int ig, k, l, m, Ng, matdim, row, col, ll, phase, i, q;
	double dtg, t0, cosph, sinph;
	mat_complx *det, *mexp1, *mexp2;
	blk_mat_complx *Ud, *T;
	complx *fidptr, *z1, *z2, *z3, *z4, *z5, zz, zval;

	//LARGE_INTEGER tv1, tv2, tv3, tv4, tickpsec;
	//QueryPerformanceFrequency(&tickpsec);

	Ng = sim->ngamma;
	matdim = sim->matdim;

	/* first check on timings */
	dtg = sim->taur/(double)Ng;
	m = (int)floor(wsp->dw/dtg+1e-6);
	if ( fabs(dtg*m - wsp->dw) > TINY) {
		fprintf(stderr,"Error: (i)gcompute - bad synchronization of acquisition and gamma-averaging\n");
		fprintf(stderr,"       dwelltime(%g) can not be split in tauR(%g)/Ngamma(%d) steps\n",wsp->dw,sim->taur,sim->ngamma);
		exit(1);
	}

	/* initialize counter for storing wsp->sigma and wsp->U */
	wsp->acqblock_sto = k = ACQBLOCK_STO_INI;
	Czero(zval);

	/* evaluate pulse sequence to generate gamma-dependent sigma and props */
	//QueryPerformanceCounter(&tv1);
	for (ig=0;ig<Ng;ig++) {
		DEBUGPRINT("new_gcompute -> ig = %d\n",ig);
		wsp->tstart = t0 = ig/(double)Ng*sim->taur;
		direct_propagate(sim,wsp);
		if (wsp->acqblock_sto == ACQBLOCK_STO_INI + ig) {
			/* there was no acq_block command in pulseq, need to do the job here */
			if (fabs(wsp->t - t0 - dtg) > TINY) {
				fprintf(stderr,"Error: gcompute - pulse sequence event not synchronized with gamma-averaging\n");
				fprintf(stderr,"                  pulseq duration = %g us, need to be %g\n",wsp->t-t0,dtg);
				exit(1);
			}
			if (ig == 0) {
				wsp->STO[wsp->acqblock_sto] = blk_cm_dup(wsp->U);
				//wsp->Nacq = sim->np;
				//wsp->acqphase = 0;
			} else {
				wsp->STO[wsp->acqblock_sto] = blk_cm_mul(wsp->U,wsp->STO[wsp->acqblock_sto - 1]);
			}
			// this input file regime does not allow initial density matrix to be changed
			//wsp->matrix[wsp->acqblock_sto] = cm_change_basis_2(wsp->sigma,wsp->U->basis,sim); // initial density matrix is gamma-independent
			acqblock_sto_incr(wsp);
		}
		wsp->cryst.gamma += 360.0/(double)sim->ngamma;
		if (sim->interpolation == 0) {
			// basis of sigma and fdetect should remain unchanged and this is allowed:
			if (sim->acq_adjoint) {
				zz = cm_trace_adjoint(wsp->fdetect,wsp->sigma);
			} else {
				zz = cm_trace(wsp->fdetect,wsp->sigma);
			}
			zval.re += zz.re;
			zval.im += zz.im;
		}
	}
	//QueryPerformanceCounter(&tv2);
	//printf("timing eval_pulseq: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);

	l = wsp->acqblock_sto - 1;
	DEBUGPRINT("new_gcompute: sigmas and propagators from %d to %d\n",k,l);
	if (wsp->acqblock_sto - ACQBLOCK_STO_INI != Ng) {
		fprintf(stderr,"Error: (i)gcompute - propagator count mismatch\n");
		exit(1);
	}
	if (l+Ng >= ACQBLOCK_STO_END) {
		fprintf(stderr,"Error: (i)gcompute - overflow in matrix buffers\n");
		exit(1);
	}
	/* check fid length (better late than never) */
	if (wsp->curr_nsig + wsp->Nacq > LEN(wsp->fid)) {
		fprintf(stderr,"Error: (i)gcompute - acq overflow in fid points\n");
		exit(1);
	}
	if (sim->interpolation == 0) {
		/* store the first fid point */
		fidptr = &(wsp->fid[++(wsp->curr_nsig)]);
		if ( fabs(wsp->acqphase) < TINY ) {
			phase = 0;
			cosph = 0;
			sinph = 0;
			fidptr->re += zval.re;
			fidptr->im += zval.im;
		} else {
			phase = 1;
			cosph = cos(wsp->acqphase*DEG2RAD);
			sinph = sin(wsp->acqphase*DEG2RAD);
			fidptr->re += cosph*zval.re+sinph*zval.im;
			fidptr->im += -sinph*zval.re+cosph*zval.im;
		}
		DEBUGPRINT("new_gcompute: fid(0) = (%g, %g)\n",fidptr->re,fidptr->im);
	}
	/* transformations of dens.matrix and detect op. */
	//QueryPerformanceCounter(&tv1);
	Ud = wsp->STO[l];
	//blk_cm_print(Ud,"new_gcompute Ud");
	if (sim->acq_adjoint) {
		if (wsp->fdetect != sim->fdetect)
			cm_adjointi(wsp->fdetect);
		else
			wsp->fdetect = cm_adjoint(sim->fdetect);
	}
	det = wsp->matrix[l+1] = cm_change_basis_2(wsp->fdetect,Ud->basis,sim);
	if (wsp->matrix[k] == NULL) {
		//if (sim->EDsymmetry == 1) { // ED symmetry present
		//	if (blk_cm_isdiag(Ud)) {
		//		for (i=k; i<l; i++) {
		//			DEBUGPRINT("new_gcompute: Ud is diagonal\n");
		//			wsp->matrix[i+Ng+1] = cm_dup(det);
		//			blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
		//			free_blk_mat_complx(wsp->STO[i]);
		//			wsp->STO[i] = NULL;
		//		}
		//	} else {
		//		DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
		//		T = blk_cm_diag(Ud);
		//		for (i=k; i<l; i++) {
		//			// this should be faster:
		//			blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
		//			free_blk_mat_complx(wsp->STO[i]);
		//			wsp->STO[i] = NULL;
		//			wsp->matrix[i+Ng+1] = cm_dup(det);
		//			blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
		//			free_blk_mat_complx(dumblk);
		//		}
		//		blk_simtrans_adj(det,T,sim);
		//		free_blk_mat_complx(T);
		//	}
		//} else {  // no ED symmetry but full job needs to be done
			wsp->matrix[k] = cm_change_basis_2(wsp->fstart,Ud->basis,sim);
			if (blk_cm_isdiag(Ud)) {
				for (i=k; i<l; i++) {
					DEBUGPRINT("new_gcompute: Ud is diagonal\n");
					wsp->matrix[i+1] = cm_dup(wsp->matrix[k]);
					blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
				}
			} else {
				DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
				// need to diagonalize final propagator
				T = blk_cm_diag(Ud);
				for (i=k; i<l; i++) {
					// this should be faster:
					blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
					wsp->matrix[i+1] = cm_dup(wsp->matrix[k]);
					blk_simtrans_adj(wsp->matrix[i+1],dumblk,sim);
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
					free_blk_mat_complx(dumblk);
				}
				blk_simtrans_adj(det,T,sim);
				blk_simtrans_adj(wsp->matrix[k],T,sim);
				free_blk_mat_complx(T);
			}
			// control print out
			//for (i=0; i<Ng; i++) {
			//	cm_print(wsp->matrix[k+i],"RHO");
			//	cm_print(wsp->matrix[k+Ng+i],"DET");
			//}
		//}
	} else { // full job with gamma-dependent density matrices
		if (blk_cm_isdiag(Ud)) {
			for (i=k; i<l; i++) {
				DEBUGPRINT("new_gcompute: Ud is diagonal\n");
				blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
				wsp->matrix[i+Ng+1] = cm_dup(det);
				blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
				free_blk_mat_complx(wsp->STO[i]);
				wsp->STO[i] = NULL;
			}
		} else {
			DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
			// need to diagonalize final propagator
			T = blk_cm_diag(Ud);
			for (i=k; i<l; i++) {
				// this should be faster:
				blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
				free_blk_mat_complx(wsp->STO[i]);
				wsp->STO[i] = NULL;
				blk_simtrans_adj(wsp->matrix[i+1],dumblk,sim);
				wsp->matrix[i+Ng+1] = cm_dup(det);
				blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
				free_blk_mat_complx(dumblk);
			}
			blk_simtrans_adj(det,T,sim);
			blk_simtrans_adj(wsp->matrix[k],T,sim);
			free_blk_mat_complx(T);
		}
		// control print out
		//for (i=0; i<Ng; i++) {
		//	cm_print(wsp->matrix[k+i],"RHO");
		//	cm_print(wsp->matrix[k+Ng+i],"DET");
		//}
	}
	// END of transformations

	// control print out
	//for (i=0; i<Ng; i++) {
	//	cm_print(wsp->matrix[k+i],"RHO");
	//	cm_print(wsp->matrix[k+Ng+i],"DET");
	//}

	//DEBUGPRINT("new_gcompute: prepare mexps\n");
	//QueryPerformanceCounter(&tv1);
	mexp1 = complx_matrix(matdim,matdim,MAT_DENSE_DIAG,0,Ud->basis);
	mexp2 = complx_matrix(matdim,matdim,MAT_DENSE_DIAG,0,Ud->basis);
	z2 = mexp1->data;
	z3 = mexp2->data;
	for (i=0; i<Ud->Nblocks; i++) {
		z1 = Ud->m[i].data;
		assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
		for (ll=0; ll<Ud->blk_dims[i]; ll++) {
			*z2 = *z3 = CRpow(*z1,1.0/(double)Ng);
			z1++; z2++; z3++;
		}
	}
	//QueryPerformanceCounter(&tv2);
	//printf("timing mexp1 a mexp2: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
	//cm_print(Ud,"Ud");
	//cm_print(mexp1,"mexp1");

	/* averaging over gamma angles and mult by exp(-i w_rs ll tau) */
	//QueryPerformanceCounter(&tv1);
	DEBUGPRINT("new_gcompute: averaging over gamma\n");
	mat_complx **tmpmx = (mat_complx**)malloc(Ng*sizeof(mat_complx*));
	for (ll=0; ll<Ng; ll++) {
		tmpmx[ll] = complx_matrix(matdim,matdim,MAT_DENSE,0,Ud->basis);
		//cm_zero(wsp->STO[k+Ng+ll]);
		cm_zero(tmpmx[ll]);
		if (ll!=0) blk_simtrans_adj(wsp->matrix[k+Ng+ll-1],Ud,sim);
		gamma_ave(tmpmx[ll],&(wsp->matrix[k]),ll,Ng);
		//cm_print(tmpmx[ll],"Frs");
		if (ll!=0) {
			//simtrans(wsp->STO[k+Ng+ll],mexp1);
			simtrans(tmpmx[ll],mexp1);
			cm_multo(mexp1,mexp2);
		}
		//cm_print(tmpmx[ll],"Frs*exp(-i wrs dtg)");
	}
	//QueryPerformanceCounter(&tv2);
	//printf("timing all gamma ave: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);


/*	for (i=k; i<=l; i++) {
		printf("\n=== Fave %d ===\n",i);
		cm_print(wsp->STO[Ng+i],"Fave");
	}
*/
	if (sim->interpolation > 0) {
		// store data to file, no acquisition
		int nnztmp = 0;
		static int nnzcheck = 0;
		cm_sparse(tmpmx[0],1.0e-12); // try sparse format storage
		if (tmpmx[0]->type == MAT_SPARSE) {
			nnztmp = cm_nnz(tmpmx[0]);
			for (i=1; i<Ng; i++) {
				cm_sparse(tmpmx[i],1.0e-12);
				assert(tmpmx[i]->type == MAT_SPARSE); // make sure all matrices have the same pattern
				for (ll=0; ll<=matdim; ll++) {
					if (tmpmx[0]->irow[ll] != tmpmx[i]->irow[ll]) {
						fprintf(stderr,"Error: different structures of F_rs, please report your simulation to developer Z.T.\n");
						exit(1);
					}
				}
				if (cm_nnz(tmpmx[i]) != nnztmp) {
					fprintf(stderr,"Error: different structures of F_rs, please report your simulation to developer Z.T.\n");
					exit(1);
				}
				for (ll=0; ll<nnztmp; ll++) {
					if (tmpmx[0]->icol[ll] != tmpmx[i]->icol[ll]) {
						fprintf(stderr,"Error: different structures of F_rs, please report your simulation to developer Z.T.\n");
						exit(1);
					}
				}
			}
		} else {
			nnztmp = matdim*matdim;
		}
		long fpsize = ftell(wsp->interpol_file);
		if (fpsize == 0) {
			nnzcheck = nnztmp;
			printf("new_gcompute: fpsize = 0, nnztmp = %d, nnzcheck = %d\n",nnztmp,nnzcheck);
			// add the header part
			fwrite(&nnztmp,sizeof(int),1,wsp->interpol_file);
			//printf("%d ",nnztmp);
			if (tmpmx[0]->type == MAT_SPARSE) {
				fwrite(tmpmx[0]->irow,sizeof(int),matdim+1,wsp->interpol_file);
				//for (i=0; i<matdim+1; i++) printf("%d ",tmpmx[0]->irow[i]);
				fwrite(tmpmx[0]->icol,sizeof(int),nnztmp,wsp->interpol_file);
				//for (i=0; i<nnztmp; i++) printf("%d ",tmpmx[0]->icol[i]);
			}
		}
		if (nnzcheck != nnztmp) {
			fprintf(stderr,"Error: different nnz of F_rs from orientations, please report your simulation to developer Z.T.\n");
			printf("nnz = %d, nnzcheck = %d\n",nnztmp,nnzcheck);
			exit(1);
		}
		fwrite(&wsp->cryst_idx,sizeof(int),1,wsp->interpol_file);
		//printf("cr=%d fpos = %ld\n",wsp->cryst_idx,ftell(wsp->interpol_file));
		for (i=0; i<Ud->Nblocks; i++) {
			assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
			fwrite(Ud->m[i].data,sizeof(complx),Ud->blk_dims[i],wsp->interpol_file);
		}
		//blk_cm_print(Ud,"ORIGINAL Ud");
		for (i=0; i<Ng; i++) {
			fwrite(tmpmx[i]->data,sizeof(complx),nnztmp,wsp->interpol_file);
			//cm_print(tmpmx[i],"ORIGINAL Frs");
		}
	} else {
		/* acquisition */
		//QueryPerformanceCounter(&tv1);
		DEBUGPRINT("new_gcompute: new mexp's\n");
		const double m_Ng = (double)m/(double)Ng;
		z2 = mexp1->data;
		z3 = mexp2->data;
		for (i=0; i<Ud->Nblocks; i++) {
			z1 = Ud->m[i].data;
			assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
			for (ll=0; ll<Ud->blk_dims[i]; ll++) {
				//printf("complex(%g,%g) ",z1->re,z1->im);
				*z2 = *z3 = CRpow(*z1,m_Ng);
				z1++; z2++; z3++;
			}
		}
		//printf("\n");
		//QueryPerformanceCounter(&tv2);
		//printf("timing next mexp's: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
		//cm_print(mexp1,"mexp2");
		DEBUGPRINT("new_gcompute: m = %d\n",m);

		//QueryPerformanceCounter(&tv1);
		DEBUGPRINT("new_gcompute: acquisition\n");
		double zzre, zzim;
		for (ll=1; ll<wsp->Nacq; ll++) {
			Czero(zval);
			q = (ll*m)%Ng;
			//printf("\tll = %d, q = %d\n",ll,q);
			z3 = mexp1->data;
			z5 = mexp2->data;
			//complx *zstart = wsp->STO[k+Ng+q]->data;
			complx *zstart = tmpmx[q]->data;
			for (row=0; row<matdim; row++) {
				z1 = z2 = zstart + row*(matdim+1);
				zval.re += z1->re;
				zval.im += z1->im;
				z4 = z3;
				for (col=row+1; col<matdim; col++) {
					z1++;
					z2 += matdim;
					z4++;
					zzre = z3->re*z4->re + z3->im*z4->im;
					zzim = -z3->re*z4->im + z3->im*z4->re;
					zval.re += z1->re*zzre - z1->im*zzim;
					zval.im += z1->re*zzim + z1->im*zzre;
					zval.re += z2->re*zzre + z2->im*zzim;
					zval.im += -z2->re*zzim + z2->im*zzre;
					//if (fabs(z1->re) > TINY) {
					//	printf("r,c = %d, %d; exp = (%g, %g), F = (%g, %g)\n",row+1,col+1,zzre,zzim,z1->re,z1->im);
					//}
					//if (fabs(z2->re) > TINY) {
					//	printf("r,c = %d, %d; exp = (%g, %g), F = (%g, %g)\n",col+1,row+1,zzre,-zzim,z2->re,z2->im);
					//}

				}
				//*z3 = Cmul(*z3,*z5);
				zzre = z3->re; zzim = z3->im;
				z3->re = zzre*z5->re - zzim*z5->im;
				z3->im = zzre*z5->im + zzim*z5->re;

				z3++;
				z5++;
			}
			//printf("zval = (%g, %g)\n",zval.re,zval.im);
			fidptr++;
			//printf("\ndense gcompute: add (%f,%f) to fid (%f,%f)\n",zval.re,zval.im,fidptr->re,fidptr->im);
			if (phase) {
				fidptr->re = cosph*zval.re+sinph*zval.im;
				fidptr->im = -sinph*zval.re+cosph*zval.im;
			} else {
				fidptr->re = zval.re;
				fidptr->im = zval.im;
			}
		}
		wsp->curr_nsig += wsp->Nacq - 1;
		cv_muld(wsp->fid, 1.0/(double)Ng);
		//QueryPerformanceCounter(&tv2);
		//printf("timing acquisition: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
	} // end if: either interpolation store or acquisition

	//QueryPerformanceCounter(&tv1);
	free_complx_matrix(mexp1);
	free_complx_matrix(mexp2);
	free_blk_mat_complx(wsp->STO[l]);
	wsp->STO[l] = NULL;
	for (i=k; i<=l; i++) {
		free_complx_matrix(wsp->matrix[i]);
		wsp->matrix[i] = NULL;
		free_complx_matrix(wsp->matrix[i+Ng]);
		wsp->matrix[i+Ng] = NULL;
		//free_complx_matrix(wsp->STO[i+Ng]);
		//wsp->STO[i+Ng] = NULL;
		free_complx_matrix(tmpmx[i-k]);
	}
	free(tmpmx);
	//QueryPerformanceCounter(&tv2);
	//printf("timing cleaning: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
	DEBUGPRINT("new_gcompute: done!\n");
}



void scan_contrib(mat_complx **A, int Ng, int *irow, int *icol)
{
	int i, j, k, contrib, nnz = 0;
	int matdim = A[0]->row;
	complx z1, z2;

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);

	irow[0] = 1;
	for (i=1; i<=matdim; i++) {
		irow[i] = irow[i-1];
		for (j=1; j<=matdim; j++) {
			contrib = 0;
			for (k=0; k<Ng; k++) {
				z1 = cm_getelem(A[k],j,i);
				z2 = cm_getelem(A[k+Ng],i,j);
				//printf("[%d,%d] = (%g,%g) * (%g,%g)\n",i,j,z1.re,z1.im,z2.re,z2.im);
				if ( fabs(z1.re) > TINY && fabs(z2.re) > TINY ) {
					contrib = 1;
					break;
				}
				if ( fabs(z1.im) > TINY && fabs(z2.im) > TINY ) {
					contrib = 1;
					break;
				}
				if ( fabs(z1.re) > TINY && fabs(z2.im) > TINY ) {
					contrib = 1;
					break;
				}
				if ( fabs(z1.im) > TINY && fabs(z2.re) > TINY ) {
					contrib = 1;
					break;
				}
			}
			if (contrib != 0) {
				irow[i]++;
				icol[nnz] = j;
				nnz++;
				//printf("contrib [%d, %d]\n",i,j);
			}
		}
	}
	TIMING_TOC(tv1,tv2,"gCOMPUTE freq full scan contrib");
}

/* this was when interpolating Frs - gives artefacts
void new_gcompute_freq_fungujici(Sim_info *sim, Sim_wsp *wsp)
{
	int ig, k, l, m, Ng, matdim, r, c, i, bin;
	double dtg, t0, diff, cosph=1.0, sinph=0.0;
	mat_complx *det;
	blk_mat_complx *Ud, *T;
	complx *z1;
	double *freq, *dptr;
	int *irow, *icol, *ic;

	//LARGE_INTEGER tv1, tv2, tv3, tv4, tickpsec;
	//QueryPerformanceFrequency(&tickpsec);

	Ng = sim->ngamma;
	matdim = sim->matdim;
	dtg = sim->taur/(double)Ng;

	// initialize counter for storing wsp->sigma and wsp->U
	wsp->acqblock_sto = k = ACQBLOCK_STO_INI;

	// evaluate pulse sequence to generate gamma-dependent sigma and props
	//QueryPerformanceCounter(&tv1);
	for (ig=0;ig<Ng;ig++) {
		DEBUGPRINT("new_gcompute_freq -> ig = %d\n",ig);
		wsp->tstart = t0 = ig/(double)Ng*sim->taur;
		direct_propagate(sim,wsp);
		if (wsp->acqblock_sto == ACQBLOCK_STO_INI + ig) {
			// there was no acq_block command in pulseq, need to do the job here
			if (fabs(wsp->t - t0 - dtg) > TINY) {
				fprintf(stderr,"Error: gcompute - pulse sequence event not synchronized with gamma-averaging\n");
				fprintf(stderr,"                  pulseq duration = %g us, need to be %g\n",wsp->t-t0,dtg);
				exit(1);
			}
			if (ig == 0) {
				wsp->STO[wsp->acqblock_sto] = blk_cm_dup(wsp->U);
			} else {
				wsp->STO[wsp->acqblock_sto] = blk_cm_mul(wsp->U,wsp->STO[wsp->acqblock_sto - 1]);
			}
			wsp->matrix[wsp->acqblock_sto] = cm_change_basis_2(wsp->sigma,wsp->U->basis,sim); // initial density matrix is gamma-independent
			acqblock_sto_incr(wsp);
		}
		wsp->cryst.gamma += 360.0/(double)sim->ngamma;
	}
	//QueryPerformanceCounter(&tv2);
	//printf("timing eval_pulseq: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);

	l = wsp->acqblock_sto - 1;
	DEBUGPRINT("new_gcompute_freq: sigmas and propagators from %d to %d\n",k,l);
	if (wsp->acqblock_sto - ACQBLOCK_STO_INI != Ng) {
		fprintf(stderr,"Error: (i)gcompute - propagator count mismatch\n");
		exit(1);
	}
	if (l+Ng >= ACQBLOCK_STO_END) {
		fprintf(stderr,"Error: (i)gcompute - overflow in matrix buffers\n");
		exit(1);
	}
	if ( fabs(wsp->acqphase) > TINY ) {
		cosph = cos(wsp->acqphase*DEG2RAD);
		sinph = sin(wsp->acqphase*DEG2RAD);
	}

	// transformations of dens.matrix and detect op.
	Ud = wsp->STO[l];
	//blk_cm_print(Ud,"new_gcompute Ud");
	if (sim->acq_adjoint) {
		if (wsp->fdetect != sim->fdetect)
			cm_adjointi(wsp->fdetect);
		else
			wsp->fdetect = cm_adjoint(sim->fdetect);
	}
	det = wsp->matrix[l+1] = cm_change_basis_2(wsp->fdetect,Ud->basis,sim);
	//cm_print(det,"new_gcompute det ");
	if (blk_cm_isdiag(Ud)) {
		for (i=k; i<l; i++) {
			DEBUGPRINT("new_gcompute: Ud is diagonal\n");
			blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
			wsp->matrix[i+Ng+1] = cm_dup(det);
			blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
			free_blk_mat_complx(wsp->STO[i]);
			wsp->STO[i] = NULL;
		}
	} else {
		DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
		// need to diagonalize final propagator
		T = blk_cm_diag(Ud);
		for (i=k; i<l; i++) {
			// this should be faster:
			blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
			free_blk_mat_complx(wsp->STO[i]);
			wsp->STO[i] = NULL;
			blk_simtrans_adj(wsp->matrix[i+1],dumblk,sim);
			wsp->matrix[i+Ng+1] = cm_dup(det);
			blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
			free_blk_mat_complx(dumblk);
		}
		blk_simtrans_adj(det,T,sim);
		blk_simtrans_adj(wsp->matrix[k],T,sim);
		free_blk_mat_complx(T);
	}
	// control print out
	//for (i=0; i<Ng; i++) {
	//	cm_print(wsp->matrix[k+i],"RHO");
	//	cm_print(wsp->matrix[k+Ng+i],"DET");
	//}

	freq = dptr = (double*)malloc(matdim*sizeof(double));
	for (i=0; i<Ud->Nblocks; i++) {
		z1 = Ud->m[i].data;
		assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
		for (m=0; m<Ud->blk_dims[i]; m++) {
			*dptr = -Carg((*z1))/sim->taur*1.0e6;
			z1++;
			dptr++;
		}
	}
	// control print-out
	for (i=0; i<matdim; i++) printf("freq[%d] = %g\n",i,freq[i]);



	// construct the spectrum
	double binsize = sim->sw*2*M_PI/sim->np;
	printf("binsize = %g\n",binsize);
    fftw_complex *fftin1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
    fftw_complex *fftout1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
    fftw_plan p1 = fftw_plan_dft_1d(Ng, fftin1, fftout1, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_complex *fftin2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
    fftw_complex *fftout2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
    fftw_plan p2 = fftw_plan_dft_1d(Ng, fftin2, fftout2, FFTW_FORWARD, FFTW_ESTIMATE);

//
//    for (r=1; r<=matdim; r++) {
//   	for (c=1; c<=matdim; c++) {
//    		double znrm1 = 0.0;
//    		double znrm2 = 0.0;
//			diff = freq[r-1] - freq[c-1];
//			complx ph = Cexpi(-diff*dtg*1.0e-6);
//			complx phmul = Complx(1.0,0.0);
//    		for (i=0; i<Ng; i++) {
//    			complx zz1 = cm_getelem(wsp->matrix[k+i],c,r);
//    			complx zz2 = cm_getelem(wsp->matrix[k+Ng+i],r,c);
//    			znrm1 += fabs(zz1.re) + fabs(zz1.im);
//    			zz2 = Cmul(zz2,phmul);
//    			znrm2 += fabs(zz2.re) + fabs(zz2.im);
//    			fftin1[i][0] = zz1.re*phmul.re + zz1.im*phmul.im;
//    			fftin1[i][1] = zz1.re*phmul.im - zz1.im*phmul.re;
//    			fftin2[i][0] = zz2.re;
//    			fftin2[i][1] = zz2.im;
//    			phmul = Cmul(phmul,ph);
//    			printf("%g %g %g %g\n",zz1.re,zz1.im,zz2.re,zz2.im);
//    		}
//    		if (znrm1 < TINY || znrm2 < TINY) continue;
//    		printf("element (%d,%d) contributing\n",r,c);
//    		fftw_execute(p1);
//    		fftw_execute(p2);
//			for (i=0; i<Ng; i++) {
//				bin = (int)(1.5-(diff + sim->wr*(i-Ng/2+1))/binsize + sim->np/2);
//				printf("\tindex %d -> freq = %g, bin %d -> ",i,diff+sim->wr*(i-Ng/2+1),bin);
//				while (bin < 1) bin += sim->np;
//				while (bin > sim->np) bin -= sim->np;
//				assert(bin >= 1 && bin <= sim->np);
//				int idx = (i + Ng/2 + 1) % Ng;
//				double zzre = fftout1[idx][0]*fftout2[idx][0] + fftout1[idx][1]*fftout2[idx][1];
//				double zzim = fftout1[idx][0]*fftout2[idx][1] - fftout1[idx][1]*fftout2[idx][0];
//				wsp->fid[bin].re += zzre*cosph - zzim*sinph;
//				wsp->fid[bin].im += zzim*cosph + zzre*sinph;
//				printf("%d : %d (%g, %g) (%g, %g) (%g, %g)\n",bin,idx,fftout1[idx][0],fftout1[idx][1],fftout2[idx][0],fftout2[idx][1],zzre,zzim);
//			}
//   	}
//    }
//

    irow = (int*)malloc((matdim+1)*sizeof(int));
    icol = ic = (int*)malloc(matdim*matdim*sizeof(int)); // just large enough
    scan_contrib(&wsp->matrix[k],Ng,irow,icol);
    //printf("\n\nBUF!\n\n");
    if (sim->interpolation > 0) {
    	// now I am sure the file exists
    	if (ftell(wsp->interpol_file) == 0) {
    		// write header
    		r = irow[matdim]-1;
    		fwrite(&r,sizeof(int),1,wsp->interpol_file);
    		fwrite(irow,sizeof(int),matdim+1,wsp->interpol_file);
    		fwrite(icol,sizeof(int),r,wsp->interpol_file);
    	}
    	// write crystallite index
    	fwrite(&wsp->cryst_idx,sizeof(int),1,wsp->interpol_file);
    	// write diagonalized period propagator
    	for (i=0; i<Ud->Nblocks; i++) {
    		fwrite(Ud->m[i].data,sizeof(complx),Ud->blk_dims[i],wsp->interpol_file);
    	}
    }
    for (r=1; r<=matdim; r++) {
    	int nc = irow[r] - irow[r-1];
    	for (c=0; c<nc; c++) {
			diff = freq[r-1] - freq[*ic-1];
			complx ph = Cexpi(-diff*dtg*1.0e-6);
			complx phmul = Complx(1.0,0.0);
    		for (i=0; i<Ng; i++) {
    			complx zz1 = cm_getelem(wsp->matrix[k+i],*ic,r);
    			complx zz2 = cm_getelem(wsp->matrix[k+Ng+i],r,*ic);
    			fftin1[i][0] = zz1.re*phmul.re + zz1.im*phmul.im;
    			fftin1[i][1] = zz1.re*phmul.im - zz1.im*phmul.re;
    			fftin2[i][0] = zz2.re*phmul.re - zz2.im*phmul.im;
    			fftin2[i][1] = zz2.re*phmul.im + zz2.im*phmul.re;
    			phmul = Cmul(phmul,ph);
    		}
    		fftw_execute(p1);
    		fftw_execute(p2);
			for (i=0; i<Ng; i++) {
				int idx = (i + Ng/2 + 1) % Ng;
				double zzre = fftout1[idx][0]*fftout2[idx][0] + fftout1[idx][1]*fftout2[idx][1];
				double zzim = fftout1[idx][0]*fftout2[idx][1] - fftout1[idx][1]*fftout2[idx][0];
				if (sim->interpolation > 0) {
					fwrite(&zzre,sizeof(double),1,wsp->interpol_file);
					fwrite(&zzim,sizeof(double),1,wsp->interpol_file);
				} else {
					bin = (int)(1.5-(diff + sim->wr*(i-Ng/2+1))/binsize + sim->np/2);
					printf("\tindex %d -> freq = %g, bin %d -> ",i,diff+sim->wr*(i-Ng/2+1),bin);
					while (bin < 1) bin += sim->np;
					while (bin > sim->np) bin -= sim->np;
					assert(bin >= 1 && bin <= sim->np);
					wsp->fid[bin].re += zzre*cosph - zzim*sinph;
					wsp->fid[bin].im += zzim*cosph + zzre*sinph;
					printf("%d : %d (%g, %g) (%g, %g) (%g, %g)\n",bin,idx,fftout1[idx][0],fftout1[idx][1],fftout2[idx][0],fftout2[idx][1],zzre,zzim);
				}
			}
			ic++;
    	}
    }

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(fftin1); fftw_free(fftout1);
    fftw_free(fftin2); fftw_free(fftout2);
    free(irow);
    free(icol);
    free(freq);

	//free_complx_matrix(mexp1);
	free_blk_mat_complx(wsp->STO[l]);
	wsp->STO[l] = NULL;
	for (i=k; i<=l; i++) {
		free_complx_matrix(wsp->matrix[i]);
		wsp->matrix[i] = NULL;
		free_complx_matrix(wsp->matrix[i+Ng]);
		wsp->matrix[i+Ng] = NULL;
	}

}
*****/

/****
 * testovani, zda interpolace bez frekvenci da vysledek bez artefaktu
 *
 */
void new_gcompute_freq(Sim_info *sim, Sim_wsp *wsp)
{
	int ig, k, l, m, Ng, matdim, r, c, i, bin;
	double dtg, t0, diff, cosph=1.0, sinph=0.0;
	mat_complx *det;
	blk_mat_complx *Ud, *T;
	complx *z1;
	double *freq, *dptr;
	int *irow, *icol, *ic;

	//LARGE_INTEGER tv1, tv2, tv3, tv4, tickpsec;
	//QueryPerformanceFrequency(&tickpsec);
	TIMING_INIT_VAR(tv1);
	TIMING_INIT_VAR(tv2);
	TIMING_INIT;

	Ng = sim->ngamma;
	matdim = sim->matdim;
	dtg = sim->taur/(double)Ng;

	/* initialize counter for storing wsp->sigma and wsp->U */
	wsp->acqblock_sto = k = ACQBLOCK_STO_INI;

	// evaluate pulse sequence to generate gamma-dependent sigma and props
	//QueryPerformanceCounter(&tv1);
	TIMING_TIC(tv1);
	for (ig=0;ig<Ng;ig++) {
		DEBUGPRINT("new_gcompute_freq -> ig = %d\n",ig);
		wsp->tstart = t0 = ig/(double)Ng*sim->taur;
		direct_propagate(sim,wsp);
		if (wsp->acqblock_sto == ACQBLOCK_STO_INI + ig) {
			// there was no acq_block command in pulseq, need to do the job here
			if (fabs(wsp->t - t0 - dtg) > TINY) {
				fprintf(stderr,"Error: gcompute - pulse sequence event not synchronized with gamma-averaging\n");
				fprintf(stderr,"                  pulseq duration = %g us, need to be %g\n",wsp->t-t0,dtg);
				exit(1);
			}
			if (ig == 0) {
				wsp->STO[wsp->acqblock_sto] = blk_cm_dup(wsp->U);
			} else {
				wsp->STO[wsp->acqblock_sto] = blk_cm_mul(wsp->U,wsp->STO[wsp->acqblock_sto - 1]);
			}
			// this input file regime does not allow initial density matrix to be changed
			//wsp->matrix[wsp->acqblock_sto] = cm_change_basis_2(wsp->sigma,wsp->U->basis,sim); // initial density matrix is gamma-independent
			acqblock_sto_incr(wsp);
		}
		wsp->cryst.gamma += 360.0/(double)sim->ngamma;
	}
	//QueryPerformanceCounter(&tv2);
	//printf("timing eval_pulseq: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
	TIMING_TOC(tv1,tv2,"eval_pulseq");

	l = wsp->acqblock_sto - 1;
	DEBUGPRINT("new_gcompute_freq: sigmas and propagators from %d to %d\n",k,l);
	if (wsp->acqblock_sto - ACQBLOCK_STO_INI != Ng) {
		fprintf(stderr,"Error: (i)gcompute - propagator count mismatch\n");
		exit(1);
	}
	if (l+Ng >= ACQBLOCK_STO_END) {
		fprintf(stderr,"Error: (i)gcompute - overflow in matrix buffers\n");
		exit(1);
	}
	if ( fabs(wsp->acqphase) > TINY ) {
		cosph = cos(wsp->acqphase*DEG2RAD);
		sinph = sin(wsp->acqphase*DEG2RAD);
	}

	// transformations of dens.matrix and detect op.
	Ud = wsp->STO[l];
	//blk_cm_print(Ud,"new_gcompute Ud");
	if (sim->acq_adjoint) {
		if (wsp->fdetect != sim->fdetect)
			cm_adjointi(wsp->fdetect);
		else
			wsp->fdetect = cm_adjoint(sim->fdetect);
	}
	det = wsp->matrix[l+1] = cm_change_basis_2(wsp->fdetect,Ud->basis,sim);
	//QueryPerformanceCounter(&tv1);
	TIMING_TIC(tv1);
	if (wsp->matrix[k] == NULL) {
		if (sim->EDsymmetry == 1) { // ED symmetry present
			if (blk_cm_isdiag(Ud)) {
				for (i=k; i<l; i++) {
					DEBUGPRINT("new_gcompute: Ud is diagonal\n");
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
				}
			} else {
				DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
				T = blk_cm_diag(Ud);
				for (i=k; i<l; i++) {
					// this should be faster:
					blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
					free_blk_mat_complx(dumblk);
				}
				blk_simtrans_adj(det,T,sim);
				free_blk_mat_complx(T);
			}
		} else {  // no ED symmetry but full job needs to be done
			wsp->matrix[k] = cm_change_basis_2(wsp->fstart,Ud->basis,sim);
			if (blk_cm_isdiag(Ud)) {
				for (i=k; i<l; i++) {
					DEBUGPRINT("new_gcompute: Ud is diagonal\n");
					wsp->matrix[i+1] = cm_dup(wsp->matrix[k]);
					blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
				}
			} else {
				DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
				// need to diagonalize final propagator
				T = blk_cm_diag(Ud);
				for (i=k; i<l; i++) {
					// this should be faster:
					blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
					wsp->matrix[i+1] = cm_dup(wsp->matrix[k]);
					blk_simtrans_adj(wsp->matrix[i+1],dumblk,sim);
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
					free_blk_mat_complx(dumblk);
				}
				blk_simtrans_adj(det,T,sim);
				blk_simtrans_adj(wsp->matrix[k],T,sim);
				free_blk_mat_complx(T);
			}
			// control print out
			//for (i=0; i<Ng; i++) {
			//	cm_print(wsp->matrix[k+i],"RHO");
			//	cm_print(wsp->matrix[k+Ng+i],"DET");
			//}
		}
	} else { // full job with gamma-dependent density matrices
		if (blk_cm_isdiag(Ud)) {
			for (i=k; i<l; i++) {
				DEBUGPRINT("new_gcompute: Ud is diagonal\n");
				blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
				wsp->matrix[i+Ng+1] = cm_dup(det);
				blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
				free_blk_mat_complx(wsp->STO[i]);
				wsp->STO[i] = NULL;
			}
		} else {
			DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
			// need to diagonalize final propagator
			T = blk_cm_diag(Ud);
			for (i=k; i<l; i++) {
				// this should be faster:
				blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
				free_blk_mat_complx(wsp->STO[i]);
				wsp->STO[i] = NULL;
				blk_simtrans_adj(wsp->matrix[i+1],dumblk,sim);
				wsp->matrix[i+Ng+1] = cm_dup(det);
				blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
				free_blk_mat_complx(dumblk);
			}
			blk_simtrans_adj(det,T,sim);
			blk_simtrans_adj(wsp->matrix[k],T,sim);
			free_blk_mat_complx(T);
		}
		// control print out
		//for (i=0; i<Ng; i++) {
		//	cm_print(wsp->matrix[k+i],"RHO");
		//	cm_print(wsp->matrix[k+Ng+i],"DET");
		//}
	}
	// END of transformations
	//QueryPerformanceCounter(&tv2);
	//printf("timing transformations: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
	TIMING_TOC(tv1,tv2,"transformations");

	freq = dptr = (double*)malloc(matdim*sizeof(double));
	for (i=0; i<Ud->Nblocks; i++) {
		z1 = Ud->m[i].data;
		assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
		for (m=0; m<Ud->blk_dims[i]; m++) {
			*dptr = -Carg((*z1))/sim->taur*1.0e6;
			z1++;
			dptr++;
		}
	}
	// control print-out
	//for (i=0; i<matdim; i++) printf("freq[%d] = %g\n",i,freq[i]);

	// fork from here to ED symmetry specific calculation
	if (sim->EDsymmetry == 1) {
		TIMING_TIC(tv1);
		double binsize = sim->sw*2*M_PI/sim->np;
	    fftw_complex *fftin1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
	    fftw_complex *fftout1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
	    fftw_plan p1 = fftw_plan_dft_1d(Ng, fftin1, fftout1, FFTW_FORWARD, FFTW_ESTIMATE);
	    irow = (int*)malloc((matdim+1)*sizeof(int));
	    icol = ic = (int*)malloc(matdim*matdim*sizeof(int)); // just large enough
	    irow[0]=1;
	    for (r=1; r<=matdim; r++) {
	    	irow[r] = irow[r-1];
	    	for (c=1; c<=matdim; c++) {
	    		m = 0;
	    		for (i=0; i<Ng; i++) {
	    			complx zz = cm_getelem(wsp->matrix[k+Ng+i],r,c);
	    			if (fabs(zz.re) > TINY || fabs(zz.im) > TINY) {
	    				m = 1;
	    				break;
	    			}
	    		}
	    		if (m != 0) {
	    			irow[r]++;
	    			*ic = c;
	    			ic++;
	    		}
	    	}
	    }
	    TIMING_TOC(tv1,tv2,"scan contrib");
	    TIMING_TIC(tv1);
	    int nnz = irow[matdim] - 1;
	    complx *qdata = (complx *)malloc(nnz*Ng*sizeof(complx));
	    ic = icol;
	    m = 0;
	    for (r=1; r<=matdim; r++) {
	    	int nc = irow[r] - irow[r-1];
	    	for (c=0; c<nc; c++) {
				diff = freq[r-1] - freq[*ic-1];
				complx ph = Cexpi(-diff*dtg*1.0e-6);
				complx phmul = Complx(1.0,0.0);
	    		for (i=0; i<Ng; i++) {
	    			complx zz1 = cm_getelem(wsp->matrix[k+Ng+i],r,*ic);  //q_rs(k)
	    			fftin1[i][0] = zz1.re*phmul.re - zz1.im*phmul.im;
	    			fftin1[i][1] = zz1.re*phmul.im + zz1.im*phmul.re;
	    			phmul = Cmul(phmul,ph);
	    		}
	    		fftw_execute(p1);
	    		//printf("\nr = %d, c = %d : ",r,*ic);
	    		for (i=0; i<Ng; i++) {
	    			qdata[i*nnz+m].re = fftout1[i][0];
	    			qdata[i*nnz+m].im = fftout1[i][1];
	    			//printf("( %g, %g ) ",fftout1[i][0],fftout1[i][1]);
	    		}
	    		//printf("\n");
	    		m++;
	    		ic++;
	    	}
	    }
	    mat_complx dum;
	    dum.irow = irow;
	    dum.icol = icol;
	    dum.type = MAT_SPARSE;
	    dum.row = dum.col = matdim;
	    ic = icol;
	    m = 0;
	    for (r=1; r<=matdim; r++) {
	    	int nc = irow[r] - irow[r-1];
	    	for (c=0; c<nc; c++) {
				diff = freq[r-1] - freq[*ic-1];
				for (i=0; i<Ng; i++) {
					int idx = (i + Ng/2 + 1) % Ng;
					complx zz1 = qdata[idx*nnz+m];
					dum.data = &qdata[( ( ((Ng-i)%Ng) + Ng/2+1 )%Ng )*nnz];
					complx zz2 = cm_getelem(&dum,*ic,r); // q_sr(-k)
					double zzre = (zz1.re*zz1.re+zz1.im*zz1.im + zz1.re*zz2.re-zz1.im*zz2.im)*0.5;
					double zzim = (zz1.re*zz2.im+zz1.im*zz2.re)*0.5;
					bin = (int)(1.5-(diff + sim->wr*(i-Ng/2+1))/binsize + sim->np/2);
					//printf("i = %d, idx = %d, -i = %d, idx = %d : freq = %g, bin %d -> ",i,idx,(Ng-i)%Ng,( ((Ng-i)%Ng) + Ng/2+1 )%Ng,diff+sim->wr*(i-Ng/2+1),bin);
					while (bin < 1) bin += sim->np;
					while (bin > sim->np) bin -= sim->np;
					assert(bin >= 1 && bin <= sim->np);
					wsp->fid[bin].re += zzre*cosph - zzim*sinph;
					wsp->fid[bin].im += zzim*cosph + zzre*sinph;
					//printf(" %d : ( %g, %g ) ( %g, %g ) -> ( %g, %g )\n",bin,zz1.re,zz1.im,zz2.re,zz2.im,zzre,zzim);
				}
				ic++;
				m++;
	    	}
	    }
		free(qdata);
		free(irow);
		free(icol);
		free(freq);
		free_blk_mat_complx(wsp->STO[l]);
		wsp->STO[l] = NULL;
		for (i=k; i<=l; i++) {
			free_complx_matrix(wsp->matrix[i+Ng]);
			wsp->matrix[i+Ng] = NULL;
		}
		//QueryPerformanceCounter(&tv2);
		TIMING_TOC(tv1,tv2,"acq spectrum generation");
		return;
	}
	// END of ED symmetry calculation

	TIMING_TIC(tv1);
	// construct the spectrum
	double binsize = sim->sw*2*M_PI/sim->np;
	//printf("binsize = %g\n",binsize);
    fftw_complex *fftin1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
    fftw_complex *fftout1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
    fftw_plan p1 = fftw_plan_dft_1d(Ng, fftin1, fftout1, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_complex *fftin2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
    fftw_complex *fftout2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
    fftw_plan p2 = fftw_plan_dft_1d(Ng, fftin2, fftout2, FFTW_FORWARD, FFTW_ESTIMATE);

    irow = (int*)malloc((matdim+1)*sizeof(int));
    icol = ic = (int*)malloc(matdim*matdim*sizeof(int)); // just large enough
    scan_contrib(&wsp->matrix[k],Ng,irow,icol);
	TIMING_TOC(tv1,tv2,"scan contrib");

    if (sim->interpolation == 1) { // FWT interpolation
    	// now I am sure the file exists
    	if (ftell(wsp->interpol_file) == 0) {
    		// write header
    		r = irow[matdim]-1;
    		fwrite(&r,sizeof(int),1,wsp->interpol_file);
    		fwrite(irow,sizeof(int),matdim+1,wsp->interpol_file);
    		fwrite(icol,sizeof(int),r,wsp->interpol_file);
    		printf("Number of NNZ contributions: %d, matdim = %d\n",r,matdim);
    	}
    	// write crystallite index
    	fwrite(&wsp->cryst_idx,sizeof(int),1,wsp->interpol_file);
    	// write diagonalized period propagator
    	for (i=0; i<Ud->Nblocks; i++) {
    		fwrite(Ud->m[i].data,sizeof(complx),Ud->blk_dims[i],wsp->interpol_file);
    	}
    } else if (sim->interpolation == 2) { // ASG interpolation
    	if (ftell(wsp->interpol_file) == 0) {
    		// write header
    		r = irow[matdim]-1;
    		fwrite(&r,sizeof(int),1,wsp->interpol_file);
    	}
    	// write crystallite index
    	fwrite(&wsp->cryst_idx,sizeof(int),1,wsp->interpol_file);
    }
    TIMING_TIC(tv1);
    int kkk = -1; // DEBUGING
    for (r=1; r<=matdim; r++) {
    	int nc = irow[r] - irow[r-1];
    	for (c=0; c<nc; c++) {
    		if (sim->interpolation == 1) { // FWT interpolation
    			kkk++; // DEBUGING
        		for (i=0; i<Ng; i++) {
        			complx zz1 = cm_getelem(wsp->matrix[k+i],*ic,r);
        			complx zz2 = cm_getelem(wsp->matrix[k+Ng+i],r,*ic);
        			fwrite(&zz1,sizeof(complx),1,wsp->interpol_file);
        			fwrite(&zz2,sizeof(complx),1,wsp->interpol_file);
        			if (kkk == 63 && i==Ng/2) printf("\t %g %g\n",zz1.re,zz1.im); // DEBUGING
        		}
        		ic++;
        		continue;
    		}
			diff = freq[r-1] - freq[*ic-1];
			if (sim->interpolation == 2) fwrite(&diff,sizeof(double),1,wsp->interpol_file);
			complx ph = Cexpi(-diff*dtg*1.0e-6);
			complx phmul = Complx(1.0,0.0);
    		for (i=0; i<Ng; i++) {
    			complx zz1 = cm_getelem(wsp->matrix[k+i],*ic,r);
    			complx zz2 = cm_getelem(wsp->matrix[k+Ng+i],r,*ic);
    			fftin1[i][0] = zz1.re*phmul.re + zz1.im*phmul.im;
    			fftin1[i][1] = zz1.re*phmul.im - zz1.im*phmul.re;
    			fftin2[i][0] = zz2.re*phmul.re - zz2.im*phmul.im;
    			fftin2[i][1] = zz2.re*phmul.im + zz2.im*phmul.re;
    			phmul = Cmul(phmul,ph);
    		}
    		fftw_execute(p1);
    		fftw_execute(p2);
			for (i=0; i<Ng; i++) {
				int idx = (i + Ng/2 + 1) % Ng;
				double zzre = fftout1[idx][0]*fftout2[idx][0] + fftout1[idx][1]*fftout2[idx][1];
				double zzim = fftout1[idx][0]*fftout2[idx][1] - fftout1[idx][1]*fftout2[idx][0];
				bin = (int)(1.5-(diff + sim->wr*(i-Ng/2+1))/binsize + sim->np/2);
				//printf("i = %d, idx = %d : freq = %g, bin %d -> ",i,idx,diff+sim->wr*(i-Ng/2+1),bin);
				while (bin < 1) bin += sim->np;
				while (bin > sim->np) bin -= sim->np;
				assert(bin >= 1 && bin <= sim->np);
				if (sim->interpolation == 2) {
					complx a;
					a.re = zzre*cosph - zzim*sinph;
					a.im = zzim*cosph + zzre*sinph;
					fwrite(&a,sizeof(complx),1,wsp->interpol_file);
				} else {
					wsp->fid[bin].re += zzre*cosph - zzim*sinph;
					wsp->fid[bin].im += zzim*cosph + zzre*sinph;
				}
				//printf("%d : (%g, %g) (%g, %g) (%g, %g)\n",bin,fftout1[idx][0],fftout1[idx][1],fftout2[idx][0],fftout2[idx][1],zzre,zzim);
			}
			ic++;
    	}
    }
    TIMING_TOC(tv1,tv2,"acq spectrum");

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(fftin1); fftw_free(fftout1);
    fftw_free(fftin2); fftw_free(fftout2);
    free(irow);
    free(icol);
    free(freq);

	//free_complx_matrix(mexp1);
	free_blk_mat_complx(wsp->STO[l]);
	wsp->STO[l] = NULL;
	for (i=k; i<=l; i++) {
		free_complx_matrix(wsp->matrix[i]);
		wsp->matrix[i] = NULL;
		free_complx_matrix(wsp->matrix[i+Ng]);
		wsp->matrix[i+Ng] = NULL;
	}

}







/****
 *  this is in time domain when interpolation done on amplitudes and frequencies
 * PROBLEM: amplitudes contain folded frequencies before interpolation
 * of frequencies, and so a mismatch of +/- sim->wr can exist, spoiling the process
 ****//*
void collect_fid_interpol(acq_data *acqdata)
{
	int i, j, q, r, s, idx, *ic;
	int fnnz = acqdata->irow[acqdata->dim] - 1;
	complx *mexp = (complx *)malloc(2*acqdata->dim*sizeof(complx));
	complx z, *zf, *zfid, *z1, *z2, *z3;
	double zzre, zzim;
	const double m_Ng = (double)(acqdata->m)/(double)(acqdata->Ng);
	double weight = acqdata->weight/(double)(acqdata->Ng);
	//printf("\tweight = %g\n",weight);

	//printf("\n\nAfter interpolation\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	//printf("lams:\n");
	//for (i=0; i<acqdata->dim; i++) {
	//	printf("(%g, %g) ",acqdata->lams[i].re, acqdata->lams[i].im);
	//}
	//printf("\nFrs\n");
	//mat_complx ccc;
	//ccc.basis = 0; ccc.col = ccc.row = acqdata->dim;
	//ccc.type = MAT_SPARSE; ccc.irow = acqdata->irow; ccc.icol = acqdata->icol;
	//for (i=0; i<acqdata->Ng; i++) {
	//	ccc.data = acqdata->frs + i*fnnz;
	//	cm_print(&ccc," ");
	//}

	z1 = acqdata->lams;
	z2 = mexp;
	z3 = mexp + acqdata->dim;
	for (i=0; i<acqdata->dim; i++) {
		z2->re = 1.0; z2->im = 0.0;
		*z3 = CRpow(*z1,m_Ng);
		z1++; z2++; z3++;
	}

	zfid = acqdata->fid + 1;
	for (i=0; i<acqdata->Nacq; i++) {
		Czero(z);
		q = (i*acqdata->m)%acqdata->Ng;
		//printf("\ti = %d, q = %d\n",i,q);
		idx = 0;
		z1 = mexp;
		zf = acqdata->frs + q*fnnz;
		ic = acqdata->icol;
		for (r=0; r<acqdata->dim; r++) {
			int ns = acqdata->irow[r+1] - acqdata->irow[r];
			for (j=0; j<ns; j++ ) {
				// s = *ic - 1;
				z3 = mexp + *ic - 1;
				zzre = z1->re*z3->re + z1->im*z3->im;
				zzim = z1->re*z3->im - z1->im*z3->re;
				z.re += zf->re*zzre - zf->im*zzim;
				z.im += zf->re*zzim + zf->im*zzre;
				//printf("r,s = %d, %d; exp = (%g, %g), F = (%g, %g)\n",r+1,*ic,zzre,zzim,zf->re,zf->im);
				ic++;
				zf++;
			}
			z1++;
		}
		// update mexp
		z1 = mexp;
		z2 = mexp + acqdata->dim;
		for (r=0; r<acqdata->dim; r++) {
			zzre = z1->re;
			z1->re = zzre*z2->re - z1->im*z2->im;
			z1->im = zzre*z2->im + z1->im*z2->re;
			z1++;
			z2++;
		}
		//printf("z = (%g, %g)\n",z.re,z.im);
		// add to fid
		if (acqdata->isphase != 0) {
			zfid->re += (z.re*acqdata->cosph - z.im*acqdata->sinph)*weight;
			zfid->im += (z.im*acqdata->cosph + z.re*acqdata->sinph)*weight;
		} else {
			zfid->re += z.re*weight;
			zfid->im += z.im*weight;
		}
		zfid++;
	}
	free(mexp);

	//printf("fid[1] = (%g, %g)\n",acqdata->fid[1].re,acqdata->fid[1].im);
}
*****/

/*****
 * this is when interpolation is done separately on rho and det
 * from which the amplitudes are created with already interpolated frequencies.
 * This method preserves fully compensation of possible frequency folds
 ****/
void collect_fid_interpol(acq_data *acqdata)
{
	int i, j, r, c, *ic;
	int fnnz = acqdata->irow[acqdata->dim] - 1;
	complx zz1, zz2;
	double diff, *freq, *dptr;
	const double dtg = (double)(acqdata->taur)/(double)(acqdata->Ng);
	double weight = acqdata->weight/(double)(acqdata->Ng);
	//printf("\tweight = %g\n",weight);

	freq = dptr = (double*)malloc(acqdata->dim*sizeof(double));
	for (i=0; i<acqdata->dim; i++) {
		*dptr = -Carg(acqdata->lams[i])/acqdata->taur*1.0e6;
		dptr++;
	}
	// control print-out
	//for (i=0; i<acqdata->dim; i++) printf("freq[%d] = %g\n",i,freq[i]);


	//QueryPerformanceCounter(&tv3);
    complx *phmul = (complx*)malloc(acqdata->Ng*sizeof(complx));
    complx *Frs = complx_vector(acqdata->Ng);
    ic = acqdata->icol;
    int dataidx = 0;
    for (r=1; r<=acqdata->dim; r++) {
    	int nc = acqdata->irow[r] - acqdata->irow[r-1];
    	for (c=0; c<nc; c++) {
    		//QueryPerformanceCounter(&tv1);
			diff = freq[r-1] - freq[*ic-1];
			complx ph = Cexpi(-diff*dtg*1.0e-6);
		    phmul[0] = Cunit;
			for (i=1; i<acqdata->Ng; i++) {
				phmul[i] = Cmul(phmul[i-1],ph);
			}
			ph = Cexpi(diff*acqdata->taur*1.0e-6);
			cv_zero(Frs);
 			int ppp = acqdata->Ng;
			for (j=0; j<acqdata->Ng; j++) {
				for (i=0; i<acqdata->Ng; i++) {
					zz1 = acqdata->frs[dataidx + i*fnnz];
					zz2 = acqdata->frs[dataidx + acqdata->Ng*fnnz + ((i+j)%acqdata->Ng)*fnnz];
					//printf("\t rho = (%g, %g), det = (%g, %g), ph = (%g, %g)\n",zz1.re,zz1.im, zz2.re, zz2.im, phmul[j].re, phmul[j].im);
					zz1 = Cmul(zz1,zz2);
					zz2 = Cmul(zz1,phmul[j]);
					if (i>=ppp) zz2 = Cmul(zz2,ph);
					Frs[j+1].re += zz2.re;
					Frs[j+1].im += zz2.im;
				}
				//phmul[j] = Cmul(phmul[j],ph);
				ppp--;
			} // Frs[] vector ready to use
			//printf("elem r = %d, s = %d\n",r, *ic);
			//cv_print(Frs,"Frs");
    		//QueryPerformanceCounter(&tv2);
    		//printf("timing g-ave contrib: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);
    		//QueryPerformanceCounter(&tv1);
    		ph = Cunit;
    		zz2 = Cexpi(diff*dtg*acqdata->m*1.0e-6);
    		for (i=0; i<acqdata->Nacq; i++) {
    			zz1 = Cmul(Frs[(i*acqdata->m)%acqdata->Ng + 1],ph);
    			acqdata->fid[i+1].re += (zz1.re*acqdata->cosph - zz1.im*acqdata->sinph)*weight;
    			acqdata->fid[i+1].im += (zz1.im*acqdata->cosph + zz1.re*acqdata->sinph)*weight;
    			ph = Cmul(ph,zz2);
    		}
    		//QueryPerformanceCounter(&tv2);
    		//printf("timing fid contrib: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);

			ic++;
			dataidx++;
    	}
    }
	//QueryPerformanceCounter(&tv4);
	//printf("timing fid contrib: %.9f\n",((float)(tv4.QuadPart-tv3.QuadPart))/(float)tickpsec.QuadPart);

	free(phmul);
	free_complx_vector(Frs);
	free(freq);
}


/*****
 * this is in frequency domain
 * when interpolating amplitudes and frequencies.
 * PROBLEM: amplitudes contain folded frequencies before interpolation
 * of frequencies and so a mismatch of +/- sim->wr can exist, spoiling the process
 *//*
void collect_spc_interpol_fungujici(acq_data *acqdata)
{
	int i, r, c, bin, *ic;
	double *freq, *dptr, diff, binsize;
	complx zz, *aptr;
	int annz = acqdata->irow[acqdata->dim] - 1;
	double weight = acqdata->weight/(double)(acqdata->Ng);

	freq = dptr = (double*)malloc(acqdata->dim*sizeof(double));
	for (i=0; i<acqdata->dim; i++) {
		*dptr = -Carg(acqdata->lams[i])/acqdata->taur*1.0e6;
		dptr++;
	}
	// control print-out
	for (i=0; i<acqdata->dim; i++) printf("freq[%d] = %g\n",i,freq[i]);

	binsize = acqdata->sw*2*M_PI/acqdata->Np;
	ic = acqdata->icol;
	aptr = acqdata->frs;
    for (r=1; r<=acqdata->dim; r++) {
    	int nc = acqdata->irow[r] - acqdata->irow[r-1];
    	for (c=0; c<nc; c++) {
			diff = freq[r-1] - freq[*ic-1];
			for (i=0; i<acqdata->Ng; i++) {
				zz = aptr[i*annz];
				bin = (int)(1.5-(diff + acqdata->wr*(i-acqdata->Ng/2+1))/binsize + acqdata->Np/2);
				while (bin < 1) bin += acqdata->Np;
				while (bin > acqdata->Np) bin -= acqdata->Np;
				assert(bin >= 1 && bin <= acqdata->Np);
				acqdata->fid[bin].re += (zz.re*acqdata->cosph - zz.im*acqdata->sinph)*weight;
				acqdata->fid[bin].im += (zz.im*acqdata->cosph + zz.re*acqdata->sinph)*weight;
				printf("\telem [%d,%d] freq = %g, bin %d, (%g, %g)\n",r,*ic,diff + acqdata->wr*(i-acqdata->Ng/2+1),bin,zz.re, zz.im);
			}
			ic++;
			aptr++;
    	}
    }
    free(freq);

}
*******/

/*****
 * this is when interpolation is done separately on rho and det
 * from which the amplitudes are created with already interpolated frequencies.
 * This method preserves fully compensation of possible frequency folds
 ****/
void collect_spc_interpol(acq_data *acqdata)
{
	int i, r, c, bin, *ic;
	double *freq, *dptr, diff, binsize;
	complx *aptr;
	int annz = acqdata->irow[acqdata->dim] - 1;
	double weight = acqdata->weight/(double)(acqdata->Ng);

	freq = dptr = (double*)malloc(acqdata->dim*sizeof(double));
	for (i=0; i<acqdata->dim; i++) {
		*dptr = -Carg(acqdata->lams[i])/acqdata->taur*1.0e6;
		dptr++;
	}
	// control print-out
	//for (i=0; i<acqdata->dim; i++) printf("freq[%d] = %g\n",i,freq[i]);

    fftw_complex *fftin1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * acqdata->Ng);
    fftw_complex *fftout1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * acqdata->Ng);
    fftw_plan p1 = fftw_plan_dft_1d(acqdata->Ng, fftin1, fftout1, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_complex *fftin2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * acqdata->Ng);
    fftw_complex *fftout2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * acqdata->Ng);
    fftw_plan p2 = fftw_plan_dft_1d(acqdata->Ng, fftin2, fftout2, FFTW_FORWARD, FFTW_ESTIMATE);

	binsize = acqdata->sw*2*M_PI/acqdata->Np;
	ic = acqdata->icol;
	aptr = acqdata->frs;
    for (r=1; r<=acqdata->dim; r++) {
    	int nc = acqdata->irow[r] - acqdata->irow[r-1];
    	for (c=0; c<nc; c++) {
			diff = freq[r-1] - freq[*ic-1];
			complx ph = Cexpi(-diff*acqdata->taur*1.0e-6/acqdata->Ng);
			complx phmul = Complx(1.0,0.0);
    		for (i=0; i<acqdata->Ng; i++) {
    			complx zz1 = aptr[i*annz];
    			complx zz2 = aptr[(i+acqdata->Ng)*annz];
    			fftin1[i][0] = zz1.re*phmul.re + zz1.im*phmul.im;
    			fftin1[i][1] = zz1.re*phmul.im - zz1.im*phmul.re;
    			fftin2[i][0] = zz2.re*phmul.re - zz2.im*phmul.im;
    			fftin2[i][1] = zz2.re*phmul.im + zz2.im*phmul.re;
    			phmul = Cmul(phmul,ph);
    		}
    		fftw_execute(p1);
    		fftw_execute(p2);
			for (i=0; i<acqdata->Ng; i++) {
				int idx = (i + acqdata->Ng/2 + 1) % acqdata->Ng;
				double zzre = fftout1[idx][0]*fftout2[idx][0] + fftout1[idx][1]*fftout2[idx][1];
				double zzim = fftout1[idx][0]*fftout2[idx][1] - fftout1[idx][1]*fftout2[idx][0];
				bin = (int)(1.5-(diff + acqdata->wr*(i-acqdata->Ng/2+1))/binsize + acqdata->Np/2);
				while (bin < 1) bin += acqdata->Np;
				while (bin > acqdata->Np) bin -= acqdata->Np;
				assert(bin >= 1 && bin <= acqdata->Np);
				acqdata->fid[bin].re += (zzre*acqdata->cosph - zzim*acqdata->sinph)*weight;
				acqdata->fid[bin].im += (zzim*acqdata->cosph + zzre*acqdata->sinph)*weight;
				//printf("\telem [%d,%d] freq = %g, bin %d, (%g, %g)\n",r,*ic,diff + acqdata->wr*(i-acqdata->Ng/2+1),bin,zzre, zzim);
			}
			ic++;
			aptr++;
    	}
    }

    free(freq);
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(fftin1); fftw_free(fftout1);
    fftw_free(fftin2); fftw_free(fftout2);

}

void collect_fid_interpol_direct(acq_data *acqdata)
{
	int i, j, r, c, *ic, ig, dataidx;
	int fnnz = acqdata->irow[acqdata->dim] - 1;
	complx *Ud, Frs, ph, phfac;
	//double diff, *freq, *dptr;
	//const double dtg = (double)(acqdata->taur)/(double)(acqdata->ncp);
	double weight = acqdata->weight/(double)(acqdata->Ng);
	//printf("\tweight = %g\n",weight);

	for (ig=0; ig<acqdata->Ng; ig++) {
		Ud = &acqdata->lams[ig*acqdata->dim];
		ic = acqdata->icol;
		dataidx = 0;
		for (r=1; r<=acqdata->dim; r++) {
			int nc = acqdata->irow[r] - acqdata->irow[r-1];
			for (c=0; c<nc; c++) {
				ph.re = Ud[r-1].re*Ud[*ic-1].re + Ud[r-1].im*Ud[*ic-1].im;
				ph.im = Ud[r-1].re*Ud[*ic-1].im - Ud[r-1].im*Ud[*ic-1].re;
				j = 0;
				phfac = Cunit;
				for (i=0; i<acqdata->Nacq; i++) {
					Frs = acqdata->frs[((i%acqdata->ncp)+ig*acqdata->ncp)*fnnz+dataidx];
					Frs = Cmul(Frs,phfac);
					acqdata->fid[i+1].re += (Frs.re*acqdata->cosph - Frs.im*acqdata->sinph)*weight;
					acqdata->fid[i+1].im += (Frs.im*acqdata->cosph + Frs.re*acqdata->sinph)*weight;
					j++;
					if (j == acqdata->ncp) {
						j = 0;
						phfac = Cmul(phfac,ph);
					}
				}
				ic++;
				dataidx++;
			}
		}
	}
}

void collect_spc_interpol_direct(acq_data *acqdata)
{
	int i, r, c, bin, *ic, ig, dataidx;
	double *freq, *dptr, diff, binsize;
	int annz = acqdata->irow[acqdata->dim] - 1;
	double weight = acqdata->weight/(double)(acqdata->Ng);

	freq = (double*)malloc(acqdata->dim*sizeof(double));
    fftw_complex *fftin1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * acqdata->ncp);
    fftw_complex *fftout1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * acqdata->ncp);
    fftw_plan p1 = fftw_plan_dft_1d(acqdata->ncp, fftin1, fftout1, FFTW_FORWARD, FFTW_ESTIMATE);
	binsize = acqdata->sw*2*M_PI/acqdata->Np;

	for (ig=0; ig<acqdata->Ng;ig++) {
		dptr = freq;
		for (i=0; i<acqdata->dim; i++) {
			*dptr = -Carg(acqdata->lams[ig*acqdata->dim + i])/acqdata->taur*1.0e6;
			dptr++;
		}
		// control print-out
		//for (i=0; i<acqdata->dim; i++) printf("freq[%d] = %g\n",i,freq[i]);

		ic = acqdata->icol;
		dataidx = 0;
    	for (r=1; r<=acqdata->dim; r++) {
    		int nc = acqdata->irow[r] - acqdata->irow[r-1];
    		for (c=0; c<nc; c++) {
    			diff = freq[r-1] - freq[*ic-1];
    			complx ph = Cexpi(-diff*acqdata->taur*1.0e-6/acqdata->ncp);
    			complx phmul = Complx(1.0,0.0);
    			for (i=0; i<acqdata->ncp; i++) {
    				complx zz1 = acqdata->frs[(i+ig*acqdata->ncp)*annz+dataidx];
    				fftin1[i][0] = zz1.re*phmul.re + zz1.im*phmul.im;
    				fftin1[i][1] = zz1.re*phmul.im - zz1.im*phmul.re;
    				phmul = Cmul(phmul,ph);
    			}
    			fftw_execute(p1);
    			for (i=0; i<acqdata->ncp; i++) {
    				int idx = (i + acqdata->ncp/2 + 1) % acqdata->ncp;
    				double zzre = fftout1[idx][0];
    				double zzim = fftout1[idx][1];
    				bin = (int)(1.5-(diff + (2.0e6*M_PI/acqdata->taur)*(i-acqdata->ncp/2+1))/binsize + acqdata->Np/2);
    				while (bin < 1) bin += acqdata->Np;
    				while (bin > acqdata->Np) bin -= acqdata->Np;
    				assert(bin >= 1 && bin <= acqdata->Np);
    				acqdata->fid[bin].re += (zzre*acqdata->cosph - zzim*acqdata->sinph)*weight;
    				acqdata->fid[bin].im += (zzim*acqdata->cosph + zzre*acqdata->sinph)*weight;
    			}
    			ic++;
    			dataidx++;
    		}
    	}
	}

	free(freq);
    fftw_destroy_plan(p1);
    fftw_free(fftin1); fftw_free(fftout1);
}

/*****
 * new version of gcompute with different
 *     calculation of sideband amplitudes
 *****/
void new_gcompute_time(Sim_info *sim, Sim_wsp *wsp)
{
	int ig, k, l, m, Ng, matdim, r, c, i, j;
	double dtg, t0, diff;
	mat_complx *det;
	blk_mat_complx *Ud, *T;
	complx *z1, zz1, zz2, acqph;
	double *freq, *dptr;
	int *irow, *icol, *ic;

	//LARGE_INTEGER tv1, tv2, tv3, tv4, tickpsec;
	//QueryPerformanceFrequency(&tickpsec);
	TIMING_INIT;
	TIMING_INIT_VAR(tv1);
	TIMING_INIT_VAR(tv2);

	Ng = sim->ngamma;
	matdim = sim->matdim;
	dtg = sim->taur/(double)Ng;
	m = (int)floor(wsp->dw/dtg+1e-6);
	if ( fabs(dtg*m - wsp->dw) > TINY) {
		fprintf(stderr,"Error: (i)gcompute - bad synchronization of acquisition and gamma-averaging\n");
		fprintf(stderr,"       dwelltime(%g) can not be split in tauR(%g)/Ngamma(%d) steps\n",wsp->dw,sim->taur,sim->ngamma);
		exit(1);
	}
	//printf("new_gcompute_time: dw = %g; dtg = %g; m = %d\n",wsp->dw,dtg,m);

	/* initialize counter for storing wsp->sigma and wsp->U */
	wsp->acqblock_sto = k = ACQBLOCK_STO_INI;

	// evaluate pulse sequence to generate gamma-dependent sigma and props
	TIMING_TIC(tv1);
	for (ig=0;ig<Ng;ig++) {
		DEBUGPRINT("new_gcompute_freq -> ig = %d\n",ig);
		wsp->tstart = t0 = ig/(double)Ng*sim->taur;
		direct_propagate(sim,wsp);
		if (wsp->acqblock_sto == ACQBLOCK_STO_INI + ig) {
			// there was no acq_block command in pulseq, need to do the job here
			if (fabs(wsp->t - t0 - dtg) > TINY) {
				fprintf(stderr,"Error: gcompute - pulse sequence event not synchronized with gamma-averaging\n");
				fprintf(stderr,"                  pulseq duration = %g us, need to be %g\n",wsp->t-t0,dtg);
				exit(1);
			}
			if (ig == 0) {
				wsp->STO[wsp->acqblock_sto] = blk_cm_dup(wsp->U);
			} else {
				wsp->STO[wsp->acqblock_sto] = blk_cm_mul(wsp->U,wsp->STO[wsp->acqblock_sto - 1]);
			}
			// this input file regime does not allow initial density matrix to be changed
			//wsp->matrix[wsp->acqblock_sto] = cm_change_basis_2(wsp->sigma,wsp->U->basis,sim); // initial density matrix is gamma-independent
			acqblock_sto_incr(wsp);
		}
		wsp->cryst.gamma += 360.0/(double)sim->ngamma;
	}
	TIMING_TOC(tv1,tv2,"eval_pulseq");

	l = wsp->acqblock_sto - 1;
	DEBUGPRINT("new_gcompute_freq: sigmas and propagators from %d to %d\n",k,l);
	if (wsp->acqblock_sto - ACQBLOCK_STO_INI != Ng) {
		fprintf(stderr,"Error: (i)gcompute - propagator count mismatch\n");
		exit(1);
	}
	if (l+Ng >= ACQBLOCK_STO_END) {
		fprintf(stderr,"Error: (i)gcompute - overflow in matrix buffers\n");
		exit(1);
	}
	if ( fabs(wsp->acqphase) > TINY ) {
		acqph = Cexpi(wsp->acqphase*DEG2RAD);
	} else {
		acqph = Cunit;
	}

	// transformations of dens.matrix and detect op.
	Ud = wsp->STO[l];
	//blk_cm_print(Ud,"new_gcompute Ud");
	if (sim->acq_adjoint) {
		if (wsp->fdetect != sim->fdetect)
			cm_adjointi(wsp->fdetect);
		else
			wsp->fdetect = cm_adjoint(sim->fdetect);
	}
	det = wsp->matrix[l+1] = cm_change_basis_2(wsp->fdetect,Ud->basis,sim);
	if (wsp->matrix[k] == NULL) {
		//if (sim->EDsymmetry == 1) { // ED symmetry present
		//	if (blk_cm_isdiag(Ud)) {
		//		for (i=k; i<l; i++) {
		//			DEBUGPRINT("new_gcompute: Ud is diagonal\n");
		//			wsp->matrix[i+Ng+1] = cm_dup(det);
		//			blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
		//			free_blk_mat_complx(wsp->STO[i]);
		//			wsp->STO[i] = NULL;
		//		}
		//	} else {
		//		DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
		//		T = blk_cm_diag(Ud);
		//		for (i=k; i<l; i++) {
		//			// this should be faster:
		//			blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
		//			free_blk_mat_complx(wsp->STO[i]);
		//			wsp->STO[i] = NULL;
		//			wsp->matrix[i+Ng+1] = cm_dup(det);
		//			blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
		//			free_blk_mat_complx(dumblk);
		//		}
		//		blk_simtrans_adj(det,T,sim);
		//		free_blk_mat_complx(T);
		//	}
		//} else {  // no ED symmetry but full job needs to be done
			wsp->matrix[k] = cm_change_basis_2(wsp->fstart,Ud->basis,sim);
			if (blk_cm_isdiag(Ud)) {
				for (i=k; i<l; i++) {
					DEBUGPRINT("new_gcompute: Ud is diagonal\n");
					wsp->matrix[i+1] = cm_dup(wsp->matrix[k]);
					blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
				}
			} else {
				DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
				// need to diagonalize final propagator
				T = blk_cm_diag(Ud);
				for (i=k; i<l; i++) {
					// this should be faster:
					blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
					wsp->matrix[i+1] = cm_dup(wsp->matrix[k]);
					blk_simtrans_adj(wsp->matrix[i+1],dumblk,sim);
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
					free_blk_mat_complx(dumblk);
				}
				blk_simtrans_adj(det,T,sim);
				blk_simtrans_adj(wsp->matrix[k],T,sim);
				free_blk_mat_complx(T);
			}
			// control print out
			//for (i=0; i<Ng; i++) {
			//	cm_print(wsp->matrix[k+i],"RHO");
			//	cm_print(wsp->matrix[k+Ng+i],"DET");
			//}
		//}
	} else { // full job with gamma-dependent density matrices
		if (blk_cm_isdiag(Ud)) {
			for (i=k; i<l; i++) {
				DEBUGPRINT("new_gcompute: Ud is diagonal\n");
				blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
				wsp->matrix[i+Ng+1] = cm_dup(det);
				blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
				free_blk_mat_complx(wsp->STO[i]);
				wsp->STO[i] = NULL;
			}
		} else {
			DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
			// need to diagonalize final propagator
			T = blk_cm_diag(Ud);
			for (i=k; i<l; i++) {
				// this should be faster:
				blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
				free_blk_mat_complx(wsp->STO[i]);
				wsp->STO[i] = NULL;
				blk_simtrans_adj(wsp->matrix[i+1],dumblk,sim);
				wsp->matrix[i+Ng+1] = cm_dup(det);
				blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
				free_blk_mat_complx(dumblk);
			}
			blk_simtrans_adj(det,T,sim);
			blk_simtrans_adj(wsp->matrix[k],T,sim);
			free_blk_mat_complx(T);
		}
		// control print out
		//for (i=0; i<Ng; i++) {
		//	cm_print(wsp->matrix[k+i],"RHO");
		//	cm_print(wsp->matrix[k+Ng+i],"DET");
		//}
	}
	// END of transformations


	//QueryPerformanceCounter(&tv1);
	freq = dptr = (double*)malloc(matdim*sizeof(double));
	for (i=0; i<Ud->Nblocks; i++) {
		z1 = Ud->m[i].data;
		assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
		for (j=0; j<Ud->blk_dims[i]; j++) {
			*dptr = -Carg((*z1))/sim->taur*1.0e6;
			z1++;
			dptr++;
		}
	}
	// control print-out
	//for (i=0; i<matdim; i++) printf("freq[%d] = %g\n",i,freq[i]);
	//QueryPerformanceCounter(&tv2);
	//printf("timing frequencies: %.9f\n",((float)(tv2.QuadPart-tv1.QuadPart))/(float)tickpsec.QuadPart);

	// fork from here to ED symmetry specific calculation
	//if (sim->EDsymmetry == 1) {
	//	if (sim->interpolation != 0) {
	//		fprintf(stderr,"Error: interpolation not implemented for ED symmetry\n");
	//		exit(1);
	//	}
	//
	//
	//	free(freq);
	//	free_blk_mat_complx(wsp->STO[l]); // this is Ud
	//	wsp->STO[l] = NULL;
	//	for (i=k; i<=l; i++) {
	//		free_complx_matrix(wsp->matrix[i+Ng]);
	//		wsp->matrix[i+Ng] = NULL;
	//	}
	//	return;
	//}
	// END of ED symmeetry code


	// get relevant elements
	TIMING_TIC(tv1);
    irow = (int*)malloc((matdim+1)*sizeof(int));
    icol = ic = (int*)malloc(matdim*matdim*sizeof(int)); // just large enough
    scan_contrib(&wsp->matrix[k],Ng,irow,icol);
	TIMING_TOC(tv1,tv2,"scan_contrib");

    if (sim->interpolation == 1) { // FWT interpolation
    	// now I am sure the file exists
    	if (ftell(wsp->interpol_file) == 0) {
    		// write header
    		r = irow[matdim]-1;
    		fwrite(&r,sizeof(int),1,wsp->interpol_file);
    		fwrite(irow,sizeof(int),matdim+1,wsp->interpol_file);
    		fwrite(icol,sizeof(int),r,wsp->interpol_file);
    	}
    	// write crystallite index
    	fwrite(&wsp->cryst_idx,sizeof(int),1,wsp->interpol_file);
    	// write diagonalized period propagator
    	for (i=0; i<Ud->Nblocks; i++) {
    		fwrite(Ud->m[i].data,sizeof(complx),Ud->blk_dims[i],wsp->interpol_file);
    	}
    }

    complx *phmul = (complx*)malloc(Ng*sizeof(complx));
    complx *Frs = complx_vector(Ng);
    for (r=1; r<=matdim; r++) {
    	int nc = irow[r] - irow[r-1];
    	for (c=0; c<nc; c++) {
    		if (sim->interpolation == 1) { // FWT interpolation
        		for (i=0; i<Ng; i++) {
        			complx zz1 = cm_getelem(wsp->matrix[k+i],*ic,r);
        			complx zz2 = cm_getelem(wsp->matrix[k+Ng+i],r,*ic);
        			fwrite(&zz1,sizeof(complx),1,wsp->interpol_file);
        			fwrite(&zz2,sizeof(complx),1,wsp->interpol_file);
        		}
        		ic++;
        		continue;
    		}
			diff = freq[r-1] - freq[*ic-1];
			complx ph = Cexpi(-diff*dtg*1.0e-6);
		    phmul[0] = Cunit;
			for (i=1; i<Ng; i++) {
				phmul[i] = Cmul(phmul[i-1],ph);
				//printf("\t phmul[%d] = %g , %g\n",i,phmul[i].re,phmul[i].im);
			}
			//ph = Cmul(phmul[Ng-1],ph);
			//ph.im = -ph.im;
			ph = Cexpi(diff*sim->taur*1.0e-6);
			//printf("\t ph = %g , %g\n",ph.re,ph.im);
			// to optimize the code above
			//complx ph, *phptr = phmul;
			//phptr->re = 1.0; phptr->im = 0.0; phptr++;
			//ph.re = phptr->re = cos(diff*dtg*1.0e-6);
			//ph.im = phptr->im = -sin(diff*dtg*1.0e-6);
			//for (i=2; i<Ng; i++) {
			//	phptr->re =
			//}
 			cv_zero(Frs);
 			int ppp = Ng;
			for (j=0; j<Ng; j++) {
				for (i=0; i<Ng; i++) {
					zz1 = cm_getelem(wsp->matrix[k+i],*ic,r);
					zz2 = cm_getelem(wsp->matrix[k+Ng+(i+j)%Ng],r,*ic);
					//printf("\t rho = (%g, %g), det = (%g, %g), ph = (%g, %g)\n",zz1.re,zz1.im, zz2.re, zz2.im, phmul[j].re, phmul[j].im);
					zz1 = Cmul(zz1,zz2);
					//zz2 = Cmul(zz1,phmul[(i+j)%Ng]);
					zz2 = Cmul(zz1,phmul[j]);
					if (i>=ppp) zz2 = Cmul(zz2,ph);
					Frs[j+1].re += zz2.re;
					Frs[j+1].im += zz2.im;
				}
				//phmul[j] = Cmul(phmul[j],ph);
				ppp--;
			} // Frs[] vector ready to use
			//printf("elem r = %d, s = %d\n",r, *ic);
			//cv_print(Frs,"Frs");
    		ph = Cunit;
    		zz2 = Cexpi(diff*dtg*m*1.0e-6);
    		for (i=0; i<wsp->Nacq; i++) {
    			//printf("\t i = %d, idx = %d\n",i,(i*m)%Ng + 1);
    			zz1 = Cmul(Frs[(i*m)%Ng + 1],ph);
    			wsp->fid[i+1].re += zz1.re;
    			wsp->fid[i+1].im += zz1.im;
    			ph = Cmul(ph,zz2);
    		}
			ic++;
    	}
    }

	if ( fabs(wsp->acqphase) > TINY ) {
    	cv_mulc(wsp->fid,acqph);
    }
    wsp->curr_nsig += wsp->Nacq - 1;
	cv_muld(wsp->fid, 1.0/(double)Ng);

    free(irow);
    free(icol);
    free(freq);
    free(phmul);
    free_complx_vector(Frs);

	//free_complx_matrix(mexp1);
	free_blk_mat_complx(wsp->STO[l]);
	wsp->STO[l] = NULL;
	for (i=k; i<=l; i++) {
		free_complx_matrix(wsp->matrix[i]);
		wsp->matrix[i] = NULL;
		free_complx_matrix(wsp->matrix[i+Ng]);
		wsp->matrix[i+Ng] = NULL;
	}

}

/*************
 * New code design for gCOMPUTE methods
 */

/****
 * evaluates pulse sequence to generate gamma-dependent sigma and props
 */
void gcompute_props(Sim_info *sim, Sim_wsp *wsp) {
	double dtg, t0;
	int ig, Ng;

	TIMING_INIT_VAR2(tv1,tv2);

	Ng = sim->ngamma;
	dtg = sim->taur/(double)Ng;
	// initialize counter for storing wsp->sigma and wsp->U
	wsp->acqblock_sto = ACQBLOCK_STO_INI;

	TIMING_TIC(tv1);
	for (ig=0;ig<Ng;ig++) {
		wsp->tstart = t0 = ig/(double)Ng*sim->taur;
		direct_propagate(sim,wsp);
		if (wsp->acqblock_sto == ACQBLOCK_STO_INI + ig) {
			// there was no acq_block command in pulseq, need to do the job here
			if (fabs(wsp->t - t0 - dtg) > TINY) {
				fprintf(stderr,"Error: gcompute - pulse sequence event not synchronized with gamma-averaging\n");
				fprintf(stderr,"                  pulseq duration = %g us, need to be %g\n",wsp->t-t0,dtg);
				exit(1);
			}
			if (ig == 0) {
				wsp->STO[wsp->acqblock_sto] = blk_cm_dup(wsp->U);
			} else {
				wsp->STO[wsp->acqblock_sto] = blk_cm_mul(wsp->U,wsp->STO[wsp->acqblock_sto - 1]);
			}
			// this input file regime does not allow initial density matrix to be changed
			acqblock_sto_incr(wsp);
		}
		wsp->cryst.gamma += 360.0/(double)Ng;
	}
	TIMING_TOC(tv1,tv2,"gCOMPUTE eval_pulseq");

	// some safety checks
	DEBUGPRINT("new_gcompute_freq: sigmas and propagators from %d to %d (inclusive)\n",ACQBLOCK_STO_INI,wsp->acqblock_sto - 1);
	if (wsp->acqblock_sto - ACQBLOCK_STO_INI != Ng) {
		fprintf(stderr,"Error: (i)gcompute - propagator count mismatch\n");
		exit(1);
	}
	if (wsp->acqblock_sto+Ng > ACQBLOCK_STO_END) {
		fprintf(stderr,"Error: (i)gcompute - overflow in matrix buffers\n");
		exit(1);
	}

}

/****
 * transformations of dens.matrix and detect op. to eigenbasis of period propagator
 */
void gcompute_transforms(Sim_info *sim, Sim_wsp *wsp) {
	blk_mat_complx *Ud, *T;
	mat_complx *det;
	int i, l, Ng;

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);

	Ng = sim->ngamma;
	l = wsp->acqblock_sto-1;
	Ud = wsp->STO[l];
	if (sim->acq_adjoint) {
		if (wsp->fdetect != sim->fdetect)
			cm_adjointi(wsp->fdetect);
		else
			wsp->fdetect = cm_adjoint(sim->fdetect);
	}
	det = wsp->matrix[l+1] = cm_change_basis_2(wsp->fdetect,Ud->basis,sim);
	if (wsp->matrix[ACQBLOCK_STO_INI] == NULL) {
		if (sim->EDsymmetry == 1) { // ED symmetry present
			if (blk_cm_isdiag(Ud)) {
				for (i=ACQBLOCK_STO_INI; i<l; i++) {
					DEBUGPRINT("new_gcompute: Ud is diagonal\n");
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
				}
			} else {
				DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
				T = blk_cm_diag(Ud);
				for (i=ACQBLOCK_STO_INI; i<l; i++) {
					// this should be faster:
					blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
					free_blk_mat_complx(dumblk);
				}
				blk_simtrans_adj(det,T,sim);
				free_blk_mat_complx(T);
			}
		} else {  // no ED symmetry but full job needs to be done
			wsp->matrix[ACQBLOCK_STO_INI] = cm_change_basis_2(wsp->fstart,Ud->basis,sim);
			if (blk_cm_isdiag(Ud)) {
				for (i=ACQBLOCK_STO_INI; i<l; i++) {
					DEBUGPRINT("new_gcompute: Ud is diagonal\n");
					wsp->matrix[i+1] = cm_dup(wsp->matrix[ACQBLOCK_STO_INI]);
					blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
				}
			} else {
				DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
				// need to diagonalize final propagator
				T = blk_cm_diag(Ud);
				for (i=ACQBLOCK_STO_INI; i<l; i++) {
					// this should be faster:
					blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
					free_blk_mat_complx(wsp->STO[i]);
					wsp->STO[i] = NULL;
					wsp->matrix[i+1] = cm_dup(wsp->matrix[ACQBLOCK_STO_INI]);
					blk_simtrans_adj(wsp->matrix[i+1],dumblk,sim);
					wsp->matrix[i+Ng+1] = cm_dup(det);
					blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
					free_blk_mat_complx(dumblk);
				}
				blk_simtrans_adj(det,T,sim);
				blk_simtrans_adj(wsp->matrix[ACQBLOCK_STO_INI],T,sim);
				free_blk_mat_complx(T);
			}
		}
	} else { // full job with gamma-dependent density matrices
		if (blk_cm_isdiag(Ud)) {
			for (i=ACQBLOCK_STO_INI; i<l; i++) {
				DEBUGPRINT("new_gcompute: Ud is diagonal\n");
				blk_simtrans_adj(wsp->matrix[i+1],wsp->STO[i],sim);
				wsp->matrix[i+Ng+1] = cm_dup(det);
				blk_simtrans_adj(wsp->matrix[i+Ng+1],wsp->STO[i],sim);
				free_blk_mat_complx(wsp->STO[i]);
				wsp->STO[i] = NULL;
			}
		} else {
			DEBUGPRINT("new_gcompute: will diagonalize Ud\n");
			// need to diagonalize final propagator
			T = blk_cm_diag(Ud);
			for (i=ACQBLOCK_STO_INI; i<l; i++) {
				// this should be faster:
				blk_mat_complx *dumblk = blk_cm_mul(wsp->STO[i],T);
				free_blk_mat_complx(wsp->STO[i]);
				wsp->STO[i] = NULL;
				blk_simtrans_adj(wsp->matrix[i+1],dumblk,sim);
				wsp->matrix[i+Ng+1] = cm_dup(det);
				blk_simtrans_adj(wsp->matrix[i+Ng+1],dumblk,sim);
				free_blk_mat_complx(dumblk);
			}
			blk_simtrans_adj(det,T,sim);
			blk_simtrans_adj(wsp->matrix[ACQBLOCK_STO_INI],T,sim);
			free_blk_mat_complx(T);
		}
	}

	TIMING_TOC(tv1,tv2,"gCOMPUTE transformations");
}


/*****
 * generate time domain FID
 */
void gcompute_fid(Sim_info *sim, Sim_wsp *wsp) {
	int i, ll, m, Ng, matdim, q, row, col;
	double dtg;
	complx *z1, *z2, *z3, *z4, *z5, zval;
	mat_complx *mexp1, *mexp2;
	blk_mat_complx *Ud;

	TIMING_INIT_VAR2(tv1,tv2);

	/* first check on timings */
	Ng = sim->ngamma;
	matdim = sim->matdim;
	dtg = sim->taur/(double)Ng;
	m = (int)floor(wsp->dw/dtg+1e-6);
	if ( fabs(dtg*m - wsp->dw) > TINY) {
		fprintf(stderr,"Error: (i)gcompute - bad synchronization of acquisition and gamma-averaging\n");
		fprintf(stderr,"       dwelltime(%g) can not be split in tauR(%g)/Ngamma(%d) steps\n",wsp->dw,sim->taur,sim->ngamma);
		exit(1);
	}

	gcompute_props(sim,wsp);
	gcompute_transforms(sim,wsp);

	TIMING_TIC(tv1);
	Ud = wsp->STO[wsp->acqblock_sto - 1];
	mexp1 = complx_matrix(matdim,matdim,MAT_DENSE_DIAG,0,Ud->basis);
	mexp2 = complx_matrix(matdim,matdim,MAT_DENSE_DIAG,0,Ud->basis);
	z2 = mexp1->data;
	z3 = mexp2->data;
	for (i=0; i<Ud->Nblocks; i++) {
		z1 = Ud->m[i].data;
		assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
		for (ll=0; ll<Ud->blk_dims[i]; ll++) {
			*z2 = *z3 = CRpow(*z1,1.0/(double)Ng);
			z1++; z2++; z3++;
		}
	}
	// averaging over gamma angles and mult by exp(-i w_rs ll tau)
	mat_complx **tmpmx = (mat_complx**)malloc(Ng*sizeof(mat_complx*));
	if (tmpmx == NULL) {
		fprintf(stderr,"Error: (i)gcompute fid - no more memory\n");
		exit(1);
	}
	if (sim->EDsymmetry == 1) {
		// not very efficient here....
		for (ll=0; ll<Ng; ll++) {
			wsp->matrix[ACQBLOCK_STO_INI+ll] = cm_adjoint(wsp->matrix[ACQBLOCK_STO_INI+Ng+ll]);
			cm_addto(wsp->matrix[ACQBLOCK_STO_INI+ll],wsp->matrix[ACQBLOCK_STO_INI+Ng+ll]);
			cm_muld(wsp->matrix[ACQBLOCK_STO_INI+ll],0.5);
		}
	}
	for (ll=0; ll<Ng; ll++) {
		tmpmx[ll] = complx_matrix(matdim,matdim,MAT_DENSE,0,Ud->basis);
		//cm_zero(wsp->STO[k+Ng+ll]);
		cm_zero(tmpmx[ll]);
		if (ll!=0) blk_simtrans_adj(wsp->matrix[ACQBLOCK_STO_INI+Ng+ll-1],Ud,sim);
		gamma_ave(tmpmx[ll],&(wsp->matrix[ACQBLOCK_STO_INI]),ll,Ng);
		if (ll!=0) {
			simtrans(tmpmx[ll],mexp1);
			cm_multo(mexp1,mexp2);
		}
	}
	//cm_print(tmpmx[ll],"Frs*exp(-i wrs dtg)");
	TIMING_TOC(tv1,tv2,"all gamma ave");

	// acquisition
	TIMING_TIC(tv1);
	complx *fidptr = wsp->fid;
	if (wsp->curr_nsig != 0) {
		fprintf(stderr,"Error: (i)gcompute - some FID points has already been acquired. That is not allowed.\n");
		exit(1);
	}
	const double m_Ng = (double)m/(double)Ng;
	z2 = mexp1->data;
	z3 = mexp2->data;
	for (i=0; i<Ud->Nblocks; i++) {
		z1 = Ud->m[i].data;
		assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
		for (ll=0; ll<Ud->blk_dims[i]; ll++) {
			*z2 = *z3 = CRpow(*z1,m_Ng);
			z1++; z2++; z3++;
		}
	}
	double zzre, zzim;
	for (ll=0; ll<wsp->Nacq; ll++) {
		Czero(zval);
		q = (ll*m)%Ng;
		//printf("\tll = %d, q = %d\n",ll,q);
		z3 = mexp1->data;
		z5 = mexp2->data;
		complx *zstart = tmpmx[q]->data;
		for (row=0; row<matdim; row++) {
			z1 = z2 = zstart + row*(matdim+1);
			zval.re += z1->re;
			zval.im += z1->im;
			z4 = z3;
			for (col=row+1; col<matdim; col++) {
				z1++;
				z2 += matdim;
				z4++;
				zzre = z3->re*z4->re + z3->im*z4->im;
				zzim = -z3->re*z4->im + z3->im*z4->re;
				zval.re += z1->re*zzre - z1->im*zzim;
				zval.im += z1->re*zzim + z1->im*zzre;
				zval.re += z2->re*zzre + z2->im*zzim;
				zval.im += -z2->re*zzim + z2->im*zzre;
				//if (fabs(z1->re) > TINY) {
				//	printf("r,c = %d, %d; exp = (%g, %g), F = (%g, %g)\n",row+1,col+1,zzre,zzim,z1->re,z1->im);
				//}
				//if (fabs(z2->re) > TINY) {
				//	printf("r,c = %d, %d; exp = (%g, %g), F = (%g, %g)\n",col+1,row+1,zzre,-zzim,z2->re,z2->im);
				//}

			}
			zzre = z3->re; zzim = z3->im;
			z3->re = zzre*z5->re - zzim*z5->im;
			z3->im = zzre*z5->im + zzim*z5->re;

			z3++;
			z5++;
		}
		//printf("zval = (%g, %g)\n",zval.re,zval.im);
		fidptr++;
		//printf("\ndense gcompute: add (%f,%f) to fid (%f,%f)\n",zval.re,zval.im,fidptr->re,fidptr->im);
		fidptr->re = zval.re;
		fidptr->im = zval.im;
	}
	wsp->curr_nsig += wsp->Nacq - 1;
	if ( fabs(wsp->acqphase) > TINY ) {
		zval.re = cos(wsp->acqphase*DEG2RAD)/(double)Ng;
		zval.im = sin(wsp->acqphase*DEG2RAD)/(double)Ng;
		cv_mulc(wsp->fid, zval);
	} else {
		cv_muld(wsp->fid,1.0/(double)Ng);
	}
	TIMING_TOC(tv1,tv2,"gcompute acquisition");

	// free memory
	free_complx_matrix(mexp1);
	free_complx_matrix(mexp2);
	free_blk_mat_complx(wsp->STO[wsp->acqblock_sto - 1]);
	wsp->STO[wsp->acqblock_sto - 1] = NULL;
	for (i=ACQBLOCK_STO_INI; i<wsp->acqblock_sto; i++) {
		free_complx_matrix(wsp->matrix[i]);
		wsp->matrix[i] = NULL;
		free_complx_matrix(wsp->matrix[i+Ng]);
		wsp->matrix[i+Ng] = NULL;
		free_complx_matrix(tmpmx[i-ACQBLOCK_STO_INI]);
	}
	free(tmpmx);

}

void scan_contrib_freq_ED(mat_complx **A, int Ng, int *irow, int *icol) {
	int i, r, c, m;
	int matdim = A[0]->row;

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);
    irow[0]=1;
    for (r=1; r<=matdim; r++) {
    	irow[r] = irow[r-1];
    	for (c=1; c<=matdim; c++) {
    		m = 0;
    		for (i=0; i<Ng; i++) {
    			complx zz = cm_getelem(A[i],r,c);
    			if (fabs(zz.re) > TINY || fabs(zz.im) > TINY) {
    				m = 1;
    				break;
    			}
    		}
    		if (m != 0) {
    			irow[r]++;
    			*ic = c;
    			ic++;
    		}
    	}
    }
    TIMING_TOC(tv1,tv2,"gCOMPUTE freq EDsym scan contrib");
}

complx * qdata_EDsym(double *freq, double dtg, mat_complx **A, int Ng, int *irow, int *icol) {
	int m, r, c, i;
	int *ic;
	int matdim = A[0]->row;
	int nnz = irow[matdim] - 1;
	double diff;

	TIMING_INIT_VAR2(tv1,tv2);
	TIMING_TIC(tv1);

    fftin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
    fftout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
    plan = fftw_plan_dft_1d(Ng, fftin1, fftout1, FFTW_FORWARD, FFTW_ESTIMATE);
    if ( !fftin || !fftout || !plan ) {
    	fprintf(stderr,"Error: gCOMPUTE freq - no more memory for fftw structures");
    	exit(1);
    }

    complx *qdata = (complx *)malloc(nnz*Ng*sizeof(complx));
    if (qdata == NULL) {
    	fprintf(stderr,"Error: gCOMPUTE freq - no more memory for qdata\n");
    	exit(1);
    }
    ic = icol;
    m = 0;
    for (r=1; r<=matdim; r++) {
    	int nc = irow[r] - irow[r-1];
    	for (c=0; c<nc; c++) {
			diff = freq[r-1] - freq[*ic-1];
			complx ph = Cexpi(-diff*dtg*1.0e-6);
			complx phmul = Complx(1.0,0.0);
    		for (i=0; i<Ng; i++) {
    			complx zz1 = cm_getelem(A[i],r,*ic);  //q_rs(k)
    			fftin[i][0] = zz1.re*phmul.re - zz1.im*phmul.im;
    			fftin[i][1] = zz1.re*phmul.im + zz1.im*phmul.re;
    			phmul = Cmul(phmul,ph);
    		}
    		fftw_execute(plan);
    		//printf("\nr = %d, c = %d : ",r,*ic);
    		for (i=0; i<Ng; i++) {
    			qdata[i*nnz+m].re = fftout1[i][0];
    			qdata[i*nnz+m].im = fftout1[i][1];
    			//printf("( %g, %g ) ",fftout1[i][0],fftout1[i][1]);
    		}
    		//printf("\n");
    		m++;
    		ic++;
    	}
    }
    TIMING_TOC(tv1,tv2,"gCOMPUTE freq EDsym FFT of Q matrix");
    return qdata;
}

double * freq_from_Ud(blk_mat_complx *Ud, double period) {
	double *freq, *dptr;
	int i, m;

	freq = dptr = (double*)malloc(Ud->dim*sizeof(double));
	if (freq == NULL) {
		fprintf(stderr,"Error: gCOMPUTE freq - no more memory for freq\n");
		exit(1);
	}

	for (i=0; i<Ud->Nblocks; i++) {
		complx *z1 = Ud->m[i].data;
		assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
		for (m=0; m<Ud->blk_dims[i]; m++) {
			*dptr = -Carg((*z1))/period*1.0e6;
			z1++;
			dptr++;
		}
	}
	// control print-out
	//for (i=0; i<Ud->matdim; i++) printf("freq[%d] = %g\n",i,freq[i]);
	return freq;
}

/*****
 * generate frequency domain spectrum
 */
void gcompute_spc(Sim_info *sim, Sim_wsp *wsp) {
	int Ng, matdim, r, c, m, i, nnz;
	int *irow, *icol, *ic;
	double binsize, diff, dtg;
	blk_mat_complx *Ud;

	TIMING_INIT_VAR2(tv1,tv2);

	gcompute_props(Sim_info *sim, Sim_wsp *wsp);
	gcompute_transforms(sim,wsp);

	TIMING_TIC(tv1);

	Ud = wsp->STO[wsp->acqblock_sto - 1];
	freq = freq_from_Ud(Ud,sim->taur);

	Ng = sim->ngamma;
	matdim = sim->matdim;
	binsize = sim->sw*2*M_PI/sim->np;
	dtg = sim->taur/(double)Ng;
    irow = (int*)malloc((matdim+1)*sizeof(int));
    icol = ic = (int*)malloc(matdim*matdim*sizeof(int)); // just large enough
    if ( !irow || !icol ) {
    	fprintf(stderr,"Error: gCOMPUTE freq - no more memory for irow, icol");
    	exit(1);
    }
	if (wsp->curr_nsig != 0) {
		fprintf(stderr,"Error: gCOMPUTE freq - some FID points has already been acquired. That is not allowed.\n");
		exit(1);
	}

	if (sim->EDsymmetry == 1) {
		scan_contrib_freq_ED(&wsp->matrix[ACQBLOCK_STO_INI+Ng],Ng,irow,icol);
	    nnz = irow[matdim] - 1;
	    complx *qdata =  qdata_EDsym(freq, dtg, &wsp->matrix[ACQBLOCK_STO_INI+Ng], Ng, irow, icol);
	    mat_complx dum;
	    dum.irow = irow;
	    dum.icol = icol;
	    dum.type = MAT_SPARSE;
	    dum.row = dum.col = matdim;
	    ic = icol;
	    m = 0;
	    for (r=1; r<=matdim; r++) {
	    	int nc = irow[r] - irow[r-1];
	    	for (c=0; c<nc; c++) {
				diff = freq[r-1] - freq[*ic-1];
				for (i=0; i<Ng; i++) {
					int idx = (i + Ng/2 + 1) % Ng;
					complx zz1 = qdata[idx*nnz+m];
					dum.data = &qdata[( ( ((Ng-i)%Ng) + Ng/2+1 )%Ng )*nnz];
					complx zz2 = cm_getelem(&dum,*ic,r); // q_sr(-k)
					double zzre = (zz1.re*zz1.re+zz1.im*zz1.im + zz1.re*zz2.re-zz1.im*zz2.im)*0.5;
					double zzim = (zz1.re*zz2.im+zz1.im*zz2.re)*0.5;
					bin = (int)(1.5-(diff + sim->wr*(i-Ng/2+1))/binsize + sim->np/2);
					//printf("i = %d, idx = %d, -i = %d, idx = %d : freq = %g, bin %d -> ",i,idx,(Ng-i)%Ng,( ((Ng-i)%Ng) + Ng/2+1 )%Ng,diff+sim->wr*(i-Ng/2+1),bin);
					while (bin < 1) bin += sim->np;
					while (bin > sim->np) bin -= sim->np;
					assert(bin >= 1 && bin <= sim->np);
					wsp->fid[bin].re += zzre*cosph - zzim*sinph;
					wsp->fid[bin].im += zzim*cosph + zzre*sinph;
					//printf(" %d : ( %g, %g ) ( %g, %g ) -> ( %g, %g )\n",bin,zz1.re,zz1.im,zz2.re,zz2.im,zzre,zzim);
				}
				ic++;
				m++;
	    	}
	    }
		free(qdata);
		// END of ED symmetry calculation
	} else {
	    fftw_complex *fftin1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
	    fftw_complex *fftout1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
	    fftw_plan p1 = fftw_plan_dft_1d(Ng, fftin1, fftout1, FFTW_FORWARD, FFTW_ESTIMATE);
	    fftw_complex *fftin2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
	    fftw_complex *fftout2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
	    fftw_plan p2 = fftw_plan_dft_1d(Ng, fftin2, fftout2, FFTW_FORWARD, FFTW_ESTIMATE);
	    if (!fftin1 || !fftin2 || !fftout1 || !fftout2 || !p1 || !p2) {
	    	fprintf(stderr,"Error: gCOMPUTE freq - no more memory for FFTW structures\n");
	    	exit(1);
	    }

	    scan_contrib(&wsp->matrix[ACQBLOCK_STO_INI],Ng,irow,icol);

	    for (r=1; r<=matdim; r++) {
	    	int nc = irow[r] - irow[r-1];
	    	for (c=0; c<nc; c++) {
				diff = freq[r-1] - freq[*ic-1];
				complx ph = Cexpi(-diff*dtg*1.0e-6);
				complx phmul = Complx(1.0,0.0);
	    		for (i=0; i<Ng; i++) {
	    			complx zz1 = cm_getelem(wsp->matrix[ACQBLOCK_STO_INI+i],*ic,r);
	    			complx zz2 = cm_getelem(wsp->matrix[ACQBLOCK_STO_INI+Ng+i],r,*ic);
	    			fftin1[i][0] = zz1.re*phmul.re + zz1.im*phmul.im;
	    			fftin1[i][1] = zz1.re*phmul.im - zz1.im*phmul.re;
	    			fftin2[i][0] = zz2.re*phmul.re - zz2.im*phmul.im;
	    			fftin2[i][1] = zz2.re*phmul.im + zz2.im*phmul.re;
	    			phmul = Cmul(phmul,ph);
	    		}
	    		fftw_execute(p1);
	    		fftw_execute(p2);
				for (i=0; i<Ng; i++) {
					int idx = (i + Ng/2 + 1) % Ng;
					double zzre = fftout1[idx][0]*fftout2[idx][0] + fftout1[idx][1]*fftout2[idx][1];
					double zzim = fftout1[idx][0]*fftout2[idx][1] - fftout1[idx][1]*fftout2[idx][0];
					bin = (int)(1.5-(diff + sim->wr*(i-Ng/2+1))/binsize + sim->np/2);
					//printf("i = %d, idx = %d : freq = %g, bin %d -> ",i,idx,diff+sim->wr*(i-Ng/2+1),bin);
					while (bin < 1) bin += sim->np;
					while (bin > sim->np) bin -= sim->np;
					assert(bin >= 1 && bin <= sim->np);
					wsp->fid[bin].re += zzre*cosph - zzim*sinph;
					wsp->fid[bin].im += zzim*cosph + zzre*sinph;
					//printf("%d : (%g, %g) (%g, %g) (%g, %g)\n",bin,fftout1[idx][0],fftout1[idx][1],fftout2[idx][0],fftout2[idx][1],zzre,zzim);
				}
				ic++;
	    	}
	    }
	    fftw_destroy_plan(p1);
	    fftw_destroy_plan(p2);
	    fftw_free(fftin1); fftw_free(fftout1);
	    fftw_free(fftin2); fftw_free(fftout2);
	    // END on calculations without ED symmetry
	}

	// freeing memory
	free(irow);
	free(icol);
	free(freq);
	free_blk_mat_complx(wsp->STO[wsp->acqblock_sto-1]);
	wsp->STO[wsp->acqblock_sto-1] = NULL;
	for (i=ACQBLOCK_STO_INI; i<wsp->acqblock_sto; i++) {
		if (wsp->matrix[i] != NULL) free_complx_matrix(wsp->matrix[i]);
		free_complx_matrix(wsp->matrix[i+Ng]);
		wsp->matrix[i+Ng] = NULL;
	}

	TIMING_TOC(tv1,tv2,"gCOMPUTE spectrum generation in TOTAL");

}

/*****
 * generate frequency domain data for ASG interpolation
 */
void gcompute_ASGdata(Sim_info *sim, Sim_wsp *wsp) {
	int Ng, matdim, r, c, m, i, nnz;
	int *irow, *icol, *ic;
	double binsize, diff, dtg;
	blk_mat_complx *Ud;

	TIMING_INIT_VAR2(tv1,tv2);

	gcompute_props(Sim_info *sim, Sim_wsp *wsp);
	gcompute_transforms(sim,wsp);

	TIMING_TIC(tv1);

	Ud = wsp->STO[wsp->acqblock_sto - 1];
	freq = freq_from_Ud(Ud,sim->taur);

	Ng = sim->ngamma;
	matdim = sim->matdim;
	dtg = sim->taur/(double)Ng;
    irow = (int*)malloc((matdim+1)*sizeof(int));
    icol = ic = (int*)malloc(matdim*matdim*sizeof(int)); // just large enough
    if ( !irow || !icol ) {
    	fprintf(stderr,"Error: gCOMPUTE freq - no more memory for irow, icol");
    	exit(1);
    }

	if (sim->EDsymmetry == 1) {
		scan_contrib_freq_ED(&wsp->matrix[ACQBLOCK_STO_INI+Ng],Ng,irow,icol);
	    nnz = irow[matdim] - 1;
	    complx *qdata =  qdata_EDsym(freq, dtg, &wsp->matrix[ACQBLOCK_STO_INI+Ng], Ng, irow, icol);
	    mat_complx dum;
	    dum.irow = irow;
	    dum.icol = icol;
	    dum.type = MAT_SPARSE;
	    dum.row = dum.col = matdim;
	    ic = icol;
	    m = 0;
	    for (r=1; r<=matdim; r++) {
	    	int nc = irow[r] - irow[r-1];
	    	for (c=0; c<nc; c++) {
				sim->ASG_freq[(wsp->cryst_idx-1)*nnz + m] = freq[r-1] - freq[*ic-1];
				for (i=0; i<Ng; i++) {
					int idx = (i + Ng/2 + 1) % Ng;
					complx *zz1 = &qdata[idx*nnz+m];
					dum.data = &qdata[( ( ((Ng-i)%Ng) + Ng/2+1 )%Ng )*nnz];
					complx zz2 = cm_getelem(&dum,*ic,r); // q_sr(-k)
					double zzre = (zz1->re*zz1->re+zz1->im*zz1->im + zz1->re*zz2.re-zz1->im*zz2.im)*0.5;
					double zzim = (zz1->re*zz2.im+zz1->im*zz2.re)*0.5;
					zz1 = &sim->ASG_ampl[(wsp->cryst_idx-1)*nnz*Ng+m*Ng+i];
					zz1->re = zzre;
					zz1->im = zzim;
				}
				ic++;
				m++;
	    	}
	    }
		free(qdata);
		// END of ED symmetry calculation
	} else {
	    fftw_complex *fftin1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
	    fftw_complex *fftout1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
	    fftw_plan p1 = fftw_plan_dft_1d(Ng, fftin1, fftout1, FFTW_FORWARD, FFTW_ESTIMATE);
	    fftw_complex *fftin2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
	    fftw_complex *fftout2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ng);
	    fftw_plan p2 = fftw_plan_dft_1d(Ng, fftin2, fftout2, FFTW_FORWARD, FFTW_ESTIMATE);
	    if (!fftin1 || !fftin2 || !fftout1 || !fftout2 || !p1 || !p2) {
	    	fprintf(stderr,"Error: gCOMPUTE freq - no more memory for FFTW structures\n");
	    	exit(1);
	    }

	    scan_contrib(&wsp->matrix[ACQBLOCK_STO_INI],Ng,irow,icol);
	    nnz = irow[matdim] - 1;
	    m = 0;
	    for (r=1; r<=matdim; r++) {
	    	int nc = irow[r] - irow[r-1];
	    	for (c=0; c<nc; c++) {
				diff = sim->ASG_freq[(wsp->cryst_idx-1)*nnz + m] = freq[r-1] - freq[*ic-1];
				complx ph = Cexpi(-diff*dtg*1.0e-6);
				complx phmul = Complx(1.0,0.0);
	    		for (i=0; i<Ng; i++) {
	    			complx zz1 = cm_getelem(wsp->matrix[ACQBLOCK_STO_INI+i],*ic,r);
	    			complx zz2 = cm_getelem(wsp->matrix[ACQBLOCK_STO_INI+Ng+i],r,*ic);
	    			fftin1[i][0] = zz1.re*phmul.re + zz1.im*phmul.im;
	    			fftin1[i][1] = zz1.re*phmul.im - zz1.im*phmul.re;
	    			fftin2[i][0] = zz2.re*phmul.re - zz2.im*phmul.im;
	    			fftin2[i][1] = zz2.re*phmul.im + zz2.im*phmul.re;
	    			phmul = Cmul(phmul,ph);
	    		}
	    		fftw_execute(p1);
	    		fftw_execute(p2);
				for (i=0; i<Ng; i++) {
					int idx = (i + Ng/2 + 1) % Ng;
					complx *zz1 = &sim->ASG_ampl[(wsp->cryst_idx-1)*nnz*Ng+m*Ng+i];
					zz1->re = fftout1[idx][0]*fftout2[idx][0] + fftout1[idx][1]*fftout2[idx][1];
					zz1->im = fftout1[idx][0]*fftout2[idx][1] - fftout1[idx][1]*fftout2[idx][0];
				}
				ic++;
				m++;
	    	}
	    }
	    fftw_destroy_plan(p1);
	    fftw_destroy_plan(p2);
	    fftw_free(fftin1); fftw_free(fftout1);
	    fftw_free(fftin2); fftw_free(fftout2);
	    // END on calculations without ED symmetry
	}

	// freeing memory
	free(irow);
	free(icol);
	free(freq);
	free_blk_mat_complx(wsp->STO[wsp->acqblock_sto-1]);
	wsp->STO[wsp->acqblock_sto-1] = NULL;
	for (i=ACQBLOCK_STO_INI; i<wsp->acqblock_sto; i++) {
		if (wsp->matrix[i] != NULL) free_complx_matrix(wsp->matrix[i]);
		free_complx_matrix(wsp->matrix[i+Ng]);
		wsp->matrix[i+Ng] = NULL;
	}

	TIMING_TOC(tv1,tv2,"gCOMPUTE ASG data generation in TOTAL");
}

/*****
 * generate frequency domain data for ASG interpolation
 */
void gcompute_FWTdata(Sim_info *sim, Sim_wsp *wsp) {
	int Ng, matdim, r, c, m, i, nnz;
	int *irow, *icol, *ic;
	double dtg;
	complx *z1, *z2;
	blk_mat_complx *Ud;

	TIMING_INIT_VAR2(tv1,tv2);

	gcompute_props(Sim_info *sim, Sim_wsp *wsp);
	gcompute_transforms(sim,wsp);

	TIMING_TIC(tv1);

	Ng = sim->ngamma;
	matdim = sim->matdim;
	Ud = wsp->STO[wsp->acqblock_sto - 1];
	z2 = &sim->FWT_lam[wsp->cryst_idx-1];
	for (i=0; i<Ud->Nblocks; i++) {
		*z1 = Ud->m[i].data;
		assert(Ud->m[i].type == MAT_DENSE_DIAG || Ud->blk_dims[i] == 1);
		for (m=0; m<Ud->blk_dims[i]; m++) {
			*z2 = *z1;
			z1++;
			z2 += LEN(sim->crdata);
		}
	}
	free_blk_mat_complx(wsp->STO[wsp->acqblock_sto-1]);
	wsp->STO[wsp->acqblock_sto-1] = NULL;

    irow = (int*)malloc((matdim+1)*sizeof(int));
    icol = ic = (int*)malloc(matdim*matdim*sizeof(int)); // just large enough
    if ( !irow || !icol ) {
    	fprintf(stderr,"Error: gCOMPUTE freq - no more memory for irow, icol");
    	exit(1);
    }

	if (sim->EDsymmetry == 1) {
		scan_contrib_freq_ED(&wsp->matrix[ACQBLOCK_STO_INI+Ng],Ng,irow,icol);
	    ic = icol;
	    z2 = &sim->FWT_frs[wsp->cryst_idx-1];
	    for (r=1; r<=matdim; r++) {
	    	int nc = irow[r] - irow[r-1];
	    	for (c=0; c<nc; c++) {
				for (i=0; i<Ng; i++) {
					*z2 = cm_getelem(wsp->matrix[ACQBLOCK_STO_INI+Ng+i],r,*ic); // q_rs(k)
					z2 += LEN(sim->crdata);
				}
				ic++;
	    	}
	    }
		// END of ED symmetry calculation
	} else {
	    scan_contrib(&wsp->matrix[ACQBLOCK_STO_INI],Ng,irow,icol);
	    z1 = &sim->FWT_frs[wsp->cryst_idx-1];
	    z2 = &sim->FWT_frs[wsp->cryst_idx-1+LEN(sim->crdata)];
	    for (r=1; r<=matdim; r++) {
	    	int nc = irow[r] - irow[r-1];
	    	for (c=0; c<nc; c++) {
	    		for (i=0; i<Ng; i++) {
	    			*z1 = cm_getelem(wsp->matrix[ACQBLOCK_STO_INI+i],*ic,r);
	    			*z2 = cm_getelem(wsp->matrix[ACQBLOCK_STO_INI+Ng+i],r,*ic);
	    			z1 += 2*LEN(sim->crdata);
	    			z2 += 2*LEN(sim->crdata);
	    		}
				ic++;
				m++;
	    	}
	    }
	    // END on calculations without ED symmetry
	}

	// freeing memory
	free(irow);
	free(icol);
	for (i=ACQBLOCK_STO_INI; i<wsp->acqblock_sto; i++) {
		if (wsp->matrix[i] != NULL) free_complx_matrix(wsp->matrix[i]);
		free_complx_matrix(wsp->matrix[i+Ng]);
		wsp->matrix[i+Ng] = NULL;
	}

	TIMING_TOC(tv1,tv2,"gCOMPUTE FWT data generation in TOTAL");
}