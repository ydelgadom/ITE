/********************************************************
* Functions in MATVEC_OPERATIONS.H
*
*	ranvec ( *vec )
*	double = norm( vec )
*	complex = dot( vec1, vec2 )
* complex = dotu( vec1, vec2 )
* vec_Mvec( *inv, *outv )
* vec_Mbvec( *inv, *outv )
* vec_Mvec_Mbvec( *inv, *outMv, *outMbv )
* vec_Dwvec( *inv, *outv )
* bool = vec_invDwvec( *inv, *outv )
* gamma5( *vec )
*
********************************************************/

// Lapack and Blas routines
extern "C" 
{ 
dcomplex zdotc_(int* n,dcomplex* zx,int* incx,dcomplex *zy,int* incy);
dcomplex zdotu_(int* n,dcomplex* zx,int* incx,dcomplex *zy,int* incy);
double dznrm2_(int* n,dcomplex* x,int* incx);
// y = alpha*x + y
void zaxpy_(int* n,dcomplex* alpha,dcomplex* x,int* incx,dcomplex* y,int* incy );
// x = alpha*x
void zscal_(int* n,dcomplex* alpha,dcomplex* x,int* incx );
}

//-------------------------------------------------------------------------------------
void ranvec( cvector *v )
{
/*
		Set elements of v to {-1,+1} randomly.
*/
	double *ran = new double[(*v).size()];

	ranlxd( ran, (*v).size() );

#ifdef MKL
#pragma omp parallel for collapse(1)
#endif
	for (unsigned i=0; i<(*v).size( ); i++)
	{
		if (ran[i]<0.5)
			(*v)[i] = ONE;
		else
			(*v)[i] = mONE;	
	}

	delete[] ran;
}

//----------------------------------------------------------------------------------------
double norm( cvector inv )
{
/* 
		sqrt(v^+ v)
*/

	int size = int(inv.size());
	int inc = 1;

	return dznrm2_( &size, inv.data(), &inc );
}

//----------------------------------------------------------------------------------------------------------
dcomplex dot( cvector inv1, cvector inv2 )
{
/*
		result = sum_i { conj[inv1(i)]*inv2(i) }
		ie. dot product.  Taking the complex conjugate of inv1
*/
	int size = int(inv1.size());
	int inc = 1;

	return zdotc_( &size, inv1.data(), &inc, inv2.data(), &inc );
}

//----------------------------------------------------------------------------------------------------------
dcomplex dotu( cvector inv1, cvector inv2 )
{
/*
		result = sum_i { inv1(i)*inv2(i) }
		ie. dot product without complex conjugate
*/
	int size = int(inv1.size());
	int inc = 1;

	return zdotu_( &size, inv1.data(), &inc, inv2.data(), &inc );
}

#ifdef MM2
//----------------------------------------------------------------------------------------------------------
void vec_Mvec( cvector *inv, cvector *outv )
{
/*
		b = M*a (b=outv, a=inv)
		M = -(1-gamma4)/2 * delta(isb,isa) * delta(icb,ica) 
										  * delta(itb,Nt-1) * delta(ita,0)

	NOTE: Vector index: 
				i = i_t*12*Ns3 + i_3*12*Ns*Ns + i_2*12*Ns + i_1*12 + 4*ic + id
*/
	cvector a((*inv).size( ));
	a = (*inv);

	int itb = 12*(Nt-1)*Ns3;

	// outb = 0
	(*outv).assign( (*outv).size( ), ZERO );

#ifdef MKL
#pragma omp parallel for collapse(2)
#endif
	for ( int i=0; i<12*Ns3; i+=4)
	{
		(*outv)[itb+i+0] = -0.5*( a[i+0] - a[i+2] );
		(*outv)[itb+i+1] = -0.5*( a[i+1] - a[i+3] );
		(*outv)[itb+i+2] = -(*outv)[itb+i+0];
		(*outv)[itb+i+3] = -(*outv)[itb+i+1];
	}

}

//----------------------------------------------------------------------------------------------------------
void vec_Mbvec( cvector *inv, cvector *outv )
{
/*
		b = Mbar*a (b=outv, a=inv)
		Mbar = -(1+gamma4)/2 * delta(isb,isa) * delta(icb,ica) 
											 	 * delta(itb,0) * delta(ita,Nt-1)

	NOTE: Vector index: 
				i = i_t*12*Ns3 + i_3*12*Ns*Ns + i_2*12*Ns + i_1*12 + 4*ic + id
*/

	cvector a((*inv).size( ));
	a = (*inv);

	int ita = 12*(Nt-1)*Ns3;

	// outb = 0
	(*outv).assign( (*outv).size( ), ZERO );

#ifdef MKL
#pragma omp parallel for collapse(2)
#endif
	for ( int i=0; i<12*Ns3; i+=4)
	{
		(*outv)[i+0] = -0.5*( a[ita+i+0] + a[ita+i+2] );
		(*outv)[i+1] = -0.5*( a[ita+i+1] + a[ita+i+3] );
		(*outv)[i+2] = (*outv)[i+0];
		(*outv)[i+3] = (*outv)[i+1];
	}

}

//----------------------------------------------------------------------------------------------------------
void vec_Mvec_Mbvec( cvector *inv, cvector *vm, cvector *vmb )
{
/*
		b = M*a (b=vm, a=inv)
		M = -(1-gamma4)/2 * delta(isb,isa) * delta(icb,ica) 
										  * delta(itb,Nt-1) * delta(ita,0)

		b = Mbar*a (b=vmb, a=inv)
		Mbar = -(1+gamma4)/2 * delta(isb,isa) * delta(icb,ica) 
											 	 * delta(itb,0) * delta(ita,Nt-1)

	NOTE: Vector index: 
				i = i_t*12*Ns3 + i_3*12*Ns*Ns + i_2*12*Ns + i_1*12 + 4*ic + id
*/
	cvector a((*inv).size( ));
	a = (*inv);

	// out = 0
	(*vm).assign( (*vm).size( ), ZERO );
	(*vmb).assign( (*vmb).size( ), ZERO );

	int it = 12*(Nt-1)*Ns3;

#ifdef MKL
#pragma omp parallel for collapse(2)
#endif
	for ( int i=0; i<12*Ns3; i+=4)
	{
		// elements for M
		(*vm)[it+i+0] = -0.5*( a[i+0] - a[i+2] );
		(*vm)[it+i+1] = -0.5*( a[i+1] - a[i+3] );
		(*vm)[it+i+2] = -(*vm)[it+i+0];
		(*vm)[it+i+3] = -(*vm)[it+i+1];

		// elements for Mbar
		(*vmb)[i+0] = -0.5*( a[it+i+0] + a[it+i+2] );
		(*vmb)[i+1] = -0.5*( a[it+i+1] + a[it+i+3] );
		(*vmb)[i+2] = (*vmb)[i+0];
		(*vmb)[i+3] = (*vmb)[i+1];
	}

}
#endif

#ifdef MM1
//----------------------------------------------------------------------------------------------------------
void vec_Mvec( cvector *inv, cvector *outv )
{
/*
		b = M*a (b=outv, a=inv)
		M = (1-gamma4)/2 * U_4(i) * delta(i+4,j)

	NOTE: Vector index: 
				i = i_t*12*Ns3 + i_3*12*Ns*Ns + i_2*12*Ns + i_1*12 + 4*ic + id
*/
	cvector a( (*inv).size() );
	a = (*inv);

	dcomplex as3032, as3133;

  int isp3, ica4, icb4, is12;

	// outb = 0
	(*outv).assign( (*outv).size( ), ZERO );

	for (int is=0; is<nsite; is++ )
	{
		is12 = 12*is;
    isp3 = 12*config->neib(is,3);

		for (int icb=0; icb<3; icb++ )
		{
			icb4 = icb*4;

			for (int ica=0; ica<3; ica++ )
			{
				ica4 = ica*4;

				as3032 = 0.5*config->U(is,3,icb,ica)*
                  ( a[isp3+ica4+0] - a[isp3+ica4+2] );

				as3133 = 0.5*config->U(is,3,icb,ica)*
                  ( a[isp3+ica4+1] - a[isp3+ica4+3] );

        (*outv)[is12+icb4+0] += as3032;
        (*outv)[is12+icb4+1] += as3133;
        (*outv)[is12+icb4+2] -= as3032;
        (*outv)[is12+icb4+3] -= as3133;
  		}	// END color a
  	} // END color b
	} // END nsite
}

//----------------------------------------------------------------------------------------------------------
void vec_Mbvec( cvector *inv, cvector *outv )
{
/*
		b = Mbar*a (b=outv, a=inv)
		Mbar = (1+gamma4)/2 * U_4(i-4)^+ * delta(i-4,j)

	NOTE: Vector index: 
				i = i_t*12*Ns3 + i_3*12*Ns*Ns + i_2*12*Ns + i_1*12 + 4*ic + id
*/
	cvector a( (*inv).size() );
	a = (*inv);

	dcomplex ad3032, ad3133;

  int ism3, is3, ica4, icb4, is12;

	// outb = 0
	(*outv).assign( (*outv).size( ), ZERO );

	for (int is=0; is<nsite; is++ )
	{
		is12 = 12*is;
    is3  = config->neib(is,7);
		ism3 = 12*is3;

		for (int icb=0; icb<3; icb++ )
		{
			icb4 = icb*4;

			for (int ica=0; ica<3; ica++ )
			{
				ica4 = ica*4;

				ad3032 = 0.5*std::conj(config->U(is3,3,ica,icb))*
                  ( a[ism3+ica4+0] + a[ism3+ica4+2] );

				ad3133 = 0.5*std::conj(config->U(is3,3,ica,icb))*
                  ( a[ism3+ica4+1] + a[ism3+ica4+3] );
				

        (*outv)[is12+icb4+0] += ad3032;
        (*outv)[is12+icb4+1] += ad3133;
        (*outv)[is12+icb4+2] += ad3032;
        (*outv)[is12+icb4+3] += ad3133;
  		}	// END color a
  	} // END color b
	} // END nsite
}

//----------------------------------------------------------------------------------------------------------
void vec_Mvec_Mbvec( cvector *inv, cvector *vm, cvector *vmb )
{
/*
		b = M*a (b=vm, a=inv)
		M = (1-gamma4)/2 * U_4(i) * delta(i+4,j)

		b = Mbar*a (b=vmb, a=inv)
		Mbar = (1+gamma4)/2 * U_4(i-4)^+ * delta(i-4,j)

	NOTE: Vector index: 
				i = i_t*12*Ns3 + i_3*12*Ns*Ns + i_2*12*Ns + i_1*12 + 4*ic + id
*/
	cvector a( (*inv).size() );
	a = (*inv);

	dcomplex a3032, a3133;
  int ism3, is3, ica4, icb4, is12, isp3;

	// outb = 0
	(*vm).assign( (*vm).size( ), ZERO );
	(*vmb).assign( (*vmb).size( ), ZERO );

	for (int is=0; is<nsite; is++ )
	{
		is12 = 12*is;
    isp3 = 12*config->neib(is,3);
    is3  = config->neib(is,7);
		ism3 = 12*is3;

		for (int icb=0; icb<3; icb++ )
		{
			icb4 = icb*4;

			for (int ica=0; ica<3; ica++ )
			{
				ica4 = ica*4;

				// elements for M
				a3032 = 0.5*config->U(is,3,icb,ica)*
                  ( a[isp3+ica4+0] - a[isp3+ica4+2] );

				a3133 = 0.5*config->U(is,3,icb,ica)*
                  ( a[isp3+ica4+1] - a[isp3+ica4+3] );

        (*vm)[is12+icb4+0] += a3032;
        (*vm)[is12+icb4+1] += a3133;
        (*vm)[is12+icb4+2] -= a3032;
        (*vm)[is12+icb4+3] -= a3133;


				// elements for Mbar
				a3032 = 0.5*std::conj(config->U(is3,3,ica,icb))*
                  ( a[ism3+ica4+0] + a[ism3+ica4+2] );

				a3133 = 0.5*std::conj(config->U(is3,3,ica,icb))*
                  ( a[ism3+ica4+1] + a[ism3+ica4+3] );
				

        (*vmb)[is12+icb4+0] += a3032;
        (*vmb)[is12+icb4+1] += a3133;
        (*vmb)[is12+icb4+2] += a3032;
        (*vmb)[is12+icb4+3] += a3133;
  		}	// END color a
  	} // END color b
	} // END nsite
}
#endif


//----------------------------------------------------------------------------------------------------------
void vec_Dwvec( cvector *inv, cvector *outv )
{
/*
	   Evaluates outv = Dw inv.

	   for Dw = Wilson Dirac operator of the form
  	 Dw(x,y) = (m+4)   - (1/2) sum_mu [
  	                       (1-g_mu) U(x,mu) delta(x+mu,y)
  	                    + (1+g_mu) U(x-mu,mu)^+ delta(x-mu,y)]
  	     or    1 - kappa sum_mu [
  	                     (1-g_mu) U(x,mu) delta(x+mu,y)
  	                    +(1+g_mu) U(x-mu,mu)^+ delta(x-mu,y)]
  	 kappa = 1/( 2(m+4) )

	NOTE: Vector index: 
				i = i_t*12*Ns3 + i_3*12*Ns*Ns + i_2*12*Ns + i_1*12 + 4*ic + id
*/      

	cvector a( (*inv).size() );
	a = (*inv);

	dcomplex as0003, ad0003, as1013, ad1013;
  dcomplex as2022, ad2022, as3032, ad3032;
  dcomplex as0102, ad0102, as1112, ad1112;
  dcomplex as2123, ad2123, as3133, ad3133;

  int isp0, isp1, isp2, isp3, ism0, ism1, ism2, ism3;
	int is0, is1, is2, is3, is12, ica4, icb4;

	int size = int(a.size( ));
	int inc = 1;
	dcomplex alpha = complex<double>(4.0+mass,0);

	(*outv) = a;
	zscal_(&size,&alpha,(*outv).data(),&inc );

	for (int is=0; is<nsite; is++ )
	{
		is12 = 12*is;
    isp0 = 12*config->neib(is,0);
    isp1 = 12*config->neib(is,1);
    isp2 = 12*config->neib(is,2);
    isp3 = 12*config->neib(is,3);
    is0  = config->neib(is,4);
		ism0 = 12*is0;
    is1  = config->neib(is,5);
		ism1 = 12*is1;
    is2  = config->neib(is,6);
		ism2 = 12*is2;
    is3  = config->neib(is,7);
		ism3 = 12*is3;

		for (int icb=0; icb<3; icb++ )
		{
			icb4 = icb*4;

			for (int ica=0; ica<3; ica++ )
			{
				ica4 = ica*4;

				as0003 = config->U(is,0,icb,ica)*
                  ( a[isp0+ica4+0] + II*a[isp0+ica4+3] );

				ad0003 = std::conj(config->U(is0,0,ica,icb))*
                  ( a[ism0+ica4+0] - II*a[ism0+ica4+3] );

				as1013 = config->U(is,1,icb,ica)*
                  ( a[isp1+ica4+0] + a[isp1+ica4+3] );

				ad1013 = std::conj(config->U(is1,1,ica,icb))*
                  ( a[ism1+ica4+0] - a[ism1+ica4+3] );

				as2022 = config->U(is,2,icb,ica)*
									( a[isp2+ica4+0] + II*a[isp2+ica4+2] );

				ad2022 = std::conj(config->U(is2,2,ica,icb))*
                  ( a[ism2+ica4+0] - II*a[ism2+ica4+2] );

				as3032 = config->U(is,3,icb,ica)*
                  ( a[isp3+ica4+0] - a[isp3+ica4+2] );

				ad3032 = std::conj(config->U(is3,3,ica,icb))*
                  ( a[ism3+ica4+0] + a[ism3+ica4+2] );

				as0102 = config->U(is,0,icb,ica)*
                  ( a[isp0+ica4+1] + II*a[isp0+ica4+2] );

				ad0102 = std::conj(config->U(is0,0,ica,icb))*
                  ( a[ism0+ica4+1] - II*a[ism0+ica4+2] );

				as1112 = config->U(is,1,icb,ica)*
                  ( a[isp1+ica4+1] - a[isp1+ica4+2] );

				ad1112 = std::conj(config->U(is1,1,ica,icb))*
                  ( a[ism1+ica4+1] + a[ism1+ica4+2] );

				as2123 = config->U(is,2,icb,ica)*
                  ( a[isp2+ica4+1] - II*a[isp2+ica4+3] );

				ad2123 = std::conj(config->U(is2,2,ica,icb))*
                  ( a[ism2+ica4+1] + II*a[ism2+ica4+3] );

				as3133 = config->U(is,3,icb,ica)*
                  ( a[isp3+ica4+1] - a[isp3+ica4+3] );

				ad3133 = std::conj(config->U(is3,3,ica,icb))*
                  ( a[ism3+ica4+1] + a[ism3+ica4+3] );
				

        (*outv)[is12+icb4+0] -= 0.5*(  as0003 + ad0003 + as1013 + ad1013
                   								+ as2022 + ad2022 + as3032 + ad3032 );

        (*outv)[is12+icb4+1] -= 0.5*(  as0102 + ad0102 + as1112 + ad1112
                   								+ as2123 + ad2123 + as3133 + ad3133 );

        (*outv)[is12+icb4+2] -= 0.5*(- II*as0102 + II*ad0102 - as1112 + ad1112
                   								- II*as2022 + II*ad2022 - as3032 + ad3032 );

        (*outv)[is12+icb4+3] -= 0.5*(- II*as0003 + II*ad0003 + as1013 - ad1013
                   								+ II*as2123 - II*ad2123 - as3133 + ad3133 );
  		}	// END color a
  	} // END color b
	} // END nsite

}

//----------------------------------------------------------------------------------------------------------
void gamma5( cvector *v )
{
	int is;
	for (is=0; is<12*nsite; is+=12)
	{
		(*v)[is+4*0+2] *= mONE;
		(*v)[is+4*1+2] *= mONE;
		(*v)[is+4*2+2] *= mONE;

		(*v)[is+4*0+3] *= mONE;
		(*v)[is+4*1+3] *= mONE;
		(*v)[is+4*2+3] *= mONE;
	}
}

//----------------------------------------------------------------------------------------------------------
#define MAXITER 1000
#define PREC2 1e-20
bool vec_invDwvec( cvector *inv, cvector *outv )
{
/*
		Conjugate Gradient
		computes:  outv = A^{-1}inv with precision sqrt(PREC2)
*/

	dcomplex alpha,beta, one=ONE;
	double res2, resnew2, inv2;
	int iter;
	int size = int((*inv).size( )), inc=1;
	
	// r = d = a = inv
	cvector a((*inv).size()), r((*inv).size()), d((*inv).size()), tmp((*inv).size());
	tmp = (*inv);

	gamma5( &tmp );

	r = tmp;
	d = tmp;
	inv2 = pow( norm(tmp), 2 );
	res2 = inv2;

	if ( inv2==0.)
	{
		logfile << "CG: entry is ZERO." << endl;
		(*outv).assign( (*outv).size( ), ZERO );
		return false;		
	}

	if ( res2<PREC2 )
	{
		logfile << "CG.ERROR: |inv|^2 too small: " << res2 << endl;
		MPI::COMM_WORLD.Abort(2);
	}

	(*outv).assign( (*outv).size( ), ZERO );

	iter = 0;
	do
	{
		iter++;

		// alpha = r.r / d.(Ad)
		vec_Dwvec( &d, &a ); 
		gamma5( &a ); // a = Ad
		alpha = complex<double>(res2 / std::real(dot( a, d )), 0.0);

		// outv = outv + alpha*d
		zaxpy_( &size, &alpha, d.data(), &inc, (*outv).data(), &inc );

		// r = r - alpha*(Ad)
		alpha *= mONE;
		zaxpy_(&size, &alpha, a.data(), &inc, r.data(), &inc ); // r = r - alpha*a

		// new residuo
		resnew2 = pow( norm(r), 2 );

		if ( resnew2<=PREC2*inv2 ) break;

		// beta = a.a/r.r   // a is r_new
		beta = complex<double>(resnew2/res2,0);

		// d = r + beta*d
		zscal_(&size, &beta, d.data(), &inc ); // d = beta*d	
		zaxpy_(&size, &one, r.data(), &inc, d.data(), &inc ); // d = 1*r + d

		res2 = resnew2;
	}
	while(iter<=MAXITER);

	// compute the true residual
	vec_Dwvec( outv, &a ); 
	one = mONE;
	zaxpy_( &size, &one, (*inv).data(), &inc, a.data(), &inc );
	res2 = norm(a);

	if ( iter>=MAXITER )
	{
		logfile << "CG.ERROR: iter > MAXITER. True residual " << resnew2 << " " << res2<< endl;
		return true;
	}

	logfile << "CG: finished in " << iter << " iterations. True residual " << resnew2 << " " << res2<< endl;
	return false;
}
