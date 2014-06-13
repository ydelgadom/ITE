//---------------------------------------------------------------------------------------------------
inline void dlnZdmu( dcomplex *finalobs, dcomplex *mean )
{
/*
	Compute the derivatives of ln Z with respect to
		m  : up to 4th order
		mu : up to 4th order
	using the mean values calculated in meanValues.

  0:P		   1:Dm	 		2:Dp 	 	 3:R
  4:R2		 5:Dm2		6:PR3		 7:RDm
  8:RDp		 9:PR	 	 10:P2	 	11:R3
 12:R2Dm 	13:R2Dp  14:RDm2
 15:PR2		16:P2R	 17:R3Dm
 18:R3Dp  19:R2Dm2 20:P2R2
*/

	// 4th order:
	// dlnZdmu = Dm + RDm - R*Dm + 1/2*R2Dm - 1/2*Dm*R2 - R*RDm + R^2*Dm
	//					+ 1/6 R3Dm - 1/2*R2*RDm - 1/6*Dm*R3 - 1/2*R*R2Dm + R*Dm*R2
	//					+ R^2*RDm - R^3*Dm
	//
	/*(*finalobs) = mean[1]*(1. - mean[3] - 0.5*mean[4] - (1./6.)*mean[11] ) 
			+ mean[3]*(mean[3]*mean[1] - mean[7] - 0.5*mean[12] + mean[1]*mean[4]) 
		  - 0.5*mean[4]*mean[7] + (1./6.)*mean[17] + mean[7] + 0.5*mean[12]
			+ mean[3]*mean[3]*(mean[7] - mean[3]*mean[1]);//*/

	(*finalobs) = mean[1] + mean[7] - mean[3]*mean[1];

	(*finalobs) /= double(nsite);
#ifdef MM2
	(*finalobs) *= double(Nt);
#endif

}

//---------------------------------------------------------------------------------------------------
inline void d2lnZdmu2( dcomplex *finalobs, dcomplex *mean )
{
/*
	Compute the derivatives of ln Z with respect to
		m  : up to 4th order
		mu : up to 4th order
	using the mean values calculated in meanValues.

  0:P		   1:Dm	 		2:Dp 	 	 3:R
  4:R2		 5:Dm2		6:PR3		 7:RDm
  8:RDp		 9:PR	 	 10:P2	 	11:R3
 12:R2Dm 	13:R2Dp  14:RDm2
 15:PR2		16:P2R	 17:R3Dm
 18:R3Dp  19:R2Dm2 20:P2R2
*/

	// 4th order:
	// d2lnZdmu2 = Dp + Dm2 + RDp - Dm^2 - R*Dp + RDm2 + 1/2*R2Dp - 1/2*Dp*R2
	//						- 2*Dm*RDm -  R*Dm2 - R*RDp + 2*R*Dm^2 + R^2*Dp + 1/2*R2Dm2 + 1/6*R3Dp
	//						+ (RDm)^2 + 1/2*R2*Dm2 - 1/2*R2*RDp - 1/6*Dp*R3 - Dm*R2Dm - R*RDm2
	//						- 1/2*R*R2Dp + Dm^2*R2 + R*Dp*R2 + 4*R*Dm*RDm
	//						+ R^2*Dm2 + R^2*RDp - 3*R^2*Dm^2 - R^3*Dp
	//
	/*(*finalobs) = mean[2] + mean[5] + mean[8]
			+ mean[1]*mean[1]*(- 1. + 2.*mean[3] + mean[4])
			+ mean[1]*mean[7]*(- 2. + 4.*mean[3])
			+ mean[3]*mean[3]*(mean[2] + mean[5] + mean[8] - 3.*mean[1]*mean[1] - mean[3]*mean[2])
			+ 0.5*mean[13] + mean[14] - mean[1]*mean[12] - (1./6.)*mean[2]*mean[11] 
			+ 0.5*mean[19] + (1./6.)*mean[18] - mean[7]*mean[7] 
			+ 0.5*mean[4]*(- mean[5] - mean[8] - mean[2])
			+ mean[3]*(mean[2]*mean[4] - mean[14] - 0.5*mean[13] - mean[5] - mean[8] - mean[2]);//*/

	(*finalobs) = mean[2] + mean[5] + mean[8] - mean[1]*mean[1] - mean[3]*mean[2];

	(*finalobs) /= double(nsite);
#ifdef MM2
	(*finalobs) *= double(Nt*Nt);
#endif

}

//---------------------------------------------------------------------------------------------------
inline void dlnZdm( dcomplex *finalobs, dcomplex *mean )
{
/*
	Compute the derivatives of ln Z with respect to
		m  : up to 4th order
		mu : up to 4th order
	using the mean values calculated in meanValues.

  0:P		   1:Dm	 		2:Dp 	 	 3:R
  4:R2		 5:Dm2		6:PR3		 7:RDm
  8:RDp		 9:PR	 	 10:P2	 	11:R3
 12:R2Dm 	13:R2Dp  14:RDm2
 15:PR2		16:P2R	 17:R3Dm
 18:R3Dp  19:R2Dm2 20:P2R2
*/

	// 4th order in propagators:
	// dlnZdm = P - P*R + PR - 1/2*P*R2 + 1/2*PR2 + R^2*P - R*PR
	//					+ R2*P*R - 1/2*R2*PR - 1/2*R*PR2
	//					- 1/6*P*R3 + 1/6*PR3 - R^3*P + R^2*PR
	//
	/*(*finalobs) = mean[0]*(1. - mean[3] - 0.5*mean[4] - (1./6.)*mean[11]) 
			+ mean[3]*(- mean[9] - 0.5*mean[15] + mean[0]*mean[4])
			+ mean[9]*(- 0.5*mean[4] + 1.)
			+ (1./6.)*mean[6] + 0.5*mean[15] 
		  + mean[3]*mean[3]*(- mean[3]*mean[0] + mean[9] + mean[0]);//*/

	(*finalobs) = mean[0] - mean[0]*mean[3] + mean[9];

	(*finalobs) /= double(nsite);

}


//----------------------------------------------------------------------------------------------------
inline void d2lnZdm2( dcomplex *finalobs, dcomplex *mean )
{
/*
	Compute the derivatives of ln Z with respect to
		m  : up to 4th order
		mu : up to 4th order
	using the mean values calculated in meanValues.

  0:P		   1:Dm	 		2:Dp 	 	 3:R
  4:R2		 5:Dm2		6:PR3		 7:RDm
  8:RDp		 9:PR	 	 10:P2	 	11:R3
 12:R2Dm 	13:R2Dp  14:RDm2
 15:PR2		16:P2R	 17:R3Dm
 18:R3Dp  19:R2Dm2 20:P2R2
*/

	// 4th order:
	// d2lnZdm2 = P2 - P^2 + 2*R*P^2 - R*P2 - 2*P*PR + P2R
	//						+ P^2*R2 - 1/2*P2*R2 - P*PR2 + 1/2*P2R2
	//						- 3*P^2*R^2 + 4*P*R*PR + R^2*P2 - PR^2 - R*P2R
	//	
	/*(*finalobs) = mean[10]*(1. + mean[3]*(-1. + mean[3]) - 0.5*mean[4])
		 	+ mean[0]*mean[0]*(-1. + mean[3]*(2. - 3.*mean[3]) + mean[4])
			+ mean[0]*( 2.*mean[9]*(- 1. + 2.*mean[3]) - mean[15])
			+ 0.5*mean[20] - mean[9]*mean[9] + mean[16]*(1. - mean[3]);//*/

	(*finalobs) = mean[10] - mean[0]*mean[0];

	(*finalobs) /= double(nsite);
}

//----------------------------------------------------------------------------------------------------
void partialObservables( double mu, dcomplex *pobs, dcomplex *traces )
{
/*
	TRACES
   0:tD			    1:tMD		       2:tMbD
	 3:tMDMD	    4:tMbDMbD	     5:tMDMbD
	 6:tMDD       7:tMbDD        8:tDD
	 9:tMDMDMD   10:tMDMDMbD    11:tMDMbDMbD   12:tMbDMbDMbD
	13:tMDMDD    14:tMDMbDD     15:tMbDMbDD    16:tMDDD         17:tMbDDD
	18:tMDMDMDMD 19:tMDMDMDMbD  20:tMDMDMbDMbD 21:tMDMbDMbDMbD  22:tMbDMbDMbDMbD
	23:tMDMDDD   24:tMDMbDDD	  25:tMbDMbDDD	 26:tMDMbDMDMbD
	27:tMDMbDD	 28:tMbDMDDD		29:tMDMDMDD		 30:tMbDMbDMbDD			
	31:tMDMDMbDD 32:tMDMbDMDD		33:tMbDMDMDD
	34:tMbDMbDMDD 35:tMbDMDMbDD	36:tMDMbDMbDD
*/

	dcomplex D1, MD1, MbD1;
	dcomplex MMD2, MbMbD2, MMbD2, MD2, MbD2, D2;
	dcomplex MMMD3, MMMbD3, MMbMbD3, MbMbMbD3, MMD3, MMbD3, MbMbD3, MD3, MbD3;
	dcomplex MMMMD4, MMMMbD4, MMMbMbD4, MMbMbMbD4, MbMbMbMbD4, MMD4, MMbD4, MbMbD4;
	dcomplex MMbMMbD4, MbMD3, MbMD4, MMMD4, MbMbMbD4;
	dcomplex MMMbD4, MMbMD4, MbMMD4, MMbMbD4, MbMMbD4, MbMbMD4;

  double emu = exp(mu);
  double emmu = exp(-mu);

#ifdef MM1
  double rho = emu - 1.;
  double rhobar = emmu - 1.;
#endif

#ifdef MM2
	emu = pow(emu,Nt);
	emmu = pow(emmu,Nt);
  double rho = emu - 1.;
  double rhobar = emmu - 1.;
#endif

	cout << mu << endl;
	//cout << "rho " << rho << " " << exp(mu)-1 << endl;
	//cout << "rhobar " << rhobar << " " << exp(-mu)-1 << endl;
	//cout << "emu " << emu << " " << exp(mu) << endl;
	//cout << "emmu " << emmu << " " << exp(-mu) << endl;

	double rho2 = rho*rho, rho3 = rho2*rho;
	double rhobar2 = rhobar*rhobar, rhobar3 = rhobar2*rhobar;
	double emu2 = emu*emu;
	double emmu2 = emmu*emmu;

	dcomplex *obs, *tr;

	for (int iconf=0; iconf<nconf; iconf++)
	{

		obs = &pobs[iconf*NMEAN];
		tr  = &traces[iconf*NTR];

		// D1 = tr[D]
		D1 = tr[0];

		// MD1 = tr[MD]
		MD1 = tr[1];

		// MbD1 = tr[MbD]
		MbD1 = tr[2];

		// MMD2 = tr[MD]^2 - tr[MDMD]
		MMD2 = MD1*MD1 - tr[3];

		// MbMbD2 = tr[MbD]^2 - tr[MbDMbD]
		MbMbD2 = MbD1*MbD1 - tr[4];

		// MMbD2 = tr[MD]*tr[MbD] - tr[MDMbD]
		MMbD2 = MD1*MbD1 - tr[5];

		// MD2 = tr[D]tr[MD] - tr[MDD]
		MD2 = D1*MD1 - tr[6];

		// MbD2 = tr[D]tr[MbD] - tr[MbDD]
		MbD2 = D1*MbD1 - tr[7];

		// D2 = tr[DD]
		D2 = tr[8];

		// MMMD3 = tr[MD]^3 - 3*Tr[MD]tr[MDMD] + 2*tr[MDMDMD]
		MMMD3 = MD1*MD1*MD1 - 3.*MD1*tr[3] + 2.*tr[9];

		// MMMbD3 = tr[MbD]*tr[MD]^2 - tr[MbD]*tr[MDMD] - 2*tr[MD]*tr[MDMbD]
		//					+ 2*tr[MDMDMbD]
		MMMbD3 = MbD1*MMD2 - 2.*MD1*tr[5] + 2.*tr[10];

		// MMbMbD3 = tr[MD]*tr[MbD]^2 - tr[MD]*tr[MbDMbD] - 2.*tr[MbD]*tr[MDMbD]
		//					+ 2*tr[MDMbDMbD]
		MMbMbD3 = MD1*MbMbD2 - 2.*MbD1*tr[5] + 2.*tr[11];

		// MbMbMbD3 = tr[MbD]^3 - 3*tr[MbD]*tr[MbDMbD] + 2*tr[MbDMbDMbD]
		MbMbMbD3 = MbD1*MbD1*MbD1 - 3.*MbD1*tr[4] + 2.*tr[12];

		// MMD3 = tr[D]*tr[MD]^2 - tr[D]*tr[MDMD] - 2*tr[MD]*tr[MDD] + 2*tr[MDMDD]
		MMD3 = D1*MMD2 - 2.*MD1*tr[6] + 2.*tr[13];

		// MMbD3 = tr[D]*tr[MD]*tr[MbD] - tr[D]*tr[MDMbD] - tr[MD]*tr[MbDD]
		//				- tr[MbD]*tr[MDD] + 2*tr[MDMbDD]
		MMbD3 = D1*MMbD2 - MD1*tr[7] - MbD1*tr[6] + 2.*tr[27];

		// MbMD3 = tr[D]*tr[MD]*tr[MbD] - tr[D]*tr[MDMbD] - tr[MD]*tr[MbDD]
		//				- tr[MbD]*tr[MDD] + 2*tr[MbDMDD]
		MbMD3 = D1*MMbD2 - MD1*tr[7] - MbD1*tr[6] + 2.*tr[14];

		// MbMbD3 = tr[D]*tr[MbD]^2 - tr[D]*tr[MbDMbD] - 2*tr[MbD]*tr[MbDD]
		//				+ 2*tr[MbDMbDD]
		MbMbD3 = D1*MbMbD2 - 2.*MbD1*tr[7] + 2.*tr[15];

		// MD3 = tr[MD]*tr[D]^2 - 2*tr[D]*tr[MDD] - tr[MD]*tr[DD] + 2*tr[MDDD]
		MD3 = MD1*D1*D1 - 2.*D1*tr[6] - MD1*D2 + 2.*tr[16];

		// MbD3 = tr[MbD]*tr[D]^2 - 2*tr[D]*tr[MbDD] - tr[MbD]*tr[DD] + 2*tr[MbDDD]
		MbD3 = MbD1*D1*D1 - 2.*D1*tr[7] - MbD1*D2 + 2.*tr[17];

		// MMMMD4 = tr[MD]^4 - 6*tr[MD]^2*tr[MDMD] + 3*tr[MDMD]^2
		//				+ 8*tr[MD]*tr[MDMDMD] - 6*tr[MDMDMDMD]
		MMMMD4 = MD1*MD1*MD1*MD1 - 6.*MD1*MD1*tr[3] + 3.*tr[3]*tr[3]
						+ 8.*MD1*tr[9] - 6.*tr[18];

		// MMMMbD4 = tr[MD]^3*tr[MbD] - 3*tr[MD]*tr[MbD]*tr[MDMD] - 3*tr[MD]^2*tr[MbDMD]
		//			+ 3*tr[MDMbD]*tr[MDMD] + 2*tr[MbD]*tr[MDMDMD] + 6*tr[MD]*tr[MDMDMbD] 
		//			-	6*tr[MDMDMDMbD]
		MMMMbD4 = MD1*MD1*MD1*MbD1 - 3.*MD1*MbD1*tr[3] - 3.*MD1*MD1*tr[5]
						+ 3.*tr[5]*tr[3] + 2.*MbD1*tr[9] + 6.*MD1*tr[10] - 6.*tr[19];

		// MMMbMbD4 = tr[MD]^2*tr[MbD]^2 - tr[MbD]^2*tr[MDMD] - 4*tr[MbD]*tr[MD]*tr[MDMbD]
		//			- tr[MD]^2*tr[MbDMbD] + tr[MbDMbD]*tr[MDMD] + 2*tr[MDMbD]^2
		//			+ 4*tr[MD]*tr[MDMbDMbD] + 4*tr[MbD]*tr[MDMDMbD] - 6*tr[MDMDMbDMbD]
		MMMbMbD4 = MD1*MD1*MbD1*MbD1 - MbD1*MbD1*tr[3] - 4.*MbD1*MD1*tr[5]
						- MD1*MD1*tr[4] + tr[4]*tr[3] + 2.*tr[5]*tr[5]
						+ 4.*MD1*tr[11] + 4.*MbD1*tr[10] - 6.*tr[20];

		// MMbMMbD4 = tr[MD]^2*tr[MbD]^2 - tr[MbD]^2*tr[MDMD] - 4*tr[MbD]*tr[MD]*tr[MDMbD]
		//			- tr[MD]^2*tr[MbDMbD] + tr[MbDMbD]*tr[MDMD] + 2*tr[MDMbD]^2
		//			+ 4*tr[MD]*tr[MDMbDMbD] + 4*tr[MbD]*tr[MDMDMbD] - 6*tr[MDMbDMDMbD]
		MMbMMbD4 = MD1*MD1*MbD1*MbD1 - MbD1*MbD1*tr[3] - 4.*MbD1*MD1*tr[5]
						- MD1*MD1*tr[4] + tr[4]*tr[3] + 2.*tr[5]*tr[5]
						+ 4.*MD1*tr[11] + 4.*MbD1*tr[10] - 6.*tr[26];

		// MMbMbMbD4 = tr[MD]*tr[MbD]^3 - 3*tr[MD]*tr[MbD]*tr[MbDMbD] - 3*tr[MbD]^2*tr[MDMbD]
		//			+ 3*tr[MDMbD]*tr[MbDMbD] + 2*tr[MD]*tr[MbDMbDMbD] + 6*tr[MbD]*tr[MDMbDMbD]
		//			- 6*tr[MDMbDMbDMbD]
		MMbMbMbD4 = MD1*MbD1*MbD1*MbD1 - 3.*MD1*MbD1*tr[4] - 3.*MbD1*MbD1*tr[5]
						+ 3.*tr[5]*tr[4] + 2.*MD1*tr[12] + 6.*MbD1*tr[11] - 6.*tr[21];

		// MbMbMbMbD4 = tr[MbD]^4 - 6*tr[MbD]^2*tr[MbDMbD] + 3*tr[MbDMbD]^2
		//				+ 8*tr[MbD]*tr[MbDMbDMbD] - 6*tr[MbDMbDMbDMbD]
		MbMbMbMbD4 = MbD1*MbD1*MbD1*MbD1 - 6.*MbD1*MbD1*tr[4] + 3.*tr[4]*tr[4]
						+ 8.*MbD1*tr[12] - 6.*tr[22];

		// MMD4 = tr[MD]^2*tr[D]^2 - tr[D]^2*tr[MDMD] - 4*tr[D]*tr[MD]*tr[MDD]
		//			- tr[MD]^2*tr[DD] + tr[DD]*tr[MDMD] + 2*tr[MDD]^2
		//			+ 4*tr[MD]*tr[MDDD] + 4*tr[D]*tr[MDMDD] - 6*tr[MDMDDD]
		MMD4 = MD1*MD1*D1*D1 - D1*D1*tr[3] - 4.*D1*MD1*tr[6]
						-	MD1*MD1*D2 + D2*tr[3] + 2.*tr[6]*tr[6]
						+ 4.*MD1*tr[16] + 4.*D1*tr[13] - 6.*tr[23];

		// MMbD4 = tr[D]^2*tr[MD]*tr[MbD] - tr[D]^2*tr[MDMbD] - 2*tr[D]*tr[MD]*tr[MbDD]
		//			- 2*tr[D]*tr[MbD]*tr[MDD] - tr[MD]*tr[MbD]*tr[DD] + tr[DD]*tr[MDMbD]
		//			+ 2*tr[MDD]*tr[MbDD] + 2*tr[MD]*tr[MbDDD] + 2*tr[MbD]*tr[MDDD]
		// 			+ 4*tr[D]*tr[MDMbDD] - 6*tr[MDMbDDD]
		MMbD4 = D1*D1*MD1*MbD1 - D1*D1*tr[5] - 2.*D1*MD1*tr[7]
						- 2.*D1*MbD1*tr[6] - MD1*MbD1*D2 + D2*tr[5] + 2.*tr[6]*tr[7]
						+ 2.*MD1*tr[17] + 2.*MbD1*tr[16] + 4.*D1*tr[14] - 6.*tr[24];

		// MbMD4 = tr[D]^2*tr[MD]*tr[MbD] - tr[D]^2*tr[MDMbD] - 2*tr[D]*tr[MD]*tr[MbDD]
		//			- 2*tr[D]*tr[MbD]*tr[MDD] - tr[MD]*tr[MbD]*tr[DD] + tr[DD]*tr[MDMbD]
		//			+ 2*tr[MDD]*tr[MbDD] + 2*tr[MD]*tr[MbDDD] + 2*tr[MbD]*tr[MDDD]
		// 			+ 4*tr[D]*tr[MDMbDD] - 6*tr[MbDMDDD]
		MbMD4 = D1*D1*MD1*MbD1 - D1*D1*tr[5] - 2.*D1*MD1*tr[7]
						- 2.*D1*MbD1*tr[6] - MD1*MbD1*D2 + D2*tr[5] + 2.*tr[6]*tr[7]
						+ 2.*MD1*tr[17] + 2.*MbD1*tr[16] + 4.*D1*tr[14] - 6.*tr[28];

		// MbMbD4 = tr[MbD]^2*tr[D]^2 - tr[D]^2*tr[MbDMbD] - 4*tr[D]*tr[MbD]*tr[MbDD]
		//			- tr[MbD]^2*tr[DD] + tr[DD]*tr[MbDMbD] + 2*tr[MbDD]^2
		//			+ 4*tr[MbD]*tr[MbDDD] + 4*tr[D]*tr[MbDMbDD] - 6*tr[MbDMbDDD]
		MbMbD4 = MbD1*MbD1*D1*D1 - D1*D1*tr[4] - 4.*D1*MbD1*tr[7]
						-	MbD1*MbD1*D2 + D2*tr[4] + 2.*tr[7]*tr[7]
						+ 4.*MbD1*tr[17] + 4.*D1*tr[15] - 6.*tr[25];

		// MMMD4 = tr[D]*tr[MD]^3 - 3*tr[D]*tr[MD]*tr[MDMD]
		//			- 3*tr[MD]^2*tr[MDD] + 3*tr[MDD]*tr[MDMD] + 2*tr[D]*tr[MDMDMD]
		//			+ 6*tr[MD]*tr[MDMDD] - 6*tr[MDMDMDD]
		MMMD4 = D1*MD1*MD1*MD1 - 3.*D1*MD1*tr[3] - 3.*MD1*MD1*tr[6]
						+ 3.*tr[6]*tr[3] + 2.*D1*tr[9] + 6.*MD1*tr[13] - 6.*tr[29];

		// MbMbMbD4 = tr[D]*tr[MbD]^3 - 3*tr[D]*tr[MbD]*tr[MbDMbD]
		//			- 3*tr[MbD]^2*tr[MbDD] + 3*tr[MbDD]*tr[MbDMbD] + 2*tr[D]*tr[MbDMbDMbD]
		//			+ 6*tr[MbD]*tr[MbDMbDD] - 6*tr[MbDMbDMbDD]
		MbMbMbD4 = D1*MbD1*MbD1*MbD1 - 3.*D1*MbD1*tr[4] - 3.*MbD1*MbD1*tr[7]
						+ 3.*tr[7]*tr[4] + 2.*D1*tr[12] + 6.*MbD1*tr[15] - 6.*tr[30];

		// MMMbD4 = tr[MD]^2*tr[MbD]*tr[D] - tr[D]*tr[MbD]*tr[MDMD] - 2*tr[D]*tr[MD]*tr[MDMbD]
		// 			- 2*tr[MbD]*tr[MD]*tr[MDD] - tr[MD]^2*tr[MbDD] + tr[MbDD]*tr[MDMD]
		//			+ 2*tr[MDD]*tr[MDMbD] + 4*tr[MD]*tr[MDMbDD] + 2*tr[D]*tr[MDMDMbD]
		//			+ 2*tr[MbD]*tr[MDMDD] - 6*tr[MDMDMbDD]
		MMMbD4 = MD1*MD1*MbD1*D1 - D1*MbD1*tr[3] - 2.*D1*MD1*tr[5]
						- 2.*MbD1*MD1*tr[6] - MD1*MD1*tr[7] + tr[7]*tr[3]
						+ 2.*tr[6]*tr[5] + 4.*MD1*tr[14] + 2.*D1*tr[10]
						+ 2.*MbD1*tr[13] - 6.*tr[31];

		// MMbMD4 = tr[MD]^2*tr[MbD]*tr[D] - tr[D]*tr[MbD]*tr[MDMD] - 2*tr[D]*tr[MD]*tr[MDMbD]
		// 			- 2*tr[MbD]*tr[MD]*tr[MDD] - tr[MD]^2*tr[MbDD] + tr[MbDD]*tr[MDMD]
		//			+ 2*tr[MDD]*tr[MDMbD] + 4*tr[MD]*tr[MDMbDD] + 2*tr[D]*tr[MDMDMbD]
		//			+ 2*tr[MbD]*tr[MDMDD] - 6*tr[MDMbDMDD]
		MMbMD4 = MD1*MD1*MbD1*D1 - D1*MbD1*tr[3] - 2.*D1*MD1*tr[5]
						- 2.*MbD1*MD1*tr[6] - MD1*MD1*tr[7] + tr[7]*tr[3]
						+ 2.*tr[6]*tr[5] + 4.*MD1*tr[14] + 2.*D1*tr[10]
						+ 2.*MbD1*tr[13] - 6.*tr[32];

		// MbMMD4 = tr[MD]^2*tr[MbD]*tr[D] - tr[D]*tr[MbD]*tr[MDMD] - 2*tr[D]*tr[MD]*tr[MDMbD]
		// 			- 2*tr[MbD]*tr[MD]*tr[MDD] - tr[MD]^2*tr[MbDD] + tr[MbDD]*tr[MDMD]
		//			+ 2*tr[MDD]*tr[MDMbD] + 4*tr[MD]*tr[MDMbDD] + 2*tr[D]*tr[MDMDMbD]
		//			+ 2*tr[MbD]*tr[MDMDD] - 6*tr[MbDMDMDD]
		MbMMD4 = MD1*MD1*MbD1*D1 - D1*MbD1*tr[3] - 2.*D1*MD1*tr[5]
						- 2.*MbD1*MD1*tr[6] - MD1*MD1*tr[7] + tr[7]*tr[3]
						+ 2.*tr[6]*tr[5] + 4.*MD1*tr[14] + 2.*D1*tr[10]
						+ 2.*MbD1*tr[13] - 6.*tr[33];

		// MbMbMD4 = tr[MbD]^2*tr[MD]*tr[D] - tr[D]*tr[MD]*tr[MbDMbD] - 2*tr[D]*tr[MbD]*tr[MbDMD]
		// 			- 2*tr[MD]*tr[MbD]*tr[MbDD] - tr[MbD]^2*tr[MDD] + tr[MDD]*tr[MbDMbD]
		//			+ 2*tr[MbDD]*tr[MDMbD] + 4*tr[MbD]*tr[MDMbDD] + 2*tr[D]*tr[MDMbDMbD]
		//			+ 2*tr[MD]*tr[MbDMbDD] - 6*tr[MbDMbDMDD]
		MbMbMD4 = MbD1*MbD1*MD1*D1 - D1*MD1*tr[4] - 2.*D1*MbD1*tr[5]
						- 2.*MD1*MbD1*tr[7] - MbD1*MbD1*tr[6] + tr[6]*tr[4]
						+ 2.*tr[7]*tr[5] + 4.*MbD1*tr[14] + 2.*D1*tr[11]
						+ 2.*MD1*tr[15] - 6.*tr[34];

		// MbMMbD4 = tr[MbD]^2*tr[MD]*tr[D] - tr[D]*tr[MD]*tr[MbDMbD] - 2*tr[D]*tr[MbD]*tr[MbDMD]
		// 			- 2*tr[MD]*tr[MbD]*tr[MbDD] - tr[MbD]^2*tr[MDD] + tr[MDD]*tr[MbDMbD]
		//			+ 2*tr[MbDD]*tr[MDMbD] + 4*tr[MbD]*tr[MDMbDD] + 2*tr[D]*tr[MDMbDMbD]
		//			+ 2*tr[MD]*tr[MbDMbDD] - 6*tr[MbDMDMbDD]
		MbMMbD4 = MbD1*MbD1*MD1*D1 - D1*MD1*tr[4] - 2.*D1*MbD1*tr[5]
						- 2.*MD1*MbD1*tr[7] - MbD1*MbD1*tr[6] + tr[6]*tr[4]
						+ 2.*tr[7]*tr[5] + 4.*MbD1*tr[14] + 2.*D1*tr[11]
						+ 2.*MD1*tr[15] - 6.*tr[35];

		// MbMbMD4 = tr[MbD]^2*tr[MD]*tr[D] - tr[D]*tr[MD]*tr[MbDMbD] - 2*tr[D]*tr[MbD]*tr[MbDMD]
		// 			- 2*tr[MD]*tr[MbD]*tr[MbDD] - tr[MbD]^2*tr[MDD] + tr[MDD]*tr[MbDMbD]
		//			+ 2*tr[MbDD]*tr[MDMbD] + 4*tr[MbD]*tr[MDMbDD] + 2*tr[D]*tr[MDMbDMbD]
		//			+ 2*tr[MD]*tr[MbDMbDD] - 6*tr[MbDMbDMDD]
		MMbMbD4 = MbD1*MbD1*MD1*D1 - D1*MD1*tr[4] - 2.*D1*MbD1*tr[5]
						- 2.*MD1*MbD1*tr[7] - MbD1*MbD1*tr[6] + tr[6]*tr[4]
						+ 2.*tr[7]*tr[5] + 4.*MbD1*tr[14] + 2.*D1*tr[11]
						+ 2.*MD1*tr[15] - 6.*tr[36];
		
		// P
	 	obs[0] = D1;
		cout << "P:" << obs[0] << endl;

		// Dm
		obs[1] = - emu*MD1 + emmu*MbD1; 
		cout << "Dm:" << obs[1] << endl;

		// Dp
		obs[2] = - emu*MD1 - emmu*MbD1; 
		cout << "Dp:" << obs[2] << endl;

		// R
		obs[3] = - rho*MD1 - rhobar*MbD1; 
		cout << "R:" << obs[3] << endl;

		// R2
		obs[4] = rho2*MMD2 + rhobar2*MbMbD2 + 2.*rho*rhobar*MMbD2;
		cout << "R2:" << obs[4] << endl;
	
		// Dm2
		obs[5] = emu2*MMD2 + emmu2*MbMbD2 - 2.*MMbD2;
		cout << "Dm2:" << obs[5] << endl;

		// PR3 = - rho^3*MMMD4 - rhobar^3*MbMbMbD4 
		//		  - rho^2*rhobar*(MMMbD4 + MMbMD4 + MbMMD4)
		//			- rho*rhobar^2*(MMbMbD4 + MbMMbD4 + MbMbMD4)
		obs[6] = - rho3*MMMD4 - rhobar3*MbMbMbD4
				- rho2*rhobar*(MMMbD4 + MMbMD4 + MbMMD4)
				- rho*rhobar2*(MMbMbD4 + MbMMbD4 + MbMbMD4);

		// RDm
		obs[7] = rho*emu*MMD2 - rhobar*emmu*MbMbD2 + (-rho*emmu + rhobar*emu)*MMbD2;
		cout << "RDm:" << obs[7] << endl;

		// RDp
		obs[8] = rho*emu*MMD2 + rhobar*emmu*MbMbD2 + ( rho*emmu + rhobar*emu)*MMbD2;
		cout << "RDp:" << obs[8] << endl;

		// PR
		obs[9] = -rho*MD2 - rhobar*MbD2;
		cout << "PR:" << obs[9] << endl;

		// P2
		obs[10] = (D1*D1 - D2);
		cout << "P2:" << obs[10] << endl;

		// R3 = - rho^3*MMMD3 - rhobar^3*MbMbMbD3 
		//		  - 3.*rho^2*rhobar*MMMbD3 - 3.*rho*rhobar^2*MMbMbD3
		obs[11] = - rho3*MMMD3 - 3.*rho2*rhobar*MMMbD3
				- 3.*rho*rhobar2*MMbMbD3 - rhobar3*MbMbMbD3;

		// R2Dm = - rho^2*emu*MMMD3 + emmu*rho^2*MMMbD3 - rhobar^2*emu*MMbMbD3
		// 			  + emmu*rhobar^2*MbMbMbD3 - 2.*rho*rhobar*emu*MMMbD3 + 2.*emmu*rhobar*rho*MMbMbD3
		obs[12] = - rho2*emu*MMMD3 + (emmu*rho2 - 2.*rho*rhobar*emu)*MMMbD3
				+ (2.*emmu*rhobar*rho - rhobar2*emu)*MMbMbD3 + emmu*rhobar2*MbMbMbD3;
	
		// R2Dp = - rho^2*emu*MMMD3 - emmu*rho^2*MMMbD3 - rhobar^2*emu*MMbMbD3
		//				- emmu*rhobar^2*MbMbMbD3 - 2.*rho*rhobar*emu*MMMbD3 - 2.*emmu*rho*rhobar*MMbMbD3
		obs[13] = - rho2*emu*MMMD3 + (-emmu*rho2 - 2.*rho*rhobar*emu)*MMMbD3
				+ (-rhobar2*emu - 2.*emmu*rho*rhobar)*MMbMbD3 - emmu*rhobar2*MbMbMbD3;

		// RDm2 = - rho*emu^2*MMMD3 - rho*emmu^2*MMbMbD3 + 2.*rho*MMMbD3
		//				- emu^2*rhobar*MMMbD3 - emmu^2*rhobar*MbMbMbD3 + 2.*rhobar*MMbMbD3
		obs[14] = - rho*emu2*MMMD3 + (2.*rho - emu2*rhobar)*MMMbD3
							+ (2.*rhobar - rho*emmu2)*MMbMbD3 - emmu2*rhobar*MbMbMbD3;

		// PR2
		obs[15] = rho2*MMD3 + rho*rhobar*(MMbD3+MbMD3) + rhobar2*MbMbD3;

		// P2R
		obs[16] = -rho*MD3 - rhobar*MbD3;

		// R3Dm = rho^3*emu*MMMMD4 + rhobar^3*emu*MMbMbMbD4 
		//			 + 3.*rhobar*rho^2*emu*MMMMbD4 - 3.*rho*rhobar^2*emmu*MMbMbMbD4
		//			 - rho^3*emmu*MMMMbD4 - rhobar^3*emmu*MbMbMbMbD4
		//			 - 2.*rhobar*rho^2*emmu*MMMbMbD4 + 2.*rho*rhobar^2*emu*MMMbMbD4
		//			 -    rhobar*rho^2*emmu*MMbMMbD4 +    rho*rhobar^2*emu*MMbMMbD4
		obs[17] = rho3*emu*MMMMD4
				+ (3.*rhobar*rho2*emu - rho3*emmu)*MMMMbD4 
				+ (2.*rho*rhobar2*emu - 2.*rhobar*rho2*emmu)*MMMbMbD4
				+ (   rho*rhobar2*emu -    rhobar*rho2*emmu)*MMbMMbD4
	 			+ (rhobar3*emu - 3.*rho*rhobar2*emmu)*MMbMbMbD4
				- rhobar3*emmu*MbMbMbMbD4;

		// R3Dp = rho^3*emu*MMMMD4 + rhobar^3*emu*MMbMbMbD4 
		//			 + 3.*rhobar*rho^2*emu*MMMMbD4 + 3.*rho*rhobar^2*emmu*MMbMbMbD4
		//			 + rho^3*emmu*MMMMbD4 + rhobar^3*emmu*MbMbMbMbD4
		//			 + 2.*rhobar*rho^2*emmu*MMMbMbD4 + 2.*rho*rhobar^2*emu*MMMbMbD4
		//			 +    rhobar*rho^2*emmu*MMbMMbD4 +    rho*rhobar^2*emu*MMbMMbD4
		obs[18] = rho3*emu*MMMMD4
				+ (3.*rhobar*rho2*emu + rho3*emmu)*MMMMbD4
				+ (2.*rho*rhobar2*emu + 2.*rhobar*rho2*emmu)*MMMbMbD4
				+ (   rho*rhobar2*emu +    rhobar*rho2*emmu)*MMbMMbD4
	 			+ (rhobar3*emu + 3.*rho*rhobar2*emmu)*MMbMbMbD4
				+ rhobar3*emmu*MbMbMbMbD4;

		// R2Dm2 = rho^2*emu^2*MMMMD4 + rho^2*emmu^2*MMMbMbD4 
		//				- 2.*rho^2*MMMMbD4 + rhobar^2*emu^2*MMMbMbD4
		//				+ rhobar^2*emmu^2*MbMbMbMbD4 - 2.*rhobar^2*MMbMbMbD4 
		//				+ 2.*rho*rhobar*emu^2*MMMMbD4 + 2.*rho*rhobar*emmu^2*MMbMbMbD4 
		//				- 2.*rho*rhobar*MMMbMbD4 - 2.*rho*rhobar*MMbMMbD4
		obs[19] = rho2*emu2*MMMMD4
				+ (2.*rho*rhobar*emu2 - 2.*rho2)*MMMMbD4
				+ (rho2*emmu2 + rhobar2*emu2 - 2.*rho*rhobar)*MMMbMbD4
				+ (2.*rho*rhobar*emmu2 - 2.*rhobar2)*MMbMbMbD4
				+ rhobar2*emmu2*MbMbMbMbD4
				- 2.*rho*rhobar*MMbMMbD4;

		// P2R2 = rho^2*MMD4 + rhobar^2*MbMbD4 + rho*rhobar*(MMbD4+MbMD4)
		obs[20] = rho2*MMD4 + rhobar2*MbMbD4 + rho*rhobar*(MMbD4+MbMD4);
	}

}
