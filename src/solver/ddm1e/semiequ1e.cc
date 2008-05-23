/*****************************************************************************/
/*              8888888         88888888         88888888                    */
/*            8                8                8                            */
/*           8                 8                8                            */
/*           8                  88888888         88888888                    */
/*           8      8888                8                8                   */
/*            8       8                 8                8                   */
/*              888888         888888888        888888888                    */
/*                                                                           */
/*       A Two-Dimensional General Purpose Semiconductor Simulator.          */
/*                                                                           */
/*  GSS 0.4x                                                                 */
/*  Last update: Nov 16, 2007                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "zonedata.h"
#include "jflux1e.h"
#include "vec3.h"

#define _FLUX1_
#define _HMOB_A_
#define _HMOB_B_
#define _HMOB_C_
#define _II_
#define _EXP_R_
#define _FERMI_
#define _HURKX_BTBT_

void rotate3Scalar( PetscScalar& a, PetscScalar& b, PetscScalar& c) 
{
  PetscScalar t;
  t = a;
  a = b;
  b = c;
  c = t;
}

void SMCZone::F1E_Tri_ddm(Tri *ptri, PetscScalar *x, PetscScalar *f, vector<int> &zofs)
{
  int A, B, C;
  PetscScalar xA, xB, xC, yA, yB, yC;
  PetscScalar Sax, Say, Sbx, Sby, Scx, Scy;
  PetscScalar La, Lb, Lc;
  PetscScalar da, db, dc;
  PetscScalar cosA, cosB, cosC;
  PetscScalar sinA, sinB, sinC;
  PetscScalar tri_area, e, T, kb, Vt, Eg;
  PetscScalar VA, VB, VC, nA, nB, nC, pA, pB, pC; // potential, electron and hole densities at nodes
  PetscScalar EcA, EcB, EcC, EvA, EvB, EvC; // conduction band energy level
  PetscScalar RA, RB, RC; // Recombination rates
  PetscScalar niA, niB, niC;
  PetscScalar Ex, Ey, E;
  PetscScalar phinA, phinB, phinC; // quasi-fermi level for electrons
  PetscScalar phipA, phipB, phipC; // quasi-fermi level for holes
  PetscScalar Jna, Jnb, Jnc, Jpa, Jpb, Jpc;
  PetscScalar Jn_scale, Jp_scale;

  for (int pEdge = 0; pEdge < 3; pEdge++)
  {
    if (pEdge == 0) {
      A = ptri->node[0]; B = ptri->node[1]; C = ptri->node[2];
    } else {
      int t; t=A; A=B; B=C; C=t;
    }

    if (pEdge == 0) {
      xA = pzone->danode[ptri->node[0]].x;
      xB = pzone->danode[ptri->node[1]].x;
      xC = pzone->danode[ptri->node[2]].x;
      yA = pzone->danode[ptri->node[0]].y;
      yB = pzone->danode[ptri->node[1]].y;
      yC = pzone->danode[ptri->node[2]].y;
      Sax = (xC-xB)/ptri->edge_len[0];
      Say = (yC-yB)/ptri->edge_len[0];
      Sbx = (xA-xC)/ptri->edge_len[1];
      Sby = (yA-yC)/ptri->edge_len[1];
      Scx = (xB-xA)/ptri->edge_len[2];
      Scy = (yB-yA)/ptri->edge_len[2];
      La = ptri->edge_len[0];
      Lb = ptri->edge_len[1];
      Lc = ptri->edge_len[2];
      da = ptri->d[0];
      db = ptri->d[1];
      dc = ptri->d[2];
      cosA = cos(ptri->angle[0]);
      cosB = cos(ptri->angle[1]);
      cosC = cos(ptri->angle[2]);
      sinA = sin(ptri->angle[0]);
      sinB = sin(ptri->angle[1]);
      sinC = sin(ptri->angle[2]);
    } else {
      rotate3Scalar(xA, xB, xC);
      rotate3Scalar(yA, yB, yC);
      rotate3Scalar(Sax,Sbx,Scx);
      rotate3Scalar(Say,Sby,Scy);
      rotate3Scalar(La,Lb,Lc);
      rotate3Scalar(da,db,dc);
      rotate3Scalar(cosA,cosB,cosC);
      rotate3Scalar(sinA,sinB,sinC);
    }

    if (pEdge == 0) {
      tri_area = ptri->area;
      e  =  mt->e;
      T  = (fs[A].T+fs[B].T+fs[C].T)/3.0;
      kb = mt->kb;
      Vt = mt->kb*T/e;
      Eg = mt->band->Eg(T);
    } else {
    }

    if (pEdge == 0) {
      VA = x[zofs[zone_index]+3*A+0];     //potential of node A
      nA = x[zofs[zone_index]+3*A+1];     //electron density of node A
      pA = x[zofs[zone_index]+3*A+2];     //hole density of node A
      VB = x[zofs[zone_index]+3*B+0];     //potential of node B
      nB = x[zofs[zone_index]+3*B+1];     //electron density of node B
      pB = x[zofs[zone_index]+3*B+2];     //hole density of node B
      VC = x[zofs[zone_index]+3*C+0];     //potential of node C
      nC = x[zofs[zone_index]+3*C+1];     //electron density of node C
      pC = x[zofs[zone_index]+3*C+2];     //hole density of node C
    } else {
      rotate3Scalar(VA,VB,VC);
      rotate3Scalar(nA,nB,nC);
      rotate3Scalar(pA,pB,pC);
    }

    if (pEdge==0) {
      mt->mapping(&pzone->danode[A],&aux[A],0);
      EcA = -(e*VA + aux[A].affinity + mt->band->EgNarrowToEc(T));//conduction band energy level
      EvA = -(e*VA + aux[A].affinity - mt->band->EgNarrowToEv(T) + Eg);//valence band energy level
      RA = mt->band->Recomb(pA,nA,fs[A].T);
      niA = mt->band->nie(T);

      mt->mapping(&pzone->danode[B],&aux[B],0);
      EcB = -(e*VB + aux[B].affinity + mt->band->EgNarrowToEc(T));//conduction band energy level
      EvB = -(e*VB + aux[B].affinity - mt->band->EgNarrowToEv(T) + Eg);//valence band energy level
      RB = mt->band->Recomb(pB,nB,fs[B].T);
      niB = mt->band->nie(T);

      mt->mapping(&pzone->danode[C],&aux[C],0);
      EcC = -(e*VC + aux[C].affinity + mt->band->EgNarrowToEc(T));//conduction band energy level
      EvC = -(e*VC + aux[C].affinity - mt->band->EgNarrowToEv(T) + Eg);//valence band energy level
      RC = mt->band->Recomb(pC,nC,fs[C].T);
      niC = mt->band->nie(T);
#ifdef _FERMI_  
      if(Fermi)
      {
        PetscScalar NcA = mt->band->Nc(fs[A].T);
        PetscScalar NvA = mt->band->Nv(fs[A].T);
        EcA = EcA - e*Vt*log(gamma_f(fabs(nA)/NcA));
        EvA = EvA + e*Vt*log(gamma_f(fabs(pA)/NvA));

        PetscScalar NcB = mt->band->Nc(fs[B].T);
        PetscScalar NvB = mt->band->Nv(fs[B].T);
        EcB = EcB - e*Vt*log(gamma_f(fabs(nB)/NcB));
        EvB = EvB + e*Vt*log(gamma_f(fabs(pB)/NvB));

        PetscScalar NcC = mt->band->Nc(fs[C].T);
        PetscScalar NvC = mt->band->Nv(fs[C].T);
        EcC = EcC - e*Vt*log(gamma_f(fabs(nC)/NcC));
        EvC = EvC + e*Vt*log(gamma_f(fabs(pC)/NvC));
      }
#endif 
    } else {
      rotate3Scalar(EcA,EcB,EcC);
      rotate3Scalar(EvA,EvB,EvC);
      rotate3Scalar(RA,RB,RC);
      rotate3Scalar(niA,niB,niC);
    } 

    if(pEdge==0) {
      Ex = -((yB-yC)*VA + (yC-yA)*VB +(yA-yB)*VC)/(2*tri_area);
      Ey = -((xC-xB)*VA + (xA-xC)*VB +(xB-xA)*VC)/(2*tri_area);
      E = sqrt(Ex*Ex+Ey*Ey+1e-20);
    }
    if(pEdge==0) {
      phinA = VA - log(fabs(nA)/niA)*Vt;
      phipA = VA + log(fabs(pA)/niA)*Vt;
      phinB = VB - log(fabs(nB)/niB)*Vt;
      phipB = VB + log(fabs(pB)/niB)*Vt;
      phinC = VC - log(fabs(nC)/niC)*Vt;
      phipC = VC + log(fabs(pC)/niC)*Vt;
    } else {
      rotate3Scalar(phinA,phinB,phinC);
      rotate3Scalar(phipA,phipB,phipC);
    }

#ifdef _II_
    PetscScalar IIn, IIp;
    if(ImpactIonization) {
      if(IIType==GradQf && pEdge==0)
      {
        PetscScalar Fnx = ((yB-yC)*phinA + (yC-yA)*phinB +(yA-yB)*phinC)/(2*tri_area);
        PetscScalar Fny = ((xC-xB)*phinA + (xA-xC)*phinB +(xB-xA)*phinC)/(2*tri_area);
        PetscScalar Fpx = ((yB-yC)*phipA + (yC-yA)*phipB +(yA-yB)*phipC)/(2*tri_area);
        PetscScalar Fpy = ((xC-xB)*phipA + (xA-xC)*phipB +(xB-xA)*phipC)/(2*tri_area);
        PetscScalar Fn  =  sqrt(Fnx*Fnx+Fny*Fny);
        PetscScalar Fp  =  sqrt(Fpx*Fpx+Fpy*Fpy);
        IIn =  mt->gen->ElecGenRate(T,Fn,Eg);
        IIp =  mt->gen->HoleGenRate(T,Fp,Eg);
      }
      if(IIType==EVector && pEdge==0)
      {
        PetscScalar Fnx = Ex;
        PetscScalar Fny = Ey;
        PetscScalar Fpx = Ex;
        PetscScalar Fpy = Ey;
        PetscScalar Fn  =  sqrt(Fnx*Fnx+Fny*Fny);
        PetscScalar Fp  =  sqrt(Fpx*Fpx+Fpy*Fpy);
        IIn =  mt->gen->ElecGenRate(T,Fn,Eg);
        IIp =  mt->gen->HoleGenRate(T,Fp,Eg);
      }
    }
#endif

#ifdef _FLUX1_
    if(pEdge==0) {
      Jnc =  In(Vt,EcA/e,EcB/e,nA,nB,Lc);
      Jpc =  Ip(Vt,EvA/e,EvB/e,pA,pB,Lc);
      Jna =  In(Vt,EcB/e,EcC/e,nB,nC,La);
      Jpa =  Ip(Vt,EvB/e,EvC/e,pB,pC,La);
      Jnb =  In(Vt,EcC/e,EcA/e,nC,nA,Lb);
      Jpb =  Ip(Vt,EvC/e,EvA/e,pC,pA,Lb);
      Jn_scale = dmax(dmax(fabs(Jna),fabs(Jnb)),dmax(fabs(Jnc),niA*niA));
      Jp_scale = dmax(dmax(fabs(Jpa),fabs(Jpb)),dmax(fabs(Jpc),niA*niA));
      Jnc /=  Jn_scale;
      Jpc /=  Jp_scale;
      Jna /=  Jn_scale;
      Jpa /=  Jp_scale;
      Jnb /=  Jn_scale;
      Jpb /=  Jp_scale;
    } else {
      rotate3Scalar(Jna,Jnb,Jnc);
      rotate3Scalar(Jpa,Jpb,Jpc);
    }
#endif
#ifdef _FLUX2_
    if(pEdge==0) {
      Jnc =  In(Vt,(EcA-EcB)/e,nA,nB,Lc);
      Jpc =  Ip(Vt,(EvA-EvB)/e,pA,pB,Lc);
      Jna =  In(Vt,(EcB-EcC)/e,nB,nC,La);
      Jpa =  Ip(Vt,(EvB-EvC)/e,pC,La);
      Jnb =  In(Vt,(EcC-EcA)/e,nC,nA,Lb);
      Jpb =  Ip(Vt,(EvC-EvA)/e,pC,pA,Lb);
      Jn_scale = dmax(dmax(fabs(Jna),fabs(Jnb)),dmax(fabs(Jnc),niA*niA));
      Jp_scale = dmax(dmax(fabs(Jpa),fabs(Jpb)),dmax(fabs(Jpc),niA*niA));
      Jnc /=  Jn_scale;
      Jpc /=  Jp_scale;
      Jna /=  Jn_scale;
      Jpa /=  Jp_scale;
      Jnb /=  Jn_scale;
      Jpb /=  Jp_scale;
    } else {
      rotate3Scalar(Jna,Jnb,Jnc);
      rotate3Scalar(Jpa,Jpb,Jpc);
    }
#endif
    //flux along A-B
    PetscScalar Epnc, Etnc, Eppc, Etpc;
    PetscScalar Jnc_norm, Jpc_norm;
    PetscScalar mun, mup;
    if(HighFieldMobility)
    {
      if(EJModel || IIType==EdotJ)
      {

        PetscScalar JncxTca  = (Jnc*Say - Jna*Scy)/sinB;
        PetscScalar JncxTbc  = (Jnb*Scy - Jnc*Sby)/sinA;
        PetscScalar Jncx = (da * JncxTca + db * JncxTbc)/(da+db);

        PetscScalar JncyTca  = (-Jnc*Sax + Jna*Scx)/sinB;
        PetscScalar JncyTbc  = (-Jnb*Scx + Jnc*Sbx)/sinA;
        PetscScalar Jncy = (da * JncyTca + db * JncyTbc)/(da+db);

        PetscScalar JpcxTca  = (Jpc*Say - Jpa*Scy)/sinB;
        PetscScalar JpcxTbc  = (Jpb*Scy - Jpc*Sby)/sinA;
        PetscScalar Jpcx = (da * JpcxTca + db * JpcxTbc)/(da+db);

        PetscScalar JpcyTca  = (-Jpc*Sax + Jpa*Scx)/sinB;
        PetscScalar JpcyTbc  = (-Jpb*Scx + Jpc*Sbx)/sinA;
        PetscScalar Jpcy = (da * JpcyTca + db * JpcyTbc)/(da+db);

        Jnc_norm = sqrt(Jncx*Jncx + Jncy*Jncy + 1e-100);
        Jpc_norm = sqrt(Jpcx*Jpcx + Jpcy*Jpcy + 1e-100);
        Epnc = (Ex*Jncx + Ey*Jncy)/Jnc_norm;
        Etnc = (Ex*Jncy - Ey*Jncx)/Jnc_norm;
        Eppc = (Ex*Jpcx + Ey*Jpcy)/Jpc_norm;
        Etpc = (Ex*Jpcy - Ey*Jpcx)/Jpc_norm;
      }
      else
      {
        Epnc = fabs(EcA-EcB)/e/Lc; //parallel electrical field for electron
        Eppc = fabs(EvA-EvB)/e/Lc; //parallel electrical field for hole
        //transvers electrical field for electron and hole
        Etnc = fabs(EcA + Lb*cosA*(EcB-EcA)/Lc - EcC)/(Lb*sinA)/e;
        Etpc = fabs(EvA + Lb*cosA*(EvB-EvA)/Lc - EvC)/(Lb*sinA)/e;
      }

      mt->mapping(&pzone->danode[A],&aux[A],0);
      aux[A].mun =  mt->mob->ElecMob(pA,nA,T,dmax(0,Epnc),fabs(Etnc),T);
      aux[A].mup =  mt->mob->HoleMob(pA,nA,T,dmax(0,Eppc),fabs(Etpc),T);

      mt->mapping(&pzone->danode[B],&aux[B],0);
      aux[B].mun =  mt->mob->ElecMob(pB,nB,T,dmax(0,Epnc),fabs(Etnc),T);
      aux[B].mup =  mt->mob->HoleMob(pB,nB,T,dmax(0,Eppc),fabs(Etpc),T);
    }  
    mun = 0.5*(aux[A].mun+aux[B].mun);
    mup = 0.5*(aux[A].mup+aux[B].mup);

    PetscScalar G=0;
#ifdef _II_
    if(ImpactIonization && IIType==EdotJ)
    {
      IIn =  mt->gen->ElecGenRate(T,dmax(0,Epnc),Eg);
      IIp =  mt->gen->HoleGenRate(T,dmax(0,Eppc),Eg);
      G   =  IIn*mun*Jnc_norm*Jn_scale + IIp*mup*Jpc_norm*Jp_scale;
    }
    if( ImpactIonization && (IIType==GradQf || IIType==EVector) )
    {
      G = IIn*mun*fabs(phinA-phinB)/Lc*nmid(Vt,fs[A].P,fs[B].P,fs[A].n,fs[B].n)
        +IIp*mup*fabs(phipA-phipB)/Lc*pmid(Vt,fs[A].P,fs[B].P,fs[A].p,fs[B].p);
      G = IIn*mun*fabs(phinA-phinB)/Lc*nmid(Vt,VA,VB,nA,nB)
        +IIp*mup*fabs(phipA-phipB)/Lc*pmid(Vt,VA,VB,pA,pB);
    }
    if( ImpactIonization && IIType==ESide )
    {
      IIn =  mt->gen->ElecGenRate(T,fabs(VA-VB)/Lc,Eg);
      IIp =  mt->gen->HoleGenRate(T,fabs(VA-VB)/Lc,Eg);
      G = IIn*mun*fabs(Jnc)*Jn_scale + IIp*mup*fabs(Jpc)*Jp_scale;
    }
#endif

    f[zofs[zone_index]+3*A+0] +=   0.5*(aux[A].eps+aux[B].eps)*(VB-VA)/Lc*dc;
    f[zofs[zone_index]+3*A+1] +=   mun*Jnc*Jn_scale*dc + G*0.25*dc*Lc;
    f[zofs[zone_index]+3*A+2] += - mup*Jpc*Jp_scale*dc + G*0.25*dc*Lc;

    f[zofs[zone_index]+3*B+0] += - 0.5*(aux[A].eps+aux[B].eps)*(VB-VA)/Lc*dc;
    f[zofs[zone_index]+3*B+1] += - mun*Jnc*Jn_scale*dc + G*0.25*dc*Lc;
    f[zofs[zone_index]+3*B+2] +=   mup*Jpc*Jp_scale*dc + G*0.25*dc*Lc;

//Generation at A
#ifdef _HURKX_BTBT_
    if (BandBandTunneling)
    {
      PetscScalar S = 0.25*dc*Lc + 0.25*db*Lb;

      PetscScalar GradQfnx = -((yB-yC)*phinA + (yC-yA)*phinB +(yA-yB)*phinC)/(2*tri_area);
      PetscScalar GradQfny = -((xC-xB)*phinA + (xA-xC)*phinB +(xB-xA)*phinC)/(2*tri_area);
      PetscScalar GradQfn = sqrt(GradQfnx*GradQfnx+GradQfny*GradQfny+1e-20);

      PetscScalar GradQfpx = -((yB-yC)*phipA + (yC-yA)*phipB +(yA-yB)*phipC)/(2*tri_area);
      PetscScalar GradQfpy = -((xC-xB)*phipA + (xA-xC)*phipB +(xB-xA)*phipC)/(2*tri_area);
      PetscScalar GradQfp = sqrt(GradQfpx*GradQfpx+GradQfpy*GradQfpy+1e-20);

      PetscScalar nn = nA * pow(niA / mt->band->Nc(fs[A].T), GradQfn/E );
      PetscScalar pp = pA * pow(niA / mt->band->Nv(fs[A].T), GradQfp/E );

      PetscScalar D = (niA*niA - nn*pp)/(niA+nn)/(niA+pp);
      PetscScalar D1 = (GradQfnx*GradQfpx + GradQfny*GradQfpy)/GradQfn/GradQfp;
      PetscScalar G_BB = D*mt->band->BB_Tunneling(T,E);

      f[zofs[zone_index]+3*A+1] += G_BB*S;
      f[zofs[zone_index]+3*A+2] += G_BB*S;
    }
#endif
 
    if(pEdge==0) {
      //optical carrier generation
      f[zofs[zone_index]+3*A+1] +=  aux[A].RealOptG*ptri->s[0];
      f[zofs[zone_index]+3*A+2] +=  aux[A].RealOptG*ptri->s[0];
      f[zofs[zone_index]+3*B+1] +=  aux[B].RealOptG*ptri->s[1];
      f[zofs[zone_index]+3*B+2] +=  aux[B].RealOptG*ptri->s[1];
      f[zofs[zone_index]+3*C+1] +=  aux[C].RealOptG*ptri->s[2];
      f[zofs[zone_index]+3*C+2] +=  aux[C].RealOptG*ptri->s[2];

#ifndef _HURKX_BTBT_
      if(BandBandTunneling)
      {
        PetscScalar G_BB = mt->band->BB_Tunneling(T,E);
        f[zofs[zone_index]+3*A+1] +=  G_BB*ptri->s[0];
        f[zofs[zone_index]+3*A+2] +=  G_BB*ptri->s[0];
        f[zofs[zone_index]+3*B+1] +=  G_BB*ptri->s[1];
        f[zofs[zone_index]+3*B+2] +=  G_BB*ptri->s[1];
        f[zofs[zone_index]+3*C+1] +=  G_BB*ptri->s[2];
        f[zofs[zone_index]+3*C+2] +=  G_BB*ptri->s[2];
      }
#endif

      //---------------------------------------------------------------------------
      // Recombination item
      //---------------------------------------------------------------------------
      PetscScalar S;
      S = 0.25*dc*Lc + 0.25*db*Lb;
      f[zofs[zone_index]+3*A+1] += -RA*S;
      f[zofs[zone_index]+3*A+2] += -RA*S;

      S = 0.25*da*La + 0.25*dc*Lc;
      f[zofs[zone_index]+3*B+1] += -RB*S;
      f[zofs[zone_index]+3*B+2] += -RB*S;

      S = 0.25*db*Lb + 0.25*da*La;
      f[zofs[zone_index]+3*C+1] += -RC*S;
      f[zofs[zone_index]+3*C+2] += -RC*S;

    }

  }
}

void rotate3ADScalar( AutoDScalar& a, AutoDScalar& b, AutoDScalar& c) 
{
  AutoDScalar t=a;
  a = b;
  b = c;
  c = t;
}


void SMCZone::J1E_Tri_ddm(Tri *ptri,PetscScalar *x,Mat *jtmp, vector<int> & zofs)
{
  int A, B, C;
  PetscScalar xA, xB, xC, yA, yB, yC;
  PetscScalar Sax, Say, Sbx, Sby, Scx, Scy;
  PetscScalar La, Lb, Lc;
  PetscScalar da, db, dc;
  PetscScalar cosA, cosB, cosC;
  PetscScalar sinA, sinB, sinC;

  AutoDScalar VA, VB, VC, nA, nB, nC, pA, pB, pC; // potential, electron and hole densities at nodes
  AutoDScalar EcA, EcB, EcC, EvA, EvB, EvC; // conduction band energy level
  AutoDScalar RA, RB, RC; // Recombination rates
  PetscScalar niA, niB, niC;
  AutoDScalar Ex, Ey, E;
  AutoDScalar phinA, phinB, phinC; // quasi-fermi level for electrons
  AutoDScalar phipA, phipB, phipC; // quasi-fermi level for electrons
  AutoDScalar Jna, Jnb, Jnc, Jpa, Jpb, Jpc;
  PetscScalar Jn_scale, Jp_scale;
  AutoDScalar TD;

  Mat3       J1,J2,J3;
  Vec3       scale;
  PetscInt   index1[3],index2[3],index3[3];
  int offset1, offset2, offset3;

  for (int pEdge = 0; pEdge < 3; pEdge++)
  {
    if (pEdge == 0) {
      A = ptri->node[0]; B = ptri->node[1]; C = ptri->node[2];
      offset1=0; offset2=3; offset3=6;
    } else {
      int t; t=A; A=B; B=C; C=t;
      t=offset1; offset1=offset2; offset2=offset3; offset3=t;
    }


    if (pEdge == 0) {
      xA = pzone->danode[ptri->node[0]].x;
      xB = pzone->danode[ptri->node[1]].x;
      xC = pzone->danode[ptri->node[2]].x;
      yA = pzone->danode[ptri->node[0]].y;
      yB = pzone->danode[ptri->node[1]].y;
      yC = pzone->danode[ptri->node[2]].y;
      Sax = (xC-xB)/ptri->edge_len[0];
      Say = (yC-yB)/ptri->edge_len[0];
      Sbx = (xA-xC)/ptri->edge_len[1];
      Sby = (yA-yC)/ptri->edge_len[1];
      Scx = (xB-xA)/ptri->edge_len[2];
      Scy = (yB-yA)/ptri->edge_len[2];
      La = ptri->edge_len[0];
      Lb = ptri->edge_len[1];
      Lc = ptri->edge_len[2];
      da = ptri->d[0];
      db = ptri->d[1];
      dc = ptri->d[2];
      cosA = cos(ptri->angle[0]);
      cosB = cos(ptri->angle[1]);
      cosC = cos(ptri->angle[2]);
      sinA = sin(ptri->angle[0]);
      sinB = sin(ptri->angle[1]);
      sinC = sin(ptri->angle[2]);
    } else {
      rotate3Scalar(xA, xB, xC);
      rotate3Scalar(yA, yB, yC);
      rotate3Scalar(Sax,Sbx,Scx);
      rotate3Scalar(Say,Sby,Scy);
      rotate3Scalar(La,Lb,Lc);
      rotate3Scalar(da,db,dc);
      rotate3Scalar(cosA,cosB,cosC);
      rotate3Scalar(sinA,sinB,sinC);
    }

    PetscScalar tri_area, e, T, kb, Vt, Eg;
    if (pEdge == 0) {
      tri_area = ptri->area;
      e  =  mt->e;
      T  = (fs[A].T+fs[B].T+fs[C].T)/3.0;
      kb = mt->kb;
      Vt = mt->kb*T/e;
      Eg = mt->band->Eg(T);
    } else {
    }

    int z = zone_index;
    index1[0] = zofs[z]+3*A+0;
    index1[1] = zofs[z]+3*A+1;
    index1[2] = zofs[z]+3*A+2;
    index2[0] = zofs[z]+3*B+0;
    index2[1] = zofs[z]+3*B+1;
    index2[2] = zofs[z]+3*B+2;
    index3[0] = zofs[z]+3*C+0;
    index3[1] = zofs[z]+3*C+1;
    index3[2] = zofs[z]+3*C+2;

    if (pEdge==0) {
      //the indepedent variable number
      adtl::AutoDScalar::numdir=9;
      //synchronize with material database
      mt->set_ad_num(adtl::AutoDScalar::numdir); 
      TD = T; // dummy ad variable.

      VA = x[zofs[zone_index]+3*A+0];     //potential of node A
      nA = x[zofs[zone_index]+3*A+1];     //electron density of node A
      pA = x[zofs[zone_index]+3*A+2];     //hole density of node A
      VB = x[zofs[zone_index]+3*B+0];     //potential of node B
      nB = x[zofs[zone_index]+3*B+1];     //electron density of node B
      pB = x[zofs[zone_index]+3*B+2];     //hole density of node B
      VC = x[zofs[zone_index]+3*C+0];     //potential of node C
      nC = x[zofs[zone_index]+3*C+1];     //electron density of node C
      pC = x[zofs[zone_index]+3*C+2];     //hole density of node C
      VA.setADValue(0,1.0);
      nA.setADValue(1,1.0);
      pA.setADValue(2,1.0);
      VB.setADValue(3,1.0);
      nB.setADValue(4,1.0);
      pB.setADValue(5,1.0);
      VC.setADValue(6,1.0);
      nC.setADValue(7,1.0);
      pC.setADValue(8,1.0);
    }else{
      rotate3ADScalar(VA,VB,VC);
      rotate3ADScalar(nA,nB,nC);
      rotate3ADScalar(pA,pB,pC);
    }

    if (pEdge==0) {
      mt->mapping(&pzone->danode[A],&aux[A],0);
      EcA = -(e*VA + aux[A].affinity + mt->band->EgNarrowToEc(T));//conduction band energy level
      EvA = -(e*VA + aux[A].affinity - mt->band->EgNarrowToEv(T) + Eg);//valence band energy level
      RA = mt->band->Recomb(pA,nA,TD);
      niA = mt->band->nie(T);

      mt->mapping(&pzone->danode[B],&aux[B],0);
      EcB = -(e*VB + aux[B].affinity + mt->band->EgNarrowToEc(T));//conduction band energy level
      EvB = -(e*VB + aux[B].affinity - mt->band->EgNarrowToEv(T) + Eg);//valence band energy level
      RB = mt->band->Recomb(pB,nB,TD);
      niB = mt->band->nie(T);

      mt->mapping(&pzone->danode[C],&aux[C],0);
      EcC = -(e*VC + aux[C].affinity + mt->band->EgNarrowToEc(T));//conduction band energy level
      EvC = -(e*VC + aux[C].affinity - mt->band->EgNarrowToEv(T) + Eg);//valence band energy level
      RC = mt->band->Recomb(pC,nC,TD);
      niC = mt->band->nie(T);
#ifdef _FERMI_  
      if(Fermi)
      {
        PetscScalar NcA = mt->band->Nc(fs[A].T);
        PetscScalar NvA = mt->band->Nv(fs[A].T);
        EcA = EcA - e*Vt*log(gamma_f(fabs(nA)/NcA));
        EvA = EvA + e*Vt*log(gamma_f(fabs(pA)/NvA));

        PetscScalar NcB = mt->band->Nc(fs[B].T);
        PetscScalar NvB = mt->band->Nv(fs[B].T);
        EcB = EcB - e*Vt*log(gamma_f(fabs(nB)/NcB));
        EvB = EvB + e*Vt*log(gamma_f(fabs(pB)/NvB));

        PetscScalar NcC = mt->band->Nc(fs[C].T);
        PetscScalar NvC = mt->band->Nv(fs[C].T);
        EcC = EcC - e*Vt*log(gamma_f(fabs(nC)/NcC));
        EvC = EvC + e*Vt*log(gamma_f(fabs(pC)/NvC));
      }
#endif 
    } else {
      rotate3ADScalar(EcA,EcB,EcC);
      rotate3ADScalar(EvA,EvB,EvC);
      rotate3ADScalar(RA,RB,RC);
      rotate3Scalar(niA,niB,niC);
    } 

    if(pEdge==0) {
      Ex = -((yB-yC)*VA + (yC-yA)*VB +(yA-yB)*VC)/(2*tri_area);
      Ey = -((xC-xB)*VA + (xA-xC)*VB +(xB-xA)*VC)/(2*tri_area);
      E = sqrt(Ex*Ex+Ey*Ey+1e-20);
    }
    if(pEdge==0) {
      phinA = VA - log(fabs(nA)/niA)*Vt;
      phipA = VA + log(fabs(pA)/niA)*Vt;
      phinB = VB - log(fabs(nB)/niB)*Vt;
      phipB = VB + log(fabs(pB)/niB)*Vt;
      phinC = VC - log(fabs(nC)/niC)*Vt;
      phipC = VC + log(fabs(pC)/niC)*Vt;
    } else {
      rotate3ADScalar(phinA,phinB,phinC);
      rotate3ADScalar(phipA,phipB,phipC);
    }

#ifdef _II_
    AutoDScalar IIn, IIp;
    if(ImpactIonization) {
      if(IIType==GradQf && pEdge==0)
      {
        AutoDScalar Fnx = ((yB-yC)*phinA + (yC-yA)*phinB +(yA-yB)*phinC)/(2*tri_area);
        AutoDScalar Fny = ((xC-xB)*phinA + (xA-xC)*phinB +(xB-xA)*phinC)/(2*tri_area);
        AutoDScalar Fpx = ((yB-yC)*phipA + (yC-yA)*phipB +(yA-yB)*phipC)/(2*tri_area);
        AutoDScalar Fpy = ((xC-xB)*phipA + (xA-xC)*phipB +(xB-xA)*phipC)/(2*tri_area);
        AutoDScalar Fn  =  sqrt(Fnx*Fnx+Fny*Fny);
        AutoDScalar Fp  =  sqrt(Fpx*Fpx+Fpy*Fpy);
        IIn =  mt->gen->ElecGenRate(TD,Fn,Eg);
        IIp =  mt->gen->HoleGenRate(TD,Fp,Eg);
      }
      if(IIType==EVector && pEdge==0)
      {
        AutoDScalar Fnx = Ex;
        AutoDScalar Fny = Ey;
        AutoDScalar Fpx = Ex;
        AutoDScalar Fpy = Ey;
        AutoDScalar Fn  =  sqrt(Fnx*Fnx+Fny*Fny);
        AutoDScalar Fp  =  sqrt(Fpx*Fpx+Fpy*Fpy);
        IIn =  mt->gen->ElecGenRate(TD,Fn,Eg);
        IIp =  mt->gen->HoleGenRate(TD,Fp,Eg);
      }
    }
#endif
#ifdef _FLUX1_
    if(pEdge==0) {
      Jnc =  In(Vt,EcA/e,EcB/e,nA,nB,Lc);
      Jpc =  Ip(Vt,EvA/e,EvB/e,pA,pB,Lc);
      Jna =  In(Vt,EcB/e,EcC/e,nB,nC,La);
      Jpa =  Ip(Vt,EvB/e,EvC/e,pB,pC,La);
      Jnb =  In(Vt,EcC/e,EcA/e,nC,nA,Lb);
      Jpb =  Ip(Vt,EvC/e,EvA/e,pC,pA,Lb);
      Jn_scale = adtl::fmax(adtl::fmax(fabs(Jna.getValue()),fabs(Jnb.getValue())),
                           adtl::fmax(fabs(Jnc.getValue()),niA*niA));
      Jp_scale = adtl::fmax(adtl::fmax(fabs(Jpa.getValue()),fabs(Jpb.getValue())),
                           adtl::fmax(fabs(Jpc.getValue()),niA*niA));
      Jnc /=  Jn_scale;
      Jpc /=  Jp_scale;
      Jna /=  Jn_scale;
      Jpa /=  Jp_scale;
      Jnb /=  Jn_scale;
      Jpb /=  Jp_scale;
    } else {
      rotate3ADScalar(Jna,Jnb,Jnc);
      rotate3ADScalar(Jpa,Jpb,Jpc);
    }
#endif
#ifdef _FLUX2_
    if(pEdge==0) {
      Jnc =  In(Vt,(EcA-EcB)/e,nA,nB,Lc);
      Jpc =  Ip(Vt,(EvA-EvB)/e,pA,pB,Lc);
      Jna =  In(Vt,(EcB-EcC)/e,nB,nC,La);
      Jpa =  Ip(Vt,(EvB-EvC)/e,pC,La);
      Jnb =  In(Vt,(EcC-EcA)/e,nC,nA,Lb);
      Jpb =  Ip(Vt,(EvC-EvA)/e,pC,pA,Lb);
      Jn_scale = adtl:fmax(adtl:fmax(fabs(Jna.getValue()),fabs(Jnb.getValue())),
                           adtl::fmax(fabs(Jnc.getValue()),niA*niA));
      Jp_scale = adtl:fmax(adtl:fmax(fabs(Jpa.getValue()),fabs(Jpb.getValue())),
                           adtl::fmax(fabs(Jpc.getValue()),niA*niA));
      Jnc /=  Jn_scale;
      Jpc /=  Jp_scale;
      Jna /=  Jn_scale;
      Jpa /=  Jp_scale;
      Jnb /=  Jn_scale;
      Jpb /=  Jp_scale;
    } else {
      rotate3ADScalar(Jna,Jnb,Jnc);
      rotate3ADScalar(Jpa,Jpb,Jpc);
    }
#endif

    //flux along A-B
    AutoDScalar Epnc, Etnc, Eppc, Etpc;
    AutoDScalar Jnc_norm, Jpc_norm;
    AutoDScalar mun, mup;
    if(HighFieldMobility)
    {
      if(EJModel || IIType==EdotJ)
      {

        AutoDScalar JncxTca  = (Jnc*Say - Jna*Scy)/sinB;
        AutoDScalar JncxTbc  = (Jnb*Scy - Jnc*Sby)/sinA;
        AutoDScalar Jncx = (da * JncxTca + db * JncxTbc)/(da+db);

        AutoDScalar JncyTca  = (-Jnc*Sax + Jna*Scx)/sinB;
        AutoDScalar JncyTbc  = (-Jnb*Scx + Jnc*Sbx)/sinA;
        AutoDScalar Jncy = (da * JncyTca + db * JncyTbc)/(da+db);

        AutoDScalar JpcxTca  = (Jpc*Say - Jpa*Scy)/sinB;
        AutoDScalar JpcxTbc  = (Jpb*Scy - Jpc*Sby)/sinA;
        AutoDScalar Jpcx = (da * JpcxTca + db * JpcxTbc)/(da+db);

        AutoDScalar JpcyTca  = (-Jpc*Sax + Jpa*Scx)/sinB;
        AutoDScalar JpcyTbc  = (-Jpb*Scx + Jpc*Sbx)/sinA;
        AutoDScalar Jpcy = (da * JpcyTca + db * JpcyTbc)/(da+db);

        Jnc_norm = sqrt(Jncx*Jncx + Jncy*Jncy + 1e-100);
        Jpc_norm = sqrt(Jpcx*Jpcx + Jpcy*Jpcy + 1e-100);
        Epnc = (Ex*Jncx + Ey*Jncy)/Jnc_norm;
        Etnc = (Ex*Jncy - Ey*Jncx)/Jnc_norm;
        Eppc = (Ex*Jpcx + Ey*Jpcy)/Jpc_norm;
        Etpc = (Ex*Jpcy - Ey*Jpcx)/Jpc_norm;
      }
      else
      {
        Epnc = fabs(EcA-EcB)/e/Lc; //parallel electrical field for electron
        Eppc = fabs(EvA-EvB)/e/Lc; //parallel electrical field for hole
        //transvers electrical field for electron and hole
        Etnc = fabs(EcA + Lb*cosA*(EcB-EcA)/Lc - EcC)/(Lb*sinA)/e;
        Etpc = fabs(EvA + Lb*cosA*(EvB-EvA)/Lc - EvC)/(Lb*sinA)/e;
      }

      mt->mapping(&pzone->danode[A],&aux[A],0);
      AutoDScalar munCA =  mt->mob->ElecMob(pA,nA,TD,adtl::fmax(0,Epnc),fabs(Etnc),TD);
      AutoDScalar mupCA =  mt->mob->HoleMob(pA,nA,TD,adtl::fmax(0,Eppc),fabs(Etpc),TD);

      mt->mapping(&pzone->danode[B],&aux[B],0);
      AutoDScalar munCB =  mt->mob->ElecMob(pB,nB,TD,adtl::fmax(0,Epnc),fabs(Etnc),TD);
      AutoDScalar mupCB =  mt->mob->HoleMob(pB,nB,TD,adtl::fmax(0,Eppc),fabs(Etpc),TD);

      mun = 0.5*(munCA+munCB);
      mup = 0.5*(mupCA+mupCB);
    } else {
      mun = 0.5*(aux[A].mun+aux[B].mun);
      mup = 0.5*(aux[A].mup+aux[B].mup);
    }

    AutoDScalar G=0;
#ifdef _II_
    if(ImpactIonization && IIType==EdotJ)
    {
      IIn =  mt->gen->ElecGenRate(TD,adtl::fmax(0,Epnc),Eg);
      IIp =  mt->gen->HoleGenRate(TD,adtl::fmax(0,Eppc),Eg);
      G   =  IIn*mun*Jnc_norm*Jn_scale + IIp*mup*Jpc_norm*Jp_scale;
    }
    if( ImpactIonization && (IIType==GradQf || IIType==EVector) )
    {
      G = IIn*mun*fabs(phinA-phinB)/Lc*nmid(Vt,VA,VB,nA,nB)
        +IIp*mup*fabs(phipA-phipB)/Lc*pmid(Vt,VA,VB,pA,pB);
    }
    if( ImpactIonization && IIType==ESide )
    {
      IIn =  mt->gen->ElecGenRate(TD,fabs(VA-VB)/Lc,Eg);
      IIp =  mt->gen->HoleGenRate(TD,fabs(VA-VB)/Lc,Eg);
      G = IIn*mun*fabs(Jnc)*Jn_scale + IIp*mup*fabs(Jpc)*Jp_scale;
    }
    if (ImpactIonization)
    {
      PetscScalar S = 0.25*Lc*dc;
      J1.m[0] =  0; J1.m[1] =  0; J1.m[2] =  0;
      J1.m[3] =  G.getADValue(offset1)*S; J1.m[4] =  G.getADValue(offset1+1)*S; J1.m[5] =  G.getADValue(offset1+2)*S;
      J1.m[6] =  G.getADValue(offset1)*S; J1.m[7] =  G.getADValue(offset1+1)*S; J1.m[8] =  G.getADValue(offset1+2)*S;

      J2.m[0] =  0; J2.m[1] =  0; J2.m[2] =  0;
      J2.m[3] =  G.getADValue(offset2)*S; J2.m[4] =  G.getADValue(offset2+1)*S; J2.m[5] =  G.getADValue(offset2+2)*S;
      J2.m[6] =  G.getADValue(offset2)*S; J2.m[7] =  G.getADValue(offset2+1)*S; J2.m[8] =  G.getADValue(offset2+2)*S;

      J3.m[0] =  0; J3.m[1] =  0; J3.m[2] =  0;
      J3.m[3] =  G.getADValue(offset3)*S; J3.m[4] =  G.getADValue(offset3+1)*S; J3.m[5] =  G.getADValue(offset3+2)*S;
      J3.m[6] =  G.getADValue(offset3)*S; J3.m[7] =  G.getADValue(offset3+1)*S; J3.m[8] =  G.getADValue(offset3+2)*S;

      MatSetValues(*jtmp,3,index1,3,index1,J1.m,ADD_VALUES);
      MatSetValues(*jtmp,3,index1,3,index2,J2.m,ADD_VALUES);
      MatSetValues(*jtmp,3,index1,3,index3,J3.m,ADD_VALUES);
      MatSetValues(*jtmp,3,index2,3,index1,J1.m,ADD_VALUES);
      MatSetValues(*jtmp,3,index2,3,index2,J2.m,ADD_VALUES);
      MatSetValues(*jtmp,3,index2,3,index3,J3.m,ADD_VALUES);
    }
#endif
    //---------------------------------------------------------------------------
    AutoDScalar FA0 =   0.5*(aux[A].eps+aux[B].eps)*(VB-VA)/Lc*dc;
    AutoDScalar FA1 =   mun*Jnc*Jn_scale*dc;
    AutoDScalar FA2 = - mup*Jpc*Jp_scale*dc;

    J1.m[0] =  FA0.getADValue(offset1);  J1.m[1] =  FA0.getADValue(offset1+1);  J1.m[2] =  FA0.getADValue(offset1+2);
    J1.m[3] =  FA1.getADValue(offset1);  J1.m[4] =  FA1.getADValue(offset1+1);  J1.m[5] =  FA1.getADValue(offset1+2);
    J1.m[6] =  FA2.getADValue(offset1);  J1.m[7] =  FA2.getADValue(offset1+1);  J1.m[8] =  FA2.getADValue(offset1+2);


    J2.m[0] =  FA0.getADValue(offset2);  J2.m[1] =  FA0.getADValue(offset2+1);  J2.m[2] =  FA0.getADValue(offset2+2);
    J2.m[3] =  FA1.getADValue(offset2);  J2.m[4] =  FA1.getADValue(offset2+1);  J2.m[5] =  FA1.getADValue(offset2+2);
    J2.m[6] =  FA2.getADValue(offset2);  J2.m[7] =  FA2.getADValue(offset2+1);  J2.m[8] =  FA2.getADValue(offset2+2);


    J3.m[0] =  FA0.getADValue(offset3);  J3.m[1] =  FA0.getADValue(offset3+1);  J3.m[2] =  FA0.getADValue(offset3+2);
    J3.m[3] =  FA1.getADValue(offset3);  J3.m[4] =  FA1.getADValue(offset3+1);  J3.m[5] =  FA1.getADValue(offset3+2);
    J3.m[6] =  FA2.getADValue(offset3);  J3.m[7] =  FA2.getADValue(offset3+1);  J3.m[8] =  FA2.getADValue(offset3+2);

    MatSetValues(*jtmp,3,index1,3,index1,J1.m,ADD_VALUES);
    MatSetValues(*jtmp,3,index1,3,index2,J2.m,ADD_VALUES);
    MatSetValues(*jtmp,3,index1,3,index3,J3.m,ADD_VALUES);
    MatSetValues(*jtmp,3,index2,3,index1,(-J1).m,ADD_VALUES);
    MatSetValues(*jtmp,3,index2,3,index2,(-J2).m,ADD_VALUES);
    MatSetValues(*jtmp,3,index2,3,index3,(-J3).m,ADD_VALUES);


 //Generation at A
 #ifdef _HURKX_BTBT_
    if (BandBandTunneling)
    {
      PetscScalar S = 0.25*dc*Lc + 0.25*db*Lb;

      AutoDScalar GradQfnx = -((yB-yC)*phinA + (yC-yA)*phinB +(yA-yB)*phinC)/(2*tri_area);
      AutoDScalar GradQfny = -((xC-xB)*phinA + (xA-xC)*phinB +(xB-xA)*phinC)/(2*tri_area);
      AutoDScalar GradQfn = sqrt(GradQfnx*GradQfnx+GradQfny*GradQfny+1e-20);

      AutoDScalar GradQfpx = -((yB-yC)*phipA + (yC-yA)*phipB +(yA-yB)*phipC)/(2*tri_area);
      AutoDScalar GradQfpy = -((xC-xB)*phipA + (xA-xC)*phipB +(xB-xA)*phipC)/(2*tri_area);
      AutoDScalar GradQfp = sqrt(GradQfpx*GradQfpx+GradQfpy*GradQfpy+1e-20);

      AutoDScalar nn = nA * pow(niA / mt->band->Nc(fs[A].T), GradQfn/E );
      AutoDScalar pp = pA * pow(niA / mt->band->Nv(fs[A].T), GradQfp/E );

      AutoDScalar D = (niA*niA - nn*pp)/(niA+nn)/(niA+pp);
      AutoDScalar D1 = (GradQfnx*GradQfpx + GradQfny*GradQfpy)/GradQfn/GradQfp;
      AutoDScalar G_BB = D*mt->band->BB_Tunneling(T,E);

      J1.m[0] =  0; J1.m[1] =  0; J1.m[2] =  0;
      J1.m[3] =  G_BB.getADValue(offset1)*S; J1.m[4] =  G_BB.getADValue(offset1+1)*S; J1.m[5] =  G_BB.getADValue(offset1+2)*S;
      J1.m[6] =  G_BB.getADValue(offset1)*S; J1.m[7] =  G_BB.getADValue(offset1+1)*S; J1.m[8] =  G_BB.getADValue(offset1+2)*S;

      J2.m[0] =  0; J2.m[1] =  0; J2.m[2] =  0;
      J2.m[3] =  G_BB.getADValue(offset2)*S; J2.m[4] =  G_BB.getADValue(offset2+1)*S; J2.m[5] =  G_BB.getADValue(offset2+2)*S;
      J2.m[6] =  G_BB.getADValue(offset2)*S; J2.m[7] =  G_BB.getADValue(offset2+1)*S; J2.m[8] =  G_BB.getADValue(offset2+2)*S;

      J3.m[0] =  0; J3.m[1] =  0; J3.m[2] =  0;
      J3.m[3] =  G_BB.getADValue(offset3)*S; J3.m[4] =  G_BB.getADValue(offset3+1)*S; J3.m[5] =  G_BB.getADValue(offset3+2)*S;
      J3.m[6] =  G_BB.getADValue(offset3)*S; J3.m[7] =  G_BB.getADValue(offset3+1)*S; J3.m[8] =  G_BB.getADValue(offset3+2)*S;

      MatSetValues(*jtmp,3,index1,3,index1,J1.m,ADD_VALUES);
      MatSetValues(*jtmp,3,index1,3,index2,J2.m,ADD_VALUES);
      MatSetValues(*jtmp,3,index1,3,index3,J3.m,ADD_VALUES);
    }
#endif

    if(pEdge==0) {
#ifndef _HURKX_BTBT_
      if(BandBandTunneling)
      {
        AutoDScalar G_BB = mt->band->BB_Tunneling(T,E);

        MatSetValue(*jtmp,index1[1],index1[0],G_BB.getADValue(0)*ptri->s[0],ADD_VALUES);
        MatSetValue(*jtmp,index1[1],index2[0],G_BB.getADValue(3)*ptri->s[0],ADD_VALUES);
        MatSetValue(*jtmp,index1[1],index3[0],G_BB.getADValue(6)*ptri->s[0],ADD_VALUES);
        MatSetValue(*jtmp,index1[2],index1[0],G_BB.getADValue(0)*ptri->s[0],ADD_VALUES);
        MatSetValue(*jtmp,index1[2],index2[0],G_BB.getADValue(3)*ptri->s[0],ADD_VALUES);
        MatSetValue(*jtmp,index1[2],index3[0],G_BB.getADValue(6)*ptri->s[0],ADD_VALUES);

        MatSetValue(*jtmp,index2[1],index1[0],G_BB.getADValue(0)*ptri->s[1],ADD_VALUES);
        MatSetValue(*jtmp,index2[1],index2[0],G_BB.getADValue(3)*ptri->s[1],ADD_VALUES);
        MatSetValue(*jtmp,index2[1],index3[0],G_BB.getADValue(6)*ptri->s[1],ADD_VALUES);
        MatSetValue(*jtmp,index2[2],index1[0],G_BB.getADValue(0)*ptri->s[1],ADD_VALUES);
        MatSetValue(*jtmp,index2[2],index2[0],G_BB.getADValue(3)*ptri->s[1],ADD_VALUES);
        MatSetValue(*jtmp,index2[2],index3[0],G_BB.getADValue(6)*ptri->s[1],ADD_VALUES);

        MatSetValue(*jtmp,index3[1],index1[0],G_BB.getADValue(0)*ptri->s[2],ADD_VALUES);
        MatSetValue(*jtmp,index3[1],index2[0],G_BB.getADValue(3)*ptri->s[2],ADD_VALUES);
        MatSetValue(*jtmp,index3[1],index3[0],G_BB.getADValue(6)*ptri->s[2],ADD_VALUES);
        MatSetValue(*jtmp,index3[2],index1[0],G_BB.getADValue(0)*ptri->s[2],ADD_VALUES);
        MatSetValue(*jtmp,index3[2],index2[0],G_BB.getADValue(3)*ptri->s[2],ADD_VALUES);
        MatSetValue(*jtmp,index3[2],index3[0],G_BB.getADValue(6)*ptri->s[2],ADD_VALUES);
      }  
#endif
      //---------------------------------------------------------------------------
      Set_Mat3_zero(J1);
      PetscScalar S;
      S = 0.25*ptri->d[2]*ptri->edge_len[2] + 0.25*ptri->d[1]*ptri->edge_len[1];
      J1.m[4] =  - RA.getADValue(1)*S;   
      J1.m[5] =  - RA.getADValue(2)*S;   
      J1.m[7] =  - RA.getADValue(1)*S;  
      J1.m[8] =  - RA.getADValue(2)*S;  
      MatSetValues(*jtmp,3,index1,3,index1,J1.m,ADD_VALUES);

      //---------------------------------------------------------------------------
      Set_Mat3_zero(J2);
      S = 0.25*ptri->d[0]*ptri->edge_len[0] + 0.25*ptri->d[2]*ptri->edge_len[2];
      J2.m[4] =  - RB.getADValue(4)*S;   
      J2.m[5] =  - RB.getADValue(5)*S;   
      J2.m[7] =  - RB.getADValue(4)*S;  
      J2.m[8] =  - RB.getADValue(5)*S;  
      MatSetValues(*jtmp,3,index2,3,index2,J2.m,ADD_VALUES);

      //---------------------------------------------------------------------------
      Set_Mat3_zero(J3); 
      S = 0.25*ptri->d[1]*ptri->edge_len[1] + 0.25*ptri->d[0]*ptri->edge_len[0];
      J3.m[4] =  - RC.getADValue(7)*S;   
      J3.m[5] =  - RC.getADValue(8)*S;   
      J3.m[7] =  - RC.getADValue(7)*S;  
      J3.m[8] =  - RC.getADValue(8)*S;  
      MatSetValues(*jtmp,3,index3,3,index3,J3.m,ADD_VALUES);
    }

  }
}

//-----------------------------------------------------------------------------
// boundaries
//-----------------------------------------------------------------------------


void SMCZone::F1E_ddm_inner(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  PetscScalar ni = x[zofs[zone_index]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+3*i+2];     //hole density of node i
  
  f[zofs[zone_index]+3*i+0] = f[zofs[zone_index]+3*i+0]/pcell->area + mt->e*((pi-ni)+(aux[i].Nd-aux[i].Na));
  f[zofs[zone_index]+3*i+1] = f[zofs[zone_index]+3*i+1]/pcell->area;
  f[zofs[zone_index]+3*i+2] = f[zofs[zone_index]+3*i+2]/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+3*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+3*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+3*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+3*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F1E_ddm_ombc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int z = zone_index;
  int equ_num = 3;
  int size = pzone->davcell.size();
  PetscScalar e  =  mt->e;
  PetscScalar kb =  mt->kb;
  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+3*i+2];     //hole density of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar nie = mt->band->nie(fs[i].T);
  PetscScalar Nc  = mt->band->Nc(fs[i].T);
  PetscScalar Nv  = mt->band->Nv(fs[i].T);
  PetscScalar Vt  = kb*fs[i].T/e;
  int    om_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1) {om_equ=j;break;}
  }
    
  if(Fermi) //Fermi
  {
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) );
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) + mt->band->Eg(fs[i].T));
    PetscScalar phin = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar phip = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    f[zofs[z]+3*i+0] = Nc*fermi_half(etan) - Nv*fermi_half(etap)+ Na - Nd;
    f[zofs[z]+3*i+1] = ni - Nc*fermi_half(etan);
    f[zofs[z]+3*i+2] = pi - Nv*fermi_half(etap);
  }
  else     //Boltzmann
  {
    f[zofs[z]+3*i+0] = Vi - kb*fs[i].T/e*asinh((Nd-Na)/(2*nie)) + kb*fs[i].T/2/e*log(Nc/Nv) + mt->band->Eg(fs[i].T)/2/e 
                       + aux[i].affinity -x[zofs[z]+equ_num*size+om_equ];
    PetscScalar electron_density,hole_density;
    if(Na>Nd)   //p-type
    {
      hole_density = (-(Nd-Na)+sqrt((Nd-Na)*(Nd-Na)+4*nie*nie))/2.0;
      electron_density = nie*nie/hole_density;
    }
    else        //n-type
    {
      electron_density = ((Nd-Na)+sqrt((Nd-Na)*(Nd-Na)+4*nie*nie))/2.0;
      hole_density = nie*nie/electron_density;
    }
    f[zofs[z]+3*i+1] = ni - electron_density;  //electron density
    f[zofs[z]+3*i+2] = pi - hole_density;      //hole density
  }
}


void SMCZone::F1E_ddm_stkbc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int>& zofs,DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 3;
  int size = pzone->davcell.size();
  PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+3*i+2];     //hole density of node i
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  int    stk_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)   
    {stk_equ=j;break;}
  //Schotty Barrier Lowerring
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  //Schottky current
  PetscScalar Fn = mt->band->SchottyJsn(ni,fs[i].T,pbc->WorkFunction-aux[i].affinity-deltaVB)*L;
  PetscScalar Fp = mt->band->SchottyJsp(pi,fs[i].T,pbc->WorkFunction-aux[i].affinity+deltaVB)*L;

  f[zofs[z]+3*i+0] = x[zofs[z]+3*i+0] + pbc->WorkFunction - deltaVB - x[zofs[z]+equ_num*size+stk_equ];
  f[zofs[zone_index]+3*i+1] = (f[zofs[zone_index]+3*i+1]+Fn)/pcell->area;
  f[zofs[zone_index]+3*i+2] = (f[zofs[zone_index]+3*i+2]-Fp)/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+3*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+3*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+3*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+3*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F1E_ddm_insulator_gate(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int equ_num = 3;
  int size = pzone->davcell.size();
  PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+3*i+2];     //hole density of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar grad_P = 0;
  int    ins_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
     {ins_equ=j;break;}
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
    //the poisson's equation on third boundary type
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar vgate = x[zofs[zone_index]+equ_num*size+ins_equ] - pbc->WorkFunction;
      PetscScalar q = mt->e*pbc->QF; //sigma is the surface change density
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar r=q + eps_ox/Thick*vgate;
      PetscScalar s=eps_ox/Thick;
      grad_P += 0.5*pcell->ilen[j]*(r-0.25*s*(3*Vi+Vj));
    }
  }
  f[zofs[zone_index]+3*i+0] = (f[zofs[zone_index]+3*i+0]+grad_P)/pcell->area
                              +  mt->e*((pi-ni)+(aux[i].Nd-aux[i].Na));
  f[zofs[zone_index]+3*i+1] = f[zofs[zone_index]+3*i+1]/pcell->area;
  f[zofs[zone_index]+3*i+2] = f[zofs[zone_index]+3*i+2]/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+3*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+3*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+3*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+3*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F1E_ddm_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                               ISZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorInterfaceBC *pbc = dynamic_cast<InsulatorInterfaceBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+3*i+2];     //hole density of node i
  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar grad_P = 0;
  
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
    grad_P += aux[i].eps*pcell->elen[j]/pcell->ilen[j]*(Vj-Vi);
  }
  //poisson's equation at interface
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int    nb = ncell->nb_array[j];
    PetscScalar Vj_n = x[zofs[pz->pzone->zone_index]+nb];     //potential of nb node
    grad_P += pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j]*(Vj_n-Vi);
  }

  f[zofs[zone_index]+3*i+0] =   grad_P + mt->e*((pi-ni)+(Nd-Na))*pcell->area + pbc->QF*L;
  f[zofs[zone_index]+3*i+1] =   f[zofs[zone_index]+3*i+1]/pcell->area;
  f[zofs[zone_index]+3*i+2] =   f[zofs[zone_index]+3*i+2]/pcell->area;
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+3*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+3*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+3*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+3*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F1E_ddm_homojunction(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                               SMCZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  HomoInterfaceBC *pbc = dynamic_cast<HomoInterfaceBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+3*i+2];     //hole density of node i

  if(zone_index > pz->pzone->zone_index)
  {
        f[zofs[zone_index]+3*i+0] =  Vi - x[zofs[pz->zone_index]+3*n+0];
        f[zofs[zone_index]+3*i+1] =  ni - x[zofs[pz->zone_index]+3*n+1];
        f[zofs[zone_index]+3*i+2] =  pi - x[zofs[pz->zone_index]+3*n+2];
        return;
  }

  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  PetscScalar grad_P=0,Fn=0,Fp=0;

  //process half cell in local zone
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
    grad_P += aux[i].eps*pcell->elen[j]/pcell->ilen[j]*(Vj-Vi);
  }
  Fn +=  f[zofs[zone_index]+3*i+1];
  Fp += -f[zofs[zone_index]+3*i+2];

  //process another half cell in dornor zone
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int  nb = ncell->nb_array[j];
    PetscScalar Vj  = x[zofs[pz->zone_index]+3*nb+0];     //potential of nb node
    grad_P += pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j]*(Vj-Vi);
  }
  Fn +=  f[zofs[pz->zone_index]+3*n+1];
  Fp += -f[zofs[pz->zone_index]+3*n+2];

  PetscScalar area = (pcell->area+ncell->area);
  f[zofs[zone_index]+3*i+0] =  grad_P + mt->e*((pi-ni)+(Nd-Na))*area;
  f[zofs[zone_index]+3*i+1] =   Fn;
  f[zofs[zone_index]+3*i+2] = - Fp;
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+3*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt)*area;
      f[zofs[zone_index]+3*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt)*area;
    }
    else //first order
    {
      f[zofs[zone_index]+3*i+1] += -(ni-fs[i].n)/ODE_F.dt*area;
      f[zofs[zone_index]+3*i+2] += -(pi-fs[i].p)/ODE_F.dt*area;
    }
  }
}


void SMCZone::F1E_ddm_heterojunction(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                               SMCZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  HeteroInterfaceBC *pbc = dynamic_cast<HeteroInterfaceBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar k = mt->kb;
  PetscScalar e  =  mt->e;
  PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+3*i+2];     //hole density of node i
  PetscScalar nd = x[zofs[pz->pzone->zone_index]+3*n+1];
  PetscScalar pd = x[zofs[pz->pzone->zone_index]+3*n+2];
  PetscScalar Ti = fs[i].T;
  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar Eci =  -e*(Vi + aux[i].affinity + mt->kb*fs[i].T/e*log(aux[i].Nc));
  PetscScalar Evi =  -e*(Vi + aux[i].affinity + aux[i].Eg/e - mt->kb*fs[i].T/e*log(aux[i].Nv));
  PetscScalar grad_P = 0;
  PetscScalar Fn=0,Fp=0;

  PetscScalar L = (pcell->ilen[0]+pcell->ilen[pcell->nb_num-1])/2.0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
    grad_P += aux[i].eps*pcell->elen[j]/pcell->ilen[j]*(Vj-Vi);
    //Thermal emit current
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Ecj = -e*(Vi + pz->aux[n].affinity + k*pz->fs[n].T/e*log(pz->aux[n].Nc));
      PetscScalar Evj = -e*(Vi + pz->aux[n].affinity - k*pz->fs[n].T/e*log(pz->aux[n].Nv) + pz->aux[n].Eg/e);
      PetscScalar nj  = x[zofs[pz->pzone->zone_index]+3*n+1];
      PetscScalar pj  = x[zofs[pz->pzone->zone_index]+3*n+2];
      PetscScalar jn=0,jp=0;
      if(Ecj > Eci)
      {
        PetscScalar pm = pz->mt->band->EffecElecMass(Ti)/mt->band->EffecElecMass(Ti);
        jn = -2*e*(pz->mt->band->ThermalVn(Ti)*nj - pm*mt->band->ThermalVn(Ti)*ni*exp(-(Ecj-Eci)/(k*Ti)));
      }
      else
      {
        PetscScalar pm = mt->band->EffecElecMass(Ti)/pz->mt->band->EffecElecMass(Ti);
        jn = -2*e*(mt->band->ThermalVn(Ti)*ni - pm*pz->mt->band->ThermalVn(Ti)*nj*exp(-(Eci-Ecj)/(k*Ti)));
      }

      if(Evj > Evi)
      {
        PetscScalar pm = pz->mt->band->EffecHoleMass(Ti)/mt->band->EffecHoleMass(Ti);
        jp = 2*e*(pz->mt->band->ThermalVp(Ti)*pj - pm*mt->band->ThermalVp(Ti)*pi*exp(-(Evj-Evi)/(k*Ti)));
      }
      else
      {
        PetscScalar pm = mt->band->EffecHoleMass(Ti)/pz->mt->band->EffecHoleMass(Ti);
        jp = 2*e*(mt->band->ThermalVp(Ti)*pi - pm*pz->mt->band->ThermalVp(Ti)*pj*exp(-(Evi-Evj)/(k*Ti)));
      }

      Fn += jn*0.5*pcell->ilen[j];
      Fp += jp*0.5*pcell->ilen[j];
    }
  }

  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int  nb = ncell->nb_array[j];
    PetscScalar Vr = x[zofs[pz->pzone->zone_index]+3*nb+0];     //potential of nb node
    grad_P += pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j]*(Vr-Vi);
  }

  if(pzone->zone_index < pz->pzone->zone_index)
        f[zofs[zone_index]+3*i+0] =  grad_P + e*((pi-ni)+(Nd-Na))*pcell->area
                                        + e*((pd-nd)+(pz->aux[n].Nd-pz->aux[n].Na))*ncell->area
                                        + pbc->QF*L;
  else
        f[zofs[zone_index]+3*i+0] =  Vi - x[zofs[pz->pzone->zone_index]+3*n+0];

  f[zofs[zone_index]+3*i+1] = (f[zofs[zone_index]+3*i+1]+Fn)/pcell->area;
  f[zofs[zone_index]+3*i+2] = (f[zofs[zone_index]+3*i+2]-Fp)/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+3*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+3*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+3*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+3*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F1E_om_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(bc_index));
  PetscScalar e  =  mt->e;
  
  //calculate the total current of ohmic electrode
  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+3*node+0];     //potential of node i
    //conduction current
    current += DeviceDepth*(f[zofs[zone_index]+3*node+1]-f[zofs[zone_index]+3*node+2]);
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
      //displacement current
      current += DeviceDepth*pcell->elen[k]*aux[node].eps*((Vi-Vj)-(fs[node].P-fs[nb].P))/pcell->ilen[k]/ODE_F.dt;
    }
  }
  
  //if the electrode connect to another electrode
  if(pbc->inner_connect!=-1)
  {
        int connect_to = pbc->inner_connect;
        int connect_zone = pbc->connect_zone;
        int connect_elec=-1;
        SMCZone * pz = dynamic_cast <SMCZone * > (pbc->pzonedata);
        for(int j=0;j<pz->electrode.size();j++)
        {
                if(bc[pz->electrode[j]].BCType==OhmicContact)
                {
                  OhmicBC *om_bc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pz->electrode[j]));
                  if(om_bc->inner_connect==bc_index)              connect_elec=j;
                }  
        }
        int connect_index=zofs[connect_zone]+equ_num*pz->pzone->davcell.size()+connect_elec;

        if(bc_index < connect_to)
        {
          f[zofs[zone_index]+equ_num*size+i]+=x[zofs[zone_index]+equ_num*size+i];          //voltage equation
          f[connect_index]+=current;   //current equation 
        }
        else
        {
          f[zofs[zone_index]+equ_num*size+i]+=current; //current equation
          f[connect_index]+=-x[zofs[zone_index]+equ_num*size+i];//voltage equation
        }
        pbc->Set_Current_new(current);
  }     
  //the electrode has an external voltage circult
  else if(pbc->electrode_type==VoltageBC)
  {
    PetscScalar V = pbc->Vapp;
    PetscScalar R = pbc->R;
    PetscScalar C = pbc->C;
    PetscScalar L = pbc->L;
    PetscScalar In = pbc->current;
    PetscScalar Icn = pbc->cap_current;
    PetscScalar Pn = pbc->potential;
    PetscScalar P = x[zofs[zone_index]+equ_num*size+i];
#ifdef    _EXP_R_
    if(R<1.0) 
      f[zofs[zone_index]+equ_num*size+i]=(L/ODE_F.dt+R)*current + (P-V) + (L/ODE_F.dt+R)*C/ODE_F.dt*P
                                         -(L/ODE_F.dt+R)*C/ODE_F.dt*Pn - L/ODE_F.dt*(In+Icn);
    else //when R is large, more like a current cource (R->+inf) 
      f[zofs[zone_index]+equ_num*size+i]=current +(P-V)/(L/ODE_F.dt+R)  + C/ODE_F.dt*P
                                         -C/ODE_F.dt*Pn - L/ODE_F.dt*(In+Icn)/(L/ODE_F.dt+R);    
#else
   f[zofs[zone_index]+equ_num*size+i]=(L/ODE_F.dt+R)*current-V+(1+(L/ODE_F.dt+R)*C/ODE_F.dt)*P
                                       -(L/ODE_F.dt+R)*C/ODE_F.dt*Pn-L/ODE_F.dt*(In+Icn);
#endif                                   
    pbc->Set_Current_new(current);
  }
  //the electrode has an external current circult
  else if(pbc->electrode_type==CurrentBC)
  {
    PetscScalar I = pbc->Iapp;
    PetscScalar C = pbc->C;
    f[zofs[zone_index]+equ_num*size+i]=current - I;
    pbc->Set_Current_new(I);
  }
  pbc->Set_Potential_new(x[zofs[zone_index]+equ_num*size+i]);
}


void SMCZone::F1E_stk_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc.Get_pointer(bc_index));

  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+3*node+0];     //potential of node i
    PetscScalar ni = x[zofs[zone_index]+3*node+1];     //electron density of node i
    PetscScalar pi = x[zofs[zone_index]+3*node+2];     //hole density of node i
    mt->mapping(&pzone->danode[node],&aux[node],ODE_F.clock);
    PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[node].eps,sqrt(aux[node].Ex*aux[node].Ex+aux[node].Ey*aux[node].Ey));
    PetscScalar Jsn = mt->band->SchottyJsn(ni,fs[node].T,pbc->WorkFunction-aux[node].affinity-deltaVB);
    PetscScalar Jsp = mt->band->SchottyJsp(pi,fs[node].T,pbc->WorkFunction-aux[node].affinity+deltaVB);
    current += -Jsn*0.5*pcell->ilen[0]*DeviceDepth;
    current += -Jsp*0.5*pcell->ilen[0]*DeviceDepth;
    current += -Jsn*0.5*pcell->ilen[pcell->nb_num-1]*DeviceDepth;
    current += -Jsp*0.5*pcell->ilen[pcell->nb_num-1]*DeviceDepth;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
      //displacement current
      current += DeviceDepth*pcell->elen[k]*aux[node].eps*((Vi-Vj)-(fs[node].P-fs[nb].P))/pcell->ilen[k]/ODE_F.dt;
    }
  }
  if(pbc->electrode_type==VoltageBC)
  {
    PetscScalar R = pbc->R;
    PetscScalar C = pbc->C;
    PetscScalar L = pbc->L;
    PetscScalar V = pbc->Vapp;
    PetscScalar In = pbc->current;
    PetscScalar Icn = pbc->cap_current;
    PetscScalar Pn = pbc->potential;
    PetscScalar P = x[zofs[zone_index]+equ_num*size+i];
    f[zofs[zone_index]+equ_num*size+i]=(L/ODE_F.dt+R)*current-V+(1+(L/ODE_F.dt+R)*C/ODE_F.dt)*P
                                       -(L/ODE_F.dt+R)*C/ODE_F.dt*Pn-L/ODE_F.dt*(In+Icn);
    pbc->Set_Current_new(current);
  }
  else if(pbc->electrode_type==CurrentBC)
  {
    PetscScalar I = pbc->Iapp;
    PetscScalar C = pbc->C;
    f[zofs[zone_index]+equ_num*size+i]=current - I;
    pbc->Set_Current_new(I);
  }
  pbc->Set_Potential_new(x[zofs[zone_index]+equ_num*size+i]);
}


void SMCZone::F1E_ins_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  InsulatorContactBC *pbc = dynamic_cast <InsulatorContactBC * > (bc.Get_pointer(bc_index));

  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+3*node+0];     //potential of node i
    for(int k=0;k<pcell->nb_num;k++)
    {
      int  nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
      //displacement current
      current += DeviceDepth*pcell->elen[k]*aux[node].eps*((Vi-Vj)-(fs[node].P-fs[nb].P))/pcell->ilen[k]/ODE_F.dt;
    }
  }
  if(pbc->electrode_type==VoltageBC)
  {
    PetscScalar R = pbc->R;
    PetscScalar C = pbc->C;
    PetscScalar L = pbc->L;
    PetscScalar V = pbc->Vapp;
    PetscScalar In = pbc->current;
    PetscScalar Icn = pbc->cap_current;
    PetscScalar Pn = pbc->potential;
    PetscScalar P = x[zofs[zone_index]+equ_num*size+i];
    f[zofs[zone_index]+equ_num*size+i]=(L/ODE_F.dt+R)*current-V+(1+(L/ODE_F.dt+R)*C/ODE_F.dt)*P
                                       -(L/ODE_F.dt+R)*C/ODE_F.dt*Pn-L/ODE_F.dt*(In+Icn);
    pbc->Set_Current_new(current);
  }
  pbc->Set_Potential_new(x[zofs[zone_index]+equ_num*size+i]);
}


//-----------------------------------------------------------------------------
// Jacobian matrix for boundaries
//-----------------------------------------------------------------------------


void SMCZone::J1E_ddm_inner(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs)
{
  Mat3      A;
  PetscInt  index[3],col[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int   z = zone_index;

  PetscScalar Vi = x[zofs[z]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i
  PetscScalar area = pcell->area;

  //--------------------------------
  index[0] = zofs[z]+3*i+0;
  index[1] = zofs[z]+3*i+1;
  index[2] = zofs[z]+3*i+2;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    col[0] = zofs[z]+3*nb+0;
    col[1] = zofs[z]+3*nb+1;
    col[2] = zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,index,3,col,A.m);
    A=A/area;
    MatSetValues(*jac,3,index,3,col,A.m,INSERT_VALUES);
  }
  MatGetValues(*jtmp,3,index,3,index,A.m);
  A=A/area;

  A.m[1] +=  -mt->e;   //dfun(0)/dn(i)
  A.m[2] +=   mt->e;   //dfun(0)/dp(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[4] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[8] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[4] += -1/ODE_F.dt;
      A.m[8] += -1/ODE_F.dt;
    }
  }

  MatSetValues(*jac,3,index,3,index,A.m,INSERT_VALUES);
  
}


void SMCZone::J1E_ddm_ombc(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  PetscScalar    A[9];
  PetscInt       index[3];
  int om_equ;
  int equ_num = 3;
  int size = pzone->davcell.size();
  PetscScalar e  =  mt->e;
  PetscScalar kb =  mt->kb;
  PetscScalar Vt =  kb*fs[i].T/e;
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  for(int j=0;j<electrode.size();j++)
  if(electrode[j]==pcell->bc_index-1)   {om_equ=j;break;}
  index[0] = zofs[zone_index]+3*i+0;
  index[1] = zofs[zone_index]+3*i+1;
  index[2] = zofs[zone_index]+3*i+2;
  if(Fermi)  //Fermi
  {
    mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
    PetscScalar Nc  = mt->band->Nc(fs[i].T);
    PetscScalar Nv  = mt->band->Nv(fs[i].T);
    PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) );
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) + mt->band->Eg(fs[i].T));
    PetscScalar phin = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar phip = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    PetscScalar detan_dVi =  1.0/Vt;
    PetscScalar detap_dVi = -1.0/Vt;
    A[0] = Nc*fermi_mhalf(etan)*detan_dVi - Nv*fermi_mhalf(etap)*detap_dVi;  
    A[1] = 0;  
    A[2] = 0;
    A[3] = -Nc*fermi_mhalf(etan)*detan_dVi;  
    A[4] = 1;  
    A[5] = 0;
    A[6] = -Nv*fermi_mhalf(etap)*detap_dVi;  
    A[7] = 0;  
    A[8] = 1;
    MatSetValues(*jac,3,index,3,index,A,INSERT_VALUES);
    //dfermi_half(eta)_dVapp = - dfermi_half(eta)_dVi
    MatSetValue(*jac,zofs[zone_index]+3*i+0,zofs[zone_index]+equ_num*size+om_equ,-A[0],INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+3*i+1,zofs[zone_index]+equ_num*size+om_equ,-A[3],INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+3*i+2,zofs[zone_index]+equ_num*size+om_equ,-A[6],INSERT_VALUES);
  }
  else //Boltzmann
  {
    A[0] = 1;  A[1] = 0;  A[2] = 0;
    A[3] = 0;  A[4] = 1;  A[5] = 0;
    A[6] = 0;  A[7] = 0;  A[8] = 1;
    MatSetValues(*jac,3,index,3,index,A,INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+3*i,zofs[zone_index]+equ_num*size+om_equ,-1,INSERT_VALUES);
  }
}


void SMCZone::J1E_ddm_stkbc(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  Mat3    A;
  PetscInt       index[3],col[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int  z = zone_index;
  int  equ_num = 3;
  int  size = pzone->davcell.size();
  PetscScalar VB = pbc->WorkFunction-aux[i].affinity;
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  PetscScalar Vi = x[zofs[z]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i

  PetscScalar area = pcell->area;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fp_dpi = 0;
  //--------------------------------
  index[0] = zofs[z]+3*i+0;
  index[1] = zofs[z]+3*i+1;
  index[2] = zofs[z]+3*i+2;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    col[0] = zofs[z]+3*nb+0;
    col[1] = zofs[z]+3*nb+1;
    col[2] = zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,index,3,col,A.m);
    A.m[0] = 0.0;
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      d_Fn_dni += mt->band->pdSchottyJsn_pdn(ni,fs[i].T,VB-deltaVB)*0.5*pcell->ilen[j];
      d_Fp_dpi += mt->band->pdSchottyJsp_pdp(pi,fs[i].T,VB+deltaVB)*0.5*pcell->ilen[j];
    }
    MatSetValues(*jac,3,index,3,col,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,3,index,3,index,A.m);
  A=A/area;
  A.m[0] =  1.0;                      //dfun(0)/dP(i)
  A.m[1] =  0.0;                      //dfun(0)/dn(i)
  A.m[2] =  0.0;                      //dfun(0)/dp(i)

  A.m[4] +=   d_Fn_dni/area;          //dfun(1)/dn(i)
  A.m[8] +=  -d_Fp_dpi/area;          //dfun(2)/dp(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[4] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[8] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[4] += -1/ODE_F.dt;
      A.m[8] += -1/ODE_F.dt;
    }
  }
  MatSetValues(*jac,3,index,3,index,A.m,INSERT_VALUES);

  int stk_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
      {stk_equ=j;break;}
  MatSetValue(*jac,zofs[zone_index]+3*i,zofs[zone_index]+equ_num*size+stk_equ,-1,INSERT_VALUES);
}


void SMCZone::J1E_ddm_insulator_gate(int i,PetscScalar *x, Mat *jac, Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  Mat3    A;
  PetscInt       index[3],col[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int  equ_num = 3;
  int  size = pzone->davcell.size();

  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i

  PetscScalar area = pcell->area;
  PetscScalar d_grad_P_dVi = 0;
  PetscScalar d_grad_P_dVapp = 0;
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  int ins_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
      {ins_equ=j;break;}
  //--------------------------------
  index[0] = zofs[z]+3*i+0;
  index[1] = zofs[z]+3*i+1;
  index[2] = zofs[z]+3*i+2;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    col[0] = zofs[z]+3*nb+0;
    col[1] = zofs[z]+3*nb+1;
    col[2] = zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,index,3,col,A.m);
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar s=eps_ox/Thick;
      d_grad_P_dVi   += -0.5*pcell->ilen[j]*0.25*s*3/area;
      d_grad_P_dVapp +=  0.5*pcell->ilen[j]*s/area;
      A.m[0]         += -0.5*pcell->ilen[j]*0.25*s/area;
    }
    MatSetValues(*jac,3,index,3,col,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,3,index,3,index,A.m);
  A=A/area;
  A.m[0] +=  d_grad_P_dVi;             //dfun(0)/dP(i)
  A.m[1] +=  -mt->e;   //dfun(0)/dn(i)
  A.m[2] +=   mt->e;   //dfun(0)/dp(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[4] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[8] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[4] += -1/ODE_F.dt;
      A.m[8] += -1/ODE_F.dt;
    }
  }
  MatSetValues(*jac,3,index,3,index,A.m,INSERT_VALUES);
  MatSetValue(*jac,zofs[z]+3*i,zofs[z]+equ_num*size+ins_equ,d_grad_P_dVapp,INSERT_VALUES);
}


void SMCZone::J1E_ddm_interface(int i,PetscScalar *x,Mat *jac, Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                               ISZone *pz, int n)
{
  Mat3    A;
  PetscInt       index[3],col[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int z = zone_index;
  PetscScalar Vi = x[zofs[z]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i
  PetscScalar area = pcell->area;
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar d_grad_P_dVi = 0;
  //--------------------------------
  index[0] = zofs[z]+3*i+0;
  index[1] = zofs[z]+3*i+1;
  index[2] = zofs[z]+3*i+2;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];

    PetscScalar Vj = x[zofs[z]+3*nb+0];     //potential of nb node
    col[0] = zofs[z]+3*nb+0;
    col[1] = zofs[z]+3*nb+1;
    col[2] = zofs[z]+3*nb+2;
    //-------------------------------------
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    MatGetValues(*jtmp,3,index,3,col,A.m);
    A=A/area;
    //-------------------------------------
    //df(x)i/dx(r)
    A.m[0] =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    A.m[1] =  0;
    A.m[2] =  0;
    MatSetValues(*jac,3,index,3,col,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,3,index,3,index,A.m);
  A=A/area;
  
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  int column = zofs[pz->pzone->zone_index]+n;
  for(int j=0;j<ncell->nb_num;j++)
  {
    int    nb = ncell->nb_array[j];
    PetscScalar Vr_n = x[zofs[pz->pzone->zone_index]+nb];     //potential of nb node
    d_grad_P_dVi += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    PetscScalar value = pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[z]+3*i+0,zofs[pz->pzone->zone_index]+nb,value,INSERT_VALUES);
  }
  
  A.m[0] =  d_grad_P_dVi;  //dfun(0)/dP(i)
  A.m[1] =  -mt->e*area;   //dfun(0)/dn(i)
  A.m[2] =   mt->e*area;   //dfun(0)/dp(i)
  
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[4] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[8] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[4] += -1/ODE_F.dt;
      A.m[8] += -1/ODE_F.dt;
    }
  }

  MatSetValues(*jac,3,index,3,index,A.m,INSERT_VALUES);

}


void SMCZone::J1E_ddm_homojunction(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                               SMCZone *pz, int n)
{
  Mat3           A1,A2,AJ;
  PetscInt       index[3],indexn[3],col[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int z = zone_index;
  PetscScalar e  =  mt->e;
  PetscScalar Vi = x[zofs[z]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i

  PetscScalar d_grad_P_dVi = 0;
  //--------------------------------
  index[0] = zofs[z]+3*i+0;
  index[1] = zofs[z]+3*i+1;
  index[2] = zofs[z]+3*i+2;
  //--------------------------------

  if(zone_index > pz->pzone->zone_index)
  {
        Set_Mat3_I(AJ);
        MatSetValues(*jac,3,index,3,index,AJ.m,INSERT_VALUES);
        AJ = -1.0*AJ;
        col[0] = zofs[pz->zone_index]+3*n+0;
        col[1] = zofs[pz->zone_index]+3*n+1;
        col[2] = zofs[pz->zone_index]+3*n+2;
        MatSetValues(*jac,3,index,3,col,AJ.m,INSERT_VALUES);
        return;
  }

  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[z]+3*nb+0];     //potential of nb node
    //-------------------------------------
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    col[0] = zofs[z]+3*nb+0;
    col[1] = zofs[z]+3*nb+1;
    col[2] = zofs[z]+3*nb+2;

    MatGetValues(*jtmp,3,index,3,col,A1.m);
    //-------------------------------------
    //df(x)i/dx(r)
    A1.m[0] =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];   //dfun(0)/dP(r)
    A1.m[1] =  0;                                       //dfun(0)/dn(r)
    A1.m[2] =  0;                                       //dfun(0)/dp(r)
    MatSetValues(*jac,3,index,3,col,A1.m,INSERT_VALUES);
  }

  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  indexn[0] = zofs[pz->zone_index]+3*n+0;
  indexn[1] = zofs[pz->zone_index]+3*n+1;
  indexn[2] = zofs[pz->zone_index]+3*n+2;
  for(int j=0;j<ncell->nb_num;j++)
  {
    int  nb = ncell->nb_array[j];
    PetscScalar Vj = x[zofs[pz->zone_index]+3*nb+0];     //potential of nb node
    //-------------------------------------
    d_grad_P_dVi += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    col[0] = zofs[pz->zone_index]+3*nb+0;
    col[1] = zofs[pz->zone_index]+3*nb+1;
    col[2] = zofs[pz->zone_index]+3*nb+2;
    MatGetValues(*jtmp,3,indexn,3,col,A2.m);
    //-------------------------------------
    //df(x)i/dx(r)
     A2.m[0] =  pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];     //dfun(0)/dP(r)
     A2.m[1] =  0;                                                     //dfun(0)/dn(r)
     A2.m[2] =  0;                                                     //dfun(0)/dp(r)
    MatSetValues(*jac,3,index,3,col,A2.m,INSERT_VALUES);
  }

  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  PetscScalar area = pcell->area + ncell->area;
  MatGetValues(*jtmp,3,index,3,index,A1.m);
  MatGetValues(*jtmp,3,indexn,3,indexn,A2.m);
  AJ=A1+A2;

  AJ.m[0] =  d_grad_P_dVi;    //dfun(0)/dP(i)
  AJ.m[1] =  -e*area;         //dfun(0)/dn(i)
  AJ.m[2] =   e*area;         //dfun(0)/dp(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      AJ.m[4] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*area;
      AJ.m[8] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*area;
    }
    else //first order
    {
      AJ.m[4] += -1/ODE_F.dt*area;
      AJ.m[8] += -1/ODE_F.dt*area;
    }
  }
  MatSetValues(*jac,3,index,3,index,AJ.m,INSERT_VALUES);
}


void SMCZone::J1E_ddm_heterojunction(int i,PetscScalar *x,Mat *jac, Mat *jtmp,ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                               SMCZone *pz, int n)
{
  Mat3    A;
  PetscInt       index[3],col[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  int z = zone_index;
  PetscScalar k = mt->kb;
  PetscScalar e  =  mt->e;
  PetscScalar Vi = x[zofs[z]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i
  PetscScalar Ti = fs[i].T;
  PetscScalar Eci =  -e*(Vi + aux[i].affinity + mt->kb*fs[i].T/e*log(aux[i].Nc));
  PetscScalar Evi =  -e*(Vi + aux[i].affinity + aux[i].Eg/e - mt->kb*fs[i].T/e*log(aux[i].Nv));
  PetscScalar area = pcell->area;
  PetscScalar d_grad_P_dVi = 0;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fp_dpi = 0;
  //--------------------------------
  index[0] = zofs[z]+3*i+0;
  index[1] = zofs[z]+3*i+1;
  index[2] = zofs[z]+3*i+2;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];

    PetscScalar Vj = x[zofs[z]+3*nb+0];     //potential of nb node
    PetscScalar nj = x[zofs[z]+3*nb+1];     //electron density of nb node
    PetscScalar pj = x[zofs[z]+3*nb+2];     //hole density of nb node
    //-------------------------------------
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];

    col[0] = zofs[z]+3*nb+0;
    col[1] = zofs[z]+3*nb+1;
    col[2] = zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,index,3,col,A.m);
    A=A/area;
    //-------------------------------------
    //df(x)i/dx(r)
    if(pzone->zone_index < pz->pzone->zone_index)
    {
        A.m[0] =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];//dfun(0)/dP(r)
        A.m[1] =  0;                                     //dfun(0)/dn(r)
        A.m[2] =  0;                                     //dfun(0)/dp(r)
    }
    else
    {
        A.m[0] =  0;                                     //dfun(0)/dP(r)
        A.m[1] =  0;                                     //dfun(0)/dn(r)
        A.m[2] =  0;                                     //dfun(0)/dp(r)
    }

    //Thermal emit current
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Ecj = -e*(Vi + pz->aux[n].affinity + k*pz->fs[n].T/e*log(pz->aux[n].Nc));
      PetscScalar Evj = -e*(Vi + pz->aux[n].affinity - k*pz->fs[n].T/e*log(pz->aux[n].Nv) + pz->aux[n].Eg/e);
      MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
      if(Ecj > Eci)
      {
        PetscScalar pm = pz->mt->band->EffecElecMass(Ti)/mt->band->EffecElecMass(Ti);
        d_Fn_dni += 2*e*pm*mt->band->ThermalVn(Ti)*exp(-(Ecj-Eci)/(k*Ti))*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fn_dnd = -2*e*pz->mt->band->ThermalVn(Ti)*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+3*i+1,zofs[pz->pzone->zone_index]+3*n+1,d_Fn_dnd,ADD_VALUES);
      }
      else
      {
        PetscScalar pm = mt->band->EffecElecMass(Ti)/pz->mt->band->EffecElecMass(Ti);
        d_Fn_dni += -2*e*mt->band->ThermalVn(Ti)*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fn_dnd = 2*e*pm*pz->mt->band->ThermalVn(Ti)*exp(-(Eci-Ecj)/(k*Ti))*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+3*i+1,zofs[pz->pzone->zone_index]+3*n+1,d_Fn_dnd,ADD_VALUES);
      }

      if(Evj > Evi)
      {
        PetscScalar pm = pz->mt->band->EffecHoleMass(Ti)/mt->band->EffecHoleMass(Ti);
        d_Fp_dpi += -2*e*pm*mt->band->ThermalVp(Ti)*exp(-(Evj-Evi)/(k*Ti))*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fp_dpd = 2*e*pz->mt->band->ThermalVp(Ti)*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+3*i+2,zofs[pz->pzone->zone_index]+3*n+2,-d_Fp_dpd,ADD_VALUES);
      }
      else
      {
        PetscScalar pm = mt->band->EffecHoleMass(Ti)/pz->mt->band->EffecHoleMass(Ti);
        d_Fp_dpi += 2*e*mt->band->ThermalVp(Ti)*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fp_dpd = -2*e*pm*pz->mt->band->ThermalVp(Ti)*exp(-(Evi-Evj)/(k*Ti))*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+3*i+2,zofs[pz->pzone->zone_index]+3*n+2,-d_Fp_dpd,ADD_VALUES);
      }
      MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
    }
    MatSetValues(*jac,3,index,3,col,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,3,index,3,index,A.m);
  A=A/area;
  if(pzone->zone_index < pz->pzone->zone_index)
  {
    A.m[0] =  d_grad_P_dVi;    //dfun(0)/dP(i)
    A.m[1] =  -e*pcell->area;  //dfun(0)/dn(i)
    A.m[2] =   e*pcell->area;  //dfun(0)/dp(i)
    for(int j=0;j<ncell->nb_num;j++)
    {
      int    nb = ncell->nb_array[j];
      A.m[0] += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
      PetscScalar value = pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
      MatSetValue(*jac,zofs[z]+3*i+0,zofs[pz->pzone->zone_index]+3*nb+0,value,INSERT_VALUES);
    }
    MatSetValue(*jac,zofs[z]+3*i+0,zofs[pz->pzone->zone_index]+3*n+1,-e*ncell->area,INSERT_VALUES);
    MatSetValue(*jac,zofs[z]+3*i+0,zofs[pz->pzone->zone_index]+3*n+2,e*ncell->area,INSERT_VALUES);
  }
  else
  {
    A.m[0] =   1;    //dfun(0)/dP(i)
    A.m[1] =   0;    //dfun(0)/dn(i)
    A.m[2] =   0;    //dfun(0)/dp(i)
    MatSetValue(*jac,zofs[z]+3*i+0,zofs[pz->pzone->zone_index]+3*n+0,-1,INSERT_VALUES);
  }

  A.m[4] +=   d_Fn_dni;  //dfun(1)/dn(i)
  A.m[8] +=  -d_Fp_dpi;  //dfun(2)/dp(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[4] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[8] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[4] += -1/ODE_F.dt;
      A.m[8] += -1/ODE_F.dt;
    }
  }

  MatSetValues(*jac,3,index,3,index,A.m,INSERT_VALUES);

}


void SMCZone::J1E_om_electrode(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(bc_index));
  PetscScalar    A1[3],A2[3],J[3];
  PetscInt       index[3],col[3];
  int matrix_row = zofs[zone_index]+equ_num*size+i;
  int connect_index;
  if(pbc->inner_connect!=-1)
  {
        int connect_to = pbc->inner_connect;
        int connect_zone = pbc->connect_zone;
        int connect_elec=-1;
        SMCZone * pz = dynamic_cast <SMCZone * > (pbc->pzonedata);
        for(int j=0;j<pz->electrode.size();j++)
        {
                if(bc[pz->electrode[j]].BCType==OhmicContact)
                {
                  OhmicBC *om_bc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pz->electrode[j]));
                  if(om_bc->inner_connect==bc_index)              connect_elec=j;
                }  
        }
        connect_index=zofs[connect_zone]+equ_num*pz->pzone->davcell.size()+connect_elec;
  }      
  PetscScalar R = pbc->R;
  PetscScalar C = pbc->C;
  PetscScalar L = pbc->L;
  PetscScalar e = mt->e;
  PetscScalar dEc_dV = -e;
  PetscScalar dEv_dV = -e;
  
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar dJdisp_dVi=0,dJdisp_dVr=0;
    index[0] = zofs[zone_index]+3*node+0;
    index[1] = zofs[zone_index]+3*node+1;
    index[2] = zofs[zone_index]+3*node+2;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      col[0] = zofs[zone_index]+3*nb+0;
      col[1] = zofs[zone_index]+3*nb+1;
      col[2] = zofs[zone_index]+3*nb+2;

      MatGetValues(*jtmp,1,&index[1],3,col,A1);
      MatGetValues(*jtmp,1,&index[2],3,col,A2);

      //for displacement current
      dJdisp_dVi += aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];
      dJdisp_dVr = -aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];
      
      //if the electrode connect to another electrode
      if(pbc->inner_connect!=-1)
      {
        int connect_to = pbc->inner_connect;
        if(bc_index < connect_to)
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
          MatSetValues(*jac,1,&connect_index,3,col,J,ADD_VALUES);
        }
        else
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
          MatSetValues(*jac,1,&matrix_row,3,col,J,ADD_VALUES);
        }                                      
      } 
      //the electrode has an external voltage circuit
      else if(pbc->electrode_type==VoltageBC)
      {
#ifdef        _EXP_R_  
        if(R<1.0)
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth*(L/ODE_F.dt+R);
          J[1]=(A1[1]-A2[1])*DeviceDepth*(L/ODE_F.dt+R);
          J[2]=(A1[2]-A2[2])*DeviceDepth*(L/ODE_F.dt+R);
        }
        else
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
        }  
#else
        J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth*(L/ODE_F.dt+R);
        J[1]=(A1[1]-A2[1])*DeviceDepth*(L/ODE_F.dt+R);
        J[2]=(A1[2]-A2[2])*DeviceDepth*(L/ODE_F.dt+R);
#endif        
        MatSetValues(*jac,1,&matrix_row,3,col,J,ADD_VALUES);
      }
      //the electrode has an external current circuit
      else if(pbc->electrode_type==CurrentBC)
      {
        J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
        J[1]=(A1[1]-A2[1])*DeviceDepth;
        J[2]=(A1[2]-A2[2])*DeviceDepth;
        MatSetValues(*jac,1,&matrix_row,3,col,J,ADD_VALUES);
      }
    }

    MatGetValues(*jtmp,1,&index[1],3,index,A1);
    MatGetValues(*jtmp,1,&index[2],3,index,A2);
    
    if(pbc->inner_connect!=-1)
    {
        int connect_to = pbc->inner_connect;
        if(bc_index < connect_to)
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
          MatSetValues(*jac,1,&connect_index,3,index,J,ADD_VALUES);
        }
        else
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
          MatSetValues(*jac,1,&matrix_row,3,index,J,ADD_VALUES);
        }                                      
    }   
    else if(pbc->electrode_type==VoltageBC)
    {
#ifdef    _EXP_R_
      if(R<1.0)
      {
        J[0]=(A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth*(L/ODE_F.dt+R);
        J[1]=(A1[1]-A2[1])*DeviceDepth*(L/ODE_F.dt+R);
        J[2]=(A1[2]-A2[2])*DeviceDepth*(L/ODE_F.dt+R);
      }
      else
      {
        J[0]=(A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth;
        J[1]=(A1[1]-A2[1])*DeviceDepth;
        J[2]=(A1[2]-A2[2])*DeviceDepth;
      } 
#else
       J[0]=(A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth*(L/ODE_F.dt+R);
       J[1]=(A1[1]-A2[1])*DeviceDepth*(L/ODE_F.dt+R);
       J[2]=(A1[2]-A2[2])*DeviceDepth*(L/ODE_F.dt+R);
#endif          
      MatSetValues(*jac,1,&matrix_row,3,index,J,ADD_VALUES);
    }
    else if(pbc->electrode_type==CurrentBC)
    {
      J[0]=(A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth;
      J[1]=(A1[1]-A2[1])*DeviceDepth;
      J[2]=(A1[2]-A2[2])*DeviceDepth;
      MatSetValues(*jac,1,&matrix_row,3,index,J,ADD_VALUES);
    }

  }

  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  
  if(pbc->inner_connect!=-1)
  {
        int connect_to = pbc->inner_connect;
        if(bc_index < connect_to)
        {
          MatSetValue(*jac,matrix_row,matrix_row,1.0,ADD_VALUES);
        }
        else
        {
          MatSetValue(*jac,connect_index,matrix_row,-1.0,ADD_VALUES);
        }                                      
  }
  else if(pbc->electrode_type==VoltageBC)
  {
#ifdef _EXP_R_  
    if(R<1.0)
      MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP
    else
      MatSetValue(*jac,matrix_row,matrix_row,1/(L/ODE_F.dt+R)+C/ODE_F.dt,INSERT_VALUES); //dJ/dP   
#else
    MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP
#endif      
  }  
}


void SMCZone::J1E_stk_electrode(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc.Get_pointer(bc_index));
  int matrix_row = zofs[zone_index]+equ_num*size+i;
  PetscScalar R = pbc->R;
  PetscScalar C = pbc->C;
  PetscScalar L = pbc->L;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
      int node = bc[bc_index].psegment->node_array[j];
      const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
      PetscScalar VB = pbc->WorkFunction-aux[node].affinity;
      PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[node].eps,sqrt(aux[node].Ex*aux[node].Ex+aux[node].Ey*aux[node].Ey));
      PetscScalar ni = x[zofs[zone_index]+3*node+1];     //electron density of node i
      PetscScalar pi = x[zofs[zone_index]+3*node+2];     //hole density of node i
      PetscScalar dJ_dni = -mt->band->pdSchottyJsn_pdn(ni,fs[node].T,VB-deltaVB)*0.5*pcell->ilen[0]
                           -mt->band->pdSchottyJsn_pdn(ni,fs[node].T,VB-deltaVB)*0.5*pcell->ilen[pcell->nb_num-1];
      PetscScalar dJ_dpi = -mt->band->pdSchottyJsp_pdp(pi,fs[node].T,VB+deltaVB)*0.5*pcell->ilen[0]
                           -mt->band->pdSchottyJsp_pdp(pi,fs[node].T,VB+deltaVB)*0.5*pcell->ilen[pcell->nb_num-1];

      //for displacement current
      MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
      for(int k=0;k<pcell->nb_num;k++)
      {
        int    nb = pcell->nb_array[k];
        PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
        PetscScalar dI_dVi =  DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
        PetscScalar dI_dVr = -DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
        if(pbc->electrode_type==VoltageBC)
        {
          MatSetValue(*jac,matrix_row,zofs[zone_index]+3*node+0,dI_dVi*(L/ODE_F.dt+R),ADD_VALUES);
          MatSetValue(*jac,matrix_row,zofs[zone_index]+3*nb+0,dI_dVr*(L/ODE_F.dt+R),ADD_VALUES);
        }
        else if(pbc->electrode_type==CurrentBC)
        {
          MatSetValue(*jac,matrix_row,zofs[zone_index]+3*node+0,dI_dVi,ADD_VALUES);
          MatSetValue(*jac,matrix_row,zofs[zone_index]+3*nb+0,dI_dVr,ADD_VALUES);
        }
      }
      MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

      if(pbc->electrode_type==VoltageBC)
      {
        MatSetValue(*jac,matrix_row,zofs[zone_index]+3*node+1,DeviceDepth*dJ_dni*(L/ODE_F.dt+R),INSERT_VALUES);
        MatSetValue(*jac,matrix_row,zofs[zone_index]+3*node+2,DeviceDepth*dJ_dpi*(L/ODE_F.dt+R),INSERT_VALUES);
      }
      else if(pbc->electrode_type==CurrentBC)
      {
         MatSetValue(*jac,matrix_row,zofs[zone_index]+3*node+1,DeviceDepth*dJ_dni,INSERT_VALUES);
         MatSetValue(*jac,matrix_row,zofs[zone_index]+3*node+2,DeviceDepth*dJ_dpi,INSERT_VALUES);
      }
  }
  if(pbc->electrode_type==VoltageBC)
    MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP
}

void SMCZone::J1E_ins_electrode(int i,PetscScalar *x,Mat *jac, Mat *jtmp,ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  InsulatorContactBC *pbc = dynamic_cast <InsulatorContactBC * > (bc.Get_pointer(bc_index));
  int matrix_row = zofs[zone_index]+equ_num*size+i;
  PetscScalar R = pbc->R;
  PetscScalar C = pbc->C;
  PetscScalar L = pbc->L;

  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
      int node = bc[bc_index].psegment->node_array[j];
      const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
      //for displacement current
      for(int k=0;k<pcell->nb_num;k++)
      {
        int    nb = pcell->nb_array[k];
        PetscScalar dI_dVi =  DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
        PetscScalar dI_dVr = -DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
        if(pbc->electrode_type==VoltageBC)
        {
          MatSetValue(*jac,matrix_row,zofs[zone_index]+3*node+0,dI_dVi*(L/ODE_F.dt+R),ADD_VALUES);
          MatSetValue(*jac,matrix_row,zofs[zone_index]+3*nb+0,dI_dVr*(L/ODE_F.dt+R),ADD_VALUES);
        }
      }
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  if(pbc->electrode_type==VoltageBC)
    MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP

}


void SMCZone::F1E_om_electrode_Trace(int electrode_index, PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,
         Vec & x, Mat *jtmp, Mat *jac,ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
{
  int equ_num = 3;
  int bc_index = electrode[electrode_index];
  int size = pzone->davcell.size();
  int extern_equ = zofs[zone_index]+equ_num*size+electrode_index;
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(bc_index));
  pdI_pdV = 0;
  PetscScalar e = mt->e;
  VecZeroEntries(pdI_pdw);
  VecZeroEntries(pdF_pdV);
  PetscScalar * xx;
  PetscScalar * apdI_pdw;
  PetscScalar * apdF_pdV;
  VecGetArray(x,&xx);
  VecGetArray(pdI_pdw,&apdI_pdw);
  VecGetArray(pdF_pdV,&apdF_pdV);
  PetscScalar    A1[3],A2[3];
  PetscInt       index[3],col[3];
  //delete electrode current equation, omit the effect of external resistance 
  MatSetValuesRow(*jac,extern_equ,apdF_pdV);
  MatSetValue(*jac,extern_equ,extern_equ,1.0,INSERT_VALUES);
  
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = xx[zofs[zone_index]+3*node+0];     //potential of node i
    PetscScalar dJdisp_dVi=0,dJdisp_dVr=0;
    index[0] = zofs[zone_index]+3*node+0;
    index[1] = zofs[zone_index]+3*node+1;
    index[2] = zofs[zone_index]+3*node+2;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = xx[zofs[zone_index]+3*nb+0];     //potential of nb node
      col[0] = zofs[zone_index]+3*nb+0;
      col[1] = zofs[zone_index]+3*nb+1;
      col[2] = zofs[zone_index]+3*nb+2;

      MatGetValues(*jtmp,1,&index[1],3,col,A1);
      MatGetValues(*jtmp,1,&index[2],3,col,A2);

      //for displacement current
      dJdisp_dVi += aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];
      dJdisp_dVr = -aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];
  
      //
      apdI_pdw[col[0]] += (A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
      apdI_pdw[col[1]] += (A1[1]-A2[1])*DeviceDepth;
      apdI_pdw[col[2]] += (A1[2]-A2[2])*DeviceDepth;
    }
    MatGetValues(*jtmp,1,&index[1],3,index,A1);
    MatGetValues(*jtmp,1,&index[2],3,index,A2);
    
    apdI_pdw[index[0]] += (A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth;
    apdI_pdw[index[1]] += (A1[1]-A2[1])*DeviceDepth;
    apdI_pdw[index[2]] += (A1[2]-A2[2])*DeviceDepth;
    
    apdF_pdV[index[0]] = 1.0;
    pdI_pdV = 0;
  }
  
  MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);
  VecRestoreArray(x,&xx);
  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}


void SMCZone::F1E_stk_electrode_Trace(int electrode_index,  PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,
         Vec & x, Mat *jtmp, Mat *jac,ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
{
  int equ_num = 3;
  int bc_index = electrode[electrode_index];
  int size = pzone->davcell.size();
  int extern_equ = zofs[zone_index]+equ_num*size+electrode_index;
  SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc.Get_pointer(bc_index));
  pdI_pdV = 0;
  PetscScalar e = mt->e;
  VecZeroEntries(pdI_pdw);
  VecZeroEntries(pdF_pdV);
  PetscScalar * xx;
  PetscScalar * apdI_pdw;
  PetscScalar * apdF_pdV;
  VecGetArray(x,&xx);
  VecGetArray(pdI_pdw,&apdI_pdw);
  VecGetArray(pdF_pdV,&apdF_pdV);
  //delete electrode current equation, omit the effect of external resistance 
  MatSetValuesRow(*jac,extern_equ,apdF_pdV);
  MatSetValue(*jac,extern_equ,extern_equ,1.0,INSERT_VALUES);
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
      int node = bc[bc_index].psegment->node_array[j];
      const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
      PetscScalar VB = pbc->WorkFunction-aux[node].affinity;
      PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[node].eps,sqrt(aux[node].Ex*aux[node].Ex+aux[node].Ey*aux[node].Ey));
      PetscScalar ni = xx[zofs[zone_index]+3*node+1];     //electron density of node i
      PetscScalar pi = xx[zofs[zone_index]+3*node+2];     //hole density of node i
      PetscScalar dJ_dni = -mt->band->pdSchottyJsn_pdn(ni,fs[node].T,VB-deltaVB)*0.5*pcell->ilen[0]
                           -mt->band->pdSchottyJsn_pdn(ni,fs[node].T,VB-deltaVB)*0.5*pcell->ilen[pcell->nb_num-1];
      PetscScalar dJ_dpi = -mt->band->pdSchottyJsp_pdp(pi,fs[node].T,VB+deltaVB)*0.5*pcell->ilen[0]
                           -mt->band->pdSchottyJsp_pdp(pi,fs[node].T,VB+deltaVB)*0.5*pcell->ilen[pcell->nb_num-1];

      //for displacement current
      for(int k=0;k<pcell->nb_num;k++)
      {
        int    nb = pcell->nb_array[k];
        PetscScalar dI_dVi =  DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
        PetscScalar dI_dVr = -DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
        apdI_pdw[zofs[zone_index]+3*node+0] += dI_dVi;
        apdI_pdw[zofs[zone_index]+3*nb+0]   += dI_dVr;
      }
      apdI_pdw[zofs[zone_index]+3*node+1] = DeviceDepth*dJ_dni;
      apdI_pdw[zofs[zone_index]+3*node+2] = DeviceDepth*dJ_dpi;
      
      apdF_pdV[zofs[zone_index]+3*node+0] = 1.0;
      pdI_pdV = 0;
 }
  MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);
  VecRestoreArray(x,&xx);
  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}


void SMCZone::F1E_efield_update(PetscScalar *x,vector<int> & zofs, DABC &bc, vector<BZoneData *> zonedata)
{
  PetscScalar P_x=0,P_y=0,w=0;
  //calculate Ex Ey with least-squares gradient construction
  for(int i=0; i<pzone->davcell.size();i++)
  {
    P_x=0,P_y=0,w=0;
    VoronoiCell *pcell = pzone->davcell.GetPointer(i);
    for(int k=0;k<pcell->nb_num;k++)
    {
      const VoronoiCell *ncell = pzone->davcell.GetPointer(pcell->nb_array[k]);
      PetscScalar dx = ncell->x - pcell->x;
      PetscScalar dy = ncell->y - pcell->y;
      PetscScalar dP = x[zofs[zone_index]+3*pcell->nb_array[k]+0] - x[zofs[zone_index]+3*i+0];
      w=1.0/sqrt(dx*dx+dy*dy);
      P_x+=w*w*dP*dx;
      P_y+=w*w*dP*dy;
    }
    if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
    {
      InsulatorInterfaceBC *pbc;
      pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
      int n_zone = pbc->pinterface->Find_neighbor_zone_index(zone_index);
      int n_node = pbc->pinterface->Find_neighbor_node_index(zone_index,i);
      ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
      const VoronoiCell* dcell = pz->pzone->davcell.GetPointer(n_node);
      for(int k=1;k<dcell->nb_num-1;k++)
      {
        const VoronoiCell *ncell = pz->pzone->davcell.GetPointer(dcell->nb_array[k]);
        PetscScalar dx = ncell->x - dcell->x;
        PetscScalar dy = ncell->y - dcell->y;
        PetscScalar dP = x[zofs[n_zone]+dcell->nb_array[k]]- x[zofs[n_zone]+n_node];
        w=1.0/sqrt(dx*dx+dy*dy);
        P_x+=w*w*dP*dx;
        P_y+=w*w*dP*dy;
      }
    }
    aux[i].Ex = -(pcell->sc*P_x-pcell->sb*P_y);
    aux[i].Ey = -(pcell->sa*P_y-pcell->sb*P_x);
  }
}

