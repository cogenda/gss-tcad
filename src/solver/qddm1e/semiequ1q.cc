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
/*  Last update: July 20, 2007                                               */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "zonedata.h"
#include "jflux1q.h"
#include "vec5.h"
#define _FLUX1_
#define _HMOB_
#define _EXP_R_


void SMCZone::F1Q_Tri_qddm(Tri *ptri,PetscScalar *x,PetscScalar *f, vector<int> &zofs)
{
  int A = ptri->node[0];
  int B = ptri->node[1];
  int C = ptri->node[2];
  PetscScalar xa = pzone->danode[ptri->node[0]].x;
  PetscScalar xb = pzone->danode[ptri->node[1]].x;
  PetscScalar xc = pzone->danode[ptri->node[2]].x;
  PetscScalar ya = pzone->danode[ptri->node[0]].y;
  PetscScalar yb = pzone->danode[ptri->node[1]].y;
  PetscScalar yc = pzone->danode[ptri->node[2]].y;
  PetscScalar tri_area = ptri->area;
  PetscScalar Sax = (xc-xb)/ptri->edge_len[0];
  PetscScalar Say = (yc-yb)/ptri->edge_len[0];
  PetscScalar Sbx = (xa-xc)/ptri->edge_len[1];
  PetscScalar Sby = (ya-yc)/ptri->edge_len[1];
  PetscScalar Scx = (xb-xa)/ptri->edge_len[2];
  PetscScalar Scy = (yb-ya)/ptri->edge_len[2];
  PetscScalar Wab = ptri->d[1]/((ptri->d[1]+ptri->d[2])*sin(ptri->angle[2])*sin(ptri->angle[2]));
  PetscScalar Wac = ptri->d[2]/((ptri->d[1]+ptri->d[2])*sin(ptri->angle[1])*sin(ptri->angle[1]));
  PetscScalar Wbc = ptri->d[2]/((ptri->d[0]+ptri->d[2])*sin(ptri->angle[0])*sin(ptri->angle[0]));
  PetscScalar Wba = ptri->d[0]/((ptri->d[0]+ptri->d[2])*sin(ptri->angle[2])*sin(ptri->angle[2]));
  PetscScalar Wca = ptri->d[0]/((ptri->d[0]+ptri->d[1])*sin(ptri->angle[1])*sin(ptri->angle[1]));
  PetscScalar Wcb = ptri->d[1]/((ptri->d[0]+ptri->d[1])*sin(ptri->angle[0])*sin(ptri->angle[0]));
  PetscScalar Ma  = -cos(ptri->angle[0]);
  PetscScalar Mb  = -cos(ptri->angle[1]);
  PetscScalar Mc  = -cos(ptri->angle[2]);

  PetscScalar e  =  mt->e;
  PetscScalar hbar = mt->hbar;
  PetscScalar T  = (fs[A].T+fs[B].T+fs[C].T)/3.0;
  PetscScalar kb = mt->kb;
  PetscScalar Vt = mt->kb*T/e;
  PetscScalar Eg = mt->band->Eg(T);
  PetscScalar me =  mt->band->EffecElecMass(T);
  PetscScalar mh =  mt->band->EffecHoleMass(T);
  PetscScalar gamman = mt->band->Gamman(); 
  PetscScalar gammap = mt->band->Gammap();
  PetscScalar bn = gamman*hbar*hbar/(6*e*me);
  PetscScalar bp = gammap*hbar*hbar/(6*e*mh);
  PetscScalar S;
  
  PetscScalar Va   = x[zofs[zone_index]+5*A+0];     //potential of node A
  PetscScalar na   = x[zofs[zone_index]+5*A+1];     //electron density of node A
  PetscScalar pa   = x[zofs[zone_index]+5*A+2];     //hole density of node A
  PetscScalar Eqca = x[zofs[zone_index]+5*A+3];     //quantum conduction band of node A
  PetscScalar Eqva = x[zofs[zone_index]+5*A+4];     //quantum valence band of node A
  mt->mapping(&pzone->danode[A],&aux[A],0);
  PetscScalar Eca = e*(Eqca/e + kb*fs[A].T/e*(log(aux[A].Nc)-1.5*log(fs[A].T)));
  PetscScalar Eva = e*(Eqva/e - kb*fs[A].T/e*(log(aux[A].Nv)-1.5*log(fs[A].T)));
  PetscScalar Ra = mt->band->Recomb(pa,na,fs[A].T);
  PetscScalar nia = mt->band->nie(T);
  PetscScalar phina = Va - log(fabs(na)/nia)*Vt;
  PetscScalar phipa = Va + log(fabs(pa)/nia)*Vt;
  PetscScalar qca = -(-e*phina-Eqca)/(kb*T);
  PetscScalar qva =  (-e*phipa-Eqva)/(kb*T);
  
  PetscScalar Vb   = x[zofs[zone_index]+5*B+0];     //potential of node B
  PetscScalar nb   = x[zofs[zone_index]+5*B+1];     //electron density of node B
  PetscScalar pb   = x[zofs[zone_index]+5*B+2];     //hole density of node B
  PetscScalar Eqcb = x[zofs[zone_index]+5*B+3];     //quantum conduction band of node B
  PetscScalar Eqvb = x[zofs[zone_index]+5*B+4];     //quantum valence band of node B
  mt->mapping(&pzone->danode[B],&aux[B],0);
  PetscScalar Ecb = e*(Eqcb/e + kb*fs[B].T/e*(log(aux[B].Nc)-1.5*log(fs[B].T)));
  PetscScalar Evb = e*(Eqvb/e - kb*fs[B].T/e*(log(aux[B].Nv)-1.5*log(fs[B].T)));
  PetscScalar Rb = mt->band->Recomb(pb,nb,fs[B].T);
  PetscScalar nib = mt->band->nie(T);
  PetscScalar phinb = Vb - log(fabs(nb)/nib)*Vt;
  PetscScalar phipb = Vb + log(fabs(pb)/nib)*Vt;
  PetscScalar qcb = -(-e*phinb-Eqcb)/(kb*T);
  PetscScalar qvb =  (-e*phipb-Eqvb)/(kb*T);
  
  PetscScalar Vc   = x[zofs[zone_index]+5*C+0];     //potential of node C
  PetscScalar nc   = x[zofs[zone_index]+5*C+1];     //electron density of node C
  PetscScalar pc   = x[zofs[zone_index]+5*C+2];     //hole density of node C
  PetscScalar Eqcc = x[zofs[zone_index]+5*C+3];     //quantum conduction band of node C
  PetscScalar Eqvc = x[zofs[zone_index]+5*C+4];     //quantum valence band of node C
  mt->mapping(&pzone->danode[C],&aux[C],0);
  PetscScalar Ecc = e*(Eqcc/e + kb*fs[C].T/e*(log(aux[C].Nc)-1.5*log(fs[C].T)));
  PetscScalar Evc = e*(Eqvc/e - kb*fs[C].T/e*(log(aux[C].Nv)-1.5*log(fs[C].T)));
  PetscScalar Rc = mt->band->Recomb(pc,nc,fs[C].T);
  PetscScalar nic = mt->band->nie(T);
  PetscScalar phinc = Vc - log(fabs(nc)/nic)*Vt;
  PetscScalar phipc = Vc + log(fabs(pc)/nic)*Vt;
  PetscScalar qcc = -(-e*phinc-Eqcc)/(kb*T);
  PetscScalar qvc =  (-e*phipc-Eqvc)/(kb*T);
  
  PetscScalar Ex = -((yb-yc)*Va + (yc-ya)*Vb +(ya-yb)*Vc)/(2*tri_area);
  PetscScalar Ey = -((xc-xb)*Va + (xa-xc)*Vb +(xb-xa)*Vc)/(2*tri_area);
  PetscScalar E = sqrt(Ex*Ex+Ey*Ey+1e-20);
         
  PetscScalar G = 0;
  PetscScalar IIn=0,IIp=0;
  PetscScalar mun,mup;
  PetscScalar Epnc=0,Etnc=0;
  PetscScalar Eppc=0,Etpc=0;
  PetscScalar Epna=0,Etna=0;
  PetscScalar Eppa=0,Etpa=0;
  PetscScalar Epnb=0,Etnb=0;
  PetscScalar Eppb=0,Etpb=0;
  PetscScalar Jna_norm,Jpa_norm;
  PetscScalar Jnb_norm,Jpb_norm;
  PetscScalar Jnc_norm,Jpc_norm;

#ifdef _FLUX1_
  PetscScalar Jnc =  In(Vt,Eca/e,Ecb/e,na,nb,ptri->edge_len[2]);
  PetscScalar Jpc =  Ip(Vt,Eva/e,Evb/e,pa,pb,ptri->edge_len[2]);
  PetscScalar Jna =  In(Vt,Ecb/e,Ecc/e,nb,nc,ptri->edge_len[0]);
  PetscScalar Jpa =  Ip(Vt,Evb/e,Evc/e,pb,pc,ptri->edge_len[0]);
  PetscScalar Jnb =  In(Vt,Ecc/e,Eca/e,nc,na,ptri->edge_len[1]);
  PetscScalar Jpb =  Ip(Vt,Evc/e,Eva/e,pc,pa,ptri->edge_len[1]);
  PetscScalar Jn_scale = dmax(dmax(fabs(Jna),fabs(Jnb)),dmax(fabs(Jnc),nia*nia));
  PetscScalar Jp_scale = dmax(dmax(fabs(Jpa),fabs(Jpb)),dmax(fabs(Jpc),nia*nia));
  //cout<<Jn_scale<<" "<<Jnc<<" "<<Jna<<" "<<Jnb<<endl;
  //cout<<Jp_scale<<" "<<Jpc<<" "<<Jpa<<" "<<Jpb<<endl;
  Jnc /=  Jn_scale;
  Jpc /=  Jp_scale;
  Jna /=  Jn_scale;
  Jpa /=  Jp_scale;
  Jnb /=  Jn_scale;
  Jpb /=  Jp_scale;
#endif
#ifdef _FLUX2_
  PetscScalar Jnc =  In(Vt,(Ecb-Eca)/e,na,nb,ptri->edge_len[2]);
  PetscScalar Jpc =  Ip(Vt,(Evb-Eva)/e,pa,pb,ptri->edge_len[2]);
  PetscScalar Jna =  In(Vt,(Ecc-Ecb)/e,nb,nc,ptri->edge_len[0]);
  PetscScalar Jpa =  Ip(Vt,(Evc-Evb)/e,pb,pc,ptri->edge_len[0]);
  PetscScalar Jnb =  In(Vt,(Eca-Ecc)/e,nc,na,ptri->edge_len[1]);
  PetscScalar Jpb =  Ip(Vt,(Eva-Evc)/e,pc,pa,ptri->edge_len[1]);
  PetscScalar Jn_scale = dmax(dmax(fabs(Jna),fabs(Jnb)),dmax(fabs(Jnc),nia*nia));
  PetscScalar Jp_scale = dmax(dmax(fabs(Jpa),fabs(Jpb)),dmax(fabs(Jpc),nia*nia));
  Jnc /=  Jn_scale;
  Jpc /=  Jp_scale;
  Jna /=  Jn_scale;
  Jpa /=  Jp_scale;
  Jnb /=  Jn_scale;
  Jpb /=  Jp_scale;
#endif
 
  //flux along A-B
  if(HighFieldMobility)
  {
    if(EJModel || IIType==EdotJ)
    {
      PetscScalar Jncx = ((Wca+Wcb)*Scx - Wca*Mb*Sax - Wcb*Ma*Sbx)*Jnc
                       + (Wca*Sax - Wca*Mb*Scx)*Jna
                       + (Wcb*Sbx - Wcb*Ma*Scx)*Jnb;
      PetscScalar Jncy = ((Wca+Wcb)*Scy - Wca*Mb*Say - Wcb*Ma*Sby)*Jnc
                       + (Wca*Say - Wca*Mb*Scy)*Jna
                       + (Wcb*Sby - Wcb*Ma*Scy)*Jnb;
      PetscScalar Jpcx = ((Wca+Wcb)*Scx - Wca*Mb*Sax - Wcb*Ma*Sbx)*Jpc
                       + (Wca*Sax - Wca*Mb*Scx)*Jpa
                       + (Wcb*Sbx - Wcb*Ma*Scx)*Jpb;
      PetscScalar Jpcy = ((Wca+Wcb)*Scy - Wca*Mb*Say - Wcb*Ma*Sby)*Jpc
                       + (Wca*Say - Wca*Mb*Scy)*Jpa
                       + (Wcb*Sby - Wcb*Ma*Scy)*Jpb;
      Jnc_norm = sqrt(Jncx*Jncx + Jncy*Jncy + 1e-100);
      Jpc_norm = sqrt(Jpcx*Jpcx + Jpcy*Jpcy + 1e-100);
      Epnc = (Ex*Jncx + Ey*Jncy)/Jnc_norm;
      Etnc = (Ex*Jncy - Ey*Jncx)/Jnc_norm;
      Eppc = (Ex*Jpcx + Ey*Jpcy)/Jpc_norm;
      Etpc = (Ex*Jpcy - Ey*Jpcx)/Jpc_norm;
    }
    else
    {
      Epnc = fabs(Eca-Ecb)/e/ptri->edge_len[2]; //parallel electrical field for electron
      Eppc = fabs(Eva-Evb)/e/ptri->edge_len[2]; //parallel electrical field for hole
      //transvers electrical field for electron and hole
      Etnc = fabs(Eca + ptri->edge_len[1]*cos(ptri->angle[0])*(Ecb-Eca)/ptri->edge_len[2] - Ecc)
           /(ptri->edge_len[1]*sin(ptri->angle[0]))/e;
      Etpc = fabs(Eva + ptri->edge_len[1]*cos(ptri->angle[0])*(Evb-Eva)/ptri->edge_len[2] - Evc)
           /(ptri->edge_len[1]*sin(ptri->angle[0]))/e;
    }

    mt->mapping(&pzone->danode[A],&aux[A],0);
    aux[A].mun =  mt->mob->ElecMob(pa,na,T,dmax(0,Epnc),fabs(Etnc),T);
    aux[A].mup =  mt->mob->HoleMob(pa,na,T,dmax(0,Eppc),fabs(Etpc),T);

    mt->mapping(&pzone->danode[B],&aux[B],0);
    aux[B].mun =  mt->mob->ElecMob(pb,nb,T,dmax(0,Epnc),fabs(Etnc),T);
    aux[B].mup =  mt->mob->HoleMob(pb,nb,T,dmax(0,Eppc),fabs(Etpc),T);
  }
  
  mun = 0.5*(aux[A].mun+aux[B].mun);
  mup = 0.5*(aux[A].mup+aux[B].mup);
  
  
  f[zofs[zone_index]+5*A+0] +=   0.5*(aux[A].eps+aux[B].eps)*(Vb-Va)/ptri->edge_len[2]*ptri->d[2];
  f[zofs[zone_index]+5*A+1] +=   mun*Jnc*Jn_scale*ptri->d[2];
  f[zofs[zone_index]+5*A+2] += - mup*Jpc*Jp_scale*ptri->d[2];
  f[zofs[zone_index]+5*A+3] += - QNFactor*ptri->d[2]*bn*qV(qca,qcb,ptri->edge_len[2]);
  f[zofs[zone_index]+5*A+4] +=   QPFactor*ptri->d[2]*bp*qV(qva,qvb,ptri->edge_len[2]);

  f[zofs[zone_index]+5*B+0] += - 0.5*(aux[A].eps+aux[B].eps)*(Vb-Va)/ptri->edge_len[2]*ptri->d[2];
  f[zofs[zone_index]+5*B+1] += - mun*Jnc*Jn_scale*ptri->d[2];
  f[zofs[zone_index]+5*B+2] +=   mup*Jpc*Jp_scale*ptri->d[2];
  f[zofs[zone_index]+5*B+3] += - QNFactor*ptri->d[2]*bn*qV(qcb,qca,ptri->edge_len[2]);
  f[zofs[zone_index]+5*B+4] +=   QPFactor*ptri->d[2]*bp*qV(qvb,qva,ptri->edge_len[2]);
  
  //flux along B-C
  if(HighFieldMobility)
  {
    if(EJModel || IIType==EdotJ)
    {
      PetscScalar Jnax = ((Wab+Wac)*Sax - Wab*Mc*Sbx - Wac*Mb*Scx)*Jna
                       + (Wab*Sbx - Wab*Mc*Sax)*Jnb
                       + (Wac*Scx - Wac*Mb*Sax)*Jnc;
      PetscScalar Jnay = ((Wab+Wac)*Say - Wab*Mc*Sby - Wac*Mb*Scy)*Jna
                       + (Wab*Sby - Wab*Mc*Say)*Jnb
                       + (Wac*Scy - Wac*Mb*Say)*Jnc;
      PetscScalar Jpax = ((Wab+Wac)*Sax - Wab*Mc*Sbx - Wac*Mb*Scx)*Jpa
                       + (Wab*Sbx - Wab*Mc*Sax)*Jpb
                       + (Wac*Scx - Wac*Mb*Sax)*Jpc;
      PetscScalar Jpay = ((Wab+Wac)*Say - Wab*Mc*Sby - Wac*Mb*Scy)*Jpa
                       + (Wab*Sby - Wab*Mc*Say)*Jpb
                       + (Wac*Scy - Wac*Mb*Say)*Jpc;
      Jna_norm = sqrt(Jnax*Jnax + Jnay*Jnay + 1e-100);
      Jpa_norm = sqrt(Jpax*Jpax + Jpay*Jpay + 1e-100);
      Epna = (Ex*Jnax + Ey*Jnay)/Jna_norm;
      Etna = (Ex*Jnay - Ey*Jnax)/Jna_norm;
      Eppa = (Ex*Jpax + Ey*Jpay)/Jpa_norm;
      Etpa = (Ex*Jpay - Ey*Jpax)/Jpa_norm;
    }
    else
    {
      Epna = fabs(Ecb-Ecc)/e/ptri->edge_len[0]; //parallel electrical field for electron
      Eppa = fabs(Evb-Evc)/e/ptri->edge_len[0]; //parallel electrical field for hole
      //transvers electrical field for electron and hole
      Etna = fabs(Ecb + ptri->edge_len[2]*cos(ptri->angle[1])*(Ecc-Ecb)/ptri->edge_len[0] - Eca)
           /(ptri->edge_len[2]*sin(ptri->angle[1]))/e;
      Etpa = fabs(Evb + ptri->edge_len[2]*cos(ptri->angle[1])*(Evc-Evb)/ptri->edge_len[0] - Eva)
           /(ptri->edge_len[2]*sin(ptri->angle[1]))/e;
    }
    mt->mapping(&pzone->danode[B],&aux[B],0);
    aux[B].mun =  mt->mob->ElecMob(pb,nb,T,dmax(0,Epna),fabs(Etna),T);
    aux[B].mup =  mt->mob->HoleMob(pb,nb,T,dmax(0,Eppa),fabs(Etpa),T);

    mt->mapping(&pzone->danode[C],&aux[C],0);
    aux[C].mun =  mt->mob->ElecMob(pc,nc,T,dmax(0,Epna),fabs(Etna),T);
    aux[C].mup =  mt->mob->HoleMob(pc,nc,T,dmax(0,Eppa),fabs(Etpa),T);
  }
  mun = 0.5*(aux[B].mun+aux[C].mun);
  mup = 0.5*(aux[B].mup+aux[C].mup);


  f[zofs[zone_index]+5*B+0] +=   0.5*(aux[B].eps+aux[C].eps)*(Vc-Vb)/ptri->edge_len[0]*ptri->d[0];
  f[zofs[zone_index]+5*B+1] +=   mun*Jna*Jn_scale*ptri->d[0];
  f[zofs[zone_index]+5*B+2] += - mup*Jpa*Jp_scale*ptri->d[0];
  f[zofs[zone_index]+5*B+3] += - QNFactor*ptri->d[0]*bn*qV(qcb,qcc,ptri->edge_len[0]);
  f[zofs[zone_index]+5*B+4] +=   QPFactor*ptri->d[0]*bp*qV(qvb,qvc,ptri->edge_len[0]);
  
  f[zofs[zone_index]+5*C+0] += - 0.5*(aux[B].eps+aux[C].eps)*(Vc-Vb)/ptri->edge_len[0]*ptri->d[0];
  f[zofs[zone_index]+5*C+1] += - mun*Jna*Jn_scale*ptri->d[0];
  f[zofs[zone_index]+5*C+2] +=   mup*Jpa*Jp_scale*ptri->d[0];
  f[zofs[zone_index]+5*C+3] += - QNFactor*ptri->d[0]*bn*qV(qcc,qcb,ptri->edge_len[0]);
  f[zofs[zone_index]+5*C+4] +=   QPFactor*ptri->d[0]*bp*qV(qvc,qvb,ptri->edge_len[0]);
   
  //flux along C-A
  if(HighFieldMobility)
  {
    if(EJModel || IIType==EdotJ)
    {
      PetscScalar Jnbx = ((Wbc+Wba)*Sbx - Wbc*Ma*Scx - Wba*Mc*Sax)*Jnb
                       + (Wbc*Scx - Wbc*Ma*Sbx)*Jnc
                       + (Wba*Sax - Wba*Mc*Sbx)*Jna;
      PetscScalar Jnby = ((Wbc+Wba)*Sby - Wbc*Ma*Scy - Wba*Mc*Say)*Jnb
                       + (Wbc*Scy - Wbc*Ma*Sby)*Jnc
                       + (Wba*Say - Wba*Mc*Sby)*Jna;
      PetscScalar Jpbx = ((Wbc+Wba)*Sbx - Wbc*Ma*Scx - Wba*Mc*Sax)*Jpb
                       + (Wbc*Scx - Wbc*Ma*Sbx)*Jpc
                       + (Wba*Sax - Wba*Mc*Sbx)*Jpa;
      PetscScalar Jpby = ((Wbc+Wba)*Sby - Wbc*Ma*Scy - Wba*Mc*Say)*Jpb
                       + (Wbc*Scy - Wbc*Ma*Sby)*Jpc
                       + (Wba*Say - Wba*Mc*Sby)*Jpa;
      Jnb_norm = sqrt(Jnbx*Jnbx + Jnby*Jnby + 1e-100);
      Jpb_norm = sqrt(Jpbx*Jpbx + Jpby*Jpby + 1e-100);
      Epnb = (Ex*Jnbx + Ey*Jnby)/Jnb_norm;
      Etnb = (Ex*Jnby - Ey*Jnbx)/Jnb_norm;
      Eppb = (Ex*Jpbx + Ey*Jpby)/Jpb_norm;
      Etpb = (Ex*Jpby - Ey*Jpbx)/Jpb_norm;
  
    }
    else
    {
      Epnb = fabs(Ecc-Eca)/e/ptri->edge_len[1]; //parallel electrical field for electron
      Eppb = fabs(Evc-Eva)/e/ptri->edge_len[1]; //parallel electrical field for hole
      //transvers electrical field for electron and hole
      Etnb = fabs(Ecc + ptri->edge_len[0]*cos(ptri->angle[2])*(Eca-Ecc)/ptri->edge_len[1] - Ecb)
           /(ptri->edge_len[0]*sin(ptri->angle[2]))/e;
      Etpb = fabs(Evc + ptri->edge_len[0]*cos(ptri->angle[2])*(Eva-Evc)/ptri->edge_len[1] - Evb)
           /(ptri->edge_len[0]*sin(ptri->angle[2]))/e;
    }
    mt->mapping(&pzone->danode[C],&aux[C],0);
    aux[C].mun =  mt->mob->ElecMob(pc,nc,T,dmax(0,Epnb),fabs(Etnb),T);
    aux[C].mup =  mt->mob->HoleMob(pc,nc,T,dmax(0,Eppb),fabs(Etpb),T);

    mt->mapping(&pzone->danode[A],&aux[A],0);
    aux[A].mun =  mt->mob->ElecMob(pa,na,T,dmax(0,Epnb),fabs(Etnb),T);
    aux[A].mup =  mt->mob->HoleMob(pa,na,T,dmax(0,Eppb),fabs(Etpb),T);
  }
  mun = 0.5*(aux[C].mun+aux[A].mun);
  mup = 0.5*(aux[C].mup+aux[A].mup);


  f[zofs[zone_index]+5*C+0] +=   0.5*(aux[C].eps+aux[A].eps)*(Va-Vc)/ptri->edge_len[1]*ptri->d[1];
  f[zofs[zone_index]+5*C+1] +=   mun*Jnb*Jn_scale*ptri->d[1];
  f[zofs[zone_index]+5*C+2] += - mup*Jpb*Jp_scale*ptri->d[1];
  f[zofs[zone_index]+5*C+3] += - QNFactor*ptri->d[1]*bn*qV(qcc,qca,ptri->edge_len[1]);
  f[zofs[zone_index]+5*C+4] +=   QPFactor*ptri->d[1]*bp*qV(qvc,qva,ptri->edge_len[1]);
  
  f[zofs[zone_index]+5*A+0] += - 0.5*(aux[C].eps+aux[A].eps)*(Va-Vc)/ptri->edge_len[1]*ptri->d[1];
  f[zofs[zone_index]+5*A+1] += - mun*Jnb*Jn_scale*ptri->d[1];
  f[zofs[zone_index]+5*A+2] +=   mup*Jpb*Jp_scale*ptri->d[1];
  f[zofs[zone_index]+5*A+3] += - QNFactor*ptri->d[1]*bn*qV(qca,qcc,ptri->edge_len[1]);
  f[zofs[zone_index]+5*A+4] +=   QPFactor*ptri->d[1]*bp*qV(qva,qvc,ptri->edge_len[1]);
  
  //---------------------------------------------------------------------------
  // Recombination item
  //---------------------------------------------------------------------------
  
  S = 0.25*ptri->d[2]*ptri->edge_len[2] + 0.25*ptri->d[1]*ptri->edge_len[1];
  f[zofs[zone_index]+5*A+1] += -Ra*S;
  f[zofs[zone_index]+5*A+2] += -Ra*S;
  
  S = 0.25*ptri->d[0]*ptri->edge_len[0] + 0.25*ptri->d[2]*ptri->edge_len[2];
  f[zofs[zone_index]+5*B+1] += -Rb*S;
  f[zofs[zone_index]+5*B+2] += -Rb*S;
  
  S = 0.25*ptri->d[1]*ptri->edge_len[1] + 0.25*ptri->d[0]*ptri->edge_len[0];
  f[zofs[zone_index]+5*C+1] += -Rc*S;
  f[zofs[zone_index]+5*C+2] += -Rc*S;
}



void SMCZone::J1Q_Tri_qddm(Tri *ptri,PetscScalar *x,Mat *jtmp, vector<int> & zofs)
{
  int z = zone_index;
  int A = ptri->node[0];
  int B = ptri->node[1];
  int C = ptri->node[2];
  
  PetscScalar xa = pzone->danode[ptri->node[0]].x;
  PetscScalar xb = pzone->danode[ptri->node[1]].x;
  PetscScalar xc = pzone->danode[ptri->node[2]].x;
  PetscScalar ya = pzone->danode[ptri->node[0]].y;
  PetscScalar yb = pzone->danode[ptri->node[1]].y;
  PetscScalar yc = pzone->danode[ptri->node[2]].y;
  PetscScalar tri_area = ptri->area;
  PetscScalar Sax = (xc-xb)/ptri->edge_len[0];
  PetscScalar Say = (yc-yb)/ptri->edge_len[0];
  PetscScalar Sbx = (xa-xc)/ptri->edge_len[1];
  PetscScalar Sby = (ya-yc)/ptri->edge_len[1];
  PetscScalar Scx = (xb-xa)/ptri->edge_len[2];
  PetscScalar Scy = (yb-ya)/ptri->edge_len[2];
  PetscScalar Wab = ptri->d[1]/((ptri->d[1]+ptri->d[2])*sin(ptri->angle[2])*sin(ptri->angle[2]));
  PetscScalar Wac = ptri->d[2]/((ptri->d[1]+ptri->d[2])*sin(ptri->angle[1])*sin(ptri->angle[1]));
  PetscScalar Wbc = ptri->d[2]/((ptri->d[0]+ptri->d[2])*sin(ptri->angle[0])*sin(ptri->angle[0]));
  PetscScalar Wba = ptri->d[0]/((ptri->d[0]+ptri->d[2])*sin(ptri->angle[2])*sin(ptri->angle[2]));
  PetscScalar Wca = ptri->d[0]/((ptri->d[0]+ptri->d[1])*sin(ptri->angle[1])*sin(ptri->angle[1]));
  PetscScalar Wcb = ptri->d[1]/((ptri->d[0]+ptri->d[1])*sin(ptri->angle[0])*sin(ptri->angle[0]));
  PetscScalar Ma  = -cos(ptri->angle[0]);
  PetscScalar Mb  = -cos(ptri->angle[1]);
  PetscScalar Mc  = -cos(ptri->angle[2]);
  
  Mat5       J1,J2,J3;
  Vec5       scale;
  PetscInt   index1[5],index2[5],index3[5];
  index1[0] = zofs[z]+5*A+0;
  index1[1] = zofs[z]+5*A+1;
  index1[2] = zofs[z]+5*A+2;
  index1[3] = zofs[z]+5*A+3;
  index1[4] = zofs[z]+5*A+4;
  
  index2[0] = zofs[z]+5*B+0;
  index2[1] = zofs[z]+5*B+1;
  index2[2] = zofs[z]+5*B+2;
  index2[3] = zofs[z]+5*B+3;
  index2[4] = zofs[z]+5*B+4;
  
  index3[0] = zofs[z]+5*C+0;
  index3[1] = zofs[z]+5*C+1;
  index3[2] = zofs[z]+5*C+2;
  index3[3] = zofs[z]+5*C+3;
  index3[4] = zofs[z]+5*C+4;
  
  //---------------------------------------------------------------------------

  PetscScalar kb = mt->kb;
  PetscScalar hbar = mt->hbar;
  PetscScalar e  =  mt->e;
  PetscScalar T  = (fs[A].T+fs[B].T+fs[C].T)/3.0;
  PetscScalar Eg = mt->band->Eg(T);
  PetscScalar Vt = kb*T/e;
  PetscScalar me =  mt->band->EffecElecMass(T);
  PetscScalar mh =  mt->band->EffecHoleMass(T);
  PetscScalar gamman = mt->band->Gamman(); 
  PetscScalar gammap = mt->band->Gammap();
  PetscScalar bn = gamman*hbar*hbar/(6*e*me);
  PetscScalar bp = gammap*hbar*hbar/(6*e*mh);
  PetscScalar S;
   
  //the indepedent variable number
  adtl::AutoDScalar::numdir=15;
  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir); 
  AutoDScalar TD = T; // dummy ad variable.
  AutoDScalar F[5];
  
  AutoDScalar Va   = x[zofs[zone_index]+5*A+0];     //potential of node A
  AutoDScalar na   = x[zofs[zone_index]+5*A+1];     //electron density of node A
  AutoDScalar pa   = x[zofs[zone_index]+5*A+2];     //hole density of node A
  AutoDScalar Eqca = x[zofs[zone_index]+5*A+3];     //quantum conduction band of node A
  AutoDScalar Eqva = x[zofs[zone_index]+5*A+4];     //quantum valence band of node A
  Va.setADValue(0,1.0);
  na.setADValue(1,1.0);
  pa.setADValue(2,1.0);
  Eqca.setADValue(3,1.0);
  Eqva.setADValue(4,1.0);
  mt->mapping(&pzone->danode[A],&aux[A],0);
  AutoDScalar Eca = e*(Eqca/e + kb*fs[A].T/e*(log(aux[A].Nc)-1.5*log(fs[A].T)));
  AutoDScalar Eva = e*(Eqva/e - kb*fs[A].T/e*(log(aux[A].Nv)-1.5*log(fs[A].T)));
  AutoDScalar Ra = mt->band->Recomb(pa,na,fs[A].T);
  PetscScalar nia = mt->band->nie(T);
  AutoDScalar phina = Va - log(fabs(na)/nia)*Vt;
  AutoDScalar phipa = Va + log(fabs(pa)/nia)*Vt;
  AutoDScalar qca = -(-e*phina-Eqca)/(kb*T);
  AutoDScalar qva =  (-e*phipa-Eqva)/(kb*T);
  
  AutoDScalar Vb   = x[zofs[zone_index]+5*B+0];     //potential of node B
  AutoDScalar nb   = x[zofs[zone_index]+5*B+1];     //electron density of node B
  AutoDScalar pb   = x[zofs[zone_index]+5*B+2];     //hole density of node B
  AutoDScalar Eqcb = x[zofs[zone_index]+5*B+3];     //quantum conduction band of node B
  AutoDScalar Eqvb = x[zofs[zone_index]+5*B+4];     //quantum valence band of node B
  Vb.setADValue(5,1.0);
  nb.setADValue(6,1.0);
  pb.setADValue(7,1.0);
  Eqcb.setADValue(8,1.0);
  Eqvb.setADValue(9,1.0);
  mt->mapping(&pzone->danode[B],&aux[B],0);
  AutoDScalar Ecb = e*(Eqcb/e + kb*fs[B].T/e*(log(aux[B].Nc)-1.5*log(fs[B].T)));
  AutoDScalar Evb = e*(Eqvb/e - kb*fs[B].T/e*(log(aux[B].Nv)-1.5*log(fs[B].T)));
  AutoDScalar Rb = mt->band->Recomb(pb,nb,fs[B].T);
  PetscScalar nib = mt->band->nie(T);
  AutoDScalar phinb = Vb - log(fabs(nb)/nib)*Vt;
  AutoDScalar phipb = Vb + log(fabs(pb)/nib)*Vt;
  AutoDScalar qcb = -(-e*phinb-Eqcb)/(kb*T);
  AutoDScalar qvb =  (-e*phipb-Eqvb)/(kb*T);
  
  AutoDScalar Vc   = x[zofs[zone_index]+5*C+0];     //potential of node C
  AutoDScalar nc   = x[zofs[zone_index]+5*C+1];     //electron density of node C
  AutoDScalar pc   = x[zofs[zone_index]+5*C+2];     //hole density of node C
  AutoDScalar Eqcc = x[zofs[zone_index]+5*C+3];     //quantum conduction band of node C
  AutoDScalar Eqvc = x[zofs[zone_index]+5*C+4];     //quantum valence band of node C
  Vc.setADValue(10,1.0);
  nc.setADValue(11,1.0);
  pc.setADValue(12,1.0);
  Eqcc.setADValue(13,1.0);
  Eqvc.setADValue(14,1.0);
  mt->mapping(&pzone->danode[C],&aux[C],0);
  AutoDScalar Ecc = e*(Eqcc/e + kb*fs[C].T/e*(log(aux[C].Nc)-1.5*log(fs[C].T)));
  AutoDScalar Evc = e*(Eqvc/e - kb*fs[C].T/e*(log(aux[C].Nv)-1.5*log(fs[C].T)));
  AutoDScalar Rc = mt->band->Recomb(pc,nc,fs[C].T);
  PetscScalar nic = mt->band->nie(T);
  AutoDScalar phinc = Vc - log(fabs(nc)/nic)*Vt;
  AutoDScalar phipc = Vc + log(fabs(pc)/nic)*Vt;
  AutoDScalar qcc = -(-e*phinc-Eqcc)/(kb*T);
  AutoDScalar qvc =  (-e*phipc-Eqvc)/(kb*T);
  
  AutoDScalar Ex = -((yb-yc)*Va + (yc-ya)*Vb +(ya-yb)*Vc)/(2*tri_area);
  AutoDScalar Ey = -((xc-xb)*Va + (xa-xc)*Vb +(xb-xa)*Vc)/(2*tri_area);
  AutoDScalar E = sqrt(Ex*Ex+Ey*Ey+1e-20);
         
  AutoDScalar G = 0;
  AutoDScalar IIn=0,IIp=0;
  AutoDScalar mun,mup;
  AutoDScalar Epnc=0,Etnc=0;
  AutoDScalar Eppc=0,Etpc=0;
  AutoDScalar Epna=0,Etna=0;
  AutoDScalar Eppa=0,Etpa=0;
  AutoDScalar Epnb=0,Etnb=0;
  AutoDScalar Eppb=0,Etpb=0;
  AutoDScalar Jna_norm=0,Jpa_norm=0;
  AutoDScalar Jnb_norm=0,Jpb_norm=0;
  AutoDScalar Jnc_norm=0,Jpc_norm=0;
  
  #ifdef _FLUX1_
  AutoDScalar Jnc =  In(Vt,Eca/e,Ecb/e,na,nb,ptri->edge_len[2]);
  AutoDScalar Jpc =  Ip(Vt,Eva/e,Evb/e,pa,pb,ptri->edge_len[2]);
  AutoDScalar Jna =  In(Vt,Ecb/e,Ecc/e,nb,nc,ptri->edge_len[0]);
  AutoDScalar Jpa =  Ip(Vt,Evb/e,Evc/e,pb,pc,ptri->edge_len[0]);
  AutoDScalar Jnb =  In(Vt,Ecc/e,Eca/e,nc,na,ptri->edge_len[1]);
  AutoDScalar Jpb =  Ip(Vt,Evc/e,Eva/e,pc,pa,ptri->edge_len[1]);
  PetscScalar Jn_scale = adtl::fmax(adtl::fmax(fabs(Jna.getValue()),fabs(Jnb.getValue())),
                                    adtl::fmax(fabs(Jnc.getValue()),nia*nia));
  PetscScalar Jp_scale = adtl::fmax(adtl::fmax(fabs(Jpa.getValue()),fabs(Jpb.getValue())),
                                    adtl::fmax(fabs(Jpc.getValue()),nia*nia));
  
  Jnc =  Jnc/Jn_scale;
  Jpc =  Jpc/Jp_scale;
  Jna =  Jna/Jn_scale;
  Jpa =  Jpa/Jp_scale;
  Jnb =  Jnb/Jn_scale;
  Jpb =  Jpb/Jp_scale;

#endif
#ifdef _FLUX2_
  AutoDScalar Jnc =  In(Vt,(Ecb-Eca)/e,na,nb,ptri->edge_len[2]);
  AutoDScalar Jpc =  Ip(Vt,(Evb-Eva)/e,pa,pb,ptri->edge_len[2]);
  AutoDScalar Jna =  In(Vt,(Ecc-Ecb)/e,nb,nc,ptri->edge_len[0]);
  AutoDScalar Jpa =  Ip(Vt,(Evc-Evb)/e,pb,pc,ptri->edge_len[0]);
  AutoDScalar Jnb =  In(Vt,(Eca-Ecc)/e,nc,na,ptri->edge_len[1]);
  AutoDScalar Jpb =  Ip(Vt,(Eva-Evc)/e,pc,pa,ptri->edge_len[1]);
  PetscScalar Jn_scale = adtl::fmax(adtl::fmax(fabs(Jna.getValue()),fabs(Jnb.getValue())),
                                    adtl::fmax(fabs(Jnc.getValue()),nia*nia));
  PetscScalar Jp_scale = adtl::fmax(adtl::fmax(fabs(Jpa.getValue()),fabs(Jpb.getValue())),
                                    adtl::fmax(fabs(Jpc.getValue()),nia*nia));
  Jnc /=  Jn_scale;
  Jpc /=  Jp_scale;
  Jna /=  Jn_scale;
  Jpa /=  Jp_scale;
  Jnb /=  Jn_scale;
  Jpb /=  Jp_scale;
#endif

  //---------------------------------------------------------------------------
  //flux along A-B
  //---------------------------------------------------------------------------
  if(HighFieldMobility)
  {
    if(EJModel || IIType==EdotJ)
    {
      AutoDScalar Jncx = ((Wca+Wcb)*Scx - Wca*Mb*Sax - Wcb*Ma*Sbx)*Jnc
                       + (Wca*Sax - Wca*Mb*Scx)*Jna
                       + (Wcb*Sbx - Wcb*Ma*Scx)*Jnb;
      AutoDScalar Jncy = ((Wca+Wcb)*Scy - Wca*Mb*Say - Wcb*Ma*Sby)*Jnc
                       + (Wca*Say - Wca*Mb*Scy)*Jna
                       + (Wcb*Sby - Wcb*Ma*Scy)*Jnb;
      AutoDScalar Jpcx = ((Wca+Wcb)*Scx - Wca*Mb*Sax - Wcb*Ma*Sbx)*Jpc
                       + (Wca*Sax - Wca*Mb*Scx)*Jpa
                       + (Wcb*Sbx - Wcb*Ma*Scx)*Jpb;
      AutoDScalar Jpcy = ((Wca+Wcb)*Scy - Wca*Mb*Say - Wcb*Ma*Sby)*Jpc
                       + (Wca*Say - Wca*Mb*Scy)*Jpa
                       + (Wcb*Sby - Wcb*Ma*Scy)*Jpb;
      Jnc_norm = sqrt(Jncx*Jncx + Jncy*Jncy + 1e-100);
      Jpc_norm = sqrt(Jpcx*Jpcx + Jpcy*Jpcy + 1e-100);
      Epnc = Ex*(Jncx/Jnc_norm) + Ey*(Jncy/Jnc_norm);
      Etnc = Ex*(Jncy/Jnc_norm) - Ey*(Jncx/Jnc_norm);
      Eppc = Ex*(Jpcx/Jpc_norm) + Ey*(Jpcy/Jpc_norm);
      Etpc = Ex*(Jpcy/Jpc_norm) - Ey*(Jpcx/Jpc_norm);
    }
    else
    {
      Epnc = fabs(Eca-Ecb)/e/ptri->edge_len[2]; //parallel electrical field for electron
      Eppc = fabs(Eva-Evb)/e/ptri->edge_len[2]; //parallel electrical field for hole
      //transvers electrical field for electron and hole
      Etnc = fabs(Eca + ptri->edge_len[1]*cos(ptri->angle[0])*(Ecb-Eca)/ptri->edge_len[2] - Ecc)
           /(ptri->edge_len[1]*sin(ptri->angle[0]))/e;
      Etpc = fabs(Eva + ptri->edge_len[1]*cos(ptri->angle[0])*(Evb-Eva)/ptri->edge_len[2] - Evc)
           /(ptri->edge_len[1]*sin(ptri->angle[0]))/e;
    }
    mt->mapping(&pzone->danode[A],&aux[A],0);
    AutoDScalar munCA =  mt->mob->ElecMob(pa,na,TD,adtl::fmax(0,Epnc),fabs(Etnc),TD);
    AutoDScalar mupCA =  mt->mob->HoleMob(pa,na,TD,adtl::fmax(0,Eppc),fabs(Etpc),TD);

    mt->mapping(&pzone->danode[B],&aux[B],0);
    AutoDScalar munCB =  mt->mob->ElecMob(pb,nb,TD,adtl::fmax(0,Epnc),fabs(Etnc),TD);
    AutoDScalar mupCB =  mt->mob->HoleMob(pb,nb,TD,adtl::fmax(0,Eppc),fabs(Etpc),TD);
  
    mun = 0.5*(munCA+munCB);
    mup = 0.5*(mupCA+mupCB);
  }
  else
  {
    mun = 0.5*(aux[A].mun+aux[B].mun);
    mup = 0.5*(aux[A].mup+aux[B].mup);
  }
  
  //---------------------------------------------------------------------------
  F[0] =   0.5*(aux[A].eps+aux[B].eps)*(Vb-Va)/ptri->edge_len[2]*ptri->d[2];
  F[1] =   mun*Jnc*Jn_scale*ptri->d[2];
  F[2] = - mup*Jpc*Jp_scale*ptri->d[2];
  F[3] = - QNFactor*ptri->d[2]*bn*qV(qca,qcb,ptri->edge_len[2]);
  F[4] =   QPFactor*ptri->d[2]*bp*qV(qva,qvb,ptri->edge_len[2]); 

  for(int i=0;i<5;i++)
    for(int j=0;j<5;j++)
    {
      J1.m[5*i+j] = F[i].getADValue(0+j);
      J2.m[5*i+j] = F[i].getADValue(5+j);
      J3.m[5*i+j] = F[i].getADValue(10+j);
    }
  MatSetValues(*jtmp,5,index1,5,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index1,5,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index1,5,index3,J3.m,ADD_VALUES);
  
  F[0] = - 0.5*(aux[A].eps+aux[B].eps)*(Vb-Va)/ptri->edge_len[2]*ptri->d[2];
  F[1] = - mun*Jnc*Jn_scale*ptri->d[2];
  F[2] =   mup*Jpc*Jp_scale*ptri->d[2];
  F[3] = - QNFactor*ptri->d[2]*bn*qV(qcb,qca,ptri->edge_len[2]);
  F[4] =   QPFactor*ptri->d[2]*bp*qV(qvb,qva,ptri->edge_len[2]);  
  for(int i=0;i<5;i++)
    for(int j=0;j<5;j++)
    {
      J1.m[5*i+j] = F[i].getADValue(0+j);
      J2.m[5*i+j] = F[i].getADValue(5+j);
      J3.m[5*i+j] = F[i].getADValue(10+j);
    }
  MatSetValues(*jtmp,5,index2,5,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index2,5,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index2,5,index3,J3.m,ADD_VALUES);
  
  //---------------------------------------------------------------------------
  //flux along B-C
  //---------------------------------------------------------------------------
  if(HighFieldMobility)
  {
    if(EJModel || IIType==EdotJ)
    {
      AutoDScalar Jnax = ((Wab+Wac)*Sax - Wab*Mc*Sbx - Wac*Mb*Scx)*Jna
                       + (Wab*Sbx - Wab*Mc*Sax)*Jnb
                       + (Wac*Scx - Wac*Mb*Sax)*Jnc;
      AutoDScalar Jnay = ((Wab+Wac)*Say - Wab*Mc*Sby - Wac*Mb*Scy)*Jna
                       + (Wab*Sby - Wab*Mc*Say)*Jnb
                       + (Wac*Scy - Wac*Mb*Say)*Jnc;
      AutoDScalar Jpax = ((Wab+Wac)*Sax - Wab*Mc*Sbx - Wac*Mb*Scx)*Jpa
                       + (Wab*Sbx - Wab*Mc*Sax)*Jpb
                       + (Wac*Scx - Wac*Mb*Sax)*Jpc;
      AutoDScalar Jpay = ((Wab+Wac)*Say - Wab*Mc*Sby - Wac*Mb*Scy)*Jpa
                       + (Wab*Sby - Wab*Mc*Say)*Jpb
                       + (Wac*Scy - Wac*Mb*Say)*Jpc;
      Jna_norm = sqrt(Jnax*Jnax + Jnay*Jnay + 1e-100);
      Jpa_norm = sqrt(Jpax*Jpax + Jpay*Jpay + 1e-100);
      Epna = Ex*(Jnax/Jna_norm) + Ey*(Jnay/Jna_norm);
      Etna = Ex*(Jnay/Jna_norm) - Ey*(Jnax/Jna_norm);
      Eppa = Ex*(Jpax/Jpa_norm) + Ey*(Jpay/Jpa_norm);
      Etpa = Ex*(Jpay/Jpa_norm) - Ey*(Jpax/Jpa_norm);
    }
    else
    {
      Epna = fabs(Ecb-Ecc)/e/ptri->edge_len[0]; //parallel electrical field for electron
      Eppa = fabs(Evb-Evc)/e/ptri->edge_len[0]; //parallel electrical field for hole
      //transvers electrical field for electron and hole
      Etna = fabs(Ecb + ptri->edge_len[2]*cos(ptri->angle[1])*(Ecc-Ecb)/ptri->edge_len[0] - Eca)
           /(ptri->edge_len[2]*sin(ptri->angle[1]))/e;
      Etpa = fabs(Evb + ptri->edge_len[2]*cos(ptri->angle[1])*(Evc-Evb)/ptri->edge_len[0] - Eva)
           /(ptri->edge_len[2]*sin(ptri->angle[1]))/e;
    }
    mt->mapping(&pzone->danode[B],&aux[B],0);
    AutoDScalar munAB =  mt->mob->ElecMob(pb,nb,TD,adtl::fmax(0,Epna),fabs(Etna),TD);
    AutoDScalar mupAB =  mt->mob->HoleMob(pb,nb,TD,adtl::fmax(0,Eppa),fabs(Etpa),TD);

    mt->mapping(&pzone->danode[C],&aux[C],0);
    AutoDScalar munAC =  mt->mob->ElecMob(pc,nc,TD,adtl::fmax(0,Epna),fabs(Etna),TD);
    AutoDScalar mupAC =  mt->mob->HoleMob(pc,nc,TD,adtl::fmax(0,Eppa),fabs(Etpa),TD);
  
    mun = 0.5*(munAB+munAC);
    mup = 0.5*(mupAB+mupAC);
  }
  else 
  {  
    mun = 0.5*(aux[B].mun+aux[C].mun);
    mup = 0.5*(aux[B].mup+aux[C].mup);
  }
  F[0] =   0.5*(aux[B].eps+aux[C].eps)*(Vc-Vb)/ptri->edge_len[0]*ptri->d[0];
  F[1] =   mun*Jna*Jn_scale*ptri->d[0];
  F[2] = - mup*Jpa*Jp_scale*ptri->d[0];
  F[3] = - QNFactor*ptri->d[0]*bn*qV(qcb,qcc,ptri->edge_len[0]);
  F[4] =   QPFactor*ptri->d[0]*bp*qV(qvb,qvc,ptri->edge_len[0]);
  for(int i=0;i<5;i++)
    for(int j=0;j<5;j++)
    {
      J1.m[5*i+j] = F[i].getADValue(0+j);
      J2.m[5*i+j] = F[i].getADValue(5+j);
      J3.m[5*i+j] = F[i].getADValue(10+j);
    }
  MatSetValues(*jtmp,5,index2,5,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index2,5,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index2,5,index3,J3.m,ADD_VALUES);
  
  F[0] = - 0.5*(aux[B].eps+aux[C].eps)*(Vc-Vb)/ptri->edge_len[0]*ptri->d[0];
  F[1] = - mun*Jna*Jn_scale*ptri->d[0];
  F[2] =   mup*Jpa*Jp_scale*ptri->d[0];
  F[3] = - QNFactor*ptri->d[0]*bn*qV(qcc,qcb,ptri->edge_len[0]);
  F[4] =   QPFactor*ptri->d[0]*bp*qV(qvc,qvb,ptri->edge_len[0]);
  for(int i=0;i<5;i++)
    for(int j=0;j<5;j++)
    {
      J1.m[5*i+j] = F[i].getADValue(0+j);
      J2.m[5*i+j] = F[i].getADValue(5+j);
      J3.m[5*i+j] = F[i].getADValue(10+j);
    }
  MatSetValues(*jtmp,5,index3,5,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index3,5,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index3,5,index3,J3.m,ADD_VALUES);
  
  //---------------------------------------------------------------------------
  //flux along C-A
  //---------------------------------------------------------------------------
  if(HighFieldMobility)
  {
    if(EJModel || IIType==EdotJ)
    {
      AutoDScalar Jnbx = ((Wbc+Wba)*Sbx - Wbc*Ma*Scx - Wba*Mc*Sax)*Jnb
                       + (Wbc*Scx - Wbc*Ma*Sbx)*Jnc
                       + (Wba*Sax - Wba*Mc*Sbx)*Jna;
      AutoDScalar Jnby = ((Wbc+Wba)*Sby - Wbc*Ma*Scy - Wba*Mc*Say)*Jnb
                       + (Wbc*Scy - Wbc*Ma*Sby)*Jnc
                       + (Wba*Say - Wba*Mc*Sby)*Jna;
      AutoDScalar Jpbx = ((Wbc+Wba)*Sbx - Wbc*Ma*Scx - Wba*Mc*Sax)*Jpb
                       + (Wbc*Scx - Wbc*Ma*Sbx)*Jpc
                       + (Wba*Sax - Wba*Mc*Sbx)*Jpa;
      AutoDScalar Jpby = ((Wbc+Wba)*Sby - Wbc*Ma*Scy - Wba*Mc*Say)*Jpb
                       + (Wbc*Scy - Wbc*Ma*Sby)*Jpc
                       + (Wba*Say - Wba*Mc*Sby)*Jpa;
      Jnb_norm = sqrt(Jnbx*Jnbx + Jnby*Jnby + 1e-100);
      Jpb_norm = sqrt(Jpbx*Jpbx + Jpby*Jpby + 1e-100);
      Epnb = Ex*(Jnbx/Jnb_norm) + Ey*(Jnby/Jnb_norm);
      Etnb = Ex*(Jnby/Jnb_norm) - Ey*(Jnbx/Jnb_norm);
      Eppb = Ex*(Jpbx/Jpb_norm) + Ey*(Jpby/Jpb_norm);
      Etpb = Ex*(Jpby/Jpb_norm) - Ey*(Jpbx/Jpb_norm);
    }
    else
    {
      Epnb = fabs(Ecc-Eca)/e/ptri->edge_len[1]; //parallel electrical field for electron
      Eppb = fabs(Evc-Eva)/e/ptri->edge_len[1]; //parallel electrical field for hole
      //transvers electrical field for electron and hole
      Etnb = fabs(Ecc + ptri->edge_len[0]*cos(ptri->angle[2])*(Eca-Ecc)/ptri->edge_len[1] - Ecb)
           /(ptri->edge_len[0]*sin(ptri->angle[2]))/e;
      Etpb = fabs(Evc + ptri->edge_len[0]*cos(ptri->angle[2])*(Eva-Evc)/ptri->edge_len[1] - Evb)
           /(ptri->edge_len[0]*sin(ptri->angle[2]))/e;
    }
    mt->mapping(&pzone->danode[C],&aux[C],0);
    AutoDScalar munBC =  mt->mob->ElecMob(pc,nc,TD,adtl::fmax(0,Epnb),fabs(Etnb),TD);
    AutoDScalar mupBC =  mt->mob->HoleMob(pc,nc,TD,adtl::fmax(0,Eppb),fabs(Etpb),TD);

    mt->mapping(&pzone->danode[A],&aux[A],0);
    AutoDScalar munBA =  mt->mob->ElecMob(pa,na,TD,adtl::fmax(0,Epnb),fabs(Etnb),TD);
    AutoDScalar mupBA =  mt->mob->HoleMob(pa,na,TD,adtl::fmax(0,Eppb),fabs(Etpb),TD);
  
    mun = 0.5*(munBC+munBA);
    mup = 0.5*(mupBC+mupBA);
  }
  else
  {
    mun = 0.5*(aux[C].mun+aux[A].mun);
    mup = 0.5*(aux[C].mup+aux[A].mup);
  }  
  F[0] =   0.5*(aux[C].eps+aux[A].eps)*(Va-Vc)/ptri->edge_len[1]*ptri->d[1];
  F[1] =   mun*Jnb*Jn_scale*ptri->d[1];
  F[2] = - mup*Jpb*Jp_scale*ptri->d[1];
  F[3] = - QNFactor*ptri->d[1]*bn*qV(qcc,qca,ptri->edge_len[1]);
  F[4] =   QPFactor*ptri->d[1]*bp*qV(qvc,qva,ptri->edge_len[1]);
  for(int i=0;i<5;i++)
    for(int j=0;j<5;j++)
    {
      J1.m[5*i+j] = F[i].getADValue(0+j);
      J2.m[5*i+j] = F[i].getADValue(5+j);
      J3.m[5*i+j] = F[i].getADValue(10+j);
    }
  MatSetValues(*jtmp,5,index3,5,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index3,5,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index3,5,index3,J3.m,ADD_VALUES);
  
  F[0] = - 0.5*(aux[C].eps+aux[A].eps)*(Va-Vc)/ptri->edge_len[1]*ptri->d[1];
  F[1] = - mun*Jnb*Jn_scale*ptri->d[1];
  F[2] =   mup*Jpb*Jp_scale*ptri->d[1];
  F[3] = - QNFactor*ptri->d[1]*bn*qV(qca,qcc,ptri->edge_len[1]);
  F[4] =   QPFactor*ptri->d[1]*bp*qV(qva,qvc,ptri->edge_len[1]);
  for(int i=0;i<5;i++)
    for(int j=0;j<5;j++)
    {
      J1.m[5*i+j] = F[i].getADValue(0+j);
      J2.m[5*i+j] = F[i].getADValue(5+j);
      J3.m[5*i+j] = F[i].getADValue(10+j);
    }
  MatSetValues(*jtmp,5,index1,5,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index1,5,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,5,index1,5,index3,J3.m,ADD_VALUES);
  
  //---------------------------------------------------------------------------
  Set_Mat5_zero(J1);
  S = 0.25*ptri->d[2]*ptri->edge_len[2] + 0.25*ptri->d[1]*ptri->edge_len[1];
  J1.m[6] =  - Ra.getADValue(1)*S;   
  J1.m[7] =  - Ra.getADValue(2)*S;   
  J1.m[11] =  - Ra.getADValue(1)*S;  
  J1.m[12] =  - Ra.getADValue(2)*S;  
  MatSetValues(*jtmp,5,index1,5,index1,J1.m,ADD_VALUES);
  
  //---------------------------------------------------------------------------
  Set_Mat5_zero(J2);
  S = 0.25*ptri->d[0]*ptri->edge_len[0] + 0.25*ptri->d[2]*ptri->edge_len[2];
  J2.m[6] =  - Rb.getADValue(4)*S;   
  J2.m[7] =  - Rb.getADValue(5)*S;   
  J2.m[11] =  - Rb.getADValue(4)*S;  
  J2.m[12] =  - Rb.getADValue(5)*S;  
  MatSetValues(*jtmp,5,index2,5,index2,J2.m,ADD_VALUES);
  
  //---------------------------------------------------------------------------
  Set_Mat5_zero(J3); 
  S = 0.25*ptri->d[1]*ptri->edge_len[1] + 0.25*ptri->d[0]*ptri->edge_len[0];
  J3.m[6] =  - Rc.getADValue(7)*S;   
  J3.m[7] =  - Rc.getADValue(8)*S;   
  J3.m[11] =  - Rc.getADValue(7)*S;  
  J3.m[12] =  - Rc.getADValue(8)*S;  
  MatSetValues(*jtmp,5,index3,5,index3,J3.m,ADD_VALUES);
}

//-----------------------------------------------------------------------------
// boundaries
//-----------------------------------------------------------------------------


void SMCZone::F1Q_qddm_inner(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  PetscScalar e    =  mt->e;
  PetscScalar hbar =  mt->hbar;
  PetscScalar Vi  = x[zofs[zone_index]+5*i+0];     //potential of node i
  PetscScalar ni  = x[zofs[zone_index]+5*i+1];     //electron density of node i
  PetscScalar pi  = x[zofs[zone_index]+5*i+2];     //hole density of node i
  PetscScalar Eqc = x[zofs[zone_index]+5*i+3];
  PetscScalar Eqv = x[zofs[zone_index]+5*i+4];
   
  f[zofs[zone_index]+5*i+0] = f[zofs[zone_index]+5*i+0]/pcell->area + mt->e*((pi-ni)+aux[i].Net_doping());
  f[zofs[zone_index]+5*i+1] = f[zofs[zone_index]+5*i+1]/pcell->area;
  f[zofs[zone_index]+5*i+2] = f[zofs[zone_index]+5*i+2]/pcell->area;
  f[zofs[zone_index]+5*i+3] = f[zofs[zone_index]+5*i+3]/pcell->area + (Eqc+Vi+aux[i].affinity);
  f[zofs[zone_index]+5*i+4] = f[zofs[zone_index]+5*i+4]/pcell->area + (Eqv+Vi+aux[i].affinity+aux[i].Eg);

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+5*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+5*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+5*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+5*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F1Q_qddm_ombc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int z = zone_index;
  int equ_num = 5;
  int size = pzone->davcell.size();
  PetscScalar e  =  mt->e;
  PetscScalar hbar =  mt->hbar;
  PetscScalar Na = aux[i].Total_Na();
  PetscScalar Nd = aux[i].Total_Nd();
  PetscScalar Vi  = x[zofs[zone_index]+5*i+0];     //potential of node i
  PetscScalar Eqc = x[zofs[zone_index]+5*i+3];
  PetscScalar Eqv = x[zofs[zone_index]+5*i+4];
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar nie = mt->band->nie(fs[i].T);
  PetscScalar Nc  = mt->band->Nc(fs[i].T);
  PetscScalar electron_density,hole_density;
  
  int    om_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1) {om_equ=j;break;}
  }
  f[zofs[z]+5*i+0] = Vi- mt->kb*fs[i].T/e*asinh((Nd-Na)/(2*nie)) + aux[i].affinity + mt->kb*fs[i].T/e*log(Nc/nie)
                     -x[zofs[z]+equ_num*size+om_equ];

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
  f[zofs[z]+5*i+1] = x[zofs[z]+5*i+1] - electron_density;  //electron density
  f[zofs[z]+5*i+2] = x[zofs[z]+5*i+2] - hole_density;      //hole density
  f[zofs[z]+5*i+3] = f[zofs[z]+5*i+3]/pcell->area + (Eqc+Vi+aux[i].affinity);
  f[zofs[z]+5*i+4] = f[zofs[z]+5*i+4]/pcell->area + (Eqv+Vi+aux[i].affinity+aux[i].Eg);
  //f[zofs[z]+5*i+3] = x[zofs[z]+5*i+3] + x[zofs[z]+5*i+0]*e + aux[i].affinity*e;
  //f[zofs[z]+5*i+4] = x[zofs[z]+5*i+4] + x[zofs[z]+5*i+0]*e + aux[i].affinity*e + aux[i].Eg;

}


void SMCZone::F1Q_qddm_stkbc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int>& zofs,DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 5;
  int size = pzone->davcell.size();
  PetscScalar Vi = x[zofs[zone_index]+5*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+5*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+5*i+2];     //hole density of node i
  PetscScalar Eqc = x[zofs[zone_index]+5*i+3];
  PetscScalar Eqv = x[zofs[zone_index]+5*i+4];
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  
  int    stk_equ;
  for(int j=0;j<electrode.size();j++)
  if(electrode[j]==pcell->bc_index-1)   {stk_equ=j;break;}
  //Schotty Barrier Lowerring
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  //Schottky current
  PetscScalar Fn = mt->band->SchottyJsn(ni,fs[i].T,pbc->WorkFunction-aux[i].affinity-deltaVB)*0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar Fp = mt->band->SchottyJsp(pi,fs[i].T,pbc->WorkFunction-aux[i].affinity+deltaVB)*0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);

  f[zofs[z]+5*i+0] = x[zofs[z]+5*i+0] + pbc->WorkFunction - deltaVB - x[zofs[z]+equ_num*size+stk_equ];
  f[zofs[z]+5*i+1] = (f[zofs[z]+5*i+1]+Fn)/pcell->area;
  f[zofs[z]+5*i+2] = (f[zofs[z]+5*i+2]-Fp)/pcell->area;
  f[zofs[z]+5*i+3] = f[zofs[z]+5*i+3]/pcell->area + (Eqc+Vi+aux[i].affinity);
  f[zofs[z]+5*i+4] = f[zofs[z]+5*i+4]/pcell->area + (Eqv+Vi+aux[i].affinity+aux[i].Eg);
  
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+5*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+5*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+5*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+5*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F1Q_qddm_insulator_gate(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int equ_num = 5;
  int size = pzone->davcell.size();
  PetscScalar e  =  mt->e;
  PetscScalar eV  =  mt->pscale->s_eV;
  PetscScalar hbar =  mt->hbar;
  PetscScalar m0 =  mt->me;
  PetscScalar Vi = x[zofs[zone_index]+5*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+5*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+5*i+2];     //hole density of node i
  PetscScalar Eqc = x[zofs[zone_index]+5*i+3];
  PetscScalar Eqv = x[zofs[zone_index]+5*i+4];
  PetscScalar grad_P = 0, grad_qpn=0, grad_qpp=0;
  int    ins_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
     {ins_equ=j;break;}
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+5*nb+0];     //potential of nb node
    
    if(j==0||j==pcell->nb_num-1)
    {
      //the poisson's equation on third boundary type
      PetscScalar vgate = x[zofs[zone_index]+equ_num*size+ins_equ] - pbc->WorkFunction;
      PetscScalar q = mt->e*pbc->QF; //sigma is the surface change density
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar r=q + eps_ox/Thick*vgate;
      PetscScalar s=eps_ox/Thick;
      grad_P += 0.5*pcell->ilen[j]*(r-0.25*s*(3*Vi+Vj));
      //the WBK approx of Si/SiO2 interface
      PetscScalar b_nox = hbar*hbar/(6*e*0.14*m0);
      PetscScalar b_pox = hbar*hbar/(6*e*1.0*m0);
      PetscScalar x_np = hbar/sqrt(2*0.4*m0*3.15*eV);
      PetscScalar x_pp = hbar/sqrt(2*0.4*m0*4.50*eV);
      grad_qpn += -b_nox/x_np*0.5*pcell->ilen[j];
      grad_qpp += -b_pox/x_pp*0.5*pcell->ilen[j];
    }
  }
  
  f[zofs[zone_index]+5*i+0] = (f[zofs[zone_index]+5*i+0]+grad_P)/pcell->area
                              + mt->e*((pi-ni)+aux[i].Net_doping());
  f[zofs[zone_index]+5*i+1] = f[zofs[zone_index]+5*i+1]/pcell->area;
  f[zofs[zone_index]+5*i+2] = f[zofs[zone_index]+5*i+2]/pcell->area;
  f[zofs[zone_index]+5*i+3] = (f[zofs[zone_index]+5*i+3]+QNFactor*grad_qpn)/pcell->area + (Eqc+Vi+aux[i].affinity);
  f[zofs[zone_index]+5*i+4] = (f[zofs[zone_index]+5*i+4]-QPFactor*grad_qpp)/pcell->area + (Eqv+Vi+aux[i].affinity+aux[i].Eg);
  
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+5*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+5*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+5*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+5*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F1Q_qddm_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                               ISZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorInterfaceBC *pbc = dynamic_cast<InsulatorInterfaceBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar e  =  mt->e;
  PetscScalar eV  =  mt->pscale->s_eV;
  PetscScalar hbar =  mt->hbar;
  PetscScalar m0 =  mt->me;
  
  PetscScalar Vi = x[zofs[zone_index]+5*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+5*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+5*i+2];     //hole density of node i
  PetscScalar Eqc = x[zofs[zone_index]+5*i+3];
  PetscScalar Eqv = x[zofs[zone_index]+5*i+4];
  PetscScalar Na = aux[i].Total_Na();
  PetscScalar Nd = aux[i].Total_Nd();
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  
  PetscScalar grad_P = 0, grad_qpn=0, grad_qpp=0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+5*nb+0];     //potential of nb node
    grad_P += aux[i].eps*pcell->elen[j]/pcell->ilen[j]*(Vj-Vi);
    if(j==0||j==pcell->nb_num-1)
    {
      //the WBK approx of Si/SiO2 interface
      PetscScalar b_nox = hbar*hbar/(6*e*0.14*m0);
      PetscScalar b_pox = hbar*hbar/(6*e*1.0*m0);
      PetscScalar x_np = hbar/sqrt(2*0.4*m0*3.15*eV);
      PetscScalar x_pp = hbar/sqrt(2*0.4*m0*4.50*eV);
      grad_qpn += -b_nox/x_np*0.5*pcell->ilen[j];
      grad_qpp += -b_pox/x_pp*0.5*pcell->ilen[j];
    }
  }
  //poisson's equation at interface
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int    nb = ncell->nb_array[j];
    PetscScalar Vj_n = x[zofs[pz->pzone->zone_index]+nb];     //potential of nb node
    grad_P += pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j]*(Vj_n-Vi);
  }

  f[zofs[zone_index]+5*i+0] = grad_P + mt->e*((pi-ni)+(Nd-Na))*pcell->area + pbc->QF*L;
  f[zofs[zone_index]+5*i+1] = f[zofs[zone_index]+5*i+1]/pcell->area;
  f[zofs[zone_index]+5*i+2] = f[zofs[zone_index]+5*i+2]/pcell->area;
  f[zofs[zone_index]+5*i+3] = (f[zofs[zone_index]+5*i+3]+QNFactor*grad_qpn)/pcell->area + (Eqc+Vi+aux[i].affinity);
  f[zofs[zone_index]+5*i+4] = (f[zofs[zone_index]+5*i+4]-QPFactor*grad_qpp)/pcell->area + (Eqv+Vi+aux[i].affinity+aux[i].Eg);
  
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+5*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+5*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+5*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+5*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F1Q_om_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 5;
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
    PetscScalar Vi = x[zofs[zone_index]+5*node+0];     //potential of node i
    //conduction current
    current += DeviceDepth*(f[zofs[zone_index]+5*node+1]-f[zofs[zone_index]+5*node+2]);
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+5*nb+0];     //potential of nb node
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


void SMCZone::F1Q_stk_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 5;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc.Get_pointer(bc_index));

  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+5*node+0];     //potential of node i
    PetscScalar ni = x[zofs[zone_index]+5*node+1];     //electron density of node i
    PetscScalar pi = x[zofs[zone_index]+5*node+2];     //hole density of node i
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
      PetscScalar Vj = x[zofs[zone_index]+5*nb+0];     //potential of nb node
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


void SMCZone::F1Q_ins_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 5;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  InsulatorContactBC *pbc = dynamic_cast <InsulatorContactBC * > (bc.Get_pointer(bc_index));

  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+5*node+0];     //potential of node i
    for(int k=0;k<pcell->nb_num;k++)
    {
      int  nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+5*nb+0];     //potential of nb node
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



void SMCZone::J1Q_qddm_inner(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs)
{
  Mat5      A;
  PetscInt  index[5],col[5];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int   z = zone_index;

  PetscScalar Vi = x[zofs[z]+5*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+5*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+5*i+2];     //hole density of node i
  PetscScalar area = pcell->area;
  PetscScalar e  =  mt->e;

  //--------------------------------
  index[0] = zofs[z]+5*i+0;
  index[1] = zofs[z]+5*i+1;
  index[2] = zofs[z]+5*i+2;
  index[3] = zofs[z]+5*i+3;
  index[4] = zofs[z]+5*i+4;
  
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    col[0] = zofs[z]+5*nb+0;
    col[1] = zofs[z]+5*nb+1;
    col[2] = zofs[z]+5*nb+2;
    col[3] = zofs[z]+5*nb+3;
    col[4] = zofs[z]+5*nb+4;
    MatGetValues(*jtmp,5,index,5,col,A.m);
    A=A/area;
    MatSetValues(*jac,5,index,5,col,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,5,index,5,index,A.m);
  A=A/area;

  A.m[1] +=  -mt->e;   //dfun(0)/dn(i)
  A.m[2] +=   mt->e;   //dfun(0)/dp(i)

  A.m[15] += 1;
  A.m[18] += 1;
  
  A.m[20] += 1;
  A.m[24] += 1;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[6] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[12] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[6]  += -1/ODE_F.dt;
      A.m[12] += -1/ODE_F.dt;
    }
  }

  MatSetValues(*jac,5,index,5,index,A.m,INSERT_VALUES);
}


void SMCZone::J1Q_qddm_ombc(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  Mat5      A;
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  PetscScalar e  =  mt->e;
  PetscScalar area = pcell->area;
  PetscInt   index[5],col[5];
  index[0] = zofs[zone_index]+5*i+0;
  index[1] = zofs[zone_index]+5*i+1;
  index[2] = zofs[zone_index]+5*i+2;
  index[3] = zofs[zone_index]+5*i+3;
  index[4] = zofs[zone_index]+5*i+4;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    col[0] = zofs[zone_index]+5*nb+0;
    col[1] = zofs[zone_index]+5*nb+1;
    col[2] = zofs[zone_index]+5*nb+2;
    col[3] = zofs[zone_index]+5*nb+3;
    col[4] = zofs[zone_index]+5*nb+4;
    MatGetValues(*jtmp,5,index,5,col,A.m);
    A=A/area;
    A.m[0]=0.0;A.m[1]=0.0;A.m[2]=0.0;A.m[3]=0.0;A.m[4]=0.0;
    A.m[5]=0.0;A.m[6]=0.0;A.m[7]=0.0;A.m[8]=0.0;A.m[9]=0.0;
    A.m[10]=0.0;A.m[11]=0.0;A.m[12]=0.0;A.m[13]=0.0;A.m[14]=0.0;
    MatSetValues(*jac,5,index,5,col,A.m,INSERT_VALUES);
  }
    
  MatGetValues(*jtmp,5,index,5,index,A.m);
  A=A/area;
  
  A.m[0] =  1.0;   
  A.m[1] =  1.0;  
  A.m[2] =  1.0;  
  
  A.m[5] =  0.0;   
  A.m[6] =  1.0;   
  A.m[7] =  0.0; 
  A.m[8] =  0.0; 
  A.m[9] =  0.0; 

  A.m[10] =  0.0;   
  A.m[11] =  0.0;   
  A.m[12] =  1.0; 
  A.m[13] =  0.0; 
  A.m[14] =  0.0; 

  A.m[15] += 1;
  A.m[18] += 1;
  
  A.m[20] += 1;
  A.m[24] += 1;
  
  MatSetValues(*jac,5,index,5,index,A.m,INSERT_VALUES);
  int om_equ;
  int equ_num = 5;
  int size = pzone->davcell.size();
  for(int j=0;j<electrode.size();j++)
  if(electrode[j]==pcell->bc_index-1)   {om_equ=j;break;}
  MatSetValue(*jac,zofs[zone_index]+5*i,zofs[zone_index]+equ_num*size+om_equ,-1,INSERT_VALUES);
}


void SMCZone::J1Q_qddm_stkbc(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  Mat5    A;
  PetscInt       index[5],col[5];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int  z = zone_index;
  int  equ_num = 5;
  int  size = pzone->davcell.size();
  PetscScalar VB = pbc->WorkFunction-aux[i].affinity;
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  PetscScalar Vi = x[zofs[z]+5*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+5*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+5*i+2];     //hole density of node i

  PetscScalar area = pcell->area;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fp_dpi = 0;
  //--------------------------------
  index[0] = zofs[z]+5*i+0;
  index[1] = zofs[z]+5*i+1;
  index[2] = zofs[z]+5*i+2;
  index[3] = zofs[z]+5*i+3;
  index[4] = zofs[z]+5*i+4;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    col[0] = zofs[z]+5*nb+0;
    col[1] = zofs[z]+5*nb+1;
    col[2] = zofs[z]+5*nb+2;
    col[3] = zofs[z]+5*nb+3;
    col[4] = zofs[z]+5*nb+4;
    MatGetValues(*jtmp,5,index,5,col,A.m);
    A.m[0] = 0.0;
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      d_Fn_dni += mt->band->pdSchottyJsn_pdn(ni,fs[i].T,VB-deltaVB)*0.5*pcell->ilen[j];
      d_Fp_dpi += mt->band->pdSchottyJsp_pdp(pi,fs[i].T,VB+deltaVB)*0.5*pcell->ilen[j];
    }
    MatSetValues(*jac,5,index,5,col,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,5,index,5,index,A.m);
  A=A/area;
  A.m[0] =  1.0;                      //dfun(0)/dP(i)
  A.m[1] =  0.0;                      //dfun(0)/dn(i)
  A.m[2] =  0.0;                      //dfun(0)/dp(i)

  A.m[6]  +=   d_Fn_dni/area;          //dfun(1)/dn(i)
  A.m[12] +=  -d_Fp_dpi/area;          //dfun(2)/dp(i)
  
  A.m[15] += 1;
  A.m[18] += 1;
  
  A.m[20] += 1;
  A.m[24] += 1;
  
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[6] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[12] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[6] += -1/ODE_F.dt;
      A.m[12] += -1/ODE_F.dt;
    }
  }
  MatSetValues(*jac,5,index,5,index,A.m,INSERT_VALUES);

  int stk_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
      {stk_equ=j;break;}
  MatSetValue(*jac,zofs[zone_index]+5*i,zofs[zone_index]+equ_num*size+stk_equ,-1,INSERT_VALUES);
}


void SMCZone::J1Q_qddm_insulator_gate(int i,PetscScalar *x, Mat *jac, Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  Mat5    A;
  PetscInt       index[5],col[5];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int  equ_num = 5;
  int  size = pzone->davcell.size();

  PetscScalar ni = x[zofs[z]+5*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+5*i+2];     //hole density of node i
  PetscScalar e  =  mt->e;
  PetscScalar area = pcell->area;
  PetscScalar d_grad_P_dVi = 0;
  PetscScalar d_grad_P_dVapp = 0;
  int ins_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
      {ins_equ=j;break;}
  //--------------------------------
  index[0] = zofs[z]+5*i+0;
  index[1] = zofs[z]+5*i+1;
  index[2] = zofs[z]+5*i+2;
  index[3] = zofs[z]+5*i+3;
  index[4] = zofs[z]+5*i+4;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    col[0] = zofs[z]+5*nb+0;
    col[1] = zofs[z]+5*nb+1;
    col[2] = zofs[z]+5*nb+2;
    col[3] = zofs[z]+5*nb+3;
    col[4] = zofs[z]+5*nb+4;
    MatGetValues(*jtmp,5,index,5,col,A.m);
    A=A/area;
   
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar s=eps_ox/Thick;
      d_grad_P_dVi += -0.5*pcell->ilen[j]*0.25*s*3;
      d_grad_P_dVapp += 0.5*pcell->ilen[j]*s/area;
      A.m[0]+= -0.5*pcell->ilen[j]*0.25*s/area;
    }
    MatSetValues(*jac,5,index,5,col,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,5,index,5,index,A.m);
  A=A/area;
  A.m[0] +=  d_grad_P_dVi/area;             //dfun(0)/dP(i)
  A.m[1] +=  -mt->e;                        //dfun(0)/dn(i)
  A.m[2] +=   mt->e;                        //dfun(0)/dp(i)

  A.m[15] += 1;
  A.m[18] += 1;
  
  A.m[20] += 1;
  A.m[24] += 1;
 
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[6] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[12] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[6] += -1/ODE_F.dt;
      A.m[12] += -1/ODE_F.dt;
    }
  }
  MatSetValues(*jac,5,index,5,index,A.m,INSERT_VALUES);
  MatSetValue(*jac,zofs[z]+5*i,zofs[z]+equ_num*size+ins_equ,d_grad_P_dVapp,INSERT_VALUES);
}


void SMCZone::J1Q_qddm_interface(int i,PetscScalar *x,Mat *jac, Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                               ISZone *pz, int n)
{
  Mat5    A;
  PetscInt       index[5],col[5];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int z = zone_index;
  PetscScalar Vi = x[zofs[z]+5*i+0];     //potential of node i
  PetscScalar area = pcell->area;
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar d_grad_P_dVi = 0;
  //--------------------------------
  index[0] = zofs[z]+5*i+0;
  index[1] = zofs[z]+5*i+1;
  index[2] = zofs[z]+5*i+2;
  index[3] = zofs[z]+5*i+3;
  index[4] = zofs[z]+5*i+4;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[z]+5*nb+0];     //potential of nb node
    col[0] = zofs[z]+5*nb+0;
    col[1] = zofs[z]+5*nb+1;
    col[2] = zofs[z]+5*nb+2;
    col[3] = zofs[z]+5*nb+3;
    col[4] = zofs[z]+5*nb+4;
    //-------------------------------------
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    MatGetValues(*jtmp,5,index,5,col,A.m);
    A=A/area;
    //-------------------------------------
    //df(x)i/dx(r)
    A.m[0] =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    A.m[1] =  0;
    A.m[2] =  0;
    MatSetValues(*jac,5,index,5,col,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,5,index,5,index,A.m);
  A=A/area;
  
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  int column = zofs[pz->pzone->zone_index]+n;
  for(int j=0;j<ncell->nb_num;j++)
  {
    int    nb = ncell->nb_array[j];
    PetscScalar Vr_n = x[zofs[pz->pzone->zone_index]+nb];     //potential of nb node
    d_grad_P_dVi += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    PetscScalar value = pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[z]+5*i+0,zofs[pz->pzone->zone_index]+nb,value,INSERT_VALUES);
  }
  
  A.m[0] =  d_grad_P_dVi;  //dfun(0)/dP(i)
  A.m[1] =  -mt->e*area;   //dfun(0)/dn(i)
  A.m[2] =   mt->e*area;   //dfun(0)/dp(i)
  
  A.m[15] += 1;
  A.m[18] += 1;
  
  A.m[20] += 1;
  A.m[24] += 1;
  
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[6] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[12] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[6] += -1/ODE_F.dt;
      A.m[12] += -1/ODE_F.dt;
    }
  }

  MatSetValues(*jac,5,index,5,index,A.m,INSERT_VALUES);

}


void SMCZone::J1Q_om_electrode(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 5;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(bc_index));
  PetscScalar    A1[5],A2[5],J[5];
  PetscInt       index[5],col[5];
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
    index[0] = zofs[zone_index]+5*node+0;
    index[1] = zofs[zone_index]+5*node+1;
    index[2] = zofs[zone_index]+5*node+2;
    index[3] = zofs[zone_index]+5*node+3;
    index[4] = zofs[zone_index]+5*node+4;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      col[0] = zofs[zone_index]+5*nb+0;
      col[1] = zofs[zone_index]+5*nb+1;
      col[2] = zofs[zone_index]+5*nb+2;
      col[3] = zofs[zone_index]+5*nb+3;
      col[4] = zofs[zone_index]+5*nb+4;

      MatGetValues(*jtmp,1,&index[1],5,col,A1);
      MatGetValues(*jtmp,1,&index[2],5,col,A2);

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
          J[3]=(A1[3]-A2[3])*DeviceDepth;
          J[4]=(A1[4]-A2[4])*DeviceDepth;
          MatSetValues(*jac,1,&connect_index,5,col,J,ADD_VALUES);
        }
        else
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
          J[3]=(A1[3]-A2[3])*DeviceDepth;
          J[4]=(A1[4]-A2[4])*DeviceDepth;
          MatSetValues(*jac,1,&matrix_row,5,col,J,ADD_VALUES);
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
          J[3]=(A1[3]-A2[3])*DeviceDepth*(L/ODE_F.dt+R);
          J[4]=(A1[4]-A2[4])*DeviceDepth*(L/ODE_F.dt+R);
        }
        else
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
          J[3]=(A1[3]-A2[3])*DeviceDepth;
          J[4]=(A1[4]-A2[4])*DeviceDepth;
        }  
#else
        J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth*(L/ODE_F.dt+R);
        J[1]=(A1[1]-A2[1])*DeviceDepth*(L/ODE_F.dt+R);
        J[2]=(A1[2]-A2[2])*DeviceDepth*(L/ODE_F.dt+R);
        J[3]=(A1[3]-A2[3])*DeviceDepth*(L/ODE_F.dt+R);
        J[4]=(A1[4]-A2[4])*DeviceDepth*(L/ODE_F.dt+R);
#endif        
        MatSetValues(*jac,1,&matrix_row,5,col,J,ADD_VALUES);
      }
      //the electrode has an external current circuit
      else if(pbc->electrode_type==CurrentBC)
      {
        J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
        J[1]=(A1[1]-A2[1])*DeviceDepth;
        J[2]=(A1[2]-A2[2])*DeviceDepth;
        J[3]=(A1[3]-A2[3])*DeviceDepth;
        J[4]=(A1[4]-A2[4])*DeviceDepth;
        MatSetValues(*jac,1,&matrix_row,5,col,J,ADD_VALUES);
      }
    }

    MatGetValues(*jtmp,1,&index[1],5,index,A1);
    MatGetValues(*jtmp,1,&index[2],5,index,A2);
    
    if(pbc->inner_connect!=-1)
    {
        int connect_to = pbc->inner_connect;
        if(bc_index < connect_to)
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
          J[3]=(A1[3]-A2[3])*DeviceDepth;
          J[4]=(A1[4]-A2[4])*DeviceDepth;
          MatSetValues(*jac,1,&connect_index,5,index,J,ADD_VALUES);
        }
        else
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
          J[3]=(A1[3]-A2[3])*DeviceDepth;
          J[4]=(A1[4]-A2[4])*DeviceDepth;
          MatSetValues(*jac,1,&matrix_row,5,index,J,ADD_VALUES);
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
        J[3]=(A1[3]-A2[3])*DeviceDepth*(L/ODE_F.dt+R);
        J[4]=(A1[4]-A2[4])*DeviceDepth*(L/ODE_F.dt+R);
      }
      else
      {
        J[0]=(A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth;
        J[1]=(A1[1]-A2[1])*DeviceDepth;
        J[2]=(A1[2]-A2[2])*DeviceDepth;
        J[3]=(A1[3]-A2[3])*DeviceDepth;
        J[4]=(A1[4]-A2[4])*DeviceDepth;
      } 
#else
       J[0]=(A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth*(L/ODE_F.dt+R);
       J[1]=(A1[1]-A2[1])*DeviceDepth*(L/ODE_F.dt+R);
       J[2]=(A1[2]-A2[2])*DeviceDepth*(L/ODE_F.dt+R);
       J[3]=(A1[3]-A2[3])*DeviceDepth*(L/ODE_F.dt+R);
       J[4]=(A1[4]-A2[4])*DeviceDepth*(L/ODE_F.dt+R);
#endif          
      MatSetValues(*jac,1,&matrix_row,5,index,J,ADD_VALUES);
    }
    else if(pbc->electrode_type==CurrentBC)
    {
      J[0]=(A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth;
      J[1]=(A1[1]-A2[1])*DeviceDepth;
      J[2]=(A1[2]-A2[2])*DeviceDepth;
      J[3]=(A1[3]-A2[3])*DeviceDepth;
      J[4]=(A1[4]-A2[4])*DeviceDepth;
      MatSetValues(*jac,1,&matrix_row,5,index,J,ADD_VALUES);
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


void SMCZone::J1Q_stk_electrode(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 5;
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
      PetscScalar ni = x[zofs[zone_index]+5*node+1];     //electron density of node i
      PetscScalar pi = x[zofs[zone_index]+5*node+2];     //hole density of node i
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
        PetscScalar Vj = x[zofs[zone_index]+5*nb+0];     //potential of nb node
        PetscScalar dI_dVi =  DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
        PetscScalar dI_dVr = -DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
        if(pbc->electrode_type==VoltageBC)
        {
          MatSetValue(*jac,matrix_row,zofs[zone_index]+5*node+0,dI_dVi*(L/ODE_F.dt+R),ADD_VALUES);
          MatSetValue(*jac,matrix_row,zofs[zone_index]+5*nb+0,dI_dVr*(L/ODE_F.dt+R),ADD_VALUES);
        }
        else if(pbc->electrode_type==CurrentBC)
        {
          MatSetValue(*jac,matrix_row,zofs[zone_index]+5*node+0,dI_dVi,ADD_VALUES);
          MatSetValue(*jac,matrix_row,zofs[zone_index]+5*nb+0,dI_dVr,ADD_VALUES);
        }
      }
      MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

      if(pbc->electrode_type==VoltageBC)
      {
        MatSetValue(*jac,matrix_row,zofs[zone_index]+5*node+1,DeviceDepth*dJ_dni*(L/ODE_F.dt+R),INSERT_VALUES);
        MatSetValue(*jac,matrix_row,zofs[zone_index]+5*node+2,DeviceDepth*dJ_dpi*(L/ODE_F.dt+R),INSERT_VALUES);
      }
      else if(pbc->electrode_type==CurrentBC)
      {
         MatSetValue(*jac,matrix_row,zofs[zone_index]+5*node+1,DeviceDepth*dJ_dni,INSERT_VALUES);
         MatSetValue(*jac,matrix_row,zofs[zone_index]+5*node+2,DeviceDepth*dJ_dpi,INSERT_VALUES);
      }
  }
  if(pbc->electrode_type==VoltageBC)
    MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP
}


void SMCZone::J1Q_ins_electrode(int i,PetscScalar *x,Mat *jac, Mat *jtmp,ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 5;
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
          MatSetValue(*jac,matrix_row,zofs[zone_index]+5*node+0,dI_dVi*(L/ODE_F.dt+R),ADD_VALUES);
          MatSetValue(*jac,matrix_row,zofs[zone_index]+5*nb+0,dI_dVr*(L/ODE_F.dt+R),ADD_VALUES);
        }
      }
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  if(pbc->electrode_type==VoltageBC)
    MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP

}

void SMCZone::F1Q_efield_update(PetscScalar *x,vector<int> & zofs, DABC &bc, vector<BZoneData *>zonedata)
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
      PetscScalar dP = x[zofs[zone_index]+5*pcell->nb_array[k]+0] - x[zofs[zone_index]+5*i+0];
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

