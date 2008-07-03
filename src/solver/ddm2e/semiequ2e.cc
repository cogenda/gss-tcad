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
/*  Last update: Jan 15, 2007                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "mathfunc.h"
#include "jflux2e.h"
#include "zonedata.h"
#include "vec4.h"

#define _HEAT_
#define _HMOB_
#define _II_
#define _IIA_
#define _IIB_
#define _IIC_
#define _FERMI_
#define _APPROXIMATE_


void SMCZone::F2E_Tri_ddm(Tri *ptri,PetscScalar *x,PetscScalar *f, vector<int> &zofs)
{
  PetscScalar grad_P,kapa,grad_T,H=0;
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
  PetscScalar kb =  mt->kb;
  PetscScalar S;
  
  PetscScalar Va = x[zofs[zone_index]+4*A+0];     //potential of node A
  PetscScalar na = x[zofs[zone_index]+4*A+1];     //electron density of node A
  PetscScalar pa = x[zofs[zone_index]+4*A+2];     //hole density of node A
  PetscScalar Ta = x[zofs[zone_index]+4*A+3];     //Temperature of node A
  mt->mapping(&pzone->danode[A],&aux[A],0);
  PetscScalar nia = mt->band->nie(Ta);
  PetscScalar Ra = mt->band->Recomb(pa,na,Ta);
#ifdef _APPROXIMATE_
  PetscScalar Nca = mt->band->Nc(fs[A].T);
  PetscScalar Nva = mt->band->Nv(fs[A].T);
  PetscScalar Ega = mt->band->Eg(fs[A].T);
  PetscScalar Eca = -(e*Va + aux[A].affinity + mt->band->EgNarrowToEc(fs[A].T) + kb*fs[A].T*(log(Nca)-1.5*log(fs[A].T)));
  PetscScalar Eva = -(e*Va + aux[A].affinity - mt->band->EgNarrowToEv(fs[A].T) - kb*fs[A].T*(log(Nva)-1.5*log(fs[A].T)) + Ega);
#else  
  PetscScalar Nca = mt->band->Nc(Ta);
  PetscScalar Nva = mt->band->Nv(Ta);
  PetscScalar Ega = mt->band->Eg(Ta);
  PetscScalar Eca = -(e*Va + aux[A].affinity + mt->band->EgNarrowToEc(Ta) + kb*Ta*(log(Nca)-1.5*log(Ta)));
  PetscScalar Eva = -(e*Va + aux[A].affinity - mt->band->EgNarrowToEv(Ta) - kb*Ta*(log(Nva)-1.5*log(Ta)) + Ega);
#endif  
#ifdef _FERMI_
  if(Fermi)
  {
    Eca = Eca - kb*Ta*log(gamma_f(fabs(na)/Nca));
    Eva = Eva + kb*Ta*log(gamma_f(fabs(pa)/Nva));
  }
#endif
  PetscScalar phina = Va - log(fabs(na)/nia)*kb*Ta/e;
  PetscScalar phipa = Va + log(fabs(pa)/nia)*kb*Ta/e;

  PetscScalar Vb = x[zofs[zone_index]+4*B+0];     //potential of node B
  PetscScalar nb = x[zofs[zone_index]+4*B+1];     //electron density of node B
  PetscScalar pb = x[zofs[zone_index]+4*B+2];     //hole density of node B
  PetscScalar Tb = x[zofs[zone_index]+4*B+3];     //Temperature of node B
  mt->mapping(&pzone->danode[B],&aux[B],0);
  PetscScalar nib = mt->band->nie(Tb);
  PetscScalar Rb = mt->band->Recomb(pb,nb,Tb);
#ifdef _APPROXIMATE_
  PetscScalar Ncb = mt->band->Nc(fs[B].T);
  PetscScalar Nvb = mt->band->Nv(fs[B].T);
  PetscScalar Egb = mt->band->Eg(fs[B].T);
  PetscScalar Ecb =  -(e*Vb + aux[B].affinity + mt->band->EgNarrowToEc(fs[B].T) + kb*fs[B].T*(log(Ncb)-1.5*log(fs[B].T)));
  PetscScalar Evb =  -(e*Vb + aux[B].affinity - mt->band->EgNarrowToEv(fs[B].T) - kb*fs[B].T*(log(Nvb)-1.5*log(fs[B].T))+ Egb);
#else    
  PetscScalar Ncb = mt->band->Nc(Tb);
  PetscScalar Nvb = mt->band->Nv(Tb);
  PetscScalar Egb = mt->band->Eg(Tb);
  PetscScalar Ecb =  -(e*Vb + aux[B].affinity + mt->band->EgNarrowToEc(Tb) + kb*Tb*(log(Ncb)-1.5*log(Tb)));
  PetscScalar Evb =  -(e*Vb + aux[B].affinity - mt->band->EgNarrowToEv(Tb) - kb*Tb*(log(Nvb)-1.5*log(Tb))+ Egb);

#endif
#ifdef _FERMI_
  if(Fermi)
  {
    Ecb = Ecb - kb*Tb*log(gamma_f(fabs(nb)/Ncb));
    Evb = Evb + kb*Tb*log(gamma_f(fabs(pb)/Nvb));
  }
#endif
  PetscScalar phinb = Vb - log(fabs(nb)/nib)*kb*Tb/e;
  PetscScalar phipb = Vb + log(fabs(pb)/nib)*kb*Tb/e;

  PetscScalar Vc = x[zofs[zone_index]+4*C+0];     //potential of node C
  PetscScalar nc = x[zofs[zone_index]+4*C+1];     //electron density of node C
  PetscScalar pc = x[zofs[zone_index]+4*C+2];     //hole density of node C
  PetscScalar Tc = x[zofs[zone_index]+4*C+3];     //Temperature of node C
  mt->mapping(&pzone->danode[C],&aux[C],0);
  PetscScalar nic = mt->band->nie(Tc);
  PetscScalar Rc = mt->band->Recomb(pc,nc,Tc);
#ifdef _APPROXIMATE_  
  PetscScalar Ncc = mt->band->Nc(fs[C].T);
  PetscScalar Nvc = mt->band->Nv(fs[C].T);
  PetscScalar Egc = mt->band->Eg(fs[C].T);
  PetscScalar Ecc =  -(e*Vc + aux[C].affinity + mt->band->EgNarrowToEc(fs[C].T) + kb*fs[C].T*(log(Ncc)-1.5*log(fs[C].T)));
  PetscScalar Evc =  -(e*Vc + aux[C].affinity - mt->band->EgNarrowToEv(fs[C].T) - kb*fs[C].T*(log(Nvc)-1.5*log(fs[C].T))+ Egc);
#else
  PetscScalar Ncc = mt->band->Nc(Tc);
  PetscScalar Nvc = mt->band->Nv(Tc);
  PetscScalar Egc = mt->band->Eg(Tc);
  PetscScalar Ecc =  -(e*Vc + aux[C].affinity + mt->band->EgNarrowToEc(Tc) + kb*Tc*(log(Ncc)-1.5*log(Tc)));
  PetscScalar Evc =  -(e*Vc + aux[C].affinity - mt->band->EgNarrowToEv(Tc) - kb*Tc*(log(Nvc)-1.5*log(Tc))+ Egc);
#endif  
#ifdef _FERMI_
  if(Fermi)
  {
    Ecc = Ecc - kb*Tc*log(gamma_f(fabs(nc)/Ncc));
    Evc = Evc + kb*Tc*log(gamma_f(fabs(pc)/Nvc));
  }
#endif
  PetscScalar phinc = Vc - log(fabs(nc)/nic)*kb*Tc/e;
  PetscScalar phipc = Vc + log(fabs(pc)/nic)*kb*Tc/e;

  PetscScalar Ex = -((yb-yc)*Va + (yc-ya)*Vb +(ya-yb)*Vc)/(2*tri_area);
  PetscScalar Ey = -((xc-xb)*Va + (xa-xc)*Vb +(xb-xa)*Vc)/(2*tri_area);
  PetscScalar E = sqrt(Ex*Ex+Ey*Ey+1e-20);
  PetscScalar Eg =  mt->band->Eg((Ta+Tb+Tc)/3.0);
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
#ifdef _II_
  if(ImpactIonization && IIType==GradQf)
  {
    PetscScalar Fnx = ((yb-yc)*phina + (yc-ya)*phinb +(ya-yb)*phinc)/(2*tri_area);
    PetscScalar Fny = ((xc-xb)*phina + (xa-xc)*phinb +(xb-xa)*phinc)/(2*tri_area);
    PetscScalar Fpx = ((yb-yc)*phipa + (yc-ya)*phipb +(ya-yb)*phipc)/(2*tri_area);
    PetscScalar Fpy = ((xc-xb)*phipa + (xa-xc)*phipb +(xb-xa)*phipc)/(2*tri_area);
    PetscScalar Fn  =  sqrt(Fnx*Fnx+Fny*Fny+1e-20);
    PetscScalar Fp  =  sqrt(Fpx*Fpx+Fpy*Fpy+1e-20);
    IIn =  mt->gen->ElecGenRate((Ta+Tb+Tc)/3.0,Fn,Eg);
    IIp =  mt->gen->HoleGenRate((Ta+Tb+Tc)/3.0,Fp,Eg);
  }

  if(ImpactIonization && IIType==EVector)
  {
    PetscScalar Fnx = Ex;
    PetscScalar Fny = Ey;
    PetscScalar Fpx = Ex;
    PetscScalar Fpy = Ey;
    PetscScalar Fn  =  sqrt(Fnx*Fnx+Fny*Fny+1e-20);
    PetscScalar Fp  =  sqrt(Fpx*Fpx+Fpy*Fpy+1e-20);
    IIn =  mt->gen->ElecGenRate((Ta+Tb+Tc)/3.0,Fn,Eg);
    IIp =  mt->gen->HoleGenRate((Ta+Tb+Tc)/3.0,Fp,Eg);
  }
#endif
  PetscScalar Jnc =  In(kb,e,(Eca-Ecb)/e,na,nb,0.5*(Ta+Tb),Tb-Ta,ptri->edge_len[2]);
  PetscScalar Jpc =  Ip(kb,e,(Eva-Evb)/e,pa,pb,0.5*(Ta+Tb),Tb-Ta,ptri->edge_len[2]);
  PetscScalar Jna =  In(kb,e,(Ecb-Ecc)/e,nb,nc,0.5*(Tb+Tc),Tc-Tb,ptri->edge_len[0]);
  PetscScalar Jpa =  Ip(kb,e,(Evb-Evc)/e,pb,pc,0.5*(Tb+Tc),Tc-Tb,ptri->edge_len[0]);
  PetscScalar Jnb =  In(kb,e,(Ecc-Eca)/e,nc,na,0.5*(Tc+Ta),Ta-Tc,ptri->edge_len[1]);
  PetscScalar Jpb =  Ip(kb,e,(Evc-Eva)/e,pc,pa,0.5*(Tc+Ta),Ta-Tc,ptri->edge_len[1]);
  PetscScalar Jn_scale = dmax(dmax(fabs(Jna),fabs(Jnb)),dmax(fabs(Jnc),nia*nia));
  PetscScalar Jp_scale = dmax(dmax(fabs(Jpa),fabs(Jpb)),dmax(fabs(Jpc),nia*nia));
  Jnc /=  Jn_scale;
  Jpc /=  Jp_scale;
  Jna /=  Jn_scale;
  Jpa /=  Jp_scale;
  Jnb /=  Jn_scale;
  Jpb /=  Jp_scale;


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
    aux[A].mun =  mt->mob->ElecMob(pa,na,Ta,dmax(0,Epnc),fabs(Etnc),Ta);
    aux[A].mup =  mt->mob->HoleMob(pa,na,Ta,dmax(0,Eppc),fabs(Etpc),Ta);

    mt->mapping(&pzone->danode[B],&aux[B],0);
    aux[B].mun =  mt->mob->ElecMob(pb,nb,Tb,dmax(0,Epnc),fabs(Etnc),Tb);
    aux[B].mup =  mt->mob->HoleMob(pb,nb,Tb,dmax(0,Eppc),fabs(Etpc),Tb);
  }

  mun = 0.5*(aux[A].mun+aux[B].mun);
  mup = 0.5*(aux[A].mup+aux[B].mup);

  grad_P = (Vb-Va)/ptri->edge_len[2];
  kapa =  mt->thermal->HeatConduction(0.5*(Ta+Tb));
  grad_T = kapa*(Tb-Ta)/ptri->edge_len[2];
#ifdef _IIC_
  if(ImpactIonization && IIType==EdotJ)
  {
    IIn =  mt->gen->ElecGenRate(0.5*(Ta+Tb),dmax(0,Epnc),Eg);
    IIp =  mt->gen->HoleGenRate(0.5*(Ta+Tb),dmax(0,Eppc),Eg);
    G   =  IIn*mun*Jnc_norm*Jn_scale + IIp*mup*Jpc_norm*Jp_scale;
  }
  if( ImpactIonization && (IIType==GradQf || IIType==EVector) )
  {
    G = IIn*mun*fabs(phina-phinb)/ptri->edge_len[2]*nmid(kb,e,(Eca-Ecb)/e,na,nb,0.5*(Ta+Tb),Tb-Ta)
        +IIp*mup*fabs(phipa-phipb)/ptri->edge_len[2]*pmid(kb,e,(Eva-Evb)/e,pa,pb,0.5*(Ta+Tb),Tb-Ta);
    //G = IIn*mun*fabs(Va-Vb)/ptri->edge_len[2]*nmid(kb,e,(Eca-Ecb)/e,na,nb,T_mid,Tb-Ta)/e
    //   +IIp*mup*fabs(Va-Vb)/ptri->edge_len[2]*pmid(kb,e,(Eva-Evb)/e,pa,pb,T_mid,Tb-Ta)/e;
  }
  if( ImpactIonization && IIType==ESide )
  {
    IIn =  mt->gen->ElecGenRate(0.5*(Ta+Tb),fabs(Va-Vb)/ptri->edge_len[2],Eg);
    IIp =  mt->gen->HoleGenRate(0.5*(Ta+Tb),fabs(Va-Vb)/ptri->edge_len[2],Eg);
    G = IIn*mun*fabs(Jnc)*Jn_scale + IIp*mup*fabs(Jpc)*Jp_scale;
  }
#endif


#ifdef _HEAT_
  H = 0.5*(Va-Vb)*(mun*Jnc*Jn_scale+mup*Jpc*Jp_scale);
#endif

  f[zofs[zone_index]+4*A+0] +=   grad_P*ptri->d[2];
  f[zofs[zone_index]+4*A+1] +=   mun*Jnc*Jn_scale*ptri->d[2] + (G-Ra)*0.25*ptri->d[2]*ptri->edge_len[2];
  f[zofs[zone_index]+4*A+2] += - mup*Jpc*Jp_scale*ptri->d[2] + (G-Ra)*0.25*ptri->d[2]*ptri->edge_len[2];
  f[zofs[zone_index]+4*A+3] +=   ( grad_T + H -  G*0.25*ptri->edge_len[2]*(Eg+3*kb*0.5*(Ta+Tb)))*ptri->d[2];


  f[zofs[zone_index]+4*B+0] += - grad_P*ptri->d[2];
  f[zofs[zone_index]+4*B+1] += - mun*Jnc*Jn_scale*ptri->d[2] + (G-Rb)*0.25*ptri->d[2]*ptri->edge_len[2];
  f[zofs[zone_index]+4*B+2] +=   mup*Jpc*Jp_scale*ptri->d[2] + (G-Rb)*0.25*ptri->d[2]*ptri->edge_len[2];
  f[zofs[zone_index]+4*B+3] +=   (-grad_T + H -  G*0.25*ptri->edge_len[2]*(Eg+3*kb*0.5*(Ta+Tb)))*ptri->d[2];


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
    aux[B].mun =  mt->mob->ElecMob(pb,nb,Tb,dmax(0,Epna),fabs(Etna),Tb);
    aux[B].mup =  mt->mob->HoleMob(pb,nb,Tb,dmax(0,Eppa),fabs(Etpa),Tb);

    mt->mapping(&pzone->danode[C],&aux[C],0);
    aux[C].mun =  mt->mob->ElecMob(pc,nc,Tc,dmax(0,Epna),fabs(Etna),Tc);
    aux[C].mup =  mt->mob->HoleMob(pc,nc,Tc,dmax(0,Eppa),fabs(Etpa),Tc);
  }

  mun = 0.5*(aux[B].mun+aux[C].mun);
  mup = 0.5*(aux[B].mup+aux[C].mup);
  grad_P = (Vc-Vb)/ptri->edge_len[0];
  kapa =  mt->thermal->HeatConduction(0.5*(Tb+Tc));
  grad_T = kapa*(Tc-Tb)/ptri->edge_len[0];
#ifdef _IIA_
  if(ImpactIonization && IIType==EdotJ)
  {
    IIn =  mt->gen->ElecGenRate(0.5*(Tb+Tc),dmax(0,Epna),Eg);
    IIp =  mt->gen->HoleGenRate(0.5*(Tb+Tc),dmax(0,Eppa),Eg);
    G   =  IIn*mun*Jna_norm*Jn_scale + IIp*mup*Jpa_norm*Jp_scale;
  }
  if( ImpactIonization && (IIType==GradQf || IIType==EVector) )
  {
    G = IIn*mun*fabs(phinb-phinc)/ptri->edge_len[0]*nmid(kb,e,(Ecb-Ecc)/e,nb,nc,0.5*(Tb+Tc),Tc-Tb)
        +IIp*mup*fabs(phipb-phipc)/ptri->edge_len[0]*pmid(kb,e,(Evb-Evc)/e,pb,pc,0.5*(Tb+Tc),Tc-Tb);
    //G = IIn*mun*fabs(Vb-Vc)*nmid(kb,e,(Ecb-Ecc)/e,nb,nc,T_mid,Tc-Tb)/e
    //   +IIp*mup*fabs(Vb-Vc)*pmid(kb,e,(Evb-Evc)/e,pb,pc,T_mid,Tc-Tb)/e;
  }
  if( ImpactIonization && IIType==ESide )
  {
    IIn =  mt->gen->ElecGenRate(0.5*(Tb+Tc),fabs(Vb-Vc)/ptri->edge_len[0],Eg);
    IIp =  mt->gen->HoleGenRate(0.5*(Tb+Tc),fabs(Vb-Vc)/ptri->edge_len[0],Eg);
    G = IIn*mun*fabs(Jna)*Jn_scale + IIp*mup*fabs(Jpa)*Jp_scale;
  }
#endif

#ifdef _HEAT_
  H = 0.5*(Vb-Vc)*(mun*Jna*Jn_scale+mup*Jpa*Jp_scale);
#endif

  f[zofs[zone_index]+4*B+0] +=   grad_P*ptri->d[0];
  f[zofs[zone_index]+4*B+1] +=   mun*Jna*Jn_scale*ptri->d[0] + (G-Rb)*0.25*ptri->d[0]*ptri->edge_len[0];
  f[zofs[zone_index]+4*B+2] += - mup*Jpa*Jp_scale*ptri->d[0] + (G-Rb)*0.25*ptri->d[0]*ptri->edge_len[0];
  f[zofs[zone_index]+4*B+3] +=   ( grad_T + H - G*0.25*ptri->edge_len[0]*(Eg+3*kb*0.5*(Tb+Tc)))*ptri->d[0];


  f[zofs[zone_index]+4*C+0] += - grad_P*ptri->d[0];
  f[zofs[zone_index]+4*C+1] += - mun*Jna*Jn_scale*ptri->d[0] + (G-Rc)*0.25*ptri->d[0]*ptri->edge_len[0];
  f[zofs[zone_index]+4*C+2] +=   mup*Jpa*Jp_scale*ptri->d[0] + (G-Rc)*0.25*ptri->d[0]*ptri->edge_len[0];
  f[zofs[zone_index]+4*C+3] +=   (-grad_T + H - G*0.25*ptri->edge_len[0]*(Eg+3*kb*0.5*(Tb+Tc)))*ptri->d[0];


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
    aux[C].mun =  mt->mob->ElecMob(pc,nc,Tc,dmax(0,Epnb),fabs(Etnb),Tc);
    aux[C].mup =  mt->mob->HoleMob(pc,nc,Tc,dmax(0,Eppb),fabs(Etpb),Tc);

    mt->mapping(&pzone->danode[A],&aux[A],0);
    aux[A].mun =  mt->mob->ElecMob(pa,na,Ta,dmax(0,Epnb),fabs(Etnb),Ta);
    aux[A].mup =  mt->mob->HoleMob(pa,na,Ta,dmax(0,Eppb),fabs(Etpb),Ta);
  }

  mun = 0.5*(aux[C].mun+aux[A].mun);
  mup = 0.5*(aux[C].mup+aux[A].mup);

  grad_P = (Va-Vc)/ptri->edge_len[1];
  kapa =  mt->thermal->HeatConduction(0.5*(Tc+Ta));
  grad_T = kapa*(Ta-Tc)/ptri->edge_len[1];
#ifdef _IIB_
  if(ImpactIonization && IIType==EdotJ)
  {
    IIn =  mt->gen->ElecGenRate(0.5*(Tc+Ta),dmax(0,Epnb),Eg);
    IIp =  mt->gen->HoleGenRate(0.5*(Tc+Ta),dmax(0,Eppb),Eg);
    G   =  IIn*mun*Jnb_norm*Jn_scale + IIp*mup*Jpb_norm*Jp_scale;
  }
  if( ImpactIonization && (IIType==GradQf || IIType==EVector) )
  {
    G = IIn*mun*fabs(phinc-phina)/ptri->edge_len[1]*nmid(kb,e,(Ecc-Eca)/e,nc,na,0.5*(Tc+Ta),Ta-Tc)
        +IIp*mup*fabs(phipc-phipa)/ptri->edge_len[1]*pmid(kb,e,(Evc-Eva)/e,pc,pa,0.5*(Tc+Ta),Ta-Tc);
    //G = IIn*mun*fabs(Vc-Va)/ptri->edge_len[1]*nmid(kb,e,(Ecc-Eca)/e,nc,na,T_mid,Ta-Tc)/e
    //   +IIp*mup*fabs(Vc-Va)/ptri->edge_len[1]*pmid(kb,e,(Evc-Eva)/e,pc,pa,T_mid,Ta-Tc)/e;
  }
  if( ImpactIonization && IIType==ESide )
  {
    IIn =  mt->gen->ElecGenRate(0.5*(Tc+Ta),fabs(Vc-Va)/ptri->edge_len[1],Eg);
    IIp =  mt->gen->HoleGenRate(0.5*(Tc+Ta),fabs(Vc-Va)/ptri->edge_len[1],Eg);
    G = IIn*mun*fabs(Jnb)*Jn_scale + IIp*mup*fabs(Jpb)*Jp_scale;
  }
#endif

#ifdef _HEAT_
  H = 0.5*(Vc-Va)*(mun*Jnb*Jn_scale+mup*Jpb*Jp_scale);
#endif
  f[zofs[zone_index]+4*C+0] +=   grad_P*ptri->d[1];
  f[zofs[zone_index]+4*C+1] +=   mun*Jnb*Jn_scale*ptri->d[1] + (G-Rc)*0.25*ptri->d[1]*ptri->edge_len[1];
  f[zofs[zone_index]+4*C+2] += - mup*Jpb*Jp_scale*ptri->d[1] + (G-Rc)*0.25*ptri->d[1]*ptri->edge_len[1];
  f[zofs[zone_index]+4*C+3] +=   ( grad_T + H - G*0.25*ptri->edge_len[1]*(Eg+3*kb*0.5*(Tc+Ta)))*ptri->d[1];
 

  f[zofs[zone_index]+4*A+0] +=  -grad_P*ptri->d[1];
  f[zofs[zone_index]+4*A+1] +=  -mun*Jnb*Jn_scale*ptri->d[1] + (G-Ra)*0.25*ptri->d[1]*ptri->edge_len[1];
  f[zofs[zone_index]+4*A+2] +=   mup*Jpb*Jp_scale*ptri->d[1] + (G-Ra)*0.25*ptri->d[1]*ptri->edge_len[1];
  f[zofs[zone_index]+4*A+3] +=   (-grad_T + H - G*0.25*ptri->edge_len[1]*(Eg+3*kb*0.5*(Tc+Ta)))*ptri->d[1];

  
  if(BandBandTunneling)
  {
    PetscScalar G_BB = mt->band->BB_Tunneling((Ta+Tb+Tc)/3.0,E);
    f[zofs[zone_index]+4*A+1] +=  G_BB*ptri->s[0];
    f[zofs[zone_index]+4*A+2] +=  G_BB*ptri->s[0];
    f[zofs[zone_index]+4*B+1] +=  G_BB*ptri->s[1];
    f[zofs[zone_index]+4*B+2] +=  G_BB*ptri->s[1];
    f[zofs[zone_index]+4*C+1] +=  G_BB*ptri->s[2];
    f[zofs[zone_index]+4*C+2] +=  G_BB*ptri->s[2];
  }
  
  //---------------------------------------------------------------------------
  // Recomb item and heat source for each node
  //---------------------------------------------------------------------------
  
  S = 0.25*ptri->d[2]*ptri->edge_len[2] + 0.25*ptri->d[1]*ptri->edge_len[1];
  PetscScalar Ha  = Ra*(Ega+3*kb*Ta);
  f[zofs[zone_index]+4*A+1] += -Ra*S;
  f[zofs[zone_index]+4*A+2] += -Ra*S;
  f[zofs[zone_index]+4*A+3] +=  Ha*S;
  
  S = 0.25*ptri->d[0]*ptri->edge_len[0] + 0.25*ptri->d[2]*ptri->edge_len[2];
  PetscScalar Hb  = Rb*(Egb+3*kb*Tb);
  f[zofs[zone_index]+4*B+1] += -Rb*S;
  f[zofs[zone_index]+4*B+2] += -Rb*S;
  f[zofs[zone_index]+4*B+3] +=  Hb*S;

  
  S = 0.25*ptri->d[1]*ptri->edge_len[1] + 0.25*ptri->d[0]*ptri->edge_len[0];
  PetscScalar Hc  = Rc*(Egc+3*kb*Tc);
  f[zofs[zone_index]+4*C+1] += -Rc*S;
  f[zofs[zone_index]+4*C+2] += -Rc*S;
  f[zofs[zone_index]+4*C+3] +=  Hc*S;


}


void SMCZone::J2E_Tri_ddm(Tri *ptri,PetscScalar *x,Mat *jtmp, vector<int> & zofs)
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

  PetscScalar e  =  mt->e;
  PetscScalar kb =  mt->kb;
  PetscScalar S;
  //the indepedent variable number
  adtl::AutoDScalar::numdir=12;
  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  Mat4        J1,J2,J3;
  AutoDScalar F[4];
  PetscInt    index1[4],index2[4],index3[4];
  index1[0] = zofs[z]+4*A+0;
  index1[1] = zofs[z]+4*A+1;
  index1[2] = zofs[z]+4*A+2;
  index1[3] = zofs[z]+4*A+3;

  index2[0] = zofs[z]+4*B+0;
  index2[1] = zofs[z]+4*B+1;
  index2[2] = zofs[z]+4*B+2;
  index2[3] = zofs[z]+4*B+3;

  index3[0] = zofs[z]+4*C+0;
  index3[1] = zofs[z]+4*C+1;
  index3[2] = zofs[z]+4*C+2;
  index3[3] = zofs[z]+4*C+3;
  //---------------------------------------------------------------------------

  AutoDScalar Va = x[zofs[zone_index]+4*A+0];     //potential of node A
  AutoDScalar na = x[zofs[zone_index]+4*A+1];     //electron density of node A
  AutoDScalar pa = x[zofs[zone_index]+4*A+2];     //hole density of node A
  AutoDScalar Ta = x[zofs[zone_index]+4*A+3];     //Temperature of node A
  Va.setADValue(0,1.0);
  na.setADValue(1,1.0);
  pa.setADValue(2,1.0);
  Ta.setADValue(3,1.0);
  mt->mapping(&pzone->danode[A],&aux[A],0);
  AutoDScalar nia = mt->band->nie(Ta);
  AutoDScalar Ra = mt->band->Recomb(pa,na,Ta);
#ifdef _APPROXIMATE_
  PetscScalar Nca = mt->band->Nc(fs[A].T);
  PetscScalar Nva = mt->band->Nv(fs[A].T);
  PetscScalar Ega = mt->band->Eg(fs[A].T);
  AutoDScalar Eca = -(e*Va + aux[A].affinity + mt->band->EgNarrowToEc(fs[A].T) + kb*fs[A].T*(log(Nca)-1.5*log(fs[A].T)));
  AutoDScalar Eva = -(e*Va + aux[A].affinity - mt->band->EgNarrowToEv(fs[A].T) - kb*fs[A].T*(log(Nva)-1.5*log(fs[A].T)) + Ega);
#else  
  AutoDScalar Nca = mt->band->Nc(Ta);
  AutoDScalar Nva = mt->band->Nv(Ta);
  AutoDScalar Ega = mt->band->Eg(Ta);
  AutoDScalar Eca = -(e*Va + aux[A].affinity + mt->band->EgNarrowToEc(Ta) + kb*Ta*(log(Nca)-1.5*log(Ta)));
  AutoDScalar Eva = -(e*Va + aux[A].affinity - mt->band->EgNarrowToEv(Ta) - kb*Ta*(log(Nva)-1.5*log(Ta)) + Ega);
#endif  
#ifdef _FERMI_
  if(Fermi)
  {
    Eca = Eca - kb*Ta*log(gamma_f(fabs(na)/Nca));
    Eva = Eva + kb*Ta*log(gamma_f(fabs(pa)/Nva));
  }
#endif
  
  AutoDScalar Vb = x[zofs[zone_index]+4*B+0];     //potential of node B
  AutoDScalar nb = x[zofs[zone_index]+4*B+1];     //electron density of node B
  AutoDScalar pb = x[zofs[zone_index]+4*B+2];     //hole density of node B
  AutoDScalar Tb = x[zofs[zone_index]+4*B+3];     //Temperature of node B
  Vb.setADValue(4,1.0);
  nb.setADValue(5,1.0);
  pb.setADValue(6,1.0);
  Tb.setADValue(7,1.0);
  mt->mapping(&pzone->danode[B],&aux[B],0);
  AutoDScalar nib = mt->band->nie(Tb);
  AutoDScalar Rb = mt->band->Recomb(pb,nb,Tb);
#ifdef _APPROXIMATE_
  PetscScalar Ncb = mt->band->Nc(fs[B].T);
  PetscScalar Nvb = mt->band->Nv(fs[B].T);
  PetscScalar Egb = mt->band->Eg(fs[B].T);
  AutoDScalar Ecb =  -(e*Vb + aux[B].affinity + mt->band->EgNarrowToEc(fs[B].T)+ kb*fs[B].T*(log(Ncb)-1.5*log(fs[B].T)));
  AutoDScalar Evb =  -(e*Vb + aux[B].affinity - mt->band->EgNarrowToEv(fs[B].T)- kb*fs[B].T*(log(Nvb)-1.5*log(fs[B].T))+ Egb);
#else  
  AutoDScalar Ncb = mt->band->Nc(Tb);
  AutoDScalar Nvb = mt->band->Nv(Tb);
  AutoDScalar Egb = mt->band->Eg(Tb);
  AutoDScalar Ecb =  -(e*Vb + aux[B].affinity + mt->band->EgNarrowToEc(Tb)+ kb*Tb*(log(Ncb)-1.5*log(Tb)));
  AutoDScalar Evb =  -(e*Vb + aux[B].affinity - mt->band->EgNarrowToEv(Tb)- kb*Tb*(log(Nvb)-1.5*log(Tb))+ Egb);
#endif
#ifdef _FERMI_
  if(Fermi)
  {
    Ecb = Ecb - kb*Tb*log(gamma_f(fabs(nb)/Ncb));
    Evb = Evb + kb*Tb*log(gamma_f(fabs(pb)/Nvb));
  }
#endif
 
  AutoDScalar Vc = x[zofs[zone_index]+4*C+0];     //potential of node C
  AutoDScalar nc = x[zofs[zone_index]+4*C+1];     //electron density of node C
  AutoDScalar pc = x[zofs[zone_index]+4*C+2];     //hole density of node C
  AutoDScalar Tc = x[zofs[zone_index]+4*C+3];     //Temperature of node C
  Vc.setADValue(8,1.0);
  nc.setADValue(9,1.0);
  pc.setADValue(10,1.0);
  Tc.setADValue(11,1.0);
  mt->mapping(&pzone->danode[C],&aux[C],0);
  AutoDScalar nic = mt->band->nie(fs[C].T);
  AutoDScalar Rc = mt->band->Recomb(pc,nc,Tc);
#ifdef _APPROXIMATE_  
  PetscScalar Ncc = mt->band->Nc(fs[C].T);
  PetscScalar Nvc = mt->band->Nv(fs[C].T);
  PetscScalar Egc = mt->band->Eg(fs[C].T);
  AutoDScalar Ecc =  -(e*Vc + aux[C].affinity + mt->band->EgNarrowToEc(fs[C].T)+ kb*fs[C].T*(log(Ncc)-1.5*log(fs[C].T)));
  AutoDScalar Evc =  -(e*Vc + aux[C].affinity - mt->band->EgNarrowToEv(fs[C].T)- kb*fs[C].T*(log(Nvc)-1.5*log(fs[C].T))+ Egc);
#else
  AutoDScalar Ncc = mt->band->Nc(Tc);
  AutoDScalar Nvc = mt->band->Nv(Tc);
  AutoDScalar Egc = mt->band->Eg(Tc);
  AutoDScalar Ecc =  -(e*Vc + aux[C].affinity + mt->band->EgNarrowToEc(Tc)+ kb*Tc*(log(Ncc)-1.5*log(Tc)));
  AutoDScalar Evc =  -(e*Vc + aux[C].affinity - mt->band->EgNarrowToEv(Tc)- kb*Tc*(log(Nvc)-1.5*log(Tc))+ Egc);
#endif  
#ifdef _FERMI_
  if(Fermi)
  {
    Ecc = Ecc - kb*Tc*log(gamma_f(fabs(nc)/Ncc));
    Evc = Evc + kb*Tc*log(gamma_f(fabs(pc)/Nvc));
  }
#endif
  
  AutoDScalar phina;
  AutoDScalar phipa;
  AutoDScalar phinb;
  AutoDScalar phipb;
  AutoDScalar phinc;
  AutoDScalar phipc;

  AutoDScalar Ex = -((yb-yc)*Va + (yc-ya)*Vb +(ya-yb)*Vc)/(2*tri_area);
  AutoDScalar Ey = -((xc-xb)*Va + (xa-xc)*Vb +(xb-xa)*Vc)/(2*tri_area);
  AutoDScalar Eg =  mt->band->Eg((Ta+Tb+Tc)/3.0);
  AutoDScalar G = 0;
  AutoDScalar IIn=0,IIp=0;
  AutoDScalar mun,mup;
  AutoDScalar Epnc=0,Etnc=0;
  AutoDScalar Eppc=0,Etpc=0;
  AutoDScalar Epna=0,Etna=0;
  AutoDScalar Eppa=0,Etpa=0;
  AutoDScalar Epnb=0,Etnb=0;
  AutoDScalar Eppb=0,Etpb=0;
  AutoDScalar Jna_norm,Jpa_norm;
  AutoDScalar Jnb_norm,Jpb_norm;
  AutoDScalar Jnc_norm,Jpc_norm;
  AutoDScalar grad_P,kapa,grad_T,T_mid,H;

#ifdef _II_
  if(ImpactIonization && IIType==GradQf)
  {
    phina = Va - log(fabs(na)/nia)*kb*Ta/e;
    phipa = Va + log(fabs(pa)/nia)*kb*Ta/e;
    phinb = Vb - log(fabs(nb)/nib)*kb*Tb/e;
    phipb = Vb + log(fabs(pb)/nib)*kb*Tb/e;
    phinc = Vc - log(fabs(nc)/nic)*kb*Tc/e;
    phipc = Vc + log(fabs(pc)/nic)*kb*Tc/e;
    AutoDScalar Fnx = ((yb-yc)*phina + (yc-ya)*phinb +(ya-yb)*phinc)/(2*tri_area);
    AutoDScalar Fny = ((xc-xb)*phina + (xa-xc)*phinb +(xb-xa)*phinc)/(2*tri_area);
    AutoDScalar Fpx = ((yb-yc)*phipa + (yc-ya)*phipb +(ya-yb)*phipc)/(2*tri_area);
    AutoDScalar Fpy = ((xc-xb)*phipa + (xa-xc)*phipb +(xb-xa)*phipc)/(2*tri_area);
    AutoDScalar Fn  =  sqrt(Fnx*Fnx+Fny*Fny+1e-20);
    AutoDScalar Fp  =  sqrt(Fpx*Fpx+Fpy*Fpy+1e-20);
    IIn =  mt->gen->ElecGenRate((Ta+Tb+Tc)/3.0,Fn,Eg);
    IIp =  mt->gen->HoleGenRate((Ta+Tb+Tc)/3.0,Fp,Eg);
  }

  if(ImpactIonization && IIType==EVector)
  {
    phina = Va - log(fabs(na)/nia)*kb*Ta/e;
    phipa = Va + log(fabs(pa)/nia)*kb*Ta/e;
    phinb = Vb - log(fabs(nb)/nib)*kb*Tb/e;
    phipb = Vb + log(fabs(pb)/nib)*kb*Tb/e;
    phinc = Vc - log(fabs(nc)/nic)*kb*Tc/e;
    phipc = Vc + log(fabs(pc)/nic)*kb*Tc/e;
    AutoDScalar F  =  sqrt(Ex*Ex+Ey*Ey+1e-20);
    IIn =  mt->gen->ElecGenRate((Ta+Tb+Tc)/3.0,F,Eg);
    IIp =  mt->gen->HoleGenRate((Ta+Tb+Tc)/3.0,F,Eg);
  }
#endif
  AutoDScalar Jnc =  In(kb,e,(Eca-Ecb)/e,na,nb,0.5*(Ta+Tb),Tb-Ta,ptri->edge_len[2]);
  AutoDScalar Jpc =  Ip(kb,e,(Eva-Evb)/e,pa,pb,0.5*(Ta+Tb),Tb-Ta,ptri->edge_len[2]);
  AutoDScalar Jna =  In(kb,e,(Ecb-Ecc)/e,nb,nc,0.5*(Tb+Tc),Tc-Tb,ptri->edge_len[0]);
  AutoDScalar Jpa =  Ip(kb,e,(Evb-Evc)/e,pb,pc,0.5*(Tb+Tc),Tc-Tb,ptri->edge_len[0]);
  AutoDScalar Jnb =  In(kb,e,(Ecc-Eca)/e,nc,na,0.5*(Tc+Ta),Ta-Tc,ptri->edge_len[1]);
  AutoDScalar Jpb =  Ip(kb,e,(Evc-Eva)/e,pc,pa,0.5*(Tc+Ta),Ta-Tc,ptri->edge_len[1]);
  PetscScalar Jn_scale = dmax(dmax(fabs(Jna.getValue()),fabs(Jnb.getValue())),
                              dmax(fabs(Jnc.getValue()),nia.getValue()*nia.getValue()));
  PetscScalar Jp_scale = dmax(dmax(fabs(Jpa.getValue()),fabs(Jpb.getValue())),
                              dmax(fabs(Jpc.getValue()),nia.getValue()*nia.getValue()));

  Jnc /=  Jn_scale;
  Jpc /=  Jp_scale;
  Jna /=  Jn_scale;
  Jpa /=  Jp_scale;
  Jnb /=  Jn_scale;
  Jpb /=  Jp_scale;

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
    AutoDScalar munCA =  mt->mob->ElecMob(pa,na,Ta,adtl::fmax(0,Epnc),fabs(Etnc),Ta);
    AutoDScalar mupCA =  mt->mob->HoleMob(pa,na,Ta,adtl::fmax(0,Eppc),fabs(Etpc),Ta);

    mt->mapping(&pzone->danode[B],&aux[B],0);
    AutoDScalar munCB =  mt->mob->ElecMob(pb,nb,Tb,adtl::fmax(0,Epnc),fabs(Etnc),Tb);
    AutoDScalar mupCB =  mt->mob->HoleMob(pb,nb,Tb,adtl::fmax(0,Eppc),fabs(Etpc),Tb);

    mun = 0.5*(munCA+munCB);
    mup = 0.5*(mupCA+mupCB);
  }
  else
  {
    mun = 0.5*(aux[A].mun+aux[B].mun);
    mup = 0.5*(aux[A].mup+aux[B].mup);
  }
  T_mid = 0.5*(Ta+Tb);
  grad_P = (Vb-Va)/ptri->edge_len[2];
  kapa =  mt->thermal->HeatConduction(T_mid);
  grad_T = kapa*(Tb-Ta)/ptri->edge_len[2];

#ifdef _IIC_
  if(ImpactIonization)
  {
    if(IIType==EdotJ)
    {
      IIn =  mt->gen->ElecGenRate(T_mid,adtl::fmax(0,Epnc),Eg);
      IIp =  mt->gen->HoleGenRate(T_mid,adtl::fmax(0,Eppc),Eg);
      G   =  IIn*mun*Jnc_norm*Jn_scale + IIp*mup*Jpc_norm*Jp_scale;
    }
    if((IIType==GradQf || IIType==EVector) )
    {
      G = IIn*mun*fabs(phina-phinb)/ptri->edge_len[2]*nmid(kb,e,(Eca-Ecb)/e,na,nb,T_mid,Tb-Ta)
         +IIp*mup*fabs(phipa-phipb)/ptri->edge_len[2]*pmid(kb,e,(Eva-Evb)/e,pa,pb,T_mid,Tb-Ta);
      //G = IIn*mun*fabs(Va-Vb)/ptri->edge_len[2]*nmid(kb,e,(Eca-Ecb)/e,na,nb,T_mid,Tb-Ta)/e
      //   +IIp*mup*fabs(Va-Vb)/ptri->edge_len[2]*pmid(kb,e,(Eva-Evb)/e,pa,pb,T_mid,Tb-Ta)/e;
    }
    if(IIType==ESide )
    {
      IIn =  mt->gen->ElecGenRate(T_mid,fabs(Va-Vb)/ptri->edge_len[2],Eg);
      IIp =  mt->gen->HoleGenRate(T_mid,fabs(Va-Vb)/ptri->edge_len[2],Eg);
      G = IIn*mun*fabs(Jnc)*Jn_scale + IIp*mup*fabs(Jpc)*Jp_scale;
    }
    AutoDScalar W = G*(Eg + 3*kb*T_mid);
    J1.m[0] =  0;
    J1.m[1] =  0;
    J1.m[2] =  0;
    J1.m[3] =  0;
    J1.m[4] =  G.getADValue(0)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[5] =  G.getADValue(1)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[6] =  G.getADValue(2)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[7] =  G.getADValue(3)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[8] =  G.getADValue(0)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[9] =  G.getADValue(1)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[10]=  G.getADValue(2)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[11]=  G.getADValue(3)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[12]=  -W.getADValue(0)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[13]=  -W.getADValue(1)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[14]=  -W.getADValue(2)*0.25*ptri->d[2]*ptri->edge_len[2];
    J1.m[15]=  -W.getADValue(3)*0.25*ptri->d[2]*ptri->edge_len[2];

    J2.m[0] =  0;
    J2.m[1] =  0;
    J2.m[2] =  0;
    J2.m[3] =  0;
    J2.m[4] =  G.getADValue(4)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[5] =  G.getADValue(5)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[6] =  G.getADValue(6)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[7] =  G.getADValue(7)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[8] =  G.getADValue(4)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[9] =  G.getADValue(5)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[10]=  G.getADValue(6)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[11]=  G.getADValue(7)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[12]=  -W.getADValue(4)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[13]=  -W.getADValue(5)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[14]=  -W.getADValue(6)*0.25*ptri->d[2]*ptri->edge_len[2];
    J2.m[15]=  -W.getADValue(7)*0.25*ptri->d[2]*ptri->edge_len[2];

    J3.m[0] =  0;
    J3.m[1] =  0;
    J3.m[2] =  0;
    J3.m[3] =  0;
    J3.m[4] =  G.getADValue(8)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[5] =  G.getADValue(9)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[6] =  G.getADValue(10)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[7] =  G.getADValue(11)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[8] =  G.getADValue(8)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[9] =  G.getADValue(9)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[10]=  G.getADValue(10)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[11]=  G.getADValue(11)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[12]=  -W.getADValue(8)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[13]=  -W.getADValue(9)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[14]=  -W.getADValue(10)*0.25*ptri->d[2]*ptri->edge_len[2];
    J3.m[15]=  -W.getADValue(11)*0.25*ptri->d[2]*ptri->edge_len[2];

    MatSetValues(*jtmp,4,index1,4,index1,J1.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index1,4,index2,J2.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index1,4,index3,J3.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index2,4,index1,J1.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index2,4,index2,J2.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index2,4,index3,J3.m,ADD_VALUES);
    
  }
#endif


#ifdef _HEAT_
  H = 0.5*(Va-Vb)*(mun*Jnc*Jn_scale+mup*Jpc*Jp_scale)*ptri->d[2];
  Set_Mat4_zero(J1);
  Set_Mat4_zero(J2);
  Set_Mat4_zero(J3);
  J1.m[12] =  H.getADValue(0);          //dH/dP(A)
  J1.m[13] =  H.getADValue(1);          //dH/dn(A)
  J1.m[14] =  H.getADValue(2);          //dH/dp(A)
  J1.m[15] =  H.getADValue(3);          //dH/dT(A)

  J2.m[12] =  H.getADValue(4);          //dH/dP(B)
  J2.m[13] =  H.getADValue(5);          //dH/dn(B)
  J2.m[14] =  H.getADValue(6);          //dH/dp(B)
  J2.m[15] =  H.getADValue(7);          //dH/dT(B)

  J3.m[12] =  H.getADValue(8);          //dH/dP(C)
  J3.m[13] =  H.getADValue(9);          //dH/dn(C)
  J3.m[14] =  H.getADValue(10);         //dH/dp(C)
  J3.m[15] =  H.getADValue(11);         //dH/dT(C)
  MatSetValues(*jtmp,4,index1,4,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index1,4,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index1,4,index3,J3.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index2,4,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index2,4,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index2,4,index3,J3.m,ADD_VALUES);
#endif
  //---------------------------------------------------------------------------
  F[0] =   grad_P*ptri->d[2];
  F[1] =   mun*Jnc*Jn_scale*ptri->d[2];
  F[2] = - mup*Jpc*Jp_scale*ptri->d[2];
  F[3] =   grad_T*ptri->d[2];
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
    {
      J1.m[4*i+j] = F[i].getADValue(0+j);
      J2.m[4*i+j] = F[i].getADValue(4+j);
      J3.m[4*i+j] = F[i].getADValue(8+j);
    }
  MatSetValues(*jtmp,4,index1,4,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index1,4,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index1,4,index3,J3.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index2,4,index1,(-J1).m,ADD_VALUES);
  MatSetValues(*jtmp,4,index2,4,index2,(-J2).m,ADD_VALUES);
  MatSetValues(*jtmp,4,index2,4,index3,(-J3).m,ADD_VALUES);
    
  
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
    AutoDScalar munAB =  mt->mob->ElecMob(pb,nb,Tb,adtl::fmax(0,Epna),fabs(Etna),Tb);
    AutoDScalar mupAB =  mt->mob->HoleMob(pb,nb,Tb,adtl::fmax(0,Eppa),fabs(Etpa),Tb);

    mt->mapping(&pzone->danode[C],&aux[C],0);
    AutoDScalar munAC =  mt->mob->ElecMob(pc,nc,Tc,adtl::fmax(0,Epna),fabs(Etna),Tc);
    AutoDScalar mupAC =  mt->mob->HoleMob(pc,nc,Tc,adtl::fmax(0,Eppa),fabs(Etpa),Tc);
  
    mun = 0.5*(munAB+munAC);
    mup = 0.5*(mupAB+mupAC);
  }
  else 
  {  
    mun = 0.5*(aux[B].mun+aux[C].mun);
    mup = 0.5*(aux[B].mup+aux[C].mup);
  }
  T_mid = 0.5*(Tb+Tc);
  grad_P = (Vc-Vb)/ptri->edge_len[0];
  kapa =  mt->thermal->HeatConduction(T_mid);
  grad_T = kapa*(Tc-Tb)/ptri->edge_len[0];
#ifdef _IIA_
  if(ImpactIonization)
  {
    if(IIType==EdotJ)
    {
      IIn =  mt->gen->ElecGenRate(T_mid,adtl::fmax(0,Epna),Eg);
      IIp =  mt->gen->HoleGenRate(T_mid,adtl::fmax(0,Eppa),Eg);
      G   =  IIn*mun*Jna_norm*Jn_scale + IIp*mup*Jpa_norm*Jp_scale;
    }
    if(IIType==GradQf || IIType==EVector )
    {
      G = IIn*mun*fabs(phinb-phinc)/ptri->edge_len[0]*nmid(kb,e,(Ecb-Ecc)/e,nb,nc,T_mid,Tc-Tb)
         +IIp*mup*fabs(phipb-phipc)/ptri->edge_len[0]*pmid(kb,e,(Evb-Evc)/e,pb,pc,T_mid,Tc-Tb);
      //G = IIn*mun*fabs(Vb-Vc)*nmid(kb,e,(Ecb-Ecc)/e,nb,nc,T_mid,Tc-Tb)/e
      //   +IIp*mup*fabs(Vb-Vc)*pmid(kb,e,(Evb-Evc)/e,pb,pc,T_mid,Tc-Tb)/e;
    }
    if(IIType==ESide )
    {
      IIn =  mt->gen->ElecGenRate(T_mid,fabs(Vb-Vc)/ptri->edge_len[0],Eg);
      IIp =  mt->gen->HoleGenRate(T_mid,fabs(Vb-Vc)/ptri->edge_len[0],Eg);
      G = IIn*mun*fabs(Jna)*Jn_scale + IIp*mup*fabs(Jpa)*Jp_scale;
    }
    AutoDScalar W = G*(Eg + 3*kb*T_mid);
    
    J1.m[0] =  0;
    J1.m[1] =  0;
    J1.m[0] =  0;
    J1.m[3] =  0;
    J1.m[4] =  G.getADValue(0)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[5] =  G.getADValue(1)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[6] =  G.getADValue(2)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[7] =  G.getADValue(3)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[8] =  G.getADValue(0)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[9] =  G.getADValue(1)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[10]=  G.getADValue(2)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[11]=  G.getADValue(3)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[12]=  -W.getADValue(0)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[13]=  -W.getADValue(1)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[14]=  -W.getADValue(2)*0.25*ptri->d[0]*ptri->edge_len[0];
    J1.m[15]=  -W.getADValue(3)*0.25*ptri->d[0]*ptri->edge_len[0];

    J2.m[0] =  0;
    J2.m[1] =  0;
    J2.m[0] =  0;
    J2.m[3] =  0;
    J2.m[4] =  G.getADValue(4)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[5] =  G.getADValue(5)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[6] =  G.getADValue(6)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[7] =  G.getADValue(7)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[8] =  G.getADValue(4)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[9] =  G.getADValue(5)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[10]=  G.getADValue(6)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[11]=  G.getADValue(7)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[12]=  -W.getADValue(4)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[13]=  -W.getADValue(5)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[14]=  -W.getADValue(6)*0.25*ptri->d[0]*ptri->edge_len[0];
    J2.m[15]=  -W.getADValue(7)*0.25*ptri->d[0]*ptri->edge_len[0];

    J3.m[0] =  0;
    J3.m[1] =  0;
    J3.m[0] =  0;
    J3.m[3] =  0;
    J3.m[4] =  G.getADValue(8)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[5] =  G.getADValue(9)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[6] =  G.getADValue(10)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[7] =  G.getADValue(11)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[8] =  G.getADValue(8)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[9] =  G.getADValue(9)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[10]=  G.getADValue(10)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[11]=  G.getADValue(11)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[12]=  -W.getADValue(8)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[13]=  -W.getADValue(9)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[14]=  -W.getADValue(10)*0.25*ptri->d[0]*ptri->edge_len[0];
    J3.m[15]=  -W.getADValue(11)*0.25*ptri->d[0]*ptri->edge_len[0];

    MatSetValues(*jtmp,4,index2,4,index2,J2.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index2,4,index3,J3.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index2,4,index1,J1.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index3,4,index2,J2.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index3,4,index3,J3.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index3,4,index1,J1.m,ADD_VALUES);
  }
#endif    

#ifdef _HEAT_
  H = 0.5*(Vb-Vc)*(mun*Jna*Jn_scale+mup*Jpa*Jp_scale)*ptri->d[0];
  Set_Mat4_zero(J1);
  Set_Mat4_zero(J2);
  Set_Mat4_zero(J3);
  J1.m[12] =  H.getADValue(0);          //dH/dP(A)
  J1.m[13] =  H.getADValue(1);          //dH/dn(A)
  J1.m[14] =  H.getADValue(2);          //dH/dp(A)
  J1.m[15] =  H.getADValue(3);          //dH/dT(A)

  J2.m[12] =  H.getADValue(4);          //dH/dP(B)
  J2.m[13] =  H.getADValue(5);          //dH/dn(B)
  J2.m[14] =  H.getADValue(6);          //dH/dp(B)
  J2.m[15] =  H.getADValue(7);          //dH/dT(B)

  J3.m[12] =  H.getADValue(8);          //dH/dP(C)
  J3.m[13] =  H.getADValue(9);          //dH/dn(C)
  J3.m[14] =  H.getADValue(10);         //dH/dp(C)
  J3.m[15] =  H.getADValue(11);         //dH/dT(C)
  MatSetValues(*jtmp,4,index2,4,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index2,4,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index2,4,index3,J3.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index3,4,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index3,4,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index3,4,index3,J3.m,ADD_VALUES);
#endif
  //---------------------------------------------------------------------------
  F[0] =   grad_P*ptri->d[0];
  F[1] =   mun*Jna*Jn_scale*ptri->d[0];
  F[2] = - mup*Jpa*Jp_scale*ptri->d[0];
  F[3] =   grad_T*ptri->d[0];
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
    {
      J1.m[4*i+j] = F[i].getADValue(0+j);
      J2.m[4*i+j] = F[i].getADValue(4+j);
      J3.m[4*i+j] = F[i].getADValue(8+j);
    }
  MatSetValues(*jtmp,4,index2,4,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index2,4,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index2,4,index3,J3.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index3,4,index1,(-J1).m,ADD_VALUES);
  MatSetValues(*jtmp,4,index3,4,index2,(-J2).m,ADD_VALUES);
  MatSetValues(*jtmp,4,index3,4,index3,(-J3).m,ADD_VALUES);
  
  
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
    AutoDScalar munBC =  mt->mob->ElecMob(pc,nc,Tc,adtl::fmax(0,Epnb),fabs(Etnb),Tc);
    AutoDScalar mupBC =  mt->mob->HoleMob(pc,nc,Tc,adtl::fmax(0,Eppb),fabs(Etpb),Tc);

    mt->mapping(&pzone->danode[A],&aux[A],0);
    AutoDScalar munBA =  mt->mob->ElecMob(pa,na,Ta,adtl::fmax(0,Epnb),fabs(Etnb),Ta);
    AutoDScalar mupBA =  mt->mob->HoleMob(pa,na,Ta,adtl::fmax(0,Eppb),fabs(Etpb),Ta);
  
    mun = 0.5*(munBC+munBA);
    mup = 0.5*(mupBC+mupBA);
  }
  else
  {
    mun = 0.5*(aux[C].mun+aux[A].mun);
    mup = 0.5*(aux[C].mup+aux[A].mup);
  }
  T_mid = 0.5*(Tc+Ta);
  grad_P = (Va-Vc)/ptri->edge_len[1];
  kapa =  mt->thermal->HeatConduction(T_mid);
  grad_T = kapa*(Ta-Tc)/ptri->edge_len[1];
#ifdef _IIB_
  if(ImpactIonization)
  {
    if(IIType==EdotJ)
    {
      IIn =  mt->gen->ElecGenRate(T_mid,adtl::fmax(0,Epnb),Eg);
      IIp =  mt->gen->HoleGenRate(T_mid,adtl::fmax(0,Eppb),Eg);
      G   =  IIn*mun*Jnb_norm*Jn_scale + IIp*mup*Jpb_norm*Jp_scale;
    }
    if((IIType==GradQf || IIType==EVector) )
    {
      G = IIn*mun*fabs(phinc-phina)/ptri->edge_len[1]*nmid(kb,e,(Ecc-Eca)/e,nc,na,T_mid,Ta-Tc)
         +IIp*mup*fabs(phipc-phipa)/ptri->edge_len[1]*pmid(kb,e,(Evc-Eva)/e,pc,pa,T_mid,Ta-Tc);
    //G = IIn*mun*fabs(Vc-Va)/ptri->edge_len[1]*nmid(kb,e,(Ecc-Eca)/e,nc,na,T_mid,Ta-Tc)/e
    //   +IIp*mup*fabs(Vc-Va)/ptri->edge_len[1]*pmid(kb,e,(Evc-Eva)/e,pc,pa,T_mid,Ta-Tc)/e;
    }
    if(IIType==ESide )
    {
      IIn =  mt->gen->ElecGenRate(T_mid,fabs(Vc-Va)/ptri->edge_len[1],Eg);
      IIp =  mt->gen->HoleGenRate(T_mid,fabs(Vc-Va)/ptri->edge_len[1],Eg);
      G = IIn*mun*fabs(Jnb)*Jn_scale + IIp*mup*fabs(Jpb)*Jp_scale;
    }
  
    AutoDScalar W = G*(Eg + 3*kb*T_mid);
    J1.m[0] =  0;
    J1.m[1] =  0;
    J1.m[1] =  0;
    J1.m[3] =  0;
    J1.m[4] =  G.getADValue(0)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[5] =  G.getADValue(1)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[6] =  G.getADValue(2)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[7] =  G.getADValue(3)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[8] =  G.getADValue(0)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[9] =  G.getADValue(1)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[10]=  G.getADValue(2)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[11]=  G.getADValue(3)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[12]=  -W.getADValue(0)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[13]=  -W.getADValue(1)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[14]=  -W.getADValue(2)*0.25*ptri->d[1]*ptri->edge_len[1];
    J1.m[15]=  -W.getADValue(3)*0.25*ptri->d[1]*ptri->edge_len[1];

    J2.m[0] =  0;
    J2.m[1] =  0;
    J2.m[1] =  0;
    J2.m[3] =  0;
    J2.m[4] =  G.getADValue(4)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[5] =  G.getADValue(5)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[6] =  G.getADValue(6)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[7] =  G.getADValue(7)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[8] =  G.getADValue(4)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[9] =  G.getADValue(5)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[10]=  G.getADValue(6)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[11]=  G.getADValue(7)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[12]=  -W.getADValue(4)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[13]=  -W.getADValue(5)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[14]=  -W.getADValue(6)*0.25*ptri->d[1]*ptri->edge_len[1];
    J2.m[15]=  -W.getADValue(7)*0.25*ptri->d[1]*ptri->edge_len[1];

    J3.m[0] =  0;
    J3.m[1] =  0;
    J3.m[1] =  0;
    J3.m[3] =  0;
    J3.m[4] =  G.getADValue(8)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[5] =  G.getADValue(9)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[6] =  G.getADValue(10)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[7] =  G.getADValue(11)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[8] =  G.getADValue(8)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[9] =  G.getADValue(9)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[10]=  G.getADValue(10)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[11]=  G.getADValue(11)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[12]=  -W.getADValue(8)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[13]=  -W.getADValue(9)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[14]=  -W.getADValue(10)*0.25*ptri->d[1]*ptri->edge_len[1];
    J3.m[15]=  -W.getADValue(11)*0.25*ptri->d[1]*ptri->edge_len[1];

    MatSetValues(*jtmp,4,index3,4,index3,J3.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index3,4,index1,J1.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index3,4,index2,J2.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index1,4,index3,J3.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index1,4,index1,J1.m,ADD_VALUES);
    MatSetValues(*jtmp,4,index1,4,index2,J2.m,ADD_VALUES);
  }
#endif  
#ifdef _HEAT_
  H = 0.5*(Vc-Va)*(mun*Jnb*Jn_scale+mup*Jpb*Jp_scale)*ptri->d[1];
  Set_Mat4_zero(J1);
  Set_Mat4_zero(J2);
  Set_Mat4_zero(J3);
  J1.m[12] =  H.getADValue(0);          //dH/dP(A)
  J1.m[13] =  H.getADValue(1);          //dH/dn(A)
  J1.m[14] =  H.getADValue(2);          //dH/dp(A)
  J1.m[15] =  H.getADValue(3);          //dH/dT(A)

  J2.m[12] =  H.getADValue(4);          //dH/dP(B)
  J2.m[13] =  H.getADValue(5);          //dH/dn(B)
  J2.m[14] =  H.getADValue(6);          //dH/dp(B)
  J2.m[15] =  H.getADValue(7);          //dH/dT(B)

  J3.m[12] =  H.getADValue(8);          //dH/dP(C)
  J3.m[13] =  H.getADValue(9);          //dH/dn(C)
  J3.m[14] =  H.getADValue(10);         //dH/dp(C)
  J3.m[15] =  H.getADValue(11);         //dH/dT(C)
  MatSetValues(*jtmp,4,index3,4,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index3,4,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index3,4,index3,J3.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index1,4,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index1,4,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index1,4,index3,J3.m,ADD_VALUES);
#endif
  F[0] =   grad_P*ptri->d[1];
  F[1] =   mun*Jnb*Jn_scale*ptri->d[1];
  F[2] = - mup*Jpb*Jp_scale*ptri->d[1];
  F[3] =   grad_T*ptri->d[1];
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
    {
      J1.m[4*i+j] = F[i].getADValue(0+j);
      J2.m[4*i+j] = F[i].getADValue(4+j);
      J3.m[4*i+j] = F[i].getADValue(8+j);
    }
  MatSetValues(*jtmp,4,index3,4,index1,J1.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index3,4,index2,J2.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index3,4,index3,J3.m,ADD_VALUES);
  MatSetValues(*jtmp,4,index1,4,index1,(-J1).m,ADD_VALUES);
  MatSetValues(*jtmp,4,index1,4,index2,(-J2).m,ADD_VALUES);
  MatSetValues(*jtmp,4,index1,4,index3,(-J3).m,ADD_VALUES);    
  
  
  //---------------------------------------------------------------------------
  Set_Mat4_zero(J1);
  AutoDScalar Qra = Ra*(Ega+3*kb*Ta);
  S = 0.25*ptri->d[2]*ptri->edge_len[2] + 0.25*ptri->d[1]*ptri->edge_len[1];
  J1.m[5]  =  -Ra.getADValue(1)*S;              
  J1.m[6]  =  -Ra.getADValue(2)*S;              
  J1.m[7]  =  -Ra.getADValue(3)*S;
  J1.m[9]  =  -Ra.getADValue(1)*S;                  
  J1.m[10] =  -Ra.getADValue(2)*S;                  
  J1.m[11] =  -Ra.getADValue(3)*S;
  J1.m[13] =  Qra.getADValue(1)*S;                  
  J1.m[14] =  Qra.getADValue(2)*S;
  J1.m[15] =  Qra.getADValue(3)*S;
  MatSetValues(*jtmp,4,index1,4,index1,J1.m,ADD_VALUES);
  
  //---------------------------------------------------------------------------
  Set_Mat4_zero(J2);
  AutoDScalar Qrb = Rb*(Egb+3*kb*Tb);
  S = 0.25*ptri->d[0]*ptri->edge_len[0] + 0.25*ptri->d[2]*ptri->edge_len[2];
  J2.m[5]  =  -Rb.getADValue(5)*S;              
  J2.m[6]  =  -Rb.getADValue(6)*S;              
  J2.m[7]  =  -Rb.getADValue(7)*S;
  J2.m[9]  =  -Rb.getADValue(5)*S;                  
  J2.m[10] =  -Rb.getADValue(6)*S;                  
  J2.m[11] =  -Rb.getADValue(7)*S; 
  J2.m[13] =  Qrb.getADValue(5)*S;                  
  J2.m[14] =  Qrb.getADValue(6)*S;
  J2.m[15] =  Qrb.getADValue(7)*S;
  MatSetValues(*jtmp,4,index2,4,index2,J2.m,ADD_VALUES);
  
  //---------------------------------------------------------------------------
  Set_Mat4_zero(J3);
  AutoDScalar Qrc = Rc*(Egc+3*kb*Tc);
  S = 0.25*ptri->d[1]*ptri->edge_len[1] + 0.25*ptri->d[0]*ptri->edge_len[0];
  J3.m[5]  =  -Rc.getADValue( 9)*S;              
  J3.m[6]  =  -Rc.getADValue(10)*S;              
  J3.m[7]  =  -Rc.getADValue(11)*S;
  J3.m[9]  =  -Rc.getADValue( 9)*S;                  
  J3.m[10] =  -Rc.getADValue(10)*S;                  
  J3.m[11] =  -Rc.getADValue(11)*S;
  J3.m[13] =  Qrc.getADValue( 9)*S;                  
  J3.m[14] =  Qrc.getADValue(10)*S;
  J3.m[15] =  Qrc.getADValue(11)*S;
  MatSetValues(*jtmp,4,index3,4,index3,J3.m,ADD_VALUES); 
}



//-----------------------------------------------------------------------------
// boundaries
//-----------------------------------------------------------------------------

void SMCZone::F2E_ddm_inner(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  PetscScalar kb = mt->kb;
  PetscScalar e  = mt->e;
  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];      //lattice temperature of node i

  f[zofs[zone_index]+4*i+0] = f[zofs[zone_index]+4*i+0]/pcell->area + e/aux[i].eps*((pi-aux[i].Na)+(aux[i].Nd-ni));
  f[zofs[zone_index]+4*i+1] = f[zofs[zone_index]+4*i+1]/pcell->area;
  f[zofs[zone_index]+4*i+2] = f[zofs[zone_index]+4*i+2]/pcell->area;
  f[zofs[zone_index]+4*i+3] = f[zofs[zone_index]+4*i+3]/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }

}


void SMCZone::F2E_ddm_neumannbc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  NeumannBC *pbc = dynamic_cast <NeumannBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar kb = mt->kb;
  PetscScalar e  = mt->e;
  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];      //lattice temperature of node i

  PetscScalar grad_T=0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[zone_index]+4*nb+3];     //lattice temperature of nb node
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar h = pbc->heat_transfer;
      PetscScalar r = h*pbc->T_external;
      grad_T += 0.5*pcell->ilen[j]*(r-0.25*h*(3*Ti+Tj));
    }
  }

  f[zofs[zone_index]+4*i+0] = f[zofs[zone_index]+4*i+0]/pcell->area + e/aux[i].eps*((pi-aux[i].Na)+(aux[i].Nd-ni));
  f[zofs[zone_index]+4*i+1] = f[zofs[zone_index]+4*i+1]/pcell->area;
  f[zofs[zone_index]+4*i+2] = f[zofs[zone_index]+4*i+2]/pcell->area;
  f[zofs[zone_index]+4*i+3] = (f[zofs[zone_index]+4*i+3]+grad_T)/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }

}


void SMCZone::F2E_ddm_ombc_segment(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  OhmicBC *pbc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 4;
  int size = pzone->davcell.size();
  PetscScalar e  =  mt->e;
  PetscScalar kb = mt->kb;
  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = fabs(x[zofs[z]+4*i+3]);
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar nie = mt->band->nie(fs[i].T);
  PetscScalar Nc  = mt->band->Nc(fs[i].T);
  PetscScalar Nv  = mt->band->Nv(fs[i].T);
  PetscScalar H = 0;
  PetscScalar electron_density,hole_density;
  int    om_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1)
    {
      om_equ=j;
      break;
    }
  }
  
  if(Fermi) //Fermi
  {
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T)); 
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) + mt->band->Eg(fs[i].T));
    PetscScalar phin = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar phip = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    f[zofs[z]+4*i+0] = Nc*fermi_half(etan) - Nv*fermi_half(etap)+ Na - Nd;
    f[zofs[z]+4*i+1] = ni - Nc*fermi_half(etan);
    f[zofs[z]+4*i+2] = pi - Nv*fermi_half(etap);
  }
  else   //Boltzmann
  {
    f[zofs[z]+4*i+0] = Vi - kb*fs[i].T/mt->e*asinh((Nd-Na)/(2*nie)) + mt->kb*fs[i].T/2/e*log(Nc/Nv) 
                       + mt->band->Eg(fs[i].T)/2/mt->e 
                       + aux[i].affinity -x[zofs[z]+equ_num*size+om_equ];
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
    f[zofs[z]+4*i+1] = x[zofs[z]+4*i+1] - electron_density;  //electron density
    f[zofs[z]+4*i+2] = x[zofs[z]+4*i+2] - hole_density;      //hole density
  }
  
  PetscScalar grad_T=0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[zone_index]+4*nb+3];            //lattice temperature of nb node
    PetscScalar T_mid    = 0.5*(Tj+Ti);
    PetscScalar dTdx_mid = (Tj-Ti)/pcell->ilen[j];
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    grad_T += kapa*pcell->elen[j]*dTdx_mid/pcell->area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar h = pbc->heat_transfer;
      PetscScalar r = h*pbc->T_external;
      grad_T += 0.5*pcell->ilen[j]*(r-0.25*h*(3*Ti+Tj))/pcell->area;
    }
  }
  f[zofs[zone_index]+4*i+3] =   grad_T + H;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }

}


void SMCZone::F2E_ddm_ombc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                     ElZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  OhmicBC *pbc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 4;
  int size = pzone->davcell.size();
  PetscScalar e  =  mt->e;
  PetscScalar kb = mt->kb;
  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = fabs(x[zofs[z]+4*i+3]);
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar nie = mt->band->nie(fs[i].T);
  PetscScalar Nc  = mt->band->Nc(fs[i].T);
  PetscScalar Nv  = mt->band->Nv(fs[i].T);
  PetscScalar H = 0;
  PetscScalar electron_density,hole_density;
  int    om_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1)
    {
      om_equ=j;
      break;
    }
  }
  
  if(Fermi) //Fermi
  {
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T)); 
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) + mt->band->Eg(fs[i].T));
    PetscScalar phin = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar phip = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    f[zofs[z]+4*i+0] = Nc*fermi_half(etan) - Nv*fermi_half(etap)+ Na - Nd;
    f[zofs[z]+4*i+1] = ni - Nc*fermi_half(etan);
    f[zofs[z]+4*i+2] = pi - Nv*fermi_half(etap);
  }
  else   //Boltzmann
  {
    f[zofs[z]+4*i+0] = Vi - kb*fs[i].T/mt->e*asinh((Nd-Na)/(2*nie)) + mt->kb*fs[i].T/2/e*log(Nc/Nv) 
                       + mt->band->Eg(fs[i].T)/2/mt->e 
                       + aux[i].affinity -x[zofs[z]+equ_num*size+om_equ];

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
    f[zofs[z]+4*i+1] = x[zofs[z]+4*i+1] - electron_density;  //electron density
    f[zofs[z]+4*i+2] = x[zofs[z]+4*i+2] - hole_density;      //hole density
  }
  
  PetscScalar grad_T=0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[zone_index]+4*nb+3];            //lattice temperature of nb node
    PetscScalar T_mid    = 0.5*(Tj+Ti);
    PetscScalar dTdx_mid = (Tj-Ti)/pcell->ilen[j];
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    grad_T += kapa*pcell->elen[j]*dTdx_mid;
  }
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  pz->mt->mapping(&pz->pzone->danode[n],&pz->aux[n],ODE_F.clock);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int   nb = ncell->nb_array[j];
    PetscScalar  Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];     //potential of nb node
    PetscScalar kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    grad_T += kapa*ncell->elen[j]/ncell->ilen[j]*(Tj_n-Ti);
  }
  f[zofs[zone_index]+4*i+3] =   grad_T + H;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity1*Tt/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*Tt/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity1*(Ti-fs[i].T)/ODE_F.dt*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*(Ti-fs[i].T)/ODE_F.dt*ncell->area;
    }
  }

}


void SMCZone::F2E_ddm_stkbc_segment(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int>& zofs,DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 4;
  int size = pzone->davcell.size();
  PetscScalar kb = mt->kb;
  PetscScalar e  = mt->e;
  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];      //lattice temperature of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);

  int  stk_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
    {
      stk_equ=j;
      break;
    }
  //Schotty Barrier Lowerring
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  //potential
  f[zofs[z]+4*i+0] = x[zofs[z]+4*i+0] + pbc->WorkFunction - deltaVB - x[zofs[z]+equ_num*size+stk_equ];

  PetscScalar Fn=0,Fp=0,grad_T=0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[zone_index]+4*nb+3];     //lattice temperature of nb node
    if(j==0||j==pcell->nb_num-1)
    {
      Fn += mt->band->SchottyJsn(ni,Ti,pbc->WorkFunction-aux[i].affinity - deltaVB)*0.5*pcell->ilen[j];
      Fp += mt->band->SchottyJsp(pi,Ti,pbc->WorkFunction-aux[i].affinity + deltaVB)*0.5*pcell->ilen[j];
      PetscScalar h = pbc->heat_transfer;
      PetscScalar r = h*pbc->T_external;
      grad_T += 0.5*pcell->ilen[j]*(r-0.25*h*(3*Ti+Tj));
    }
  }

  f[zofs[zone_index]+4*i+1] = (f[zofs[zone_index]+4*i+1]+Fn)/pcell->area;
  f[zofs[zone_index]+4*i+2] = (f[zofs[zone_index]+4*i+2]-Fp)/pcell->area;
  f[zofs[zone_index]+4*i+3] = (f[zofs[zone_index]+4*i+3]+grad_T)/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }
}


void SMCZone::F2E_ddm_stkbc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int>& zofs,DABC &bc,
                                      ElZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 4;
  int size = pzone->davcell.size();
  PetscScalar kb = mt->kb;
  PetscScalar e  = mt->e;
  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];      //lattice temperature of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);

  int  stk_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
    {
      stk_equ=j;
      break;
    }
  //Schotty Barrier Lowerring
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  //potential
  f[zofs[z]+4*i+0] = x[zofs[z]+4*i+0] + pbc->WorkFunction - deltaVB - x[zofs[z]+equ_num*size+stk_equ];

  PetscScalar Fn=0,Fp=0,grad_T=0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    if(j==0||j==pcell->nb_num-1)
    {
      Fn += mt->band->SchottyJsn(ni,Ti,pbc->WorkFunction-aux[i].affinity - deltaVB)*0.5*pcell->ilen[j];
      Fp += mt->band->SchottyJsp(pi,Ti,pbc->WorkFunction-aux[i].affinity + deltaVB)*0.5*pcell->ilen[j];
    }
  }
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  pz->mt->mapping(&pz->pzone->danode[n],&pz->aux[n],ODE_F.clock);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int   nb = ncell->nb_array[j];
    PetscScalar  Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];     //potential of nb node
    PetscScalar kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    grad_T += kapa*ncell->elen[j]/ncell->ilen[j]*(Tj_n-Ti);
  }
  f[zofs[zone_index]+4*i+1] = (f[zofs[zone_index]+4*i+1]+Fn)/pcell->area;
  f[zofs[zone_index]+4*i+2] = (f[zofs[zone_index]+4*i+2]-Fp)/pcell->area;
  f[zofs[zone_index]+4*i+3] = (f[zofs[zone_index]+4*i+3]+grad_T);

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity1*Tt/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*Tt/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt;
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity1*(Ti-fs[i].T)/ODE_F.dt*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*(Ti-fs[i].T)/ODE_F.dt*ncell->area;
    }
  }
}


void SMCZone::F2E_ddm_insulator_gate(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int equ_num = 4;
  int size = pzone->davcell.size();
  PetscScalar kb = mt->kb;
  PetscScalar e  = mt->e;
  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];      //lattice temperature of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar grad_P = 0, grad_T = 0;
  int    ins_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
    {ins_equ=j;break;}
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
    PetscScalar Tj = x[zofs[zone_index]+4*nb+3];            //lattice temperature of nb node

    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar vgate = x[zofs[zone_index]+equ_num*size+ins_equ] - pbc->WorkFunction;
      PetscScalar q = e*pbc->QF; //sigma is the surface change density
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar eps = aux[i].eps;
      PetscScalar r=q/eps + eps_ox/eps/Thick*vgate;
      PetscScalar s=eps_ox/eps/Thick;
      grad_P += 0.5*pcell->ilen[j]*(r-0.25*s*(3*Vi+Vj));

      PetscScalar h = pbc->heat_transfer;
      grad_T += 0.5*pcell->ilen[j]*(h*pbc->T_external-0.25*h*(3*Ti+Tj));
    }
  }

  f[zofs[zone_index]+4*i+0] = (f[zofs[zone_index]+4*i+0]+grad_P)/pcell->area
                              + e/aux[i].eps*((pi-aux[i].Na)+(aux[i].Nd-ni));
  f[zofs[zone_index]+4*i+1] = f[zofs[zone_index]+4*i+1]/pcell->area;
  f[zofs[zone_index]+4*i+2] = f[zofs[zone_index]+4*i+2]/pcell->area;
  f[zofs[zone_index]+4*i+3] = (f[zofs[zone_index]+4*i+3]+grad_T)/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt;
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }
}


void SMCZone::F2E_ddm_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                ISZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorInterfaceBC *pbc = dynamic_cast<InsulatorInterfaceBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar kb = mt->kb;
  PetscScalar e  = mt->e;
  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];      //lattice temperature of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar grad_P = 0,grad_T=0;

  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
    grad_P += aux[i].eps*pcell->elen[j]/pcell->ilen[j]*(Vj-Vi);
  }
  //semiconductor-insulator interface
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  pz->mt->mapping(&pz->pzone->danode[n],&pz->aux[n],ODE_F.clock);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int   nb = ncell->nb_array[j];
    PetscScalar  Vj_n = x[zofs[pz->pzone->zone_index]+2*nb+0];     //potential of nb node
    PetscScalar  Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];     //potential of nb node
    grad_P += pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j]*(Vj_n-Vi);
    PetscScalar kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    grad_T += kapa*ncell->elen[j]/ncell->ilen[j]*(Tj_n-Ti);
  }

  f[zofs[zone_index]+4*i+0] = grad_P + e*((pi-aux[i].Na)+(aux[i].Nd-ni))*pcell->area + pbc->QF*L;
  f[zofs[zone_index]+4*i+1] = f[zofs[zone_index]+4*i+1]/pcell->area;
  f[zofs[zone_index]+4*i+2] = f[zofs[zone_index]+4*i+2]/pcell->area;
  f[zofs[zone_index]+4*i+3] = (f[zofs[zone_index]+4*i+3]+grad_T);

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity1*Tt/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*Tt/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt;
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity1*(Ti-fs[i].T)/ODE_F.dt*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*(Ti-fs[i].T)/ODE_F.dt*ncell->area;
    }
  }

}


void SMCZone::F2E_ddm_homojunction(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                   SMCZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  PetscScalar kb = mt->kb;
  PetscScalar e  = mt->e;
  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];      //lattice temperature of node i

  if(zone_index > pz->pzone->zone_index)
  {
    f[zofs[zone_index]+4*i+0] =  Vi - x[zofs[pz->zone_index]+4*n+0];
    f[zofs[zone_index]+4*i+1] =  ni - x[zofs[pz->zone_index]+4*n+1];
    f[zofs[zone_index]+4*i+2] =  pi - x[zofs[pz->zone_index]+4*n+2];
    f[zofs[zone_index]+4*i+3] =  Ti - x[zofs[pz->zone_index]+4*n+3];
    return;
  }

  //process half cell in local zone
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  PetscScalar grad_P=0,grad_T=0;
  PetscScalar Fn=0,Fp=0;

  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
    grad_P += aux[i].eps*pcell->elen[j]/pcell->ilen[j]*(Vj-Vi);
  }
  Fn     =  f[zofs[zone_index]+4*i+1];
  Fp     =  f[zofs[zone_index]+4*i+2];
  grad_T =  f[zofs[zone_index]+4*i+3];
  //process another half cell in dornor zone
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);

  for(int j=0;j<ncell->nb_num;j++)
  {
    int  nb = ncell->nb_array[j];
    PetscScalar Vj  = x[zofs[pz->zone_index]+4*nb+0];     //potential of nb node
    grad_P += pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j]*(Vj-Vi);
  }
  Fn     +=  f[zofs[pz->zone_index]+4*n+1];
  Fp     +=  f[zofs[pz->zone_index]+4*n+2];
  grad_T +=  f[zofs[pz->zone_index]+4*n+3];

  PetscScalar area = pcell->area+ncell->area;

  f[zofs[zone_index]+4*i+0] =   grad_P + e*(pi-ni+Nd-Na)*area;
  f[zofs[zone_index]+4*i+1] =   Fn;
  f[zofs[zone_index]+4*i+2] =   Fp;
  f[zofs[zone_index]+4*i+3] =   grad_T;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt)*area;
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt)*area;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt)*area;
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt*area;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt*area;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt*area;
    }
  }
}


void SMCZone::F2E_ddm_heterojunction(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> &                                          zofs, DABC &bc, SMCZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  HeteroInterfaceBC *pbc = dynamic_cast<HeteroInterfaceBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar kb = mt->kb;
  PetscScalar e  =  mt->e;

  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];
  PetscScalar Vt = kb*Ti/e;

  PetscScalar nd = x[zofs[pz->pzone->zone_index]+4*n+1];
  PetscScalar pd = x[zofs[pz->pzone->zone_index]+4*n+2];

  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar Eci = -e*(Vi + aux[i].affinity + kb*fs[i].T/e*(log(aux[i].Nc)-1.5*log(fs[i].T)));
  PetscScalar Evi = -e*(Vi + aux[i].affinity - kb*fs[i].T/e*(log(aux[i].Nv)-1.5*log(fs[i].T)) + aux[i].Eg/e);

  PetscScalar grad_P=0;
  PetscScalar Fn=0,Fp=0;
  PetscScalar L = (pcell->ilen[0]+pcell->ilen[pcell->nb_num-1])/2.0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
    PetscScalar nj = x[zofs[zone_index]+4*nb+1];     //electron density of nb node
    PetscScalar pj = x[zofs[zone_index]+4*nb+2];     //hole density of nb node
    PetscScalar Tj = x[zofs[zone_index]+4*nb+3];
    grad_P += aux[i].eps*pcell->elen[j]/pcell->ilen[j]*(Vj-Vi);
    //Thermal emit current
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Ecj = -e*(Vi + pz->aux[n].affinity + kb*pz->fs[n].T/e*(log(pz->aux[n].Nc)-1.5*log(pz->fs[n].T)));
      PetscScalar Evj = -e*(Vi + pz->aux[n].affinity - kb*pz->fs[n].T/e*(log(pz->aux[n].Nv)-1.5*log(pz->fs[n].T)) + pz->aux[n].Eg/e);
      nj  = x[zofs[pz->pzone->zone_index]+4*n+1];
      pj  = x[zofs[pz->pzone->zone_index]+4*n+2];
      PetscScalar jn=0,jp=0;
      if(Ecj > Eci)
      {
        PetscScalar pm = pz->mt->band->EffecElecMass(fs[i].T)/mt->band->EffecElecMass(fs[i].T);
        jn = -2*e*(pz->mt->band->ThermalVn(fs[i].T)*nj - pm*mt->band->ThermalVn(fs[i].T)*ni*exp(-(Ecj-Eci)/(kb*fs[i].T)));
      }
      else
      {
        PetscScalar pm = mt->band->EffecElecMass(fs[i].T)/pz->mt->band->EffecElecMass(fs[i].T);
        jn = -2*e*(mt->band->ThermalVn(fs[i].T)*ni - pm*pz->mt->band->ThermalVn(fs[i].T)*nj*exp(-(Eci-Ecj)/(kb*fs[i].T)));
      }

      if(Evj > Evi)
      {
        PetscScalar pm = pz->mt->band->EffecHoleMass(fs[i].T)/mt->band->EffecHoleMass(fs[i].T);
        jp = 2*e*(pz->mt->band->ThermalVp(fs[i].T)*pj - pm*mt->band->ThermalVp(fs[i].T)*pi*exp(-(Evj-Evi)/(kb*fs[i].T)));
      }
      else
      {
        PetscScalar pm = mt->band->EffecHoleMass(fs[i].T)/pz->mt->band->EffecHoleMass(fs[i].T);
        jp = 2*e*(mt->band->ThermalVp(fs[i].T)*pi - pm*pz->mt->band->ThermalVp(fs[i].T)*pj*exp(-(Evi-Evj)/(kb*fs[i].T)));
      }

      Fn += jn*0.5*pcell->ilen[j];
      Fp += jp*0.5*pcell->ilen[j];
    }

  }

  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int  nb = ncell->nb_array[j];
    PetscScalar Vr = x[zofs[pz->pzone->zone_index]+4*nb+0];     //potential of nb node
    grad_P += pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j]*(Vr-Vi);
  }

  if(pzone->zone_index < pz->pzone->zone_index)
  {
    f[zofs[zone_index]+4*i+0] =  grad_P + e*(pi-ni+Nd-Na)*pcell->area
                                 + e*(pd-nd+pz->aux[n].Nd-pz->aux[n].Na)*ncell->area
                                 + pbc->QF*L;
    f[zofs[zone_index]+4*i+3] =  f[zofs[zone_index]+4*i+3] + f[zofs[pz->zone_index]+4*n+3];
  }
  else
  {
    f[zofs[zone_index]+4*i+0] =  Vi - x[zofs[pz->pzone->zone_index]+4*n+0];
    f[zofs[zone_index]+4*i+3] =  Ti - x[zofs[pz->pzone->zone_index]+4*n+3];
  }

  f[zofs[zone_index]+4*i+1] = (f[zofs[zone_index]+4*i+1]+Fn)/pcell->area;
  f[zofs[zone_index]+4*i+2] = (f[zofs[zone_index]+4*i+2]-Fp)/pcell->area;


  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F2E_om_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 4;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(bc_index));
  PetscScalar e  = mt->e;
  PetscScalar current=0;

  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+4*node+0];     //potential of node i
    current += DeviceDepth*(f[zofs[zone_index]+4*node+1]-f[zofs[zone_index]+4*node+2]);
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
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
        if(om_bc->inner_connect==bc_index)                connect_elec=j;
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
    f[zofs[zone_index]+equ_num*size+i]=(L/ODE_F.dt+R)*current-V+(1+(L/ODE_F.dt+R)*C/ODE_F.dt)*P
                                       -(L/ODE_F.dt+R)*C/ODE_F.dt*Pn-L/ODE_F.dt*(In+Icn);
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


void SMCZone::F2E_stk_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 4;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc.Get_pointer(bc_index));

  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+4*node+0];     //electron density of node i
    PetscScalar ni = x[zofs[zone_index]+4*node+1];     //electron density of node i
    PetscScalar pi = x[zofs[zone_index]+4*node+2];     //hole density of node i
    PetscScalar Ti = x[zofs[zone_index]+4*node+3];
    mt->mapping(&pzone->danode[node],&aux[node],ODE_F.clock);
    PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[node].eps,sqrt(aux[node].Ex*aux[node].Ex+aux[node].Ey*aux[node].Ey));
    PetscScalar Jsn = mt->band->SchottyJsn(ni,Ti,pbc->WorkFunction-aux[node].affinity-deltaVB);
    PetscScalar Jsp = mt->band->SchottyJsp(pi,Ti,pbc->WorkFunction-aux[node].affinity+deltaVB);
    current += -Jsn*0.5*pcell->ilen[0]*DeviceDepth;
    current += -Jsp*0.5*pcell->ilen[0]*DeviceDepth;
    current += -Jsn*0.5*pcell->ilen[pcell->nb_num-1]*DeviceDepth;
    current += -Jsp*0.5*pcell->ilen[pcell->nb_num-1]*DeviceDepth;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
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

void SMCZone::F2E_ins_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 4;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  InsulatorContactBC *pbc = dynamic_cast <InsulatorContactBC * > (bc.Get_pointer(bc_index));

  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+4*node+0];     //electron density of node i
    for(int k=0;k<pcell->nb_num;k++)
    {
      int  nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
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

void SMCZone::J2E_ddm_inner(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs)
{
  Mat4    A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int z = zone_index;
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];     //lattice temperature of node i
  PetscScalar area = pcell->area;

  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;

    MatGetValues(*jtmp,4,index,4,col,A.m);
    A=A/area;
    MatSetValues(*jac,4,index,4,col,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,4,index,4,index,A.m);
  A=A/area;
  A.m[1] +=  -e/aux[i].eps;                                //dfun(0)/dn(i)
  A.m[2] +=   e/aux[i].eps;                                //dfun(0)/dp(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[5]  += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar        HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[5]  += -1/ODE_F.dt;
      A.m[10] += -1/ODE_F.dt;
      PetscScalar        HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
}


void SMCZone::J2E_ddm_neumannbc(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  Mat4    A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  NeumannBC *pbc = dynamic_cast <NeumannBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i

  PetscScalar area = pcell->area;
  PetscScalar d_grad_T_dTi = 0;

  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node

    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;

    MatGetValues(*jtmp,4,index,4,col,A.m);
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar h = pbc->heat_transfer;
      d_grad_T_dTi += -0.5*pcell->ilen[j]*0.25*h*3/pcell->area;
      A.m[15] += -0.5*pcell->ilen[j]*0.25*h/pcell->area;
    }
    MatSetValues(*jac,4,index,4,col,A.m,INSERT_VALUES);

  }
  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  MatGetValues(*jtmp,4,index,4,index,A.m);
  A=A/area;
  A.m[1] +=  -e/aux[i].eps;                                //dfun(0)/dn(i)
  A.m[2] +=   e/aux[i].eps;                                //dfun(0)/dp(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[5] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[5] += -1/ODE_F.dt;
      A.m[10] += -1/ODE_F.dt;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
}



void SMCZone::J2E_ddm_ombc_segment(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  Mat4    A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  OhmicBC *pbc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int om_equ;
  int equ_num = 4;
  int size = pzone->davcell.size();
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1)
    {
      om_equ=j;
      break;
    }
  }
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i
  PetscScalar Vt = kb*Ti/e;
  PetscScalar area = pcell->area;
  PetscScalar d_grad_T_dTi = 0;
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node
    PetscScalar kapa =  mt->thermal->HeatConduction(0.5*(Ti+Tj));
    //-------------------------------------
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j]/area;
    //-------------------------------------
    PetscScalar d_grad_T_dTj =  kapa*pcell->elen[j]/pcell->ilen[j]/area; //dfun(3)/dT(r)
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar h = pbc->heat_transfer;
      d_grad_T_dTi += -0.5*pcell->ilen[j]*0.25*h*3/pcell->area;
      d_grad_T_dTj += -0.5*pcell->ilen[j]*0.25*h/pcell->area;
    }
    MatSetValue(*jac,zofs[z]+4*i+3,zofs[z]+4*nb+3,d_grad_T_dTj,INSERT_VALUES);
  }
  if(Fermi)  //Fermi
  {
    PetscScalar Vt =  kb*fs[i].T/e;
    mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
    PetscScalar Nc  = mt->band->Nc(fs[i].T);
    PetscScalar Nv  = mt->band->Nv(fs[i].T);
    PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) );
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) + mt->band->Eg(fs[i].T));
    PetscScalar phin = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar phip = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    PetscScalar detan_dVi =  1.0/Vt;
    PetscScalar detap_dVi = -1.0/Vt;
    Set_Mat4_zero(A);
    A.m[0] = Nc*fermi_mhalf(etan)*detan_dVi - Nv*fermi_mhalf(etap)*detap_dVi;  
    A.m[4] = -Nc*fermi_mhalf(etan)*detan_dVi;  
    A.m[5] = 1;  
    A.m[8] = -Nv*fermi_mhalf(etap)*detap_dVi;  
    A.m[10] = 1;
    A.m[15] =  d_grad_T_dTi;
    MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
    //dfermi_half(eta)_dVapp = - dfermi_half(eta)_dVi
    MatSetValue(*jac,zofs[zone_index]+4*i+0,zofs[zone_index]+equ_num*size+om_equ,-A.m[0],INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+4*i+1,zofs[zone_index]+equ_num*size+om_equ,-A.m[4],INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+4*i+2,zofs[zone_index]+equ_num*size+om_equ,-A.m[8],INSERT_VALUES);
  }
  else //Boltzmann
  {
    Set_Mat4_I(A);
    A.m[15] =  d_grad_T_dTi;     //dfun(3)/dT(i)
  }
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);

  MatSetValue(*jac,zofs[zone_index]+4*i+0,zofs[zone_index]+equ_num*size+om_equ,-1,INSERT_VALUES);

}


void SMCZone::J2E_ddm_ombc_interface(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc,
                                     ElZone *pz, int n)
{
  Mat4    A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  OhmicBC *pbc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int om_equ;
  int equ_num = 4;
  int size = pzone->davcell.size();
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1)
    {
      om_equ=j;
      break;
    }
  }
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i
  PetscScalar Vt = kb*Ti/e;
  PetscScalar area = pcell->area;
  PetscScalar d_grad_T_dTi = 0;
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node
    PetscScalar kapa =  mt->thermal->HeatConduction(0.5*(Ti+Tj));
    //-------------------------------------
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j];
    //-------------------------------------
    PetscScalar d_grad_T_dTj =  kapa*pcell->elen[j]/pcell->ilen[j]; //dfun(3)/dT(r)
    MatSetValue(*jac,zofs[z]+4*i+3,zofs[z]+4*nb+3,d_grad_T_dTj,INSERT_VALUES);
  }
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  pz->mt->mapping(&pz->pzone->danode[n],&pz->aux[n],ODE_F.clock);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int    nb = ncell->nb_array[j];
    PetscScalar   Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];
    PetscScalar   kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    d_grad_T_dTi += -kapa*ncell->elen[j]/ncell->ilen[j];
    PetscScalar d_grad_T_dTj = kapa*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[z]+4*i+3,zofs[pz->pzone->zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
  }

  if(Fermi)  //Fermi
  {
    PetscScalar Vt =  kb*fs[i].T/e;
    mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
    PetscScalar Nc  = mt->band->Nc(fs[i].T);
    PetscScalar Nv  = mt->band->Nv(fs[i].T);
    PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T)); 
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) + mt->band->Eg(fs[i].T));
    PetscScalar phin = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar phip = x[zofs[zone_index]+equ_num*size+om_equ];
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    PetscScalar detan_dVi =  1.0/Vt;
    PetscScalar detap_dVi = -1.0/Vt;
    Set_Mat4_zero(A);
    A.m[0] = Nc*fermi_mhalf(etan)*detan_dVi - Nv*fermi_mhalf(etap)*detap_dVi;  
    A.m[4] = -Nc*fermi_mhalf(etan)*detan_dVi;  
    A.m[5] = 1;  
    A.m[8] = -Nv*fermi_mhalf(etap)*detap_dVi;  
    A.m[10] = 1;
    A.m[15] =  d_grad_T_dTi;
    MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
    //dfermi_half(eta)_dVapp = - dfermi_half(eta)_dVi
    MatSetValue(*jac,zofs[zone_index]+4*i+0,zofs[zone_index]+equ_num*size+om_equ,-A.m[0],INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+4*i+1,zofs[zone_index]+equ_num*size+om_equ,-A.m[4],INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+4*i+2,zofs[zone_index]+equ_num*size+om_equ,-A.m[8],INSERT_VALUES);
  }
  else //Boltzmann
  {
    Set_Mat4_I(A);
    A.m[15] =  d_grad_T_dTi;     //dfun(3)/dT(i)
  }
  
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity1*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*pcell->area
               -pz->aux[n].density*HeatCapacity2*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity1/ODE_F.dt*pcell->area
               -pz->aux[n].density*HeatCapacity2/ODE_F.dt*ncell->area;
    }
  }

  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);

  MatSetValue(*jac,zofs[zone_index]+4*i+0,zofs[zone_index]+equ_num*size+om_equ,-1,INSERT_VALUES);

}


void SMCZone::J2E_ddm_stkbc_segment(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  Mat4   A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 4;
  int size = pzone->davcell.size();

  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar VB = pbc->WorkFunction-aux[i].affinity;
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i

  PetscScalar area = pcell->area;
  PetscScalar d_grad_T_dTi = 0;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fn_dTi = 0;
  PetscScalar d_Fp_dpi = 0;
  PetscScalar d_Fp_dTi = 0;
  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar nj = x[zofs[z]+4*nb+1];     //electron density of nb node
    PetscScalar pj = x[zofs[z]+4*nb+2];     //hole density of nb node
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node
    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;
    MatGetValues(*jtmp,4,index,4,col,A.m);
    A.m[0] = 0.0;
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      d_Fn_dni += mt->band->pdSchottyJsn_pdn (ni,Ti,VB-deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fn_dTi += mt->band->pdSchottyJsn_pdTl(ni,Ti,VB-deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fp_dpi += mt->band->pdSchottyJsp_pdp (pi,Ti,VB+deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fp_dTi += mt->band->pdSchottyJsp_pdTl(pi,Ti,VB+deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      PetscScalar h = pbc->heat_transfer;
      d_grad_T_dTi += -0.5*pcell->ilen[j]*0.25*h*3/pcell->area;
      A.m[15] += -0.5*pcell->ilen[j]*0.25*h/pcell->area;
    }
    MatSetValues(*jac,4,index,4,col,A.m,INSERT_VALUES);
  }
  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  MatGetValues(*jtmp,4,index,4,index,A.m);
  A=A/area;

  A.m[0] =  1;

  A.m[5] +=  d_Fn_dni;      //dfun(1)/dn(i)
  A.m[7] +=  d_Fn_dTi;

  A.m[10]+=  -d_Fp_dpi;     //dfun(2)/dp(i)
  A.m[11]+=  -d_Fp_dTi;

  A.m[15]+=  d_grad_T_dTi;  //dfun(3)/dT(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[5] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[5] += -1/ODE_F.dt;
      A.m[10] += -1/ODE_F.dt;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);

  int stk_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
    {
      stk_equ=j;
      break;
    }
  MatSetValue(*jac,zofs[zone_index]+4*i,zofs[zone_index]+equ_num*size+stk_equ,-1,INSERT_VALUES);
}


void SMCZone::J2E_ddm_stkbc_interface(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc,
                                      ElZone *pz, int n)
{
  Mat4   A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 4;
  int size = pzone->davcell.size();
  PetscScalar area = pcell->area;
  Vec4 scale;
  Set_Vec4(scale,area,area,area,1.0);
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar VB = pbc->WorkFunction-aux[i].affinity;
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i

  PetscScalar d_grad_T_dTi = 0;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fn_dTi = 0;
  PetscScalar d_Fp_dpi = 0;
  PetscScalar d_Fp_dTi = 0;
  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar nj = x[zofs[z]+4*nb+1];     //electron density of nb node
    PetscScalar pj = x[zofs[z]+4*nb+2];     //hole density of nb node
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node
    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;
    MatGetValues(*jtmp,4,index,4,col,A.m);
    A.m[0] = 0.0;
    A=A/scale;
    if(j==0||j==pcell->nb_num-1)
    {
      d_Fn_dni += mt->band->pdSchottyJsn_pdn (ni,Ti,VB-deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fn_dTi += mt->band->pdSchottyJsn_pdTl(ni,Ti,VB-deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fp_dpi += mt->band->pdSchottyJsp_pdp (pi,Ti,VB+deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fp_dTi += mt->band->pdSchottyJsp_pdTl(pi,Ti,VB+deltaVB)*0.5*pcell->ilen[j]/pcell->area;
    }
    MatSetValues(*jac,4,index,4,col,A.m,INSERT_VALUES);
  }
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  pz->mt->mapping(&pz->pzone->danode[n],&pz->aux[n],ODE_F.clock);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int    nb = ncell->nb_array[j];
    PetscScalar   Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];
    PetscScalar   kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    d_grad_T_dTi += -kapa*ncell->elen[j]/ncell->ilen[j];
    PetscScalar d_grad_T_dTj = kapa*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[z]+4*i+3,zofs[pz->pzone->zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
  }
  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  MatGetValues(*jtmp,4,index,4,index,A.m);
  A=A/scale;

  A.m[0] =  1;

  A.m[5] +=  d_Fn_dni;                              //dfun(1)/dn(i)
  A.m[7] +=  d_Fn_dTi;

  A.m[10]+=  -d_Fp_dpi;                                       //dfun(2)/dp(i)
  A.m[11]+=  -d_Fp_dTi;

  A.m[15]+=  d_grad_T_dTi; //dfun(3)/dT(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[5] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity1*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                 -pz->aux[n].density*HeatCapacity2*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      A.m[5] += -1/ODE_F.dt;
      A.m[10] += -1/ODE_F.dt;
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity1/ODE_F.dt*pcell->area
                 -pz->aux[n].density*HeatCapacity2/ODE_F.dt*ncell->area;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);

  int stk_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
    {
      stk_equ=j;
      break;
    }
  MatSetValue(*jac,zofs[zone_index]+4*i,zofs[zone_index]+equ_num*size+stk_equ,-1,INSERT_VALUES);
}


void SMCZone::J2E_ddm_insulator_gate(int i,PetscScalar *x, Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  Mat4    A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int  equ_num = 4;
  int  size = pzone->davcell.size();
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i

  PetscScalar area = pcell->area;
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar d_grad_P_dVi = 0;
  PetscScalar d_grad_P_dVapp = 0;
  PetscScalar d_grad_T_dTi = 0;
  int ins_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
    {ins_equ=j;break;}
  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[z]+4*nb+0];     //potential of nb node
    PetscScalar nj = x[zofs[z]+4*nb+1];     //electron density of nb node
    PetscScalar pj = x[zofs[z]+4*nb+2];     //hole density of nb node
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node
    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;
    MatGetValues(*jtmp,4,index,4,col,A.m);
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar eps = aux[i].eps;
      PetscScalar s=eps_ox/eps/Thick;
      d_grad_P_dVi += -0.5*pcell->ilen[j]*0.25*s*3/area;
      d_grad_P_dVapp += 0.5*pcell->ilen[j]*s/area;
      A.m[0]+= -0.5*pcell->ilen[j]*0.25*s/area;
      PetscScalar h = pbc->heat_transfer;
      d_grad_T_dTi += -0.5*pcell->ilen[j]*0.25*h*3/area;
      A.m[15] += -0.5*pcell->ilen[j]*0.25*h/area;
    }
    MatSetValues(*jac,4,index,4,col,A.m,INSERT_VALUES);

  }
  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  MatGetValues(*jtmp,4,index,4,index,A.m);
  A=A/area;
  A.m[0] +=  d_grad_P_dVi;
  A.m[1] +=  -e/aux[i].eps;                                //dfun(0)/dn(i)
  A.m[2] +=   e/aux[i].eps;                                //dfun(0)/dp(i)

  A.m[15]+=  d_grad_T_dTi; //dfun(3)/dT(i)
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[5] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[5] += -1/ODE_F.dt;
      A.m[10] += -1/ODE_F.dt;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
  MatSetValue(*jac,zofs[z]+4*i,zofs[z]+equ_num*size+ins_equ,d_grad_P_dVapp,INSERT_VALUES);
}


void SMCZone::J2E_ddm_interface(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                ISZone *pz, int n)
{
  Mat4    A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorInterfaceBC *pbc = dynamic_cast<InsulatorInterfaceBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];     //lattice temperature of node i

  PetscScalar area = pcell->area;
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar d_grad_P_dVi = 0;
  PetscScalar d_grad_T_dTi = 0;

  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];

    PetscScalar Vj = x[zofs[z]+4*nb+0];     //potential of nb node
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;
    MatGetValues(*jtmp,4,index,4,col,A.m);
    Vec4 scale;
    Set_Vec4(scale,1.0,area,area,1.0);
    A=A/scale;
    A.m[0] =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];  //dfun(0)/dP(r)
    A.m[1] =  0;                                 //dfun(0)/dn(r)
    A.m[2] =  0;                                 //dfun(0)/dp(r)
    A.m[3] =  0;                                 //dfun(0)/dT(r)

    MatSetValues(*jac,4,index,4,col,A.m,INSERT_VALUES);

  }

  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  pz->mt->mapping(&pz->pzone->danode[n],&pz->aux[n],ODE_F.clock);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int    nb = ncell->nb_array[j];
    PetscScalar   Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];
    PetscScalar   kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    d_grad_P_dVi += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    PetscScalar d_grad_P_dVj = pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[z]+4*i+0,zofs[pz->pzone->zone_index]+2*nb+0,d_grad_P_dVj,INSERT_VALUES);
    d_grad_T_dTi += -kapa*ncell->elen[j]/ncell->ilen[j];
    PetscScalar d_grad_T_dTj = kapa*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[z]+4*i+3,zofs[pz->pzone->zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
  }

  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  MatGetValues(*jtmp,4,index,4,index,A.m);

  A.m[0] =   d_grad_P_dVi;                                         //dfun(0)/dP(i)
  A.m[1] =  -e*pcell->area;                                         //dfun(0)/dn(i)
  A.m[2] =   e*pcell->area;                                         //dfun(0)/dp(i)
  A.m[3] =   0;                                                  //dfun(0)/dT(i)

  A.m[4]  =  A.m[4]/area;
  A.m[5]  =  A.m[5]/area;                                       //dfun(1)/dn(i)
  A.m[6]  =  A.m[6]/area;                                       //dfun(1)/dp(i)
  A.m[7]  =  A.m[7]/area;

  A.m[8]  =  A.m[8]/area;
  A.m[9]  =  A.m[9]/area;                                       //dfun(2)/dn(i)
  A.m[10] =  A.m[10]/area;                                      //dfun(2)/dp(i)
  A.m[11] =  A.m[11]/area;

  A.m[15] += d_grad_T_dTi; //dfun(3)/dT(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[5] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity1*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                 -pz->aux[n].density*HeatCapacity2*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      A.m[5] += -1/ODE_F.dt;
      A.m[10] += -1/ODE_F.dt;
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity1/ODE_F.dt*pcell->area
                 -pz->aux[n].density*HeatCapacity2/ODE_F.dt*ncell->area;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
}


void SMCZone::J2E_ddm_homojunction(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, SMCZone *pz, int n)
{
  Mat4    A1,A2,AJ;
  PetscInt       index[4],indexn[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int z = zone_index;
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];     //lattice temperature of node i

  PetscScalar d_grad_P_dVi = 0;
  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  //--------------------------------
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  indexn[0] = zofs[pz->zone_index]+4*n+0;
  indexn[1] = zofs[pz->zone_index]+4*n+1;
  indexn[2] = zofs[pz->zone_index]+4*n+2;
  indexn[3] = zofs[pz->zone_index]+4*n+3;
  //--------------------------------
  if(zone_index > pz->pzone->zone_index)
  {
    Set_Mat4_I(AJ);
    MatSetValues(*jac,4,index,4,index,AJ.m,INSERT_VALUES);
    AJ = -1.0*AJ;
    MatSetValues(*jac,4,index,4,indexn,AJ.m,INSERT_VALUES);
    return;
  }

  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;
    //-------------------------------------
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    //-------------------------------------
    MatGetValues(*jtmp,4,index,4,col,A1.m);
    //-------------------------------------
    //df(x)i/dx(r)
    A1.m[0] =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];   //dfun(0)/dP(r)
    A1.m[1] =  0;                                       //dfun(0)/dn(r)
    A1.m[2] =  0;                                       //dfun(0)/dp(r)
    A1.m[3] =  0;
    MatSetValues(*jac,4,index,4,col,A1.m,INSERT_VALUES);
  }

  //process another half cell in dornor zone

  for(int j=0;j<ncell->nb_num;j++)
  {
    int  nb = ncell->nb_array[j];
    col[0] = zofs[pz->zone_index]+4*nb+0;
    col[1] = zofs[pz->zone_index]+4*nb+1;
    col[2] = zofs[pz->zone_index]+4*nb+2;
    col[3] = zofs[pz->zone_index]+4*nb+3;
    //-------------------------------------
    d_grad_P_dVi += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    //-------------------------------------
    //df(x)i/dx(r)
    MatGetValues(*jtmp,4,indexn,4,col,A2.m);
    A2.m[0] =  pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];       //dfun(0)/dP(r)
    A2.m[1] =  0;                                //dfun(0)/dn(r)
    A2.m[2] =  0;                                //dfun(0)/dp(r)
    A2.m[3] =  0;                                //dfun(0)/dT(r)
    MatSetValues(*jac,4,index,4,col,A2.m,INSERT_VALUES);
  }

  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  PetscScalar area = pcell->area + ncell->area;

  MatGetValues(*jtmp,4,index,4,index,A1.m);
  MatGetValues(*jtmp,4,indexn,4,indexn,A2.m);
  AJ=A1+A2;

  AJ.m[0] =  d_grad_P_dVi;    //dfun(0)/dP(i)
  AJ.m[1] =  -e*area;  //dfun(0)/dn(i)
  AJ.m[2] =   e*area;  //dfun(0)/dp(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      AJ.m[5] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*area;
      AJ.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*area;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      AJ.m[15] += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*area;
    }
    else //first order
    {
      AJ.m[5] += -1/ODE_F.dt*area;
      AJ.m[10] += -1/ODE_F.dt*area;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      AJ.m[15] += -aux[i].density*HeatCapacity/ODE_F.dt*area;
    }
  }
  MatSetValues(*jac,4,index,4,index,AJ.m,INSERT_VALUES);
}


void SMCZone::J2E_ddm_heterojunction(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F,
                                     vector<int> &zofs, DABC &bc,SMCZone *pz, int n)
{
  Mat4    A,An;
  PetscInt       index[4],indexn[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);

  int z = zone_index;
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i

  // conduct band and valance band
  PetscScalar Eci = -e*(Vi + aux[i].affinity + kb*fs[i].T/e*(log(aux[i].Nc)-1.5*log(fs[i].T)));
  PetscScalar Evi = -e*(Vi + aux[i].affinity - kb*fs[i].T/e*(log(aux[i].Nv)-1.5*log(fs[i].T)) + aux[i].Eg/e);

  PetscScalar area = pcell->area;
  PetscScalar d_grad_P_dVi = 0;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fp_dpi = 0;
  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;

  indexn[0] = zofs[pz->zone_index]+4*n+0;
  indexn[1] = zofs[pz->zone_index]+4*n+1;
  indexn[2] = zofs[pz->zone_index]+4*n+2;
  indexn[3] = zofs[pz->zone_index]+4*n+3;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[z]+4*nb+0];     //potential of nb node
    PetscScalar nj = x[zofs[z]+4*nb+1];     //electron density of nb node
    PetscScalar pj = x[zofs[z]+4*nb+2];     //hole density of nb node
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node

    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;

    PetscScalar Ecj = -e*(Vj + aux[nb].affinity + kb*fs[nb].T/e*(log(aux[nb].Nc)-1.5*log(fs[nb].T)));
    PetscScalar Evj = -e*(Vj + aux[nb].affinity - kb*fs[nb].T/e*(log(aux[nb].Nv)-1.5*log(fs[nb].T)) + aux[nb].Eg/e);

    //-------------------------------------
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    //-------------------------------------
    MatGetValues(*jtmp,4,index,4,col,A.m);
    Vec4 scale;
    Set_Vec4(scale,1.0,area,area,1.0);
    A=A/scale;
    //df(x)i/dx(r)
    if(pzone->zone_index < pz->pzone->zone_index)
    {
      A.m[0] =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];  //dfun(0)/dP(r)
      A.m[1] =  0;                                       //dfun(0)/dn(r)
      A.m[2] =  0;                                       //dfun(0)/dp(r)
      A.m[3] =  0;                                       //dfun(0)/dT(r)

      A.m[12] =  A.m[12];   //dfun(3)/dP(r)
      A.m[13] =  A.m[13];   //dfun(3)/dn(r)
      A.m[14] =  A.m[14];   //dfun(3)/dp(r)
      A.m[15] =  A.m[15];
    }
    else
    {
      A.m[0] =  0;                                 //dfun(0)/dP(r)
      A.m[1] =  0;                                       //dfun(0)/dn(r)
      A.m[2] =  0;                                       //dfun(0)/dp(r)
      A.m[3] =  0;                                       //dfun(0)/dT(r)

      A.m[12] =  0;   //dfun(3)/dP(r)
      A.m[13] =  0;   //dfun(3)/dn(r)
      A.m[14] =  0;   //dfun(3)/dp(r)
      A.m[15] =  0;
    }

    //Thermal emit current
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Ecj = -e*(Vi + pz->aux[n].affinity + kb*pz->fs[n].T/e*(log(pz->aux[n].Nc)-1.5*log(pz->fs[n].T)));
      PetscScalar Evj = -e*(Vi + pz->aux[n].affinity - kb*pz->fs[n].T/e*(log(pz->aux[n].Nv)-1.5*log(pz->fs[n].T)) + pz->aux[n].Eg/e);

      MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

      if(Ecj > Eci)
      {
        PetscScalar pm = pz->mt->band->EffecElecMass(fs[i].T)/mt->band->EffecElecMass(fs[i].T);
        d_Fn_dni += 2*e*pm*mt->band->ThermalVn(fs[i].T)*exp(-(Ecj-Eci)/(kb*fs[i].T))*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fn_dnd = -2*e*pz->mt->band->ThermalVn(fs[i].T)*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+4*i+1,zofs[pz->pzone->zone_index]+4*n+1,d_Fn_dnd,ADD_VALUES);
      }
      else
      {
        PetscScalar pm = mt->band->EffecElecMass(fs[i].T)/pz->mt->band->EffecElecMass(fs[i].T);
        d_Fn_dni += -2*e*mt->band->ThermalVn(fs[i].T)*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fn_dnd = 2*e*pm*pz->mt->band->ThermalVn(fs[i].T)*exp(-(Eci-Ecj)/(kb*fs[i].T))*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+4*i+1,zofs[pz->pzone->zone_index]+4*n+1,d_Fn_dnd,ADD_VALUES);
      }

      if(Evj > Evi)
      {
        PetscScalar pm = pz->mt->band->EffecHoleMass(fs[i].T)/mt->band->EffecHoleMass(fs[i].T);
        d_Fp_dpi += -2*e*pm*mt->band->ThermalVp(fs[i].T)*exp(-(Evj-Evi)/(kb*fs[i].T))*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fp_dpd = 2*e*pz->mt->band->ThermalVp(fs[i].T)*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+4*i+2,zofs[pz->pzone->zone_index]+4*n+2,-d_Fp_dpd,ADD_VALUES);
      }
      else
      {
        PetscScalar pm = mt->band->EffecHoleMass(fs[i].T)/pz->mt->band->EffecHoleMass(fs[i].T);
        d_Fp_dpi += 2*e*mt->band->ThermalVp(fs[i].T)*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fp_dpd = -2*e*pm*pz->mt->band->ThermalVp(fs[i].T)*exp(-(Evi-Evj)/(kb*fs[i].T))*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+4*i+2,zofs[pz->pzone->zone_index]+4*n+2,-d_Fp_dpd,ADD_VALUES);
      }
      MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
    }
    MatSetValues(*jac,4,index,4,col,A.m,INSERT_VALUES);

  }

  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  MatGetValues(*jtmp,4,index,4,index,A.m);
  MatGetValues(*jtmp,4,indexn,4,indexn,An.m);

  if(pzone->zone_index < pz->pzone->zone_index)
  {
    A.m[0] =  d_grad_P_dVi;    //dfun(0)/dP(i)
    A.m[1] =  -e*pcell->area;  //dfun(0)/dn(i)
    A.m[2] =   e*pcell->area;  //dfun(0)/dp(i)
    A.m[3] =   0;  //dfun(0)/dT(i)
    for(int j=0;j<ncell->nb_num;j++)
    {
      int    nb = ncell->nb_array[j];
      PetscScalar kapa =  pz->mt->thermal->HeatConduction(Ti);
      A.m[0] += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
      A.m[15] += -kapa*ncell->elen[j]/ncell->ilen[j];
      PetscScalar Pvar = pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
      MatSetValue(*jac,zofs[z]+4*i+0,zofs[pz->pzone->zone_index]+4*nb+0,Pvar,INSERT_VALUES);
      PetscScalar Tvar = kapa*ncell->elen[j]/ncell->ilen[j];
      MatSetValue(*jac,zofs[z]+4*i+3,zofs[pz->pzone->zone_index]+4*nb+3,Tvar,INSERT_VALUES);
    }
    MatSetValue(*jac,zofs[z]+4*i+0,zofs[pz->pzone->zone_index]+4*n+1,-e*ncell->area,INSERT_VALUES);
    MatSetValue(*jac,zofs[z]+4*i+0,zofs[pz->pzone->zone_index]+4*n+2,e*ncell->area,INSERT_VALUES);
  }
  else
  {
    A.m[0] =   1;    //dfun(0)/dP(i)
    A.m[1] =   0;  //dfun(0)/dn(i)
    A.m[2] =   0;  //dfun(0)/dp(i)
    A.m[3] =   0;  //dfun(0)/dT(i)

    A.m[12] =  0;   //dfun(3)/dP(i)
    A.m[13] =  0;   //dfun(3)/dn(i)
    A.m[14] =  0;   //dfun(3)/dp(i)
    A.m[15] =  1;   //dfun(3)/dT(i)

    MatSetValue(*jac,zofs[z]+4*i+0,zofs[pz->pzone->zone_index]+4*n+0,-1,INSERT_VALUES);
    MatSetValue(*jac,zofs[z]+4*i+3,zofs[pz->pzone->zone_index]+4*n+3,-1,INSERT_VALUES);
  }

  A.m[4]  =  A.m[4]/area;
  A.m[5]  =  A.m[5]/area + d_Fn_dni;
  A.m[6]  =  A.m[6]/area;
  A.m[7]  =  A.m[7]/area;
  A.m[8]  =  A.m[8]/area;
  A.m[9]  =  A.m[9]/area;
  A.m[10] =  A.m[10]/area - d_Fp_dpi;
  A.m[11] =  A.m[11]/area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[5] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[5] += -1/ODE_F.dt;
      A.m[10] += -1/ODE_F.dt;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);

}


void SMCZone::J2E_om_electrode(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 4;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(bc_index));
  PetscScalar    A1[4],A2[4],J[4];
  PetscInt       index[4],col[4];
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
        if(om_bc->inner_connect==bc_index)                connect_elec=j;
      }
    }
    connect_index=zofs[connect_zone]+equ_num*pz->pzone->davcell.size()+connect_elec;
  }
  PetscScalar R = pbc->R;
  PetscScalar C = pbc->C;
  PetscScalar L = pbc->L;
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;

  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

  for(int n=0;n<bc[bc_index].psegment->node_array.size();n++)
  {
    int node=bc[bc_index].psegment->node_array[n];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar dJdisp_dVi=0,dJdisp_dVr=0;
    index[0] = zofs[zone_index]+4*node+0;
    index[1] = zofs[zone_index]+4*node+1;
    index[2] = zofs[zone_index]+4*node+2;
    index[3] = zofs[zone_index]+4*node+3;
    for(int j=0;j<pcell->nb_num;j++)
    {
      int    nb = pcell->nb_array[j];
      col[0] = zofs[zone_index]+4*nb+0;
      col[1] = zofs[zone_index]+4*nb+1;
      col[2] = zofs[zone_index]+4*nb+2;
      col[3] = zofs[zone_index]+4*nb+3;

      MatGetValues(*jtmp,1,&index[1],4,col,A1);
      MatGetValues(*jtmp,1,&index[2],4,col,A2);

      //for displacement current
      dJdisp_dVi += aux[node].eps/pcell->ilen[j]/ODE_F.dt*pcell->elen[j];
      dJdisp_dVr = -aux[node].eps/pcell->ilen[j]/ODE_F.dt*pcell->elen[j];
      if(pbc->inner_connect!=-1)
      {
        int connect_to = pbc->inner_connect;
        if(bc_index < connect_to)
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
          J[3]=(A1[3]-A2[3])*DeviceDepth;
          MatSetValues(*jac,1,&connect_index,4,col,J,ADD_VALUES);
        }
        else
        {
          J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
          J[1]=(A1[1]-A2[1])*DeviceDepth;
          J[2]=(A1[2]-A2[2])*DeviceDepth;
          J[3]=(A1[3]-A2[3])*DeviceDepth;
          MatSetValues(*jac,1,&matrix_row,4,col,J,ADD_VALUES);
        }
      }
      else if(pbc->electrode_type==VoltageBC)
      {
        J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth*(L/ODE_F.dt+R);
        J[1]=(A1[1]-A2[1])*DeviceDepth*(L/ODE_F.dt+R);
        J[2]=(A1[2]-A2[2])*DeviceDepth*(L/ODE_F.dt+R);
        J[3]=(A1[3]-A2[3])*DeviceDepth*(L/ODE_F.dt+R);
        MatSetValues(*jac,1,&matrix_row,4,col,J,ADD_VALUES);
      }
      else if(pbc->electrode_type==CurrentBC)
      {
        J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
        J[1]=(A1[1]-A2[1])*DeviceDepth;
        J[2]=(A1[2]-A2[2])*DeviceDepth;
        J[3]=(A1[3]-A2[3])*DeviceDepth;
        MatSetValues(*jac,1,&matrix_row,4,col,J,ADD_VALUES);
      }
    }
    MatGetValues(*jtmp,1,&index[1],4,index,A1);
    MatGetValues(*jtmp,1,&index[2],4,index,A2);
    if(pbc->inner_connect!=-1)
    {
      int connect_to = pbc->inner_connect;
      if(bc_index < connect_to)
      {
        J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
        J[1]=(A1[1]-A2[1])*DeviceDepth;
        J[2]=(A1[2]-A2[2])*DeviceDepth;
        J[3]=(A1[3]-A2[3])*DeviceDepth;
        MatSetValues(*jac,1,&connect_index,4,index,J,ADD_VALUES);
      }
      else
      {
        J[0]=(A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
        J[1]=(A1[1]-A2[1])*DeviceDepth;
        J[2]=(A1[2]-A2[2])*DeviceDepth;
        J[3]=(A1[3]-A2[3])*DeviceDepth;
        MatSetValues(*jac,1,&matrix_row,4,index,J,ADD_VALUES);
      }
    }
    else if(pbc->electrode_type==VoltageBC)
    {
      J[0]=(A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth*(L/ODE_F.dt+R);
      J[1]=(A1[1]-A2[1])*DeviceDepth*(L/ODE_F.dt+R);
      J[2]=(A1[2]-A2[2])*DeviceDepth*(L/ODE_F.dt+R);
      J[3]=(A1[3]-A2[3])*DeviceDepth*(L/ODE_F.dt+R);
      MatSetValues(*jac,1,&matrix_row,4,index,J,ADD_VALUES);
    }
    else if(pbc->electrode_type==CurrentBC)
    {
      J[0]=(A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth;
      J[1]=(A1[1]-A2[1])*DeviceDepth;
      J[2]=(A1[2]-A2[2])*DeviceDepth;
      J[3]=(A1[3]-A2[3])*DeviceDepth;
      MatSetValues(*jac,1,&matrix_row,4,index,J,ADD_VALUES);
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
    MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP

}


void SMCZone::J2E_stk_electrode(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 4;
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
    PetscScalar ni = x[zofs[zone_index]+4*node+1];     //electron density of node i
    PetscScalar pi = x[zofs[zone_index]+4*node+2];     //hole density of node i
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
      PetscScalar dI_dVi =  DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
      PetscScalar dI_dVr = -DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
      if(pbc->electrode_type==VoltageBC)
      {
        MatSetValue(*jac,matrix_row,zofs[zone_index]+4*node+0,dI_dVi*(L/ODE_F.dt+R),ADD_VALUES);
        MatSetValue(*jac,matrix_row,zofs[zone_index]+4*nb+0,dI_dVr*(L/ODE_F.dt+R),ADD_VALUES);
      }
      else if(pbc->electrode_type==CurrentBC)
      {
        MatSetValue(*jac,matrix_row,zofs[zone_index]+4*node+0,dI_dVi,ADD_VALUES);
        MatSetValue(*jac,matrix_row,zofs[zone_index]+4*nb+0,dI_dVr,ADD_VALUES);
      }
    }
    MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

    if(pbc->electrode_type==VoltageBC)
    {
      MatSetValue(*jac,matrix_row,zofs[zone_index]+4*node+1,DeviceDepth*dJ_dni*(L/ODE_F.dt+R),INSERT_VALUES);
      MatSetValue(*jac,matrix_row,zofs[zone_index]+4*node+2,DeviceDepth*dJ_dpi*(L/ODE_F.dt+R),INSERT_VALUES);
    }
    else if(pbc->electrode_type==CurrentBC)
    {
      MatSetValue(*jac,matrix_row,zofs[zone_index]+4*node+1,DeviceDepth*dJ_dni,INSERT_VALUES);
      MatSetValue(*jac,matrix_row,zofs[zone_index]+4*node+2,DeviceDepth*dJ_dpi,INSERT_VALUES);
    }
  }
  if(pbc->electrode_type==VoltageBC)
    MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP
}


void SMCZone::J2E_ins_electrode(int i,PetscScalar *x,Mat *jac,Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 4;
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
        MatSetValue(*jac,matrix_row,zofs[zone_index]+4*node+0,dI_dVi*(L/ODE_F.dt+R),ADD_VALUES);
        MatSetValue(*jac,matrix_row,zofs[zone_index]+4*nb+0,dI_dVr*(L/ODE_F.dt+R),ADD_VALUES);
      }
    }

  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

  if(pbc->electrode_type==VoltageBC)
    MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP
}


void SMCZone::F2E_om_electrode_Trace(int electrode_index, PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,
                                     Vec & x, Mat *jtmp, Mat *jac,ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
{
  int equ_num = 4;
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
  PetscScalar    A1[4],A2[4];
  PetscInt       index[4],col[4];
  //delete electrode current equation, omit the effect of external resistance
  MatSetValuesRow(*jac,extern_equ,apdF_pdV);
  MatSetValue(*jac,extern_equ,extern_equ,1.0,INSERT_VALUES);

  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = xx[zofs[zone_index]+4*node+0];     //potential of node i
    PetscScalar dJdisp_dVi=0,dJdisp_dVr=0;
    index[0] = zofs[zone_index]+4*node+0;
    index[1] = zofs[zone_index]+4*node+1;
    index[2] = zofs[zone_index]+4*node+2;
    index[3] = zofs[zone_index]+4*node+3;

    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = xx[zofs[zone_index]+4*nb+0];     //potential of nb node
      col[0] = zofs[zone_index]+4*nb+0;
      col[1] = zofs[zone_index]+4*nb+1;
      col[2] = zofs[zone_index]+4*nb+2;
      col[3] = zofs[zone_index]+4*nb+3;

      MatGetValues(*jtmp,1,&index[1],4,col,A1);
      MatGetValues(*jtmp,1,&index[2],4,col,A2);

      //for displacement current
      dJdisp_dVi += aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];
      dJdisp_dVr = -aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];

      //
      apdI_pdw[col[0]] += (A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
      apdI_pdw[col[1]] += (A1[1]-A2[1])*DeviceDepth;
      apdI_pdw[col[2]] += (A1[2]-A2[2])*DeviceDepth;
      apdI_pdw[col[3]] += (A1[3]-A2[3])*DeviceDepth;
    }
    MatGetValues(*jtmp,1,&index[1],4,index,A1);
    MatGetValues(*jtmp,1,&index[2],4,index,A2);

    apdI_pdw[index[0]] += (A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth;
    apdI_pdw[index[1]] += (A1[1]-A2[1])*DeviceDepth;
    apdI_pdw[index[2]] += (A1[2]-A2[2])*DeviceDepth;
    apdI_pdw[index[3]] += (A1[3]-A2[3])*DeviceDepth;

    apdF_pdV[index[0]] = 1.0;
    pdI_pdV = 0;
  }

  MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);
  VecRestoreArray(x,&xx);
  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}


void SMCZone::F2E_stk_electrode_Trace(int electrode_index,  PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,
                                      Vec & x, Mat *jtmp, Mat *jac,ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
{
  int equ_num = 4;
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
    PetscScalar ni = xx[zofs[zone_index]+4*node+1];     //electron density of node i
    PetscScalar pi = xx[zofs[zone_index]+4*node+2];     //hole density of node i
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
      apdI_pdw[zofs[zone_index]+4*node+0] += dI_dVi;
      apdI_pdw[zofs[zone_index]+4*nb+0]   += dI_dVr;
    }
    apdI_pdw[zofs[zone_index]+4*node+1] = DeviceDepth*dJ_dni;
    apdI_pdw[zofs[zone_index]+4*node+2] = DeviceDepth*dJ_dpi;

    apdF_pdV[zofs[zone_index]+4*node+0] = 1.0;
    pdI_pdV = 0;
  }
  MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);
  VecRestoreArray(x,&xx);
  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}


void SMCZone::F2E_efield_update(PetscScalar *x,vector<int> & zofs, DABC &bc, vector<BZoneData *>zonedata)
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
      PetscScalar dP = x[zofs[zone_index]+4*pcell->nb_array[k]+0] - x[zofs[zone_index]+4*i+0];
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
        PetscScalar dP = x[zofs[n_zone]+2*dcell->nb_array[k]]- x[zofs[n_zone]+2*n_node];
        w=1.0/sqrt(dx*dx+dy*dy);
        P_x+=w*w*dP*dx;
        P_y+=w*w*dP*dy;
      }
    }
    aux[i].Ex = -(pcell->sc*P_x-pcell->sb*P_y);
    aux[i].Ey = -(pcell->sa*P_y-pcell->sb*P_x);
  }
}
