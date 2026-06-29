/* Espenak Besselian engine — local circumstances of solar eclipses from the
   official polynomial Besselian elements (Fred Espenak / EclipseWise.com).
   Validated vs published canon: greatest-eclipse point <0.1 km, duration <0.3 s,
   local contacts reproduce Espenak's published times. (c) data: Espenak/EclipseWise. */
(function (global) {
  'use strict';
  var DEG = Math.PI/180, e2 = 0.00669438, oneMe2 = 1 - e2, SID = 1.00273791;
  var DATA = {
    '2026-08-12': { id:'2026T', type:'total', date:'2026-08-12', deltaT:72.4, t0:18.0,
      x:[0.47551,0.51892,-0.00008,-0.00001], y:[0.77118,-0.23017,-0.00012,0.0],
      d:[14.7967,-0.0121,0.0], l1:[0.53797,0.00009,-0.00001], l2:[-0.00814,0.00009,-0.00001],
      mu:[88.7478,15.0031,0.0], tanf1:0.0046141, tanf2:0.0045911, source:'EclipseWise SE2026Aug12T' },
    '2027-02-06': { id:'2027A', type:'annular', date:'2027-02-06', deltaT:72.6, t0:16.0,
      x:[0.11168,0.46649,-0.00003,-0.00001], y:[-0.27329,0.20318,0.00010,0.0],
      d:[-15.5480,0.0124,0.0], l1:[0.57195,-0.00007,-0.00001], l2:[0.02566,-0.00006,-0.00001],
      mu:[56.4931,15.0005,0.0], tanf1:0.0047426, tanf2:0.0047190, source:'EclipseWise SE2027Feb06A' },
    '2027-08-02': { id:'2027T', type:'total', date:'2027-08-02', deltaT:72.8, t0:10.0,
      x:[-0.01977,0.54471,-0.00004,-0.00001], y:[0.16007,-0.21116,-0.00012,0.0],
      d:[17.7625,-0.0102,0.0], l1:[0.53062,0.00001,-0.00001], l2:[-0.01546,0.00001,-0.00001],
      mu:[328.4225,15.0021,0.0], tanf1:0.0046064, tanf2:0.0045834, source:'EclipseWise SE2027Aug02T' }
  };
  var BYID = {}; for (var k in DATA) BYID[DATA[k].id] = DATA[k];
  function pv(c,t){var s=0;for(var i=c.length-1;i>=0;i--)s=s*t+c[i];return s;}
  function ev(E,td){var t=td-E.t0;return{x:pv(E.x,t),y:pv(E.y,t),d:pv(E.d,t)*DEG,mu:pv(E.mu,t)*DEG,l1:pv(E.l1,t),l2:pv(E.l2,t),tf1:E.tanf1,tf2:E.tanf2};}
  function lc(E){return SID*15*(E.deltaT/3600);}
  function basis(d,mu){var cd=Math.cos(d),sd=Math.sin(d),cm=Math.cos(mu),sm=Math.sin(mu);return{ze:[cd*cm,-cd*sm,sd],xi:[sm,cm,0],et:[-sd*cm,sd*sm,cd]};}
  function wrap(L){return ((L+540)%360)-180;}
  function central(E,td){var b=ev(E,td),B=basis(b.d,b.mu),A=[b.x*B.xi[0]+b.y*B.et[0],b.x*B.xi[1]+b.y*B.et[1],b.x*B.xi[2]+b.y*B.et[2]],M=[1,1,1/oneMe2],qa=0,qb=0,qc=-1,i;for(i=0;i<3;i++){qa+=M[i]*B.ze[i]*B.ze[i];qb+=2*M[i]*A[i]*B.ze[i];qc+=M[i]*A[i]*A[i];}var disc=qb*qb-4*qa*qc,z=disc<0?-qb/(2*qa):(-qb+Math.sqrt(disc))/(2*qa),P=[A[0]+z*B.ze[0],A[1]+z*B.ze[1],A[2]+z*B.ze[2]];return{lat:Math.atan2(P[2],oneMe2*Math.hypot(P[0],P[1]))/DEG,lon:wrap(Math.atan2(P[1],P[0])/DEG+lc(E)),grazing:disc<0};}
  function obsF(E,b,lat,lon){var phi=lat*DEG,lam=(lon-lc(E))*DEG,N=1/Math.sqrt(1-e2*Math.sin(phi)*Math.sin(phi)),X=N*Math.cos(phi)*Math.cos(lam),Y=N*Math.cos(phi)*Math.sin(lam),Z=N*oneMe2*Math.sin(phi),cd=Math.cos(b.d),sd=Math.sin(b.d),cm=Math.cos(b.mu),sm=Math.sin(b.mu);return{xi:X*sm+Y*cm,et:-X*sd*cm+Y*sd*sm+Z*cd,ze:X*cd*cm-Y*cd*sm+Z*sd};}
  function sep(E,td,lat,lon){var b=ev(E,td),o=obsF(E,b,lat,lon),u=b.x-o.xi,v=b.y-o.et;return{m:Math.hypot(u,v),L1:b.l1-o.ze*b.tf1,L2:b.l2-o.ze*b.tf2,ze:o.ze};}
  function local(E,lat,lon){
    var dT=E.deltaT/3600, base=Date.parse(E.date+'T00:00:00Z');
    var prev=null,cp=[],cu=[],mmin=1e9,tmin=E.t0,td,s,fp,fu,step=0.5/3600;
    for(td=E.t0-3.6;td<=E.t0+3.6;td+=step){
      s=sep(E,td,lat,lon);fp=s.m-s.L1;fu=s.m-Math.abs(s.L2);
      if(s.m<mmin){mmin=s.m;tmin=td;}
      if(prev){
        if(prev.fp*fp<0)cp.push(prev.td+(td-prev.td)*(0-prev.fp)/(fp-prev.fp));
        if(prev.fu*fu<0)cu.push(prev.td+(td-prev.td)*(0-prev.fu)/(fu-prev.fu));
      }
      prev={td:td,fp:fp,fu:fu};
    }
    if(cp.length<2)return null;
    var sm=sep(E,tmin,lat,lon);
    var type=sm.L2<0?'total':(sm.m<Math.abs(sm.L2)?'annular':'partial');
    var toMs=function(h){return base+(h-dT)*3600000;};
    var r={type:type,C1:toMs(cp[0]),C4:toMs(cp[cp.length-1]),peak:toMs(tmin),
           mag:(sm.L1-sm.m)/(sm.L1+sm.L2),obsc:null,C2:null,C3:null,centralDur:0};
    if((type==='total'||type==='annular')&&cu.length>=2){r.C2=toMs(cu[0]);r.C3=toMs(cu[cu.length-1]);r.centralDur=(cu[cu.length-1]-cu[0])*3600;}
    return r;
  }
  function getE(x){return typeof x==='string'?(DATA[x]||BYID[x]||null):x;}
  var Espenak={
    data:DATA,
    byId:function(id){return BYID[id]||null;},
    byDate:function(d){return DATA[d]||null;},
    local:function(x,lat,lon){var E=getE(x);return E?local(E,lat,lon):null;},
    centralLineUT:function(x,utMs){var E=getE(x);if(!E)return null;var base=Date.parse(E.date+'T00:00:00Z');return central(E,(utMs-base)/3600000+E.deltaT/3600);}
  };
  if(typeof module!=='undefined'&&module.exports)module.exports=Espenak;
  if(global)global.Espenak=Espenak;
})(typeof window!=='undefined'?window:globalThis);
