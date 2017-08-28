(* ::Package:: *)

(* 
SMFpack.m
a Mathematica package for implementing genus 2 Siegel Modular Forms

The definitions used in https://arxiv.org/abs/0710.2129 are used. 

Author : Shouvik Datta
*)



(* The following SMFs are defined*)

Psi4::usage = "Siegel Eisenstein series of weight 4 : \!\(\*SubscriptBox[\(\[Psi]\), \(4\)]\). Takes the period matrix (\[CapitalOmega]) as its argument.";
Psi6::usage = "Siegel Eisenstein series of weight 6 : \!\(\*SubscriptBox[\(\[Psi]\), \(6\)]\). Takes the period matrix (\[CapitalOmega]) as its argument.";
Chi10::usage = "Siegel Cusp form of weight 10 : \!\(\*SubscriptBox[\(\[Chi]\), \(10\)]\). Takes the period matrix (\[CapitalOmega]) as its argument.";
Chi12::usage = "Siegel Cusp form of weight 12 : \!\(\*SubscriptBox[\(\[Chi]\), \(12\)]\). Takes the period matrix (\[CapitalOmega]) as its argument.";



(* ::Input::Initialization:: *)


siegelTheta[a1_,a2_,b1_,b2_,\[CapitalOmega]_]:=SiegelTheta[{{a1/2,a2/2},{b1/2,b2/2}},\[CapitalOmega],0];

(*non vanishing list*)
nvl[\[CapitalOmega]_]:={siegelTheta[0,0,0,0,\[CapitalOmega]],siegelTheta[0,0,0,1,\[CapitalOmega]],siegelTheta[0,0,1,0,\[CapitalOmega]],siegelTheta[0,0,1,1,\[CapitalOmega]],siegelTheta[0,1,0,0,\[CapitalOmega]],siegelTheta[0,1,1,0,\[CapitalOmega]],siegelTheta[1,0,0,0,\[CapitalOmega]],siegelTheta[1,0,0,1,\[CapitalOmega]],siegelTheta[1,1,0,0,\[CapitalOmega]],siegelTheta[1,1,1,1,\[CapitalOmega]]};

(*p functions*)
P0[\[CapitalOmega]_]:=(1/4) * Plus@@({siegelTheta[0,0,0,0,\[CapitalOmega]],siegelTheta[0,0,0,1,\[CapitalOmega]],siegelTheta[0,0,1,0,\[CapitalOmega]],siegelTheta[0,0,1,1,\[CapitalOmega]]}^4);
P1[\[CapitalOmega]_]:=(siegelTheta[0,1,0,0,\[CapitalOmega]]^4+siegelTheta[0,1,1,0,\[CapitalOmega]]^4)/4;
P2[\[CapitalOmega]_]:=(siegelTheta[1,0,0,0,\[CapitalOmega]]^4+siegelTheta[1,0,0,1,\[CapitalOmega]]^4)/4;
P3[\[CapitalOmega]_]:=(siegelTheta[1,1,0,0,\[CapitalOmega]]^4+siegelTheta[1,1,1,1,\[CapitalOmega]]^4)/4;
P4[\[CapitalOmega]_]:=(siegelTheta[0,1,0,0,\[CapitalOmega]]^4-siegelTheta[0,1,1,1,\[CapitalOmega]]^4)/4;

(*siegel modular forms*)
\[Psi]T4[\[CapitalOmega]_]:=(1/4) Plus@@(nvl[\[CapitalOmega]]^8);
\[CapitalDelta]10[\[CapitalOmega]_]:=(1/2^12) Times@@(nvl[\[CapitalOmega]]^2);
F12[\[CapitalOmega]_]:=(1/4) Plus@@(nvl[\[CapitalOmega]]^24);
\[Psi]T6[W_]:=P0[W]^3-9P0[W](P1[W]^2+P2[W]^2+P3[W]^2-4P4[W]^2)+54P1[W]*P2[W]*P3[W];

(*conventional normalization*)
\[Psi]4[w_] := (1/4)*\[Psi]T4[w]; 
\[Psi]6[w_] := (1/16)*\[Psi]T6[w]; 
\[Chi]10[w_] := 4*\[CapitalDelta]10[w]; 
\[Chi]12[w_] := (96/(2^13*3^4))*((9/11)*F12[w] - \[Psi]T4[w]^3 + (2/11)*\[Psi]T6[w]^2); 

Psi4[w_]:=\[Psi]4[w];
Psi6[w_]:=\[Psi]6[w];
Chi10[w_]:=\[Chi]10[w];
Chi12[w_]:=\[Chi]12[w];




Column[{
"SMFpack \n\n"<>
"A Mathematica package for implementing Siegel Modular forms\n"<>
"Author : Shouvik Datta\n"<>
"The definitions used are in Appendix A of https://arxiv.org/abs/0710.2129 \n"
}]  
