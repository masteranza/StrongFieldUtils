(* ::Package:: *)

(* ::Title:: *)
(*StrongFieldUtils*)


(* ::Author:: *)
(*author: Michal Mandrysz*)


(* ::Affiliation:: *)
(*Marian Smoluchowski Institute of Physics, Jagiellonian University, Krakow, Poland*)


(* ::Abstract:: *)
(*Collection of functions useful in Strong Field Atomic physics*)


(* ::Text:: *)
(*Version: 1.0.0*)


BeginPackage["StrongFieldUtils`"];
Begin["`Private`"];


(* ::Section:: *)
(*Essenatials*)


(* ::Text:: *)
(*Constants in Atomic Units*)


StrongFieldUtils`fsAtomic::usage="One femtosecond in atomic units";
fsAtomic=QuantityMagnitude@UnitConvert[Quantity[1,"fs"],"AtomicUnitOfTime"];

StrongFieldUtils`cAtomic::usage="Speed of light in atomic units";
cAtomic=QuantityMagnitude@UnitConvert[Quantity["ReducedPlanckConstant"]/( Quantity["FineStructureConstant"] Quantity["BohrRadius"] Quantity["ElectronMass"]),"AtomicUnitOfVelocity"];

StrongFieldUtils`IpAtomic::usage="IpAtomic[atom,n] returns n'th (def. n) ionization energy of \"atom\" in atomic units.";
IpAtomic[atom_String,n_:1]:=QuantityMagnitude[UnitConvert[ElementData[atom,"IonizationEnergies"]/Quantity["AvogadroConstant"],"Hartrees"][[n]]];


(* ::Text:: *)
(*Useful conversions*)


StrongFieldUtils`EfromIAtomic::usage="EfromIAtomic[I] takes intensity (in W/cm^2) and transforms it to electric field intensity in atomic units.";
EfromIAtomic[I_]:=QuantityMagnitude@UnitConvert[Sqrt[Quantity[I,"W/cm^2"]/(1/2*Quantity[1,"ElectricConstant"]*Quantity[cAtomic,"AtomicUnitOfVelocity"])]*Quantity[1,"ElementaryCharge"]*Quantity[1,"BohrRadius"],"Hartrees"];

StrongFieldUtils`IfromEAtomic::usage="IfromEAtomic[E] takes electric field in atomic units and transforms it to intensity in W/cm^2.";
IfromEAtomic[E_]:=QuantityMagnitude@UnitConvert[(Quantity[E,"Hartrees"]/(Quantity[1,"ElementaryCharge"]*Quantity[1,"BohrRadius"]))^2 (1/2*Quantity[1,"ElectricConstant"]*Quantity[cAtomic,"AtomicUnitOfVelocity"]),"W/cm^2"];

StrongFieldUtils`\[Omega]from\[Lambda]Atomic::usage="Given number of nm turns into \[Omega] in AtomicUnits";
\[Omega]from\[Lambda]Atomic[nm_]:=QuantityMagnitude@UnitConvert[Quantity[1,"SpeedOfLight"]Quantity[1,"PlanckConstant"]/Quantity[nm,"nm"],"Hartrees"];

StrongFieldUtils`Tfrom\[Omega]Atomic::usage="Given \[Omega] [a.u.] turns it into period [a.u.]";
Tfrom\[Omega]Atomic[\[Omega]_]:=QuantityMagnitude@UnitConvert[ Quantity[1,"PlanckConstant"]/Quantity[\[Omega],"Hartrees"],"AtomicUnitOfTime"];

StrongFieldUtils`QuiverAtomic::usage="Quiver[E,\[Omega]] returns E/\!\(\*SuperscriptBox[\(\[Omega]\), \(2\)]\).";
QuiverAtomic[E_,\[Omega]_]:=E/\[Omega]^2;

StrongFieldUtils`PonderoAtomic::usage="Pondero[E,\[Omega]] returns \!\(\*SuperscriptBox[\(E\), \(2\)]\)/(2\[Omega]\!\(\*SuperscriptBox[\()\), \(2\)]\)";
PonderoAtomic[E_,\[Omega]_]:=E^2/(2\[Omega])^2;

StrongFieldUtils`CorkumCutoff::usage="CorkumCutoff[\!\(\*SubscriptBox[\(I\), \(p\)]\),E,\[Omega]]";
CorkumCutoff[Ip_,E_,\[Omega]_]:=(Ip+3.17 PonderoAtomic[E,\[Omega]])/\[Omega];


(* ::Text:: *)
(*HHG spectrum*)


StrongFieldUtils`CalcSpectrum::usage="CalcSpectrum[{<time>,<dipole position>},window,freq,\[Omega]power] returns HHG spectrum calculated (in the desired atomic freq) from transposed dipole displacement time series with a cosine window of strength [0,1]";
CalcSpectrum[d_,window_:0,freq_:1,\[Omega]power_:0,cutoff_:None]:=Block[{t,dip0,df,half,magnitude,dt,freqAxis},
t=d[[1]];dip0=d[[2]];
dt=t[[-1]]-t[[-2]];
dip0*=Most@Table[CosineWindow[i-0.5,window],{i,0,1,1/Length@dip0}];
df=2\[Pi]/((Length[dip0]-1)*dt);
magnitude=Fourier[dip0];
freqAxis=If[cutoff===None,Range[0,Length[magnitude]-1]df/freq,Range[0,cutoff*freq/df]df/freq];
If[cutoff=!=None,magnitude=magnitude[[;;Length[freqAxis]]]];
magnitude=If[\[Omega]power>0,Power[freqAxis^\[Omega]power Abs[magnitude],2],Power[Abs[magnitude],2]];
{freqAxis, magnitude}\[Transpose]
];

stylings=Sequence[FrameStyle->Directive[FontFamily->"Latin Modern Math"],FrameTicksStyle->22,LabelStyle->Directive[FontSize->20,FontFamily->"Latin Modern Math"], Frame->True];
StrongFieldUtils`PresentSpectrum::usage="PresentSpectrum[data,maxHarmonics:{},label,minP,Z,legend]";
StrongFieldUtils`TheoCutoff::usage="Pass this to mark the theoretical cutoff position";
Options[PresentSpectrum]={TheoCutoff->{}};
PresentSpectrum[data_,opt:OptionsPattern[]]:=ListLinePlot[data,ScalingFunctions->"Log", stylings, FrameLabel->{"Harmonic order","Intensity"},GridLines->Automatic,Prolog->{Black,Thick,Line[{{#,-100},{#,100}}&/@ OptionValue[PresentSpectrum,FilterRules[{opt},Options[PresentSpectrum]],TheoCutoff]]},FilterRules[{opt},Options[ListLinePlot]]];

StrongFieldUtils`GaborSpectrum::usage="GaborSpectrum";
GaborSpectrum[d_,\[Omega]_:1,cutoff_:200,points_:260,windowWidthMultiplier_:1]:=Block[{t,dip0,ldip0,df,dt,len,res,parts,hann,\[Sigma]},
t=d[[1]];
dip0=d[[2]];
ldip0=Length[dip0];
dt=t[[-1]]-t[[-2]];
\[Sigma]=windowWidthMultiplier \[Pi]/(cutoff/8*\[Omega]*dt);
len=IntegerPart[Length[dip0]/points];
df=2\[Pi]/((ldip0-1)*dt);
parts=Table[i,{i,2len,Length[dip0]-1-2len,len}];
res=ParallelTable[Reverse@Log[Abs[Fourier[dip0*Table[N@Exp[-(j-i)^2/(2 \[Sigma]^2)],{j,1,ldip0}]][[1;;Round[cutoff*\[Omega]/df]]]]^2],{i,parts}]\[Transpose];
MatrixPlot[res,ColorFunction->"SunsetColors",AspectRatio->1/2,DataRange->{{t[[First@parts]],t[[Last@parts]]},{0,cutoff}}]
];


(* ::Section:: *)
(*Utils*)


StrongFieldUtils`LoadColumn::usage="LoadColumn[path, {list of columns}] returns all rows and the selected columns of a file at a given path";
LoadColumn[path_,pick_:{1},howmuch_:All]:=Block[{d},d=Import[path,"Data"][[howmuch,pick]]];


(* ::Section:: *)
(*End package*)


End[];
EndPackage[];
