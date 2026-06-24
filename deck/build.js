const PptxGenJS = require("pptxgenjs");
const p = new PptxGenJS();
p.defineLayout({ name: "W", width: 13.333, height: 7.5 });
p.layout = "W";

const DARK="2C5F2D", GREEN="1B9E77", ORANGE="D95F02", MOSS="6FA046", INK="1F2421", MUTE="6B6B6B", LINE="DCE4DC", PAPER="FFFFFF";
const HF="Georgia", BF="Calibri";
const DIR="C:/Users/crop/OneDrive - University of Guelph/R code/Caleb/OFE_Stats/CEMvsPSM/";
const FIG=n=>DIR+n;

// helpers
const title=(s,t,sub)=>{ s.addText(t,{x:0.6,y:0.45,w:12.1,h:0.8,fontFace:HF,fontSize:32,bold:true,color:DARK});
  if(sub) s.addText(sub,{x:0.62,y:1.18,w:12.0,h:0.5,fontFace:BF,fontSize:15,italic:true,color:MUTE}); };

// fit an image of native ratio r=w/h into a box, top-left anchored, centered
function fitImg(s,path,r,bx,by,bw,bh){ let w=bw,h=bw/r; if(h>bh){h=bh;w=bh*r;} s.addImage({path,x:bx+(bw-w)/2,y:by+(bh-h)/2,w,h}); }

/* 1 — TITLE */
let s=p.addSlide(); s.background={color:DARK};
// plot-row motif: small squares, a few orange (dropped)
const drops=new Set([3,11,17]);
for(let i=0;i<22;i++){ s.addShape(p.ShapeType.rect,{x:0.85+i*0.52,y:6.35,w:0.34,h:0.34,fill:{color:drops.has(i)?ORANGE:GREEN},line:{color:DARK,width:1}});}
s.addText("Fair comparisons in unreplicated strip trials",{x:0.85,y:1.9,w:11.6,h:1.6,fontFace:HF,fontSize:44,bold:true,color:PAPER,lineSpacingMultiple:1.0});
s.addText("Coarsened Exact Matching for on-farm nitrogen response",{x:0.87,y:3.5,w:11.0,h:0.7,fontFace:BF,fontSize:22,color:"CADCBF"});
s.addText("Home farm  ·  corn  ·  30 vs 214 lb N/ac (240 kg/ha)  ·  one strip each (24 vs 24 plots)",{x:0.87,y:4.25,w:11.0,h:0.5,fontFace:BF,fontSize:15,italic:true,color:MOSS});

/* 1b — THE HOOK: do you believe it? */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Start with the trial as it comes off the combine","Home farm · corn · one high-N strip beside one control · 24 plots each · no replication");
const hkY=1.95, hkH=2.55;
s.addShape(p.ShapeType.roundRect,{x:1.0,y:hkY,w:3.4,h:hkH,rectRadius:0.08,fill:{color:"EFEFEF"},line:{color:"9A9A9A",width:1.5}});
s.addText([{text:"CONTROL",options:{fontSize:15,bold:true,color:MUTE,breakLine:true}},{text:"30 lb N/ac",options:{fontSize:13,color:MUTE,breakLine:true}},{text:"\n220",options:{fontSize:52,bold:true,color:"5C5C5C",breakLine:true}},{text:"bu/ac",options:{fontSize:13,color:MUTE}}],{x:1.0,y:hkY+0.15,w:3.4,h:hkH-0.25,align:"center",valign:"middle",fontFace:BF,lineSpacingMultiple:0.95});
s.addShape(p.ShapeType.roundRect,{x:4.95,y:hkY,w:3.4,h:hkH,rectRadius:0.08,fill:{color:"E7F1EA"},line:{color:DARK,width:1.5}});
s.addText([{text:"HIGH N",options:{fontSize:15,bold:true,color:DARK,breakLine:true}},{text:"214 lb N/ac",options:{fontSize:13,color:GREEN,breakLine:true}},{text:"\n263",options:{fontSize:52,bold:true,color:DARK,breakLine:true}},{text:"bu/ac",options:{fontSize:13,color:MUTE}}],{x:4.95,y:hkY+0.15,w:3.4,h:hkH-0.25,align:"center",valign:"middle",fontFace:BF,lineSpacingMultiple:0.95});
s.addShape(p.ShapeType.roundRect,{x:8.9,y:hkY,w:3.4,h:hkH,rectRadius:0.08,fill:{color:DARK},line:{type:"none"}});
s.addText([{text:"RAW GAP",options:{fontSize:14,bold:true,color:"CADCBF",breakLine:true}},{text:"\n+43",options:{fontSize:58,bold:true,color:PAPER,breakLine:true}},{text:"bu/ac",options:{fontSize:15,color:"CADCBF",breakLine:true}},{text:"95% CI 34–52",options:{fontSize:12,italic:true,color:MOSS}}],{x:8.9,y:hkY+0.15,w:3.4,h:hkH-0.25,align:"center",valign:"middle",fontFace:BF,lineSpacingMultiple:0.95});
s.addText([{text:"A big, clean-looking nitrogen response — ",options:{color:INK}},{text:"but it is one strip against one strip.",options:{bold:true,color:ORANGE}}],{x:1.0,y:5.05,w:11.3,h:0.5,fontFace:BF,fontSize:19,align:"center"});
s.addText("Would you change your nitrogen program on this number?",{x:1.0,y:5.7,w:11.3,h:0.9,fontFace:HF,fontSize:30,bold:true,color:DARK,align:"center"});

/* 2 — THE PROBLEM */
s=p.addSlide(); s.background={color:PAPER};
title(s,"The on-farm trial problem","Strips are not replicated — so they were never randomized");
// left text
s.addText([
 {text:"A strip trial = ",options:{bold:false}},{text:"one",options:{bold:true,color:DARK}},{text:" N-rich strip beside ",options:{}},{text:"one",options:{bold:true,color:DARK}},{text:" control strip.",options:{breakLine:true}},
 {text:"No replication means treatment is ",options:{breakLine:false}},{text:"perfectly confounded",options:{bold:true,color:ORANGE}},{text:" with wherever the two strips happen to sit.",options:{breakLine:true}},
 {text:"\nIf the strips fall on different soil, the raw yield gap is part fertilizer, part field.",options:{}}
],{x:0.6,y:1.9,w:6.0,h:3.4,fontFace:BF,fontSize:18,color:INK,lineSpacingMultiple:1.25,valign:"top"});
s.addText("Treated minus control ≠ the N effect — until you account for the ground.",{x:0.6,y:5.5,w:6.0,h:0.9,fontFace:BF,fontSize:16,italic:true,color:DARK});
// right: two strips over a soil gradient
const gx=7.4, gw=5.2, gy=2.0, gh=4.6;
const shades=["C9B79C","BBA888","AD9A74","9F8C61","91804F"];
for(let i=0;i<5;i++) s.addShape(p.ShapeType.rect,{x:gx,y:gy+i*(gh/5),w:gw,h:gh/5,fill:{color:shades[i]},line:{type:"none"}});
s.addShape(p.ShapeType.rect,{x:gx+0.7,y:gy+0.25,w:1.7,h:gh-0.5,fill:{color:GREEN,transparency:18},line:{color:"145C44",width:1.5}});
s.addShape(p.ShapeType.rect,{x:gx+2.8,y:gy+0.25,w:1.7,h:gh-0.5,fill:{color:"8A8A8A",transparency:18},line:{color:"5C5C5C",width:1.5}});
s.addText("Treated (+N)",{x:gx+0.7,y:gy+0.3,w:1.7,h:0.4,fontFace:BF,fontSize:12,bold:true,color:PAPER,align:"center"});
s.addText("Control",{x:gx+2.8,y:gy+0.3,w:1.7,h:0.4,fontFace:BF,fontSize:12,bold:true,color:PAPER,align:"center"});
s.addText("soil / terrain gradient ↓",{x:gx,y:gy+gh+0.05,w:gw,h:0.35,fontFace:BF,fontSize:11,italic:true,color:MUTE,align:"center"});

/* 2b — THE FIELD (scene-setter) */
s=p.addSlide(); s.background={color:PAPER};
title(s,"The field the trial sits in","Yield tracks the ground — the strips cross an eroded knoll");
fitImg(s,FIG("fig9_yieldmap.png"),8.5/6.5,0.3,1.55,8.3,5.6);
s.addText([
 {text:"Two strips, side by side, the length of the field.",options:{bold:true,color:DARK,breakLine:true}},
 {text:"\nThey can’t help but sit on different ground — soil and terrain shift under them. That difference is the confounding an unreplicated trial bakes in, and exactly what matching has to undo.",options:{}}
],{x:9.0,y:2.2,w:3.9,h:4.2,fontFace:BF,fontSize:16,color:INK,lineSpacingMultiple:1.3,valign:"top"});

/* 3 — WHAT CONFOUNDS THIS FIELD */
s=p.addSlide(); s.background={color:PAPER};
title(s,"What actually confounds the comparison?","A confounder must be BOTH imbalanced across the strips AND related to yield");
fitImg(s,FIG("fig14_screen.png"),9/5.6,0.2,1.55,8.4,5.0);
s.addText([
 {text:"Only RSP & ApDepth qualify",options:{bold:true,color:DARK,fontSize:18,breakLine:true}},
 {text:"— imbalanced AND yield-related (top-right).",options:{fontSize:15,color:INK,breakLine:true}},
 {text:"\nOM is imbalanced but doesn’t track yield here; Clay tracks yield but is balanced; LS and Slope are neither.",options:{fontSize:15,color:INK,breakLine:true}},
 {text:"\nSo the bias to remove is depth & position (ApDepth, RSP) — not organic matter, and not the slope gradient.",options:{fontSize:15,bold:true,color:DARK}}
],{x:8.9,y:2.0,w:4.0,h:4.4,fontFace:BF,fontSize:15,color:INK,lineSpacingMultiple:1.25,valign:"top"});

/* 4 — THE IDEA */
s=p.addSlide(); s.background={color:PAPER};
title(s,"The fix: compare like with like","Coarsened Exact Matching (CEM)");
s.addText([
 {text:"Don’t compare a treated plot to ",options:{}},{text:"whatever",options:{italic:true}},{text:" control sits across the line — only to controls on the ",options:{}},{text:"same kind of ground",options:{bold:true,color:DARK}},{text:".",options:{breakLine:true}},
 {text:"\n“Coarsen” means: cut each variable (RSP, ApDepth, …) into a few ranges — here, thirds (low / mid / high). A plot’s",options:{}},{text:"cell",options:{bold:true,color:DARK}},{text:" is its combination of those ranges across the variables (e.g. shallow-ApDepth · high-RSP).",options:{breakLine:true}},
 {text:"\nTreated and control plots that land in the same cell are alike, so their yield gap is the N effect. Cells holding only one side are set aside.",options:{}}
],{x:0.6,y:1.95,w:7.0,h:4.1,fontFace:BF,fontSize:17,color:INK,lineSpacingMultiple:1.3,valign:"top"});
// right: a visual of coarsening one variable into bins
s.addShape(p.ShapeType.roundRect,{x:8.0,y:2.1,w:4.7,h:3.95,rectRadius:0.12,fill:{color:DARK},line:{type:"none"}});
s.addText("Coarsen, then match",{x:8.0,y:2.3,w:4.7,h:0.5,fontFace:HF,fontSize:20,bold:true,color:PAPER,align:"center"});
const binshade=["3B6B3C","6FA046","A8C68A"], binlab=["low","mid","high"];
[["RSP",3.2],["ApDepth",3.95]].forEach(r=>{ const y=r[1];
 s.addText(r[0],{x:8.2,y,w:1.0,h:0.4,fontFace:BF,fontSize:14,bold:true,color:"CADCBF",valign:"middle"});
 for(let b=0;b<3;b++){ const x=9.3+b*1.12;
  s.addShape(p.ShapeType.rect,{x,y:y+0.02,w:1.04,h:0.36,fill:{color:binshade[b]},line:{color:DARK,width:1}});
  s.addText(binlab[b],{x,y:y+0.02,w:1.04,h:0.36,fontFace:BF,fontSize:11,bold:true,color:b===2?INK:"FFFFFF",align:"center",valign:"middle"});}});
s.addText("A plot’s cell = its bins combined.\nMatch treated ↔ control inside a cell.",
 {x:8.2,y:4.75,w:4.3,h:1.1,fontFace:BF,fontSize:13,italic:true,color:"CADCBF",align:"center",lineSpacingMultiple:1.25,valign:"top"});

/* 5 — HOW CEM WORKS */
s=p.addSlide(); s.background={color:PAPER};
title(s,"How CEM works — five steps");
const steps=[["1","Coarsen","Slice each soil/terrain variable into a few bins"],
 ["2","Stratify","A plot’s stratum = its combination of bins"],
 ["3","Keep","Keep strata that hold both a treated and a control plot"],
 ["4","Drop","Set aside strata missing one side — no fair comparison"],
 ["5","Estimate","Average the within-stratum treated−control differences"]];
steps.forEach((st,i)=>{ const x=0.55+i*2.5;
 s.addShape(p.ShapeType.ellipse,{x:x+0.75,y:2.3,w:0.95,h:0.95,fill:{color:i===3?ORANGE:GREEN},line:{type:"none"}});
 s.addText(st[0],{x:x+0.75,y:2.3,w:0.95,h:0.95,fontFace:HF,fontSize:30,bold:true,color:PAPER,align:"center",valign:"middle"});
 s.addText(st[1],{x:x,y:3.45,w:2.45,h:0.5,fontFace:HF,fontSize:19,bold:true,color:DARK,align:"center"});
 s.addText(st[2],{x:x+0.1,y:3.95,w:2.25,h:1.4,fontFace:BF,fontSize:14,color:INK,align:"center",lineSpacingMultiple:1.15,valign:"top"});
 if(i<4) s.addShape(p.ShapeType.line,{x:x+1.75,y:2.77,w:0.7,h:0,line:{color:MOSS,width:2,endArrowType:"triangle"}});});
s.addText("Dropping plots is a feature, not a flaw — it refuses to compare apples to oranges.",
 {x:0.6,y:5.9,w:12.1,h:0.7,fontFace:BF,fontSize:17,italic:true,color:DARK,align:"center"});

/* 5b — BINNING WORKED EXAMPLE */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Binning, worked: a 3 x 3 grid of subclasses","The two confounders (RSP, ApDepth), 3 bins each — on the real 48 plots");
fitImg(s,FIG("fig5_binning.png"),9.7/7,0.2,1.55,8.1,5.6);
s.addText([
 {text:"Bin",options:{bold:true,color:DARK}},{text:" = a slice of one confounder (3 on RSP, 3 on ApDepth).",options:{breakLine:true}},
 {text:"\nSubclass",options:{bold:true,color:DARK}},{text:" = a cell = one RSP bin x one ApDepth bin → up to 3x3 = 9.",options:{breakLine:true}},
 {text:"\nA cell is kept only if it holds ",options:{}},{text:"both",options:{bold:true,color:DARK}},{text:" a treated and a control plot.",options:{breakLine:true}},
 {text:"\nHere 8 of 9 cells hold both (green) and 1 is empty — on the two confounders the strips overlap, so the coarse grid keeps essentially all plots. Dropping kicks in once you add variables or finer bins (ahead).",options:{}}
],{x:8.5,y:2.0,w:4.6,h:4.7,fontFace:BF,fontSize:17,color:INK,lineSpacingMultiple:1.3,valign:"top"});

/* 5b — BIN MEMBERSHIP MAP (the grid, mapped to the field) */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Bin membership, on the ground","The same RSP × ApDepth subclasses from the grid — mapped onto the real plots");
fitImg(s,FIG("fig15_binmap.png"),12/2.2,0.3,1.75,12.73,2.5);
s.addText([
 {text:"Along the field, ",options:{}},{text:"RSP runs high → low (orange → green → blue)",options:{bold:true,color:DARK}},{text:" — a terrain gradient both strips ride together.",options:{breakLine:true}},
 {text:"\nWhere the two bands share a ",options:{}},{text:"colour",options:{bold:true}},{text:", the high-N and control plots at that spot fall in the same subclass — directly comparable. Where the ",options:{}},{text:"shade differs",options:{bold:true,color:ORANGE}},{text:" (ApDepth), the strips sit on different ground at the same RSP — the confounding matching corrects.",options:{}}
],{x:0.8,y:4.55,w:11.7,h:2.2,fontFace:BF,fontSize:17,color:INK,align:"center",lineSpacingMultiple:1.3,valign:"top"});

/* 5c — MATCHING ≈ RANDOMIZATION (the thesis) */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Matching ≈ randomization for an unreplicated trial","It recovers the covariate balance you would have gotten by randomizing — after the fact");
fitImg(s,FIG("fig7_loveplot.png"),1092/520,0.15,1.7,6.45,3.4);
fitImg(s,FIG("fig8_ecdf.png"),9/4.5,6.7,1.7,6.45,3.4);
s.addText([
 {text:"An unreplicated strip trial can’t randomize — so the strips sit on different ground (RSP and ApDepth imbalanced). ",options:{}},
 {text:"Matching pulls them toward balance",options:{bold:true,color:DARK}},
 {text:" — CEM (the headline method) balances ApDepth and nudges RSP in; full matching tightens both SMDs under 0.1 (left). After matching the ApDepth distributions overlap (right) — the like-for-like comparison the design never built, and randomization’s main benefit recovered without replication.",options:{}}
],{x:0.7,y:5.35,w:11.9,h:1.5,fontFace:BF,fontSize:15,color:INK,align:"center",lineSpacingMultiple:1.25});

/* 10b — CDF DETAIL: RSP + ApDepth */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Covariate balance in full — coarse 3-bin CEM","Empirical CDFs before vs after CEM matching on the two confounders RSP + ApDepth (3 bins)");
fitImg(s,FIG("fig16_cdf_rspapdepth.png"),1300/780,0.4,1.65,12.5,4.95);
s.addText([
 {text:"Top row = before matching (curves apart = imbalance); bottom = after coarse CEM. ",options:{}},
 {text:"ApDepth’s curves close up; RSP’s only partly",options:{bold:true,color:DARK}},
 {text:" — 3 coarse bins can’t satisfy both confounders at once (46 of 48 kept). The next slides tighten the bins, then switch engines.",options:{}}
],{x:0.6,y:6.75,w:12.1,h:0.6,fontFace:BF,fontSize:13,italic:true,color:MUTE,align:"center"});

/* 10b2 — CDF DETAIL: RSP + ApDepth, finer 5-bin CEM */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Covariate balance in full — finer 5-bin CEM","Same two confounders, tighter bins: RSP + ApDepth (5 bins)");
fitImg(s,FIG("fig17_cdf_rspapdepthls.png"),1300/780,0.4,1.65,12.5,4.95);
s.addText([
 {text:"Five bins close ",options:{}},
 {text:"both",options:{bold:true,color:DARK}},
 {text:" RSP and ApDepth — but at a price: only 31 of 48 plots survive (17 dropped). Tighter bins = better balance on a smaller, more-selected sample; the N-response answer is unchanged.",options:{}}
],{x:0.6,y:6.75,w:12.1,h:0.6,fontFace:BF,fontSize:13,italic:true,color:MUTE,align:"center"});

/* 10c — CDF DETAIL: full matching (the clean balance) */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Covariate balance in full — full matching","The other engine: full matching on RSP + ApDepth, keeping every plot");
fitImg(s,FIG("fig18_cdf_fullmatch.png"),1300/780,0.4,1.65,12.5,4.95);
s.addText([
 {text:"Full matching reweights controls instead of dropping them, so it closes ",options:{}},
 {text:"both",options:{bold:true,color:DARK}},
 {text:" RSP and ApDepth (SMDs ~0.01) while keeping all 48 plots. The imbalance was modest to begin with — this is a clean tightening, not a dramatic rescue.",options:{}}
],{x:0.6,y:6.75,w:12.1,h:0.6,fontFace:BF,fontSize:13,italic:true,color:MUTE,align:"center"});

/* 5c — RANDOMIZATION ASSUMES, MATCHING SHOWS */
s=p.addSlide(); s.background={color:PAPER};
title(s,"“Doesn’t randomization already handle this?”","In one small trial, randomization balances in expectation — but verifies nothing");
const rzY=1.95, rzH=3.6;
s.addShape(p.ShapeType.roundRect,{x:0.7,y:rzY,w:5.8,h:rzH,rectRadius:0.08,fill:{color:"F4F1EC"},line:{color:"C9B79C",width:1.5}});
s.addText("Randomly placed mini-strips",{x:0.95,y:rzY+0.2,w:5.3,h:0.5,fontFace:BF,fontSize:18,bold:true,color:"8A6D3B"});
s.addText([
 {text:"Balance confounders ",options:{}},{text:"in expectation",options:{bold:true}},{text:" — averaged over many re-randomizations.",options:{breakLine:true}},
 {text:"\nIn one small layout the actual draw can still land high-N on shallower soil — and the design ",options:{}},{text:"never checks",options:{bold:true,color:ORANGE}},{text:".",options:{breakLine:true}},
 {text:"\nReal edge: it also balances what you ",options:{}},{text:"did not measure",options:{bold:true}},{text:" (in expectation).",options:{}}
],{x:0.95,y:rzY+0.8,w:5.35,h:rzH-1.0,fontFace:BF,fontSize:15,color:INK,lineSpacingMultiple:1.2,valign:"top"});
s.addShape(p.ShapeType.roundRect,{x:6.8,y:rzY,w:5.8,h:rzH,rectRadius:0.08,fill:{color:"E7F1EA"},line:{color:DARK,width:1.5}});
s.addText("Matching, after the fact",{x:7.05,y:rzY+0.2,w:5.3,h:0.5,fontFace:BF,fontSize:18,bold:true,color:DARK});
s.addText([
 {text:"Balances the ",options:{}},{text:"measured",options:{bold:true}},{text:" ground (RSP, ApDepth) and ",options:{}},{text:"shows you it did",options:{bold:true,color:DARK}},{text:" — SMDs, CDF overlap.",options:{breakLine:true}},
 {text:"\nCompares only like-for-like; drops plots with no counterpart instead of extrapolating.",options:{breakLine:true}},
 {text:"\nThe unmeasured part is ",options:{}},{text:"bounded, not assumed",options:{bold:true}},{text:" — the E-value / Γ stress test.",options:{}}
],{x:7.05,y:rzY+0.8,w:5.35,h:rzH-1.0,fontFace:BF,fontSize:15,color:INK,lineSpacingMultiple:1.2,valign:"top"});
s.addText([{text:"Randomization balances everything and verifies nothing; matching verifies what you measured and stress-tests the rest. ",options:{bold:true,color:DARK}},{text:"For one unreplicated trial, that is the one you can actually inspect.",options:{italic:true,color:INK}}],{x:0.7,y:5.85,w:11.9,h:1.0,fontFace:BF,fontSize:16,align:"center",lineSpacingMultiple:1.2});

/* 6 — DROPPING IS HONEST (3D) */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Which plots get dropped, on the real field","CEM on RSP + ApDepth at 5 bins — trial plots extruded on the LiDAR DEM");
fitImg(s,FIG("fig3_dem3d.png"),1150/860,0.3,1.6,8.6,5.7);
s.addText([{text:"Red plots have no control on matching ground, so they drop — here 17 of 48, matching the two real confounders (RSP+ApDepth) at 5 bins. ",options:{bold:true,color:DARK}},{text:"Coarse 3-bin matching keeps almost every plot but leaves RSP loose; tighten the bins and real dropping kicks in — matching refuses to compare apples to oranges, which is why the kept comparison is fair.",options:{}}],
 {x:9.1,y:2.0,w:3.9,h:3.5,fontFace:BF,fontSize:17,color:INK,lineSpacingMultiple:1.3,valign:"top"});

/* 7 — FIVE MODELS, ONE ANSWER */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Two models, one answer","RSP alone vs both confounders — CEM at 3 & 5 bins and PSM full matching");
fitImg(s,FIG("fig11_models5.png"),1170/572,0.15,1.7,7.9,5.2);
s.addText([
 {text:"~42–47 bu/ac",options:{fontSize:28,bold:true,color:DARK,breakLine:true}},
 {text:"both models, three estimators",options:{fontSize:15,color:MUTE,breakLine:true}},
 {text:"\nRSP alone and RSP+ApDepth — under CEM at 3 bins, CEM at 5 bins, and PSM full matching — all land in the same window as the naive 43.",options:{fontSize:16,color:INK,breakLine:true}},
 {text:"\nThe delta yield doesn’t hinge on which confounders you adjust for or how finely you bin. That robustness is what makes it safe to feed a delta-yield / economics calculator.",options:{fontSize:16,bold:true,color:DARK}}
],{x:8.25,y:2.0,w:4.7,h:4.7,fontFace:BF,lineSpacingMultiple:1.25,valign:"top"});

/* 7b — BINNING SENSITIVITY */
s=p.addSlide(); s.background={color:PAPER};
title(s,"More bins & more confounders shrink the sample","The estimate holds — but the population it applies to does not");
fitImg(s,FIG("fig10_binsens.png"),1430/598,0.3,1.95,12.7,3.7);
s.addText([
 {text:"Point SIZE = % of plots retained. The response holds ~42–47, but retention falls as bins rise and confounders are added — finer bins drop more plots. ",options:{}},
 {text:"The dropped ground isn’t low-yielding — its N response is simply unmeasured: no comparable control exists there, and matching refuses to guess.",options:{bold:true,color:DARK}},
 {text:" So the finer-bin number (e.g. 47.5) is the response for the well-matched subset, not the whole field. Watch retention / ESS, not just the point.",options:{}}
],{x:0.8,y:5.9,w:11.7,h:1.35,fontFace:BF,fontSize:14,color:INK,align:"center",lineSpacingMultiple:1.2});

/* 7c — CROSS-STRIP GRADIENT = APDEPTH MADE SPATIAL */
s=p.addSlide(); s.background={color:PAPER};
title(s,"The cross-strip gradient is ApDepth, made spatial","Where high-N sits upslope, it's on shallower soil than its control");
fitImg(s,FIG("fig12_gradientmap.png"),7.6/6.6,0.15,1.6,5.85,4.7);
fitImg(s,FIG("fig13_gradient_delta.png"),8.8/5.4,6.1,1.75,7.0,4.2);
s.addText("The strips aren't iso-elevation: 17 of 24 high-N plots sit upslope of their control, where high-N is on shallower soil (ApDepth ~0.4 below its control) → lower yield potential → raw delta 42 vs 45 where the strips are level. ApDepth difference correlates +0.52 with delta and explains much of it (the elevation–delta link drops from −0.40 to −0.25 once ApDepth is partialled out). Matching on ApDepth corrects it — confounding, not N transport.",
 {x:0.6,y:6.45,w:12.1,h:0.95,fontFace:BF,fontSize:13,color:INK,align:"center",lineSpacingMultiple:1.2});

/* 8 — TRANSECTS */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Walking the field: per-plot response","Raw neighbour comparison vs like-for-like matched comparison");
fitImg(s,FIG("fig2_transects.png"),9.5/7,0.4,1.55,12.5,5.2);
s.addText("At the resolution needed to match fairly on soil depth (5-bin), 12 of 24 high-N plots have no like-for-like control and drop (red ×) — including the trial’s most-eroded plot, which coarse bins instead mis-matched to deeper ground and showed as a fake negative. Every kept delta is positive.",
 {x:0.6,y:6.8,w:12.1,h:0.6,fontFace:BF,fontSize:13,italic:true,color:MUTE,align:"center"});

/* 8b — TRANSECT UNDER PSM */
s=p.addSlide(); s.background={color:PAPER};
title(s,"The same walk under the other engine (PSM)","Full matching never drops — so it can’t refuse the eroded plot’s bad comparison");
fitImg(s,FIG("fig19_transects_psm.png"),9.5/7,0.4,1.55,12.5,5.2);
s.addText([
 {text:"PSM full matching reweights controls instead of dropping, so the most-eroded plot (~208 m) is ",options:{}},
 {text:"forced onto a deeper control → a spurious −14.6",options:{bold:true,color:ORANGE}},
 {text:" — the same fake negative CEM avoided by dropping it. The negative is the bad comparison, not a real N response.",options:{}}
],{x:0.6,y:6.85,w:12.1,h:0.55,fontFace:BF,fontSize:13,italic:true,color:MUTE,align:"center"});

/* 9 — AGRONOMY */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Why ApDepth confounds — and makes agronomic sense");
const rows=[["Soil depth",GREEN,"Deeper soil (higher ApDepth) → more water & rooting → higher yield potential.","Where the high-N strip sits upslope it's on shallower soil, so its plots have lower potential — a handicap that drags the raw gap down until you match on ApDepth."],
 ["Why it confounds",ORANGE,"The strips differ in ApDepth, and ApDepth drives yield (partial cor +0.46).","Imbalanced AND yield-related = a real confounder. Matching on ApDepth removes the depth handicap, not the N effect."]];
rows.forEach((r,i)=>{ const y=2.1+i*2.3;
 s.addShape(p.ShapeType.roundRect,{x:0.6,y,w:12.1,h:2.0,rectRadius:0.1,fill:{color:"F4F7F2"},line:{color:LINE,width:1}});
 s.addShape(p.ShapeType.ellipse,{x:0.95,y:y+0.55,w:0.9,h:0.9,fill:{color:r[1]},line:{type:"none"}});
 s.addText(i===0?"▲":"▼",{x:0.95,y:y+0.55,w:0.9,h:0.9,fontFace:BF,fontSize:24,bold:true,color:PAPER,align:"center",valign:"middle"});
 s.addText(r[0],{x:2.1,y:y+0.25,w:4.2,h:0.6,fontFace:HF,fontSize:21,bold:true,color:DARK});
 s.addText(r[2],{x:2.1,y:y+0.85,w:4.4,h:0.9,fontFace:BF,fontSize:15,color:INK,lineSpacingMultiple:1.1});
 s.addText(r[3],{x:6.7,y:y+0.35,w:5.7,h:1.3,fontFace:BF,fontSize:16,color:INK,valign:"middle",lineSpacingMultiple:1.2});});
s.addText("Where the high-N strip sits on shallower (lower-ApDepth) ground, the naive ~43 understates the response; matching on ApDepth nudges it to ~44–46 bu/ac.",
 {x:0.6,y:6.85,w:12.1,h:0.5,fontFace:BF,fontSize:16,italic:true,color:DARK,align:"center"});

/* (sensitivity-to-hidden-bias slide moved to the Backup section below) */

/* 10 — TAKEAWAYS */
s=p.addSlide(); s.background={color:DARK};
s.addText("So — would you change your N program now?",{x:0.85,y:0.55,w:11.8,h:0.9,fontFace:HF,fontSize:31,bold:true,color:PAPER});
s.addText([{text:"Yes — because the +43 is no longer one strip against one strip. ",options:{bold:true,color:"CADCBF"}},{text:"It survived being turned into a fair, like-for-like comparison:",options:{color:"CADCBF"}}],{x:0.87,y:1.35,w:11.8,h:0.6,fontFace:BF,fontSize:16,italic:true});
const tk=["Unreplicated strips aren’t randomized — they sit on different ground, so the raw gap is confounded.",
 "Matching recovers the covariate balance randomization would have given — randomization’s main benefit, without replication.",
 "The delta yield (~42–47 bu/ac) holds across both confounder models, every bin count, both engines, and a hidden-bias stress test — trustworthy to feed a delta-yield / economics calculator.",
 "A negative N response is counter-theoretical — a red flag for a plot with no comparable control, not a real effect. Honest matching drops it rather than forcing a misleading comparison (as PSM and coarse bins do)."];
tk.forEach((t,i)=>{ const y=2.2+i*1.15;
 s.addShape(p.ShapeType.ellipse,{x:0.9,y:y,w:0.55,h:0.55,fill:{color:GREEN},line:{type:"none"}});
 s.addText((i+1).toString(),{x:0.9,y:y,w:0.55,h:0.55,fontFace:HF,fontSize:20,bold:true,color:PAPER,align:"center",valign:"middle"});
 s.addText(t,{x:1.75,y:y-0.1,w:10.7,h:0.9,fontFace:BF,fontSize:18,color:"E7EFE3",valign:"middle",lineSpacingMultiple:1.1});});

/* BACKUP — DIVIDER */
s=p.addSlide(); s.background={color:DARK};
s.addText("Backup — methods detail",{x:0.85,y:3.0,w:11.6,h:1.0,fontFace:HF,fontSize:36,bold:true,color:PAPER,align:"center"});
s.addText("Sensitivity to hidden bias (E-value, Rosenbaum Γ) and supporting technical detail.",{x:0.85,y:4.1,w:11.6,h:0.6,fontFace:BF,fontSize:16,italic:true,color:"CADCBF",align:"center"});

/* BACKUP — SENSITIVITY TO HIDDEN BIAS */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Backup: how robust is it to hidden bias?","Per model & engine. RSP+ApDepth (bold) targets the actual confounders; naive response = 43.");
fitImg(s,FIG("fig6_sensitivity.png"),1320/560,0.3,1.95,7.7,3.5);
s.addText([
 {text:"E-value ≈ 8",options:{bold:true,color:DARK,fontSize:19,breakLine:true}},
 {text:"a hidden confounder would need a risk-ratio link of ≈8 with BOTH the N rate and yield to explain the response away — implausibly strong.",options:{fontSize:14,color:INK,breakLine:true}},
 {text:"\nCritical Γ ≈ 6",options:{bold:true,color:DARK,fontSize:19,breakLine:true}},
 {text:"even the least-robust model (RSP alone) needs roughly a 6× shift (Γ ≈ 6.2) in the odds of high-N assignment to lose significance. CEM's small exact-matched strata push Γ much higher (to ~20), but that's only ~8 strata — read it directionally, not literally.",options:{fontSize:14,color:INK,breakLine:true}},
 {text:"\nPSM, CEM, and plain covariate-adjusted regression (no matching, 44 [36–53]) all land at ~43–45 bu/ac.",options:{fontSize:15,bold:true,color:DARK}}
],{x:8.15,y:1.7,w:4.95,h:4.8,fontFace:BF,lineSpacingMultiple:1.12,valign:"top"});
s.addText("Caveat: E-value and Γ both rise with effect size and assume the matched comparison is otherwise valid — they bound hidden bias, not the two-strip design itself. A stress test, not proof of randomization.",
 {x:0.6,y:6.5,w:12.1,h:0.7,fontFace:BF,fontSize:13,italic:true,color:MUTE,align:"center"});

/* BACKUP — SPATIAL AUTOCORRELATION */
s=p.addSlide(); s.background={color:PAPER};
title(s,"Spatial autocorrelation — the cause, not a proxy","Condition on the measured drivers (RSP, ApDepth), not a distance stand-in for them");
const saY=1.95, saH=3.6;
s.addShape(p.ShapeType.roundRect,{x:0.7,y:saY,w:5.8,h:saH,rectRadius:0.08,fill:{color:"F4F1EC"},line:{color:"C9B79C",width:1.5}});
s.addText("SAR / spatial-error model",{x:0.95,y:saY+0.2,w:5.3,h:0.5,fontFace:BF,fontSize:18,bold:true,color:"8A6D3B"});
s.addText([
 {text:"Soaks up spatial structure through a ",options:{}},{text:"distance / adjacency weight matrix",options:{bold:true}},{text:" — neighbours predict neighbours.",options:{breakLine:true}},
 {text:"\nBut distance is a ",options:{}},{text:"proxy",options:{bold:true,color:ORANGE}},{text:" for the real drivers — soil depth, slope position — the very properties you already measured.",options:{breakLine:true}},
 {text:"\nAnd at n = 48 the SAR-probit is numerically unstable — its slope coefficients come back degenerate (NaN), so it can't be relied on for inference here.",options:{italic:true}}
],{x:0.95,y:saY+0.8,w:5.35,h:saH-1.0,fontFace:BF,fontSize:15,color:INK,lineSpacingMultiple:1.2,valign:"top"});
s.addShape(p.ShapeType.roundRect,{x:6.8,y:saY,w:5.8,h:saH,rectRadius:0.08,fill:{color:"E7F1EA"},line:{color:DARK,width:1.5}});
s.addText("Matching on the measured ground",{x:7.05,y:saY+0.2,w:5.3,h:0.5,fontFace:BF,fontSize:18,bold:true,color:DARK});
s.addText([
 {text:"Conditions ",options:{}},{text:"directly on the drivers themselves",options:{bold:true}},{text:" (RSP, ApDepth) — the thing, not a stand-in.",options:{breakLine:true}},
 {text:"\nThe CI is ",options:{}},{text:"clustered by matched set",options:{bold:true,color:DARK}},{text:" — the within-subclass correlation matching induces is carried into the standard error.",options:{breakLine:true}},
 {text:"\nMost spatial autocorrelation in yield IS the terrain; remove it and little structured signal is left.",options:{}}
],{x:7.05,y:saY+0.8,w:5.35,h:saH-1.0,fontFace:BF,fontSize:15,color:INK,lineSpacingMultiple:1.2,valign:"top"});
s.addText([{text:"Honest boundary: ",options:{bold:true,color:DARK}},{text:"spatial structure beyond the measured covariates is not modelled explicitly — that residual is exactly the unmeasured bias the E-value / Γ stress test bounds.",options:{italic:true,color:INK}}],{x:0.7,y:5.85,w:11.9,h:0.9,fontFace:BF,fontSize:15,align:"center",lineSpacingMultiple:1.2});

p.writeFile({fileName:DIR+"CEM_for_Agronomists.pptx"}).then(f=>console.log("WROTE",f));
