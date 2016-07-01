

GET DATA
  /TYPE=TXT
  /FILE="I:\2 fMRI behaviour analysis\3 Analysis inputs\(28-Jun-2013) Learning results.txt"
  /DELCASE=LINE
  /DELIMITERS="\t"
  /ARRANGEMENT=DELIMITED
  /FIRSTCASE=2
  /IMPORTCASE=ALL
  /VARIABLES=
  Subject A7
  cRT_Overall F11.2
  cRT_summ.Sim F11.2
  cRT_summ.Dis F11.2
  cRT_summ.ContextR F10.2
  cRT_summ.ContextN F11.2
  cRT_cell.Sim_cR F10.2
  cRT_cell.Sim_cN F11.2
  cRT_cell.Dis_cR F10.2
  cRT_cell.Dis_cN F11.2
  eJ_Overall F8.2
  eJ_summ.Sim F8.2
  eJ_summ.Dis F8.2
  eJ_summ.ContextR F8.2
  eJ_summ.ContextN F8.2
  eJ_cell.Sim_cR F8.2
  eJ_cell.Sim_cN F8.2
  eJ_cell.Dis_cR F8.2
  eJ_cell.Dis_cN F8.2
  eRT_Overall F10.2
  eRT_summ.Sim F10.2
  eRT_summ.Dis F10.2
  eRT_summ.ContextR F10.2
  eRT_summ.ContextN F10.2
  eRT_cell.Sim_cR F10.2
  eRT_cell.Sim_cN F10.2
  eRT_cell.Dis_cR F10.2
  eRT_cell.Dis_cN F10.2.
CACHE.
EXECUTE.
DATASET NAME DataSet1 WINDOW=FRONT.

* Exclusion criteria

* Conditioning.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=MEANSE(cRT_cell.Sim_cR, 1) MEANSE(cRT_cell.Sim_cN, 1) 
    MEANSE(cRT_cell.Dis_cR, 1) MEANSE(cRT_cell.Dis_cN, 1) MISSING=LISTWISE REPORTMISSING=NO
    TRANSFORM=VARSTOCASES(SUMMARY="#SUMMARY" INDEX="#INDEX" LOW="#LOW" HIGH="#HIGH")
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: SUMMARY=col(source(s), name("#SUMMARY"))
  DATA: INDEX=col(source(s), name("#INDEX"), unit.category())
  DATA: LOW=col(source(s), name("#LOW"))
  DATA: HIGH=col(source(s), name("#HIGH"))
  GUIDE: axis(dim(2), label("Mean RT"))
  GUIDE: text.title(label("Conditioning-stage RT"))
  GUIDE: text.subtitle(label("                       Similar                                                              Dissimilar"))
  GUIDE: text.footnote(label("Error Bars: +/- 1 SE"))
  SCALE: cat(dim(1), include("0", "1", "1.5","2", "3"))
  SCALE: linear(dim(2), include(0))
  ELEMENT: interval(position(INDEX*SUMMARY), shape.interior(shape.square))
  ELEMENT: interval(position(region.spread.range(INDEX*(LOW+HIGH))), shape.interior(shape.ibeam))
END GPL.

GLM cRT_cell.Sim_cR cRT_cell.Sim_cN cRT_cell.Dis_cR cRT_cell.Dis_cN
  /WSFACTOR=SIMILARITY 2 Polynomial VALENCE 2 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=SIMILARITY VALENCE SIMILARITY*VALENCE.


* Encoding accuracy.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=MEANSE(eJ_cell.Sim_cR, 1) MEANSE(eJ_cell.Sim_cN, 1) 
    MEANSE(eJ_cell.Dis_cR, 1) MEANSE(eJ_cell.Dis_cN, 1) MISSING=LISTWISE REPORTMISSING=NO
    TRANSFORM=VARSTOCASES(SUMMARY="#SUMMARY" INDEX="#INDEX" LOW="#LOW" HIGH="#HIGH")
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: SUMMARY=col(source(s), name("#SUMMARY"))
  DATA: INDEX=col(source(s), name("#INDEX"), unit.category())
  DATA: LOW=col(source(s), name("#LOW"))
  DATA: HIGH=col(source(s), name("#HIGH"))
  GUIDE: axis(dim(2), label("Mean accuracy"))
  GUIDE: text.title(label("Encoding-stage accuracy"))
  GUIDE: text.subtitle(label("                       Similar                                                              Dissimilar"))
  GUIDE: text.footnote(label("Error Bars: +/- 1 SE"))
  SCALE: cat(dim(1), include("0", "1", "1.5","2", "3"))
  SCALE: linear(dim(2), include(0))
  ELEMENT: interval(position(INDEX*SUMMARY), shape.interior(shape.square))
  ELEMENT: interval(position(region.spread.range(INDEX*(LOW+HIGH))), shape.interior(shape.ibeam))
END GPL.

GLM eJ_cell.Sim_cR eJ_cell.Sim_cN eJ_cell.Dis_cR eJ_cell.Dis_cN
  /WSFACTOR=SIMILARITY 2 Polynomial VALENCE 2 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=SIMILARITY VALENCE SIMILARITY*VALENCE.

* Encoding RT.

GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=MEANSE(eRT_cell.Sim_cR, 1) MEANSE(eRT_cell.Sim_cN, 1) 
    MEANSE(eRT_cell.Dis_cR, 1) MEANSE(eRT_cell.Dis_cN, 1) MISSING=LISTWISE REPORTMISSING=NO
    TRANSFORM=VARSTOCASES(SUMMARY="#SUMMARY" INDEX="#INDEX" LOW="#LOW" HIGH="#HIGH")
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: SUMMARY=col(source(s), name("#SUMMARY"))
  DATA: INDEX=col(source(s), name("#INDEX"), unit.category())
  DATA: LOW=col(source(s), name("#LOW"))
  DATA: HIGH=col(source(s), name("#HIGH"))
  GUIDE: axis(dim(2), label("Mean RT"))
  GUIDE: text.title(label("Encoding-stage RT"))
  GUIDE: text.subtitle(label("                       Similar                                                              Dissimilar"))
  GUIDE: text.footnote(label("Error Bars: +/- 1 SE"))
  SCALE: cat(dim(1), include("0", "1", "1.5","2", "3"))
  SCALE: linear(dim(2), include(0))
  ELEMENT: interval(position(INDEX*SUMMARY), shape.interior(shape.square))
  ELEMENT: interval(position(region.spread.range(INDEX*(LOW+HIGH))), shape.interior(shape.ibeam))
END GPL.

GLM eRT_cell.Sim_cR eRT_cell.Sim_cN eRT_cell.Dis_cR eRT_cell.Dis_cN
  /WSFACTOR=SIMILARITY 2 Polynomial VALENCE 2 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=SIMILARITY VALENCE SIMILARITY*VALENCE.
