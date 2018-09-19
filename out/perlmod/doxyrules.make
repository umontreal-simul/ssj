DOXY_EXEC_PATH = /home/florian/eclipse-workspace/ssj/src/main/docs
DOXYFILE = /home/florian/eclipse-workspace/ssj/src/main/docs/Doxyfile
DOXYDOCS_PM = /home/florian/eclipse-workspace/ssj/out/perlmod/DoxyDocs.pm
DOXYSTRUCTURE_PM = /home/florian/eclipse-workspace/ssj/out/perlmod/DoxyStructure.pm
DOXYRULES = /home/florian/eclipse-workspace/ssj/out/perlmod/doxyrules.make
DOXYLATEX_PL = /home/florian/eclipse-workspace/ssj/out/perlmod/doxylatex.pl
DOXYLATEXSTRUCTURE_PL = /home/florian/eclipse-workspace/ssj/out/perlmod/doxylatex-structure.pl
DOXYSTRUCTURE_TEX = /home/florian/eclipse-workspace/ssj/out/perlmod/doxystructure.tex
DOXYDOCS_TEX = /home/florian/eclipse-workspace/ssj/out/perlmod/doxydocs.tex
DOXYFORMAT_TEX = /home/florian/eclipse-workspace/ssj/out/perlmod/doxyformat.tex
DOXYLATEX_TEX = /home/florian/eclipse-workspace/ssj/out/perlmod/doxylatex.tex
DOXYLATEX_DVI = /home/florian/eclipse-workspace/ssj/out/perlmod/doxylatex.dvi
DOXYLATEX_PDF = /home/florian/eclipse-workspace/ssj/out/perlmod/doxylatex.pdf

.PHONY: clean-perlmod
clean-perlmod::
	rm -f $(DOXYSTRUCTURE_PM) \
	$(DOXYDOCS_PM) \
	$(DOXYLATEX_PL) \
	$(DOXYLATEXSTRUCTURE_PL) \
	$(DOXYDOCS_TEX) \
	$(DOXYSTRUCTURE_TEX) \
	$(DOXYFORMAT_TEX) \
	$(DOXYLATEX_TEX) \
	$(DOXYLATEX_PDF) \
	$(DOXYLATEX_DVI) \
	$(addprefix $(DOXYLATEX_TEX:tex=),out aux log)

$(DOXYRULES) \
$(DOXYMAKEFILE) \
$(DOXYSTRUCTURE_PM) \
$(DOXYDOCS_PM) \
$(DOXYLATEX_PL) \
$(DOXYLATEXSTRUCTURE_PL) \
$(DOXYFORMAT_TEX) \
$(DOXYLATEX_TEX): \
	$(DOXYFILE)
	cd $(DOXY_EXEC_PATH) ; doxygen "$<"

$(DOXYDOCS_TEX): \
$(DOXYLATEX_PL) \
$(DOXYDOCS_PM)
	perl -I"$(<D)" "$<" >"$@"

$(DOXYSTRUCTURE_TEX): \
$(DOXYLATEXSTRUCTURE_PL) \
$(DOXYSTRUCTURE_PM)
	perl -I"$(<D)" "$<" >"$@"

$(DOXYLATEX_PDF) \
$(DOXYLATEX_DVI): \
$(DOXYLATEX_TEX) \
$(DOXYFORMAT_TEX) \
$(DOXYSTRUCTURE_TEX) \
$(DOXYDOCS_TEX)

$(DOXYLATEX_PDF): \
$(DOXYLATEX_TEX)
	pdflatex -interaction=nonstopmode "$<"

$(DOXYLATEX_DVI): \
$(DOXYLATEX_TEX)
	latex -interaction=nonstopmode "$<"
