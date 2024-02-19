.PHONY: proposal ppt clean

all: proposal clean

CC = xelatex

proposal: Proposal/Proposal.tex
	$(eval DIR := Proposal)
	$(CC) -output-directory=$(DIR) -include-directory=$(DIR)/ $<
	$(CC) -output-directory=$(DIR) -include-directory=$(DIR)/ $<

ppt: Presentation/Presentation.tex
	$(eval DIR := Presentation)
	$(CC) -output-directory=$(DIR) -include-directory=$(DIR)/ $<
	$(CC) -output-directory=$(DIR) -include-directory=$(DIR)/ $<


clean:
	del /Q $(DIR)\*.aux
	del /Q $(DIR)\*.lof
	del /Q $(DIR)\*.log
	del /Q $(DIR)\*.lot
	del /Q $(DIR)\*.out
	del /Q $(DIR)\*.toc
	del /Q $(DIR)\*.fdb_latexmk
	del /Q $(DIR)\*.fls
	del /Q $(DIR)\*.gz