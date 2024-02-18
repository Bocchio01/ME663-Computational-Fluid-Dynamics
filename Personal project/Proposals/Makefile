.PHONY: proposal clean

all: proposal clean

CC = xelatex

proposal: Proposal.tex
	$(CC) -output-directory= -include-directory= $<
	$(CC) -output-directory= -include-directory= $<

clean:
	del /Q *.aux
	del /Q *.lof
	del /Q *.log
	del /Q *.lot
	del /Q *.out
	del /Q *.toc
	del /Q *.fdb_latexmk
	del /Q *.fls
	del /Q *.gz
	del /Q *.acn
	del /Q *.glo
	del /Q *.ist
