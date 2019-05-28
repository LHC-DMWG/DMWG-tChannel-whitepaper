#!/bin/bash

if [ $# -ne 0 ]; then
pdflatex tChannel-whitepaper
bibtex tChannel-whitepaper
pdflatex tChannel-whitepaper
pdflatex tChannel-whitepaper
else
pdflatex tChannel-whitepaper
fi

