#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
matlab -nodisplay -nosplash -nodesktop -r "run('./experiments.m');exit;"
pdflatex -output-directory=diagrams diagrams/matmul_binary16_acc.tex
pdflatex -output-directory=diagrams diagrams/matmul_binary32_acc.tex
cd diagrams
latexmk -c matmul_binary16_acc.tex
latexmk -c matmul_binary32_acc.tex
cd ..
