#!/bin/bash

make clean;
make;
./bin/postscript automateBezier.txt 5
#./bin/postscript automateCantor.txt 5
#./bin/postscript automateCantorCarre.txt 5
#./bin/postscript automateCarreSierp.txt 5
#./bin/postscript automateAutreCantorCarre.txt 5
#./bin/postscript automateAutreCantorCarreEtCarre.txt 5