#!/bin/sh

$bin/MissStereo.sh input_0.png input_1.png

mv input_0.png input_0_original.png
mv input_1.png input_1_original.png
mv H_input_0.png                           input_0.png
mv H_input_1.png                           input_1.png


rm input_0.png_input_1.png_pairs_orsa.txt  
rm input_0.png_h.txt                       
rm input_1.png_h.txt                       
rm show_H_input_0.png                      
rm show_H_input_1.png                      
rm toto.png

