#!/bin/sh

#auto sincho in this example

p2c -m LB -p 5FNU_REC.pdb -l HIT.pdb -d 10

sincho -p 5FNU_REC.pdb -l HIT.pdb

pymol sincho-output/sincho.pse
