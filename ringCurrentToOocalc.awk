#!/usr/bin/awk -f


BEGIN {
    print "Atom", "Predicted", "Actual", "Basis", "PredShift"
}
{
    if (($2 ~ /ALA/) || ($2 ~ /MET/) || ($2 ~ /THR/)){
	print to_1_letter($2) $1, $6, $7, $4, $5
    }
    else {
	print to_1_letter($2) $1 "-" $3, $6, $7, $4, $5
    }
}

function to_1_letter(res,str) {
    str=res
    if       (res=="ALA") str="A"
    else if  (res~/^CYS/) str="C"
    else if  (res=="ASP") str="D"
    else if  (res=="GLU") str="E"
    else if  (res=="PHE") str="F"
    else if  (res=="GLY") str="G"
    else if  (res=="HIS") str="H"
    else if  (res=="ILE") str="I"
    else if  (res=="LYS") str="K"
    else if  (res=="LEU") str="L"
    else if  (res=="MET") str="M"
    else if  (res=="ASN") str="N"
    else if  (res=="PRO") str="P"
    else if  (res=="GLN") str="Q"
    else if  (res=="ARG") str="R"
    else if  (res=="SER") str="S"
    else if  (res=="THR") str="T"
    else if  (res=="VAL") str="V"
    else if  (res=="TRP") str="W"
    else if  (res=="TYR") str="Y"
    else if  (res=="LYS+") str="K"
    else if  (res=="ARG+") str="R"
    else if  (res=="HIS+") str="H"
    else if  (res=="ASP-") str="D"
    else if  (res=="GLU-") str="E"

    return str
}
