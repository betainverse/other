from pymol import cmd

def showmethyls():
    '''
    Show methyl groups on the protein as spheres, color coded by residue.
    Usage:
    use the run... item in the pymol File menu to select this script or
    run it from the command line:
    run /home/userName/path/toscript/showmethyls.py
    then type showmethyls on the command line.
    '''
    cmd.show('spheres','resn ALA and name CB')
    cmd.color('red','resn ALA and name CB')
    cmd.show('spheres','resn ILE and name CD1')
    cmd.color('orange','resn ILE and name CD1')
    cmd.show('spheres','resn LEU and name CD1+CD2')
    cmd.color('cyan','resn LEU and name CD1+CD2')
    cmd.show('spheres','resn MET and name CE')
    cmd.color('green','resn MET and name CE')
    cmd.show('spheres','resn THR and name CG2')
    cmd.color('blue','resn THR and name CG2')
    cmd.show('spheres','resn VAL and name CG1+CG2')
    cmd.color('purple','resn VAL and name CG1+CG2')


cmd.extend('showmethyls', showmethyls)
