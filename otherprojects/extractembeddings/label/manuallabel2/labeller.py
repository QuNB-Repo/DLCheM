
from label.manuallabel2.utils import utils, Hlabeller, Clabeller, Nlabeller, Olabeller, Flabeller, Clabellerldacont, Nlabellerldacont, Olabellerldacont, Flabellerldacont, ClabellerLDAnotaccount, NlabellerLDAnotaccount, HlabellerLDAnotaccount, OlabellerLDAnotaccount, FlabellerLDAnotaccount, Clabellersolubility



def label(mol_filename,number_atoms,each_atom,atomic_number,each_molecule):

    #open mol file to construct connection matrix
    mol_file = open(mol_filename,mode='r')
    mol_file_read = mol_file.read()

    xyz_filename=mol_filename.replace('.mol','.xyz')

    xyz_file = open(xyz_filename,mode='r')
    xyz_file_read = xyz_file.read()

    each_atom = each_atom+1
    connections_to_atom = utils.find_connections(number_atoms,each_atom,mol_file_read)
    
    if atomic_number == 1:
        fg_key,ldalabel,gnucolor,gnumarker,hexcolor = Hlabeller.label(connections_to_atom,xyz_file_read,mol_file_read,number_atoms,each_molecule,each_atom)
    if atomic_number == 6:
        fg_key,ldalabel,gnucolor,gnumarker,hexcolor = Clabeller.label(connections_to_atom,xyz_file_read,mol_file_read,number_atoms,each_molecule,each_atom)
    if atomic_number == 7:
        fg_key,ldalabel,gnucolor,gnumarker,hexcolor = Nlabeller.label(connections_to_atom,xyz_file_read,mol_file_read,number_atoms,each_molecule,each_atom)
    if atomic_number == 8:
        fg_key,ldalabel,gnucolor,gnumarker,hexcolor = Olabeller.label(connections_to_atom,xyz_file_read,mol_file_read,number_atoms,each_molecule,each_atom)
    if atomic_number == 9:
        fg_key,ldalabel,gnucolor,gnumarker,hexcolor = Flabeller.label(connections_to_atom,xyz_file_read,mol_file_read,number_atoms,each_molecule,each_atom)

    return fg_key,ldalabel,gnucolor,gnumarker,hexcolor









        
    


