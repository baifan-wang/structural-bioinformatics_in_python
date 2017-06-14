import numpy as np
class Atom():
    def __init__(self, atom):
        """Create Atom object.
        The Atom object stores atom serial, name, coordinate, element type
        as well as residue name and serial, chain id in which atom belong to.
        parameters, 'atom' should be a list contains the following information"""
        self._container = ''             #indicate which molecule this atom belong to, object
        self.serial =        atom[0]    #atom serial number, int
        self.name =          atom[1]    #atom name, eg. 'CA', string
        self.res_name =      atom[2]    #residues name, eg. 'ALA', sring
        self.chain_id =      atom[3]    #chain id, one letter, eg. 'A', string
        self.res_serial =    atom[4]    #residue serial number, int
        self.coord =         np.array(atom[5:8]) #x, y, z coordinates of atom, [float]
        assert self.coord.shape == (3, )
        self.x =             self.coord[0]    #x coordinate, float
        self.y =             self.coord[1]    #y coordinate, float
        self.z =             self.coord[2]    #z coordinate, float
        self.element =       atom[8]    #type of element, eg. 'C'
        self.id = self.chain_id+'.'+self.res_name+str(self.res_serial)+'.'+self.name  #identification of an atom, eg., A14.CA

    def __str__(self):
        if self._container is '':
            return('Atom object: '+self.id)
        else:
            return(str(self._container)+'.'+self.name)


    __repr__=__str__

    def __eq__(self, other):
        """
        Check whether this atom and other atom are the same
        """
        return (self._container == other._container) and (self.id == other.id) and (self.coord == other.coord)

    def distance(self, other):
        """
        distance between two point :(sum((ai - bi)**2 for ai, bi in zip(self.coord, other.coord)))**.5
        """
        assert other.coord.shape == (3, )
        d = self.coord - other.coord
        return np.sqrt(np.sum(d*d))

    def angle(self, other1, other2, other3=None):
        """
        if 4 args are provided, will calculate the angle between 
        the vector of other1->self, and 
        vector of other3->other2,
        else will compute the angle of between the point self, other1, other2.
        """
        assert other1.coord.shape == (3, )
        assert other2.coord.shape == (3, )
        v1 = other1.coord - self.coord # vector of A->B = [b.x-a.x, b.y-a.y, b.z-a.z]
        distance1 = np.sqrt(v1.dot(v1)) #the length of the vector
        v2 = other1.coord - other2.coord
        if other3 == None:
            distance2 =np.sqrt(v2.dot(v2))
            sca = np.dot(v1, v2)/(distance1*distance2)
        else:
            assert other3.coord.shape == (3, )
            v3 = other3.coord - other2.coord
            distance3 = np.sqrt(v3.dot(v3))
            sca = np.dot(v1, v3)/(distance1*distance3)
        if sca < -1.0: 
            sca = -1.0
        elif sca > 1.0: 
            sca = 1.0
        return np.arccos(sca)*180/np.pi

    def _vector_prod(self, a, b):
        """
        Compute a vector product for 3D vectors
        """
        res = np.zeros(3, 'f')
        res[0] = a[1]*b[2] - a[2]*b[1]
        res[1] = a[2]*b[0] - a[0]*b[2]
        res[2] = a[0]*b[1] - a[1]*b[0]
        return res

    def torsion(self, other2, other3, other4):
        """
        Compute the torsion angle between self, other2, other3, other4.
        Return angle in degrees.
        """
        tang=0.0
        assert other2.coord.shape == (3, )
        assert other3.coord.shape == (3, )
        assert other4.coord.shape == (3, )

        a = self.coord-other2.coord
        b = other3.coord-other2.coord
        c = self._vector_prod(a, b)

        a = other2.coord-other3.coord
        b = other4.coord-other3.coord
        d = self._vector_prod(a, b)

        dd = np.sqrt(np.sum(c*c))
        de = np.sqrt(np.sum(d*d))

        if dd<0.001 or de<0.001:
            raise ValueError ('Torsion angle undefined, degenerate points')

        vv = np.dot(c, d) / (dd*de);
        if vv<1.0: 
            tang = vv
        else: 
            tang = 1.0
        if tang < -1.0: 
            tang = -1.0
        tang = np.arccos(tang)
        tang = tang*57.296

        b = self._vector_prod(c, d)
        if np.dot(a, b) > 0.0: tang = -tang
        return tang

    def plane(self, other2, other3):
        """
        Compute the function of plane defined by self, other2, other3.
        the norm of a plane is (a,b,c)"""
        assert other2.coord.shape == (3, )
        assert other3.coord.shape == (3, )

        a = (other2.y-self.y)*(other3.z-self.z)-(other2.z-self.z)*(other3.y-self.y)
        b = (other2.z-self.z)*(other3.x-self.x)-(other2.x-self.x)*(other3.z-self.z)
        c = (other2.x-self.x)*(other3.y-self.y)-(other2.y-self.y)*(other3.x-self.x)
        d = (0-(a*self.x+b*self.y+c*self.z))
        return (a,b,c,d)

    def move(self, direction, distance):
        """move atom to the direction with distance, will change its coordinates."""
        pass

    def to_pdb(self):
        """
        Return a PDB format line.
        """
        self.character = 'ATOM'
        self.alter_local_indicater = ''
        self.code_for_insertions_of_residues = ''
        self.occupancy = 1.00
        self.temp_factor = 0.00
        self.segment_indent = ''
        self.charge = ''

        if len(self.name) <4:
            atom_name=(' '+self.name).ljust(4)
        else:
            atom_name=self.name.ljust(4)
        s = "%s%5d %s %3s %1s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s" \
                % (self.character.ljust(6) , self.serial , atom_name,  self.res_name.rjust(3) , \
                self.chain_id , self.res_serial , self.code_for_insertions_of_residues , \
                self.x , self.y , self.z , self.occupancy ,\
                self.temp_factor , self.segment_indent.ljust(4) , \
                self.element.rjust(2) , self.charge)
        return s

    def transform(self, rotation, translation):
        """
        Apply rotation and translation to the atomic coordinates.

        Example:
                >>> rotation=rotmat(pi, Vector(1, 0, 0))
                >>> translation=array((0, 0, 1), 'f')
                >>> atom.transform(rotation, translation)

        @param rot: A right multiplying rotation matrix
        @type rot: 3x3 Numeric array

        @param tran: the translation vector
        @type tran: size 3 Numeric array
        """
        self.coord = numpy.dot(self.coord, rot) + tran


class Dummy(Atom):
    """
    Creat Dummy atom class.
    Dummy atom only has a name and coordinates.
    """
    def __init__(self, name, coord):
        self.name = name             #name of dummy atom, str
        self.coord = np.array(coord) #coordinates of dummy atom, [float]
        assert self.coord.shape == (3, )
        self.x = self.coord[0]       #x,y,z coordinates of dummy atom.
        self.y = self.coord[1]
        self.z = self.coord[2]

    def __str__(self):
        return('Dummy atom object: %s' %self.name)

    def to_pdb(self, serial, res_serial, res_name='MOL', chian_id='A', element='X'):
        """res_serial
        return a PDB format line.
        extra information must be provided.
        """
        self.character = 'ATOM'
        self.alter_local_indicater = ''
        self.code_for_insertions_of_residues = ''
        self.occupancy = 1.00
        self.temp_factor = 0.00
        self.segment_indent = ''
        self.element = ''
        self.charge = ''

        if len(self.name) <4:
            atom_name=(' '+self.name).ljust(4)
        else:
            atom_name=self.name.ljust(4)
        s = "%s%5d %s %3s %1s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s" \
                % (self.character.ljust(6) , self.serial , atom_name,  self.res_name.rjust(3) , \
                self.chain_id , self.res_serial , self.code_for_insertions_of_residues , \
                self.x , self.y , self.z , self.occupancy ,\
                self.temp_factor , self.segment_indent.ljust(4) , \
                self.element.rjust(2) , self.charge)
        return s