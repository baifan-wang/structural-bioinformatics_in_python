import numpy as np

class Atom:
    def __init__(self, atom):
        '''atom should be a list contains the following information'''
        self.serial =        atom[0]    #atom serial number, must be int
        self.atom_name =     atom[1]    #atom name, string
        self.res_name =      atom[2]    #residues name(type), sring
        self.chain_id =      atom[3]    #chain id, one letter, string
        self.res_serial =    atom[4]    #residue serial number, int
        self.coord =         np.array(atom[5:8]) #x,y,z coordinates of atom
        self.x =             self.coord[0]    #x coordinate, float
        self.y =             self.coord[1]    #y coordinate, float
        self.z =             self.coord[2]    #z coordinate, float
        self.element =       atom[8]    #type of element, such as 'N', 'H', 'O'

    def __str__(self):
        return('Atom object: ' +self.chain_id+ '.'+str(self.res_serial)+'.'+self.res_name+'.'+self.atom_name)

    __repr__=__str__

    def distance(self, other):
        """
        distance between two point :(sum((ai - bi)**2 for ai, bi in zip(self.coord, other.coord)))**.5
        """
        d = self.coord - other.coord
        return np.sqrt(np.sum(d*d))

    def angle(self, other1, other2, other3=None):
        """
        if 4 args are provided, will calculate the angle between 
        the vector of other1->self, and 
        vector of other3->other2,
        else will compute the angle of between the point self, other1, other2.
        """
        v1 = other1.coord - self.coord # vector of A->B = [b.x-a.x, b.y-a.y, b.z-a.z]
        distance1 = np.sqrt(v1.dot(v1)) #the length of the vector
        v2 = other1.coord - other2.coord
        if other3 == None:
            distance2 =np.sqrt(v2.dot(v2))
            sca = np.dot(v1, v2)/(distance1*distance2)
        else:
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
        assert self.coord.shape == (3, )
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
        compute the function of plane defined by self, other2, other3.
        the norm of a plane is (a,b,c)"""
        a = (other2.y-self.y)*(other3.z-self.z)-(other2.z-self.z)*(other3.y-self.y)
        b = (other2.z-self.z)*(other3.x-self.x)-(other2.x-self.x)*(other3.z-self.z)
        c = (other2.x-self.x)*(other3.y-self.y)-(other2.y-self.y)*(other3.x-self.x)
        d = (0-(a*self.x+b*self.y+c*self.z))
        return (a,b,c,d)

    def move(self, direction, distance):
        """move atom to the direction with distance, will change its coordinates."""
        pass

    def to_pdb(self):
        """return a PDB format line"""
        self.character = 'ATOM'
        self.alter_local_indicater = ''
        self.code_for_insertions_of_residues = ''
        self.occupancy = 1.00
        self.temp_factor = 0.00
        self.segment_indent = ''
        self.element_symbol = ''
        self.charge = ''

        if len(self.atom_name) <4:
            atom_name=(' '+self.atom_name).ljust(4)
        else:
            atom_name=self.atom_name.ljust(4)
        s = "%s%5d %s %3s %1s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s" \
                % (self.character.ljust(6) , self.atom_serial , atom_name,  self.res_name.rjust(3) , \
                self.chain_id , self.res_serial , self.code_for_insertions_of_residues , \
                self.x , self.y , self.z , self.occupancy ,\
                self.temp_factor , self.segment_indent.ljust(4) , \
                self.element_symbol.rjust(2) , self.charge)
        return s


class Dummy(Atom):
    """Dummy atom class, only has a name and coordinates"""
    def __init__(self, name, coord):
        self.name = name
        self.coord = np.array(coord)
        self.x = self.coord[0]
        self.y = self.coord[1]
        self.z = self.coord[2]
    def __str__(self):
        return('Dummy atom object: %s' %self.name)