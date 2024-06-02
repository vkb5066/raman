#!/usr/bin/python3.6

#Current version: 4.0

HOW_TO_USE = \
"""
use:
./script.py -setup -io </path/to/orig_outcar> -m <modes> -d <displacement>
    --> writes "ramansetup" see below.
        <modes> is a space-seperated list of modes to calculate
        <displacement> is the offset, in angstroms, for finite differences
./script.py -gen -s </path/to/setup_file>
    --> reads ramansetup and generates displaced poscars to run vasp on
./script.py -eval -s </path/to/setup_file> -io <outcars> -t <temperature>
            -o </path/to/output_file>
    --> reads ramansetup and parses outcars to generate ramen spectra

Make sure NWRITE = 3, LEPSILON = .TRUE. for the VASP runs
IBRION = (5, 8) is only necessary for the initial mode calculation
"""


#Citations:
#very similar python code (similar to this, has no external deps!):
#   A. Fonari, S. Stauffer
#   "raman-sc"
#   https://github.com/raman-sc/VASP
#theoretical background, useful equations:
#   M. Bagheri, H. Komsa
#   "High-throughput computation of Raman spectra from first principles"
#    https://doi.org/10.1038/s41597-023-01988-5
#theoretical background, useful equations:
#   D. Porezag, M. Pederson
#   "I.R. intensities and Raman-scattering activities within DFT"
#   https://doi.org/10.1103/PhysRevB.54.7830


#some improvements (?) between this and raman-sc:
#   1.  this is written in python3
#   2.  this uses command line arguments instead of environment parameters
#   3.  this writes a ramansetup file including useful info (normal vectors of
#       atom movements not scaled by their masses, displacements for each
#       generated POSCAR, etc.)
#   4.  this has setup-generate-evaluate phases instead of a serial "all in one"
#       run mode.  This allows for more efficient jobs (100 jobs all on one core
#       each is much better than 1 job, 100 times, on 100 cores each)
#   5.  this allows for the use of IBRION=6 in the phonon portion
#   6.  this fixes a bug in the original where use of a universal scale
#       factor != 1 would make incorrect POSCARs ... this cost me over a day of
#       debugging to figure out :)
#   7.  this calculates the raman intensities in addition to activities
#now, some caveats:
#   1.  this uses OUTCAR for everything.  OUTCAR uses less precision than
#       POSCAR, so you may still want to use reman-sc for *very* accurate
#       results (or recompile VASP with some edits lol)
#   2.  this has much less error checking!  maybe I'll fix this later ;)


#Changelog:
#4.0:
# eval loop no longer crashes with missing or errored OUTCAR files
#3.0:
# changed the 'read the documentation' message to an actual help message
#2.0:
# removes the possibility for only one sample point, as using it an awful idea -
# maybe I can figure out how to improve it later, but even for simple structures
# like CdTe, the actual atom movement along the normal modes are non-standard
# enough that any approximation lower than 2-pt is awful
#1.0:
# initial code


from sys import argv

DEBUG = 0

RAMAN_SETUP_PATH = "ramansetup"
RAMAN_OUTPUT_PATH = "ramanout"

# --- helper functions ---------------------------------------------------------
def fmt(f, p=15, ff=''):
    return f'{f:{ff}.{p}f}'

def fmti(i, p=5):
    si = str(i)
    return '0'*(p-len(si)) + si

#printf character codes are insane -- I'll just do this one myself.
#fit a float to length l
def fmto(flt, l=5):
    return f'{flt: .{l}e}'

def getvol(A):
    #vol = a dot (b cross c)
    bxc = [+(A[1][1]*A[2][2] - A[1][2]*A[2][1]),
           -(A[1][0]*A[2][2] - A[1][2]*A[2][0]),
           +(A[1][0]*A[2][1] - A[1][1]*A[2][0])]
    return A[0][0]*bxc[0] + A[0][1]*bxc[1] + A[0][2]*bxc[2]

#get cartesian coordinates from direct vector x
#assumes A is stored in the usual VASP way - lattice vectors are the ROWS of A
# !!! unused !!!
#def dirtocar(A, x):
#    c = [None, None, None]
#    c[0] = A[0][0]*x[0] + A[1][0]*x[1] + A[2][0]*x[2]
#    c[1] = A[0][1]*x[0] + A[1][1]*x[1] + A[2][1]*x[2]
#    c[2] = A[0][2]*x[0] + A[1][2]*x[1] + A[2][2]*x[2]
#    return c

#TODO: very slow!
#TODO: a lot of calls to this lazily start from the beginning, but OUTCAR has
#      a format that should allow me to avoid this.  I'll change it if it
#      becomes a problem
#find first instance of 's' from file pointer 'fp'
#returns the actual string containing s and leavs fp pointing to the subsequent
#line
def ffi(s, fp, frombeginning=0):
    if(frombeginning):
        fp.seek(0)
    while(1):
        line = fp.readline()
        if(not line):
            break
        if s in line:
            return line
    return None


#returns a the comment string, a matrix of lattice vectors, atom counts, and
#atom positions (in cartesian coords) from OUTCAR
def posfromout(fp):
    comment, A, counts, positions = "", [[None, None, None] for i in range(3)],\
                                    [], []
    #get the comment
    s = ffi("POSCAR:", fp, 1).split()
    for s_ in s[1:]: comment += s_
    #get the lattice matrix A (the first instance may be symmetry-reduced - we
    #don't want that)
    s = ffi("direct lattice vectors", fp, 1)
    fp.readline()
    s = ffi("direct lattice vectors", fp)
    A[0] = [float(f) for f in fp.readline().split()[:3]]   ##a1
    A[1] = [float(f) for f in fp.readline().split()[:3]]   ##a2
    A[2] = [float(f) for f in fp.readline().split()[:3]]   ##a3
    #get the number of ions of each type
    s = ffi("ions per type", fp, 1)
    counts = [int(i) for i in s.split()[4:]]
    #get the cartesian atom positions
    positions = [None for i in range(0, sum(counts))]
    s = ffi("position of ions in cartesian coordinates", fp, 1)
    for i in range(0, sum(counts)):
        positions[i] = [float(f) for f in fp.readline().split()]

    return comment, A, counts, positions

#returns a the comment string, a matrix of lattice vectors, atom counts, and
#atom positions (in cartesian coords) from OUTCAR
def posfromsetup(fp):
    comment, A, counts, positions = "", [[None, None, None] for i in range(3)],\
                                    [], []
    ffi("BASE POSCAR", fp, 1)
    comment = fp.readline()
    fp.readline() ##1.0
    for i in range(0, 3):
        A[i] = [float(f) for f in fp.readline().split()]
    counts = [int(i) for i in fp.readline().split()]
    fp.readline() ##cartesian
    for i in range(0, sum(counts)):
        positions.append([float(f) for f in fp.readline().split()])

    return comment, A, counts, positions

#returns modes_of_interest mass-scaled eigenpairs from outcar
def modesfromout(fp, moi, A):
    #get the atomic masses from OUTCAR
    s = ffi("ions per type", fp, 1)
    elemcounts = [int(i) for i in s.split()[4:]]
    nelems = len(elemcounts)
    nions = sum(elemcounts)
    masses = [None for i in range(0, nelems)]
    fp.seek(0)
    for i in range(0, nelems):
        s = ffi("POMASS", fp)
        masses[i] = float(s.split()[2][:-1]) ##[-1] to get rid of ';'
    #get an array of inverse sqrt masses to multiply by eigenvectors
    invsqrtmass = [None for i in range(0, nions)]
    counter = 0
    for i in range(0, nelems):
        for j in range(0, elemcounts[i]):
            invsqrtmass[counter] = 1/(masses[i]**(1/2))
            counter += 1

    #now, the unscaled modes (IBRION=6 doesn't write the scaled masses, so
    #I have to do this manually)
    energies = []
    vectorsunscaled = []
    vectors = []
    s = ffi("Eigenvectors and eigenvalues of the dynamical matrix", fp)
    for i in range(0, 3): fp.readline()
    while(1):
        s = fp.readline().split()
        if(len(s) < 2): break
        if(s[1] == 'f' or s[1] == "f/i="): ##this is an actual mode
            if(int(s[0]) in moi): ##and we actually want to consider it
                energies.append(float(s[-2]))
                fp.readline() ##skip "X Y Z dx dy dz"
                vecu = [[None, None, None] for i in range(0, nions)]
                vecs = [[None, None, None] for i in range(0, nions)]
                for i in range(0, nions):
                    a = [float(f) for f in fp.readline().split()[3:]]
                    for j in range(0, 3):
                        vecu[i][j] = a[j]
                        vecs[i][j] = a[j]*invsqrtmass[i]
#                    vecu[i] = dirtocar(A, vecu[i]) VASP docs say eigvecs are in
#                    vecs[i] = dirtocar(A, vecs[i]) direct coords; this is false
                vectorsunscaled.append(vecu[:])
                vectors.append(vecs[:])
                fp.readline() ##skip blank line
            else: ##we're not interested - just skip
                for i in range(0, nions+2):
                    fp.readline()
        else: ##this is the end of the mode readouts
            break

    return energies, vectorsunscaled, vectors

def dielfromout(fp):
    D = [[None, None, None] for i in range(0, 3)]
    s = ffi("MACROSCOPIC STATIC DIELECTRIC TENSOR (", fp, 1)
    fp.readline() ##skip "---------..." line
    for i in range(0, 3):
        D[i] = [float(f) for f in fp.readline().split()]
    return D

def dielfromsetup(fp):
    D = [[None, None, None] for i in range(0, 3)]
    s = ffi("BASE DIELECTRIC TENSOR", fp, 1)
    for i in range(0, 3):
        D[i] = [float(f) for f in fp.readline().split()]
    return D

#returns a printable string representing a properly formatted POSCAR, assuming
#positions are in cartesian coords
def postostr(comment, A, counts, positions):
    #comment block, explicit universal scaling factor
    s = comment[:] + "\n" + fmt(1.0) + "\n"
    #lattice
    for i in range(0, 3):
        for j in range(0, 3): s += fmt(A[i][j]) + ' '
        s += "\n"
    #atom counts, explicit cartesian
    for i in range(0, len(counts)): s += str(counts[i]) + ' '
    s +="\nCartesian\n"
    #atom positions
    for i in range(0, sum(counts)):
        for j in range(0, 3): s += fmt(positions[i][j]) + ' '
        s += "\n"

    return s

#this will make the eval function a lot less confusing
class modeinfo:
    step = None
    diels = [None, None]

    index = None
    nrg = None
    avgpold = None   #mean polarizability derivative - alpha prime term
    anisopold = None #anisotropy of polarizability derivitive - beta prime term
    depolrat = None  #depolarization ratio
    activity = None  #response with no experimental considerations

    intensity = None #intensity with experimental considerations

    def __init__(self):
        self.diels = [None, None] ##the fact that I have to do this is retarded
        self.intensity = 0.0

#logic for adding dielectric tensor to a mode info struct
#modeinfo->diels is a 2 x (3x3) tensor, where the 0th index is the "low"
#sample point (x0-h for 3-point derivative, or x0 for 2-point derivative) and
#the 1st index is the "high" sample point (x0+h for both 2- and 3-pt derivs)
#this will have to be adjusted if we want more accurate derivatives
def adddiel(modeinfo, stepsize, D):
    if(stepsize <= 0.0 + 1e-8): ##if step <= 0, this must be the low index
        modeinfo.diels[0] = [[d[0], d[1], d[2]] for d in D]
    else:                           ##otherwise this must be the high index
        modeinfo.diels[1] = [[d[0], d[1], d[2]] for d in D]


#sets all of the physical properties for a single node
def setphysprops(modeinfo, ucvol, temp=300):
    if(modeinfo.diels[0] == None or modeinfo.diels[1] == None): return
    if(DEBUG): print("idx: " + str(modeinfo.index) +\
                     " h'=" + str(modeinfo.step))

    scale = ucvol / (4*3.141592654) / (2*modeinfo.step)
    #make the raman matrix: d(xhi)/d(normal mode coord)
    R = [[None, None, None] for i in range(0, 3)]
    for i in range(0, 3):
        for j in range(0, 3):
            R[i][j] = modeinfo.diels[1][i][j] - modeinfo.diels[0][i][j]
            R[i][j] *= scale

    #physical properties: expt. independent
    modeinfo.avgpold = 1/3*(R[0][0] + R[1][1] + R[2][2])
    modeinfo.anisopold = 1/2**(1/2) * (\
                                        (R[0][0]-R[1][1])**2 +\
                                        (R[0][0]-R[2][2])**2 +\
                                        (R[1][1]-R[2][2])**2 +\
                                        6*(R[0][1]**2+R[0][2]**2+R[1][2]**2)\
                                      )**(1/2)
    a, b = modeinfo.avgpold, modeinfo.anisopold
    modeinfo.depolrat = 3*b*b / (45*a*a + 4*b*b + 1e-9)
    modeinfo.activity = 45*a*a + 7*b*b
    #physical properties: expt. dependent
    #note that I've dropped a bunch of constants and assumed incident energy
    #is greater than mode energy to avoid dealing with scattered phonon energies
    kB = 8.617333262e-2 ##meV/K
    n = 1 / (2.718281828**(modeinfo.nrg/kB/temp) - 1) ##1/(exp[E/kT] - 1)
    modeinfo.intensity = (n+1)/modeinfo.nrg * modeinfo.activity

def writephysicalprops(modedict, name, temp):
    fp = open(name, 'w')
    #basic info
    fp.write("[energy] = meV\n")
    fp.write("[alpha,beta] = A^2 / amu^1/2\n")
    fp.write("[activity] = A^4 / amu\n")
    fp.write("intensity at " + fmt(temp, 1) + " K\n")
    fp.write("index  energy       alpha       beta        depol rat   "
             "activity     intensity\n")
    fp.write('-'*80 + "\n")
    #info of interest
    for ke, va in modedict.items():
        s  = (fmti(ke,p=4) if ke!=None else"****")+"  "
        s += (fmto(va.nrg) if va.nrg!=None else"************")+" "
        s += (fmto(va.avgpold,4) if va.avgpold!=None else"***********")+" "
        s += (fmto(va.anisopold,4) if va.anisopold!=None else"***********")+" "
        s += (fmto(va.depolrat,4) if va.depolrat!=None else"***********")+" "
        s += (fmto(va.activity) if va.activity!=None else"************")+" "
        s += (fmto(va.intensity) if va.intensity!=None else"************")+"\n"

        fp.write(s)

    fp.close()

# --- main routine functions ---------------------------------------------------
def setup():
    print("setup")

    #set reasonable defaults
    po = "OUTCAR"
    d = 0.01 ##angstroms
    m = [] ##will calculate all modes

    #parse command line args
    for i in range(1, len(argv)):
        if(argv[i] == "-io"):   ##path to original outcar
            po = argv[i+1]
        elif(argv[i] == "-m"):  ##list of modes
            for j in range(i+1, len(argv)):
                if(argv[j][0] == '-'):
                    i = j-1 ##useless optimization, but I can't help myself
                    break
                m.append(int(argv[j]))
        elif(argv[i] == "-d"):  ##displacement magnitude, in angstroms
            d = float(argv[i+1])


    #parse unperturbed OUTCAR
    print(" setup(): parsing outcar")
    fp = open(po, 'r')
    comment, A, ions, positions = posfromout(fp)
    if(m == []): m = [i for i in range(1, 3*sum(ions)+1)]
    modeenergies, modeunscaledvecs, modescaledvecs = modesfromout(fp, m, A)
    fp.close()

    #get basic strings so the below writing process isn't as messy
    poss = postostr("built from " + str(po), A, ions, positions)

    #Now all important information is known - write the clas, original poscar,
    #dielectric tensor, and eigenpairs of the modes of interest to a file
    print(" setup(): writing setup file to", RAMAN_SETUP_PATH)
    fp = open(RAMAN_SETUP_PATH, 'w')
    ##command line args
    fp.write("OUTCAR_NAME N_DISPLACEMENTS MAGNITUDE_DISPLACEMENT (A)\n")
    fp.write(po + ' ' + str(2) + ' ' + fmt(d) + "\n")
    fp.write("MODE_INDICES_TO_INCLUDE (base-1, in-line with OUTCAR)")
    s = "\n" + str(m[0])
    for i in range(1, len(m)):
        s += ' ' + str(m[i])
    fp.write(s + "\n")
    ##poscar, dielectric
    fp.write("BASE POSCAR\n" + poss)
    ##eigenpairs
    fp.write("BASE MODE EIGENPAIRS OF INTEREST\n"
             "INDEX_IN_OUTCAR EIGVAL (meV) \\n "
             "EIGVEC (not mass scaled)\n")
    for i in range(0, len(modeenergies)):
        fp.write(str(m[i]) + ' ' + fmt(modeenergies[i]) + "\n")
        for j in range(0, sum(ions)):
            v = modeunscaledvecs[i][j]
            for k in range(0, 3): fp.write(fmt(v[k], ff=' ') + ' ')
            fp.write("\n")
    ##now, write the displaced poscars and step sizes
    stepsizes = [None for i in range(0, len(m))]
    for i in range(0, len(m)):
        magvec = 0.
        for j in range(0, sum(ions)):
            for k in range(0, 3): magvec += (modescaledvecs[i][j][k])**2
        magvec = magvec**(1/2)
        stepsizes[i] = d / magvec
    fp.write("DISPLACEMENT MAPPINGS\n")
    fp.write("INDEX_IN_PROGRAM INDEX_IN_OUTCAR STEPSIZE (A*sqrt(amu)) "
             "\\n DISPLACEMENTS (A)\n")
    for i in range(0, 2*len(m)):
        step = stepsizes[i//2] * (-1 + 2*(i%2))
        fp.write(fmti(i)+' '+fmti(m[i//2])+' '+fmt(step, 15, ' ')+"\n")
        for j in range(0, sum(ions)):
            if(DEBUG):
                print("mode "+str(i//2) + "  ion "+str(j), end="")
                print("   h'="+str(step) +\
                      " e'=("+str(modescaledvecs[i//2][j][0])+','+\
                              str(modescaledvecs[i//2][j][1])+','+\
                              str(modescaledvecs[i//2][j][2])+')')
            v = [step*modescaledvecs[i//2][j][k] for k in range(0, 3)]
            for k in range(0, 3): fp.write(fmt(v[k], ff=' ') + ' ')
            fp.write("\n")
    fp.close()

    return


def generate():
    print("generate")

    #set reasonable defaults
    ps = RAMAN_SETUP_PATH

    #parse command line args
    for i in range(1, len(argv)):
        if(argv[i] == "-s"): ##path to setup file
            ps = argv[i+1]

    #parse setup file and make output files at the same time
    print(" generate(): parsing", ps)
    fpi = open(ps, 'r')
    comment, A, ions, origpositions = posfromsetup(fpi)

    print(" generate(): creating displaced POSCARs")
    ffi("DISPLACEMENT MAPPINGS", fpi, 1)
    fpi.readline()
    while(1): ##this is the last line in file, so go until file end
        s = fpi.readline()
        if(not s): break

        s = s.split()
        pid, mid = int(s[0]), int(s[1])
        fpo = open("POSCAR" + fmti(pid), 'w')
        ##create displaced atom positions
        displacedpositions = [x[:] for x in origpositions]
        for j in range(0, sum(ions)):
            d = [float(f) for f in fpi.readline().split()]
            for k in range(0, 3): displacedpositions[j][k] += d[k]
        ##write our this poscar
        poss = postostr(fmti(pid), A, ions, displacedpositions)
        fpo.write(poss)
        fpo.close()
    fpi.close()

    return




def eval():
    print("evaluate")

    #set reasonable defaults
    ps = RAMAN_SETUP_PATH
    polist = [] ##no reasonable default for list of OUTCARS
    t = 300 #kelvin
    outfile = RAMAN_OUTPUT_PATH

    #parse command line args
    for i in range(1, len(argv)):
        if(argv[i] == '-s'):    ##path to setup file
            so = argv[i+1]
        elif(argv[i] == "-io"): ##paths to OUTCARS
            for j in range(i+1, len(argv)):
                if(argv[j][0] == '-'):
                    i = j-1
                    break
                polist.append(argv[j])
        elif(argv[i] == '-t'):  ##simulation temp
            t = float(argv[i+1])
        elif(argv[i] == '-o'):  ##output file name
            outfile = argv[i+1]

    #parse ramen setup
    print(" eval(): parsing", ps)
    fp = open(ps, 'r')
    ##get all (OUTCAR-indexed) modes of interest
    fp.readline()
    fp.readline()
    fp.readline()
    outdict = {int(mode): modeinfo() for mode in fp.readline().split()}
    ##and the base information
    comment, A, ions, origpositions = posfromsetup(fp)
    ucvolume = getvol(A)
    ##fill in the modeinfo's energies
    ffi("BASE MODE EIGENPAIRS OF INTEREST", fp)
    fp.readline()
    for i in range(0, len(outdict.keys())):
        line = fp.readline().split()
        outdict[int(line[0])].nrg = float(line[1])
        for j in range(0, sum(ions)): fp.readline() ##skip eigenvectors
    #finally, create a mapping from program indices to mode indices, and a
    #mapping from program indices to displacements (overkill now, easy to
    #generalize later if necessary)
    mappitomi = [] ##mappitomi[program_index] = mode_index
    mappitoss = [] ##mappitoss[program_index] = (signed) step size
    ffi("DISPLACEMENT MAPPINGS", fp)
    fp.readline()
    while(1):  ##these are the last lines in file, so go until file end
        s = fp.readline()
        if (not s): break
        s = s.split()
        modeindex = int(s[1])
        stepsize = float(s[2])
        mappitomi.append(modeindex)
        mappitoss.append(stepsize)
        outdict[modeindex].index = modeindex
        outdict[modeindex].step = abs(stepsize)
        for j in range(0, sum(ions)): fp.readline()  ##skip eigenvectors
    fp.close()

    #now, the fun part: finish populating the output dictionary with all
    #necessary data.  OUTCARS can be in any order, so need to be careful
    print(" eval(): parsing", len(polist), "OUTCAR files")
    for infilename in polist:
        fp = open(infilename, 'r')
        s = ffi("POSCAR:", fp).split() ##POSCAR: <poscar's comment line>
        if(len(s) == 2):
            programindex = int(s[1])
            modeindex = mappitomi[programindex]
            D = dielfromout(fp)
            adddiel(modeinfo=outdict[modeindex],
                    stepsize=mappitoss[programindex],
                    D=D)
        else:
            print(" eval(): incorrectly formatted OUTCAR:", infilename)
        fp.close()

    #compute physical properties, normalize the intensities
    for val in outdict.values(): setphysprops(val, ucvolume, t)
    norm = max([val.intensity for val in outdict.values()])
    for val in outdict.values(): val.intensity /= norm

    #print result
    writephysicalprops(outdict, outfile, t)

# --- main ---------------------------------------------------------------------
runmode = 0
for i in range(1, len(argv)):
    if(argv[i] == "-setup"): runmode = 1
    if(argv[i] == "-gen"):   runmode = 2
    if(argv[i] == "-eval"):  runmode = 3
if(not runmode):
    print("main(): no runmode supplied!")
    print(HOW_TO_USE)
else: [setup, generate, eval][runmode-1]()
print("main(): finished")
