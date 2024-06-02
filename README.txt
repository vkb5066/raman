For calculating raman spectra using VASP.  Current version is v4.

Theory is from:
	M. Bagheri, H. Komsa
	"High-throughput computation of Raman spectra from first principles"
	https://doi.org/10.1038/s41597-023-01988-5
and
	D. Porezag, M. Pederson
	"I.R. intensities and Raman-scattering activities within DFT"
	https://doi.org/10.1103/PhysRevB.54.7830

Calculating raman modes is a five step process.
Step 1: Calculate phonon eigenmodes using VASP.  You may use either finite     
	differences or DFPT.  Make sure that NWRITE = 3 is used in your INCAR 
	file.
Step 2: Generate the raman setup file.  The minimal command is ...
		./raman.py -setup -io </path/to/orig_outcar> -m <modes>
	where <modes> is a space-seperated list of modes to calculate.
Step 3: Generate the displaced POSCAR files.  The minimal command is ...
		./raman.py -gen
Step 4: Use VASP to do SCF calculations on each POSCAR file (in a unique 
	directory!).  Make sure that NWRITE = 3 and LEPSILON = .TRUE. are used
	in your INCAR files.
Step 5: Generate the raman spectra files.  The minimal command is ...
		./raman.py -eval -io <outcars>
	where <outcars> is a space seperated list of OUTCAR files resulting
	from step 4.

Tip: in the bash shell, typing, for example,
	{00000..00006}
expands to:
	00000 00001 00002 00003 00004 00005 00006
Use this to avoid typing out long lists of indices.