import os
import numpy as np

__author__ = 'Reto Trappitsch, 2016, trappitsch@uchicago.edu'
__version__ = '20160804'
__copyright__ = 'Copyright 2016 Reto Trappitsch (GPLv3)'

"""
Copyright 2016 Reto Trappitsch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


class INCLGridRun:
    """
    Automatic run:
    ==============
    Let's say you want to do an automatic run and would like to calculate the cross sections for Ne-21, then
    you shoul do the following (> lines are commands you give in, e.g., ipython)

        Import this file as ig
        > import incl_gridrun as ig

        Now we want to load an instance i
        > i = ig.INCLGridRun(cosnuc='ne21', projectile='p')

        Now let us run the the INCL simulations, the INCL reader, and the cross section writer:
        > i.runner()

    Now, wait and see, everything should be pretty much done for you.

    Setup:
    ======
    If you want to setup this script there is several things you need to take care of:
        1. Install root (https://root.cern.ch), this script is tested with version 6.06/06
        1. Compile INCLxx (http://irfu.cea.fr/Sphn/Spallation/incl.html)
        2. Go through the def__init__ routine and make sure that you set correctly the following variables
           > self.pathtoincl :          This is the absolute path to the INCL folder where the executable is
           > self.pathtoscript :        This is the path were this script here lives
           > self.pathtosimulations :   This is the path where your simulations will be stored. Might need some space
           > self.executable :          This is the name of the INCL executable.
    The rest of the script should just run and work without any issues.

    Manual run:
    ===========
    This is probably for advanced users. Check out the __init__ rouine and have a look at the manual setup section.
    If you run an automatic run, this is automatially disregarded.

    The most important subroutines for a manual run are the following
        > To run INCL setups for your setup:        incl_simulations()
        > To read cross sections from root files:   incl_csreader()
        > To write cross sections to sig.10 format: cs_writer()

    In a manual run you can still use the runner() script, which executes these simulatiosn in order.

    Folder structure for this script:
    =================================
    There are two subfolders:
        > gdata :   Contains important files for this to run. Please only modify if you know what you're doing
        > cs :      Here, the cross section outputs will be written, i.e., the final product.

    Changing target elements, isobars, etc.:
    ========================================
    If you would like to change the target elements, have a look at the routines that are outside the class.
    Especially the following routines are of importance here:
        > targetdict() :    Creates a dictionary for the target elements. Set your target elements here.
        > isobardict() :    Set your isobar here, i.e., if you want to change the what isotopes are considered.
                            Creates a dictionary for it.
        > isotopesdict() :  Reads in the file (see header of that routine) for all isotopic information. Then it
                            creates a dictionary for it.

    Changing the grid of energies:
    ==============================
    If you would like to change permanently the setup of the energy grid, please have a look at the routine
        > engrid_setup()
    which is located outside the class. This subroutine should be pretty self explanatory.

    Questions:
    ==========
    Finally, if you have any questions, feel free to drop me a line. My contact information is:
        Reto Trappitsch
        trappitsch@uchicago.edu

    ToDos for the future:
    =====================
    >   cross section files, even if already read out, are currently overwritten. A toggle should be installed
        to turn this overwriting on and off.
    >   Split isobaric decays, e.g., K-40, are currently not included, i.e., you can in this version either decay
        it all (by selecting the isobar) or not decay it at all. This will come in a future update
    """

    def __init__(self, cosnuc=None, projectile=None):
        """
        Values for the init. This needs to be set for each setup.
        :parameter cosnuc:      give your cosmogenic nuclide here, if None, then manual setup, list or string
        :parameter projectile:  give incident particle, p->proton, n->neutron, d->deuteron, t-> triton, a->alpha
                                He4, ..., string
        """
        # ########################################################################################################
        # ############################################ MANUAL SETUP # ############################################
        # ########################################################################################################

        # If you choose the manual setup, then make sure that the inputs cosnuc and projectiles are both None

        # cosmogenic nuclide
        self.cosnuc = ['ne21', 'ne22']

        # if you want to set up your own energygrid
        self.energygrid = [500., 1000., 1500.]
        # if you dont want that, please comment the line above and uncomment this next line
        # self.energygrid = engrid_setup()

        # target elements and isotopes
        # self.targeteles = ['Fe']
        # self.targetisos = [[54, 56, 57, 58]]
        # self.targetisoabu = [[0.05845, 0.91754, 0.02119, 0.00282]]
        self.targeteles = ['Na', 'Si']
        self.targetisos = [[23], [28, 29, 30]]
        self.targetisoabu = [[1.], [0.92223, 0.04685, 0.03092]]

        # mapper for target to cosmogenic nuclide
        # self.d_targcosmapper = dict(zip(['Fe'], [['Fe', ['ne21']]]))
        self.d_targcosmapper = dict(zip(['Na', 'Si'], [['Na', ['ne21']], ['Si', ['ne21', 'ne22']]]))

        # #########################################################################################################
        # ######################################### END OF MANUAL SETUP # #########################################
        # #########################################################################################################

        # projectile type (string)
        #   	proton, p
        #   	neutron, n
        #   	pi+, piplus, pion+, pionplus
        #   	pi0, pizero, pion0, pionzero
        #   	pi-, piminus, pion-, pionminus
        #   	d, t, a, deuteron, triton, alpha
        #   	He-4, He4, 4He (and so on)
        self.projectile = 'p'

        # ### Stuff that is not that often changed
        # number of shots
        self.numshots = 100000

        # which de-excitation model to use:
        #   	none
        #   	ABLA07 (default)
        #   	ABLAv3p
        #   	GEMINIXX
        #   	SMM
        self.deexmod = 'abla07'

        # Pathes setup
        # ### GEPARD ###
        self.pathtoincl = '/home/reto/opt/inclxx/'
        self.pathtoscript = '/home/reto/Documents/modeling/code/incl4.6.6_grid/'
        self.pathtosimulations = self.pathtoincl + 'simulations/'
        self.executable = 'INCLCascade'
        # ### HOLZWURM ###
        # self.pathtoincl = '/home/reto/opt/inclxx/'
        # self.pathtoscript = '/Users/reto/Documents/chicago_c3/modeling/code/incl4.6.6_grid/'
        # self.pathtosimulations = '/Users/reto/Documents/chicago_c3/modeling/modeldata/INCL_raw/'
        # self.executable = 'INCLCascade'

        # load some important dictionaries
        self.d_isobar = isobardict()
        self.d_isotopes = isotopesdict()
        self.d_targets = targetdict()

        # ### AUTOSETUP ###
        if cosnuc is not None:
            self.autosetup(cosnuc, projectile)

        print 'Copyright 2016 Reto Trappitsch\n\n' \
              'This program is free software: you can redistribute it and/or modify\n' \
              'it under the terms of the GNU General Public License as published by\n' \
              'the Free Software Foundation, either version 3 of the License, or\n' \
              '(at your option) any later version.\n\n' \
              'This program is distributed in the hope that it will be useful,\n' \
              'but WITHOUT ANY WARRANTY; without even the implied warranty of\n' \
              'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n' \
              'GNU General Public License for more details.\n\n' \
              'You should have received a copy of the GNU General Public License\n' \
              'along with this program.  If not, see <http://www.gnu.org/licenses/>.'

    def runner(self):
        """
        This routine runs everything after each other, for a full simulation, this is perfect
        :return:
        """
        # run incl simulations
        self.incl_simulations()

        # run the readout
        self.incl_readroot()

        # run the readout
        self.cs_writer()

    def autosetup(self, cosnuc, projectile):
        if cosnuc is None or projectile is None:
            print 'Choose your setup more wisely, either all or nothing has to be None'
            return

        # set the easy stuff
        if type(cosnuc) == str:
            self.cosnuc = [cosnuc]
        else:
            self.cosnuc = cosnuc
        self.projectile = projectile

        # do targets
        targets = []
        for cos in self.cosnuc:
            try:
                tmp = self.d_targets[cosnucsplitter(cos)[0].lower()][1]
                for jt in tmp:
                    targets.append(jt.lower())
            except KeyError:
                print 'The chosen cosmogenic nuclide ' + cos + ' does not have a target specified.'

        # now loop through the targets and make the targeteles, targetisos, and targetisoabu
        targeteles = []
        targetisos = []
        targetisoabu = []
        finder = False
        for it in range(len(targets)):
            for jt in range(len(targeteles)):
                if targets[it] == targeteles[jt].lower():
                    finder = True
                    break
            if not finder:
                try:
                    targeteles.append(self.d_isotopes[targets[it]][0])
                    targetisos.append(self.d_isotopes[targets[it]][2])
                    targetisoabu.append(self.d_isotopes[targets[it]][3])
                except KeyError:
                    print 'Somehow, the element ' + targets[it] + ' does not exist in self.d_isotopes dictionary. ' \
                                                                  ' This is not a user problem but a software' \
                                                                  ' issue... good luck!'

        # now attach all these to self
        self.targeteles = targeteles
        self.targetisos = targetisos
        self.targetisoabu = targetisoabu

        # now create the target to cosnuc mapper

        # # mapper for target to cosmogenic nuclide
        # self.targcosmapper = ['Fe', ['ne21']]
        targcosmapper = []
        for ele in targeteles:
            cosnuctmp = []
            for cos in self.cosnuc:
                costargs = self.d_targets[cosnucsplitter(cos)[0]][1]
                try:
                    costargs.index(ele)
                    cosnuctmp.append(cos)
                except ValueError:
                    continue
            targcosmapper.append([ele, cosnuctmp])

        self.d_targcosmapper = dict(zip(targeteles, targcosmapper))

        # and finally, set up the energy grid
        self.energygrid = engrid_setup()

    def incl_simulations(self):
        """
        This routine runs the INCL simulations. Everything is defined in the __init__ setup. If files exist, this
        is not run anymore.
        :return:
        """
        # loop through elements
        for ele in range(len(self.targeteles)):
            pathtoele = self.pathtosimulations + self.targeteles[ele] + '/'
            try:
                os.mkdir(pathtoele)
            except OSError:
                print 'Folder ' + pathtoele + ' already exists'

            # loop through self.targetisos
            for iso in range(len(self.targetisos[ele])):
                # create isotope folders
                pathtoiso = pathtoele + str(self.targetisos[ele][iso]) + '/'
                try:
                    os.mkdir(pathtoiso)
                except OSError:
                    print 'Folder ' + pathtoiso + ' already exists'

                # change to isotope folder
                os.chdir(pathtoiso)

                # loop through energies
                for ene in self.energygrid:
                    # create output filename
                    runname = self.targeteles[ele] + str(self.targetisos[ele][iso]) + '-' + str(ene) + 'mev-incpart-' \
                              + self.projectile

                    # check if file exists, if not, proceed, if it exists, don't proceed
                    try:
                        os.stat(runname + '.root')
                    except OSError:
                        # create cu46.in
                        finp = open(runname + '.in', 'w')
                        # write the file
                        finp.writelines('title = ' + runname + '\n')
                        finp.writelines('output = ' + runname + '\n')
                        finp.writelines('number-shots = ' + str(int(self.numshots)) + '\n')
                        finp.writelines('target = ' + self.targeteles[ele] + str(int(self.targetisos[ele][iso])) + '\n')
                        finp.writelines('projectile = ' + self.projectile + '\n')
                        finp.writelines('energy = ' + str(ene) + '\n')
                        finp.writelines('de-excitation = ' + self.deexmod + '\n')

                        # flush and close
                        finp.flush()
                        finp.close()

                        # now run the incl calculation
                        os.system(self.executable + ' ' + runname + '.in')

                        print runname + ' done'
        os.chdir(self.pathtoscript)
        print 'All INCL simulations done'

    def incl_readroot(self):
        """
        This subroutine sets up the C files and runs root with the myLooper. In every folder, there will be text files
        afterwards that contain the cross section files, i.e., cs_si28pXne21.txt. These files will contain two numbers,
        the calculated cross section (normalized and everything) plus its statistical uncertainty.
        :return:
        """

        # tested for a sigle nuclei, and all works well! August 4, 2016

        # loop through elements
        for ele in range(len(self.targeteles)):
            pathtoele = self.pathtosimulations + self.targeteles[ele] + '/'

            # loop through self.targetisos
            for iso in range(len(self.targetisos[ele])):
                # create isotope folders
                pathtoiso = pathtoele + str(self.targetisos[ele][iso]) + '/'
                # change to isotope folder
                os.chdir(pathtoiso)

                # loop through energies
                for ene in self.energygrid:
                    # ToDo: Include a check that looks if files exist and if so, act according to toggle in input
                    # root filename
                    rootname = self.targeteles[ele] + str(self.targetisos[ele][iso]) + '-' + str(ene) + 'mev-incpart-' \
                               + self.projectile + '.root'

                    # create the header file
                    headerfname = headerwriter(rootname, self.targeteles[ele], self.targetisos[ele][iso], ene,
                                               self.pathtoscript)

                    # ### C FILE ###
                    # create the .C file
                    cfilefname = 'INCLet_' + self.targeteles[ele] + str(int(self.targetisos[ele][iso])) + '_' + \
                                 str(ene) + 'MeV.C'
                    fcout = open(cfilefname, 'w')
                    # write the header
                    fcout.writelines('#define INCLet_cxx\n' +
                                     '#include "' + headerfname + '"\n' +
                                     '#include <TH2.h>\n' +
                                     '#include <TStyle.h>\n' +
                                     '#include <TCanvas.h>\n\n' +
                                     'void INCLet::Loop()\n' +
                                     '{\n' +
                                     '   if (fChain == 0) return;\n\n' +
                                     '   // read in the config variable from the gt tree\n' +
                                     '   double normalization;\n\n' +
                                     '   // read in definitions file\n' +
                                     '   ifstream fin;\n' +
                                     '   string fname = "conftmp";\n' +
                                     '   fin.open(fname.c_str(), ios_base::in);\n' +
                                     '   if(!fin){\n' +
                                     '      cout << "cannot open conftmp file for reading" << endl;\n' +
                                     '   }\n\n' +
                                     '   fin >> normalization;\n' +
                                     '   fin.close();\n\n')
                    # write initial energy
                    fcout.writelines('   // incident energy\n' +
                                     '   string incenergy = "' + str(ene) + '";\n\n')

                    # now initialize the counter ints and cs floats
                    for cos in self.d_targcosmapper[self.targeteles[ele]][1]:

                        fcout.writelines('   // ' + cos + '\n' +
                                         '   // initialize counter and cross section float\n' +
                                         '   int ' + cos + '_count = 0;\n' +
                                         '   float ' + cos + '_cs;\n' +
                                         '   float ' + cos + '_csunc;\n\n')

                    # now start writing the looper
                    fcout.writelines('   Long64_t nentries = fChain->GetEntriesFast();\n' +
                                     '   Long64_t nbytes = 0, nb = 0;\n' +
                                     '   for (Long64_t jentry=0; jentry<nentries; jentry++) {\n' +
                                     '      Long64_t ientry = LoadTree(jentry);\n' +
                                     '      if (ientry < 0) break;\n' +
                                     '      nb = fChain->GetEntry(jentry);   nbytes += nb;\n' +
                                     '      if (Cut(ientry) < 0) continue;\n\n' +
                                     '      // START CUSTOM PROGRAM\n\n\n')

                    # loop for each cosmogenic nuclide in the list
                    for cos in self.d_targcosmapper[self.targeteles[ele]][1]:
                        fcout.writelines('      // going for nuclide ' + cos + '\n\n')
                        # now loop through the isobar
                        acurrent = self.d_isobar[cos][1]
                        for zcurrent in self.d_isobar[cos][2]:
                            fcout.writelines('      if (A[nParticles] == ' + str(acurrent) + ' && Z[nParticles]==' +
                                             str(zcurrent) + ') {\n' +
                                             '         ' + cos + '_count += nParticles;\n' +
                                             '      }\n\n')

                    # end the loop
                    fcout.writelines('   }\n')

                    # write the closing statements outside the loop, like write all cs files
                    for cos in self.d_targcosmapper[self.targeteles[ele]][1]:
                        offname = self.targeteles[ele] + str(int(self.targetisos[ele][iso])) + self.projectile + 'X' + \
                                  cos + '_' + str(ene) + 'MeV.sig'
                        fcout.writelines('   // Open file for ' + cos + '\n\n' +
                                         '// Calculate cross sections ' + cos + ' and close file\n' +
                                         '   ' + cos + '_cs = ' + cos + '_count * normalization;\n' +
                                         '   ' + cos + '_csunc = sqrt(' + cos + '_count) * normalization;\n\n' +
                                         '   // Open files for writing out\n' +
                                         '   ofstream ' + cos + '_out;\n' +
                                         '   string ' + cos + '_fname = "' + offname + '";\n' +
                                         '   // Open file and check if opened properly\n' +
                                         '   ' + cos + '_out.open(' + cos + '_fname.c_str(), ios_base::out);\n' +
                                         '   if(!' + cos + '_out){\n' +
                                         '      cout << "Error while opening ' + offname + ' file" << endl;\n' +
                                         '   }\n' +
                                         '   // write to files\n' +
                                         '   ' + cos + '_out << ' + cos + '_cs << endl << ' + cos +
                                         '_csunc << endl;\n\n' +
                                         '   // close all the files\n' +
                                         '   ' + cos + '_out.flush();\n' +
                                         '   ' + cos + '_out.close();\n\n\n')

                    # now end the routine
                    fcout.writelines('}\n')

                    # close the fcout file
                    fcout.flush()
                    fcout.close()

                    # create the myLooper file
                    mylooperpfname = mylooper(rootname, cfilefname, self.targeteles[ele], iso, ene)

                    # now run the readout!
                    os.system('root -b -q ' + mylooperpfname)

        # cleanup
        os.chdir(self.pathtoscript)
        print 'All readout done'

    def cs_writer(self):
        """
        Finally, after running all the INCL simulations and reading out the cross sections from the root files,
        put them together here!
        :return:
        """
        for cos in self.cosnuc:
            # so we want to create a file for this specific cosmogenic nuclide, with the standard name.
            csout = open('cs/' + cos + self. projectile + 'sig.10', 'w')
            for eleit in range(len(self.targeteles)):
                ele = self.targeteles[eleit]
                # first check if elements is actually used for this isotope
                try:
                    # check if the current target produces the cosmogenic nuclide we want.
                    self.d_targcosmapper[ele][1].index(cos)
                    # change to the folder for the element
                    pathtoele = self.pathtosimulations + ele + '/'

                    # make elename and eleaa for writing to header
                    elename, eleaa = cosnucsplitter(cos)

                    # write the header to the file for this element
                    csout.writelines(' ' + ele.upper() + '-  0(' + self.projectile.upper() + ',     )' +
                                     elename.upper() + '- ' + str(eleaa) + '\n')

                    # loop through energygrid
                    for ene in self.energygrid:
                        # create a cross section list for this energy and isotope
                        csisos = []   # each line has tuple of cs and uncertainty
                        for isoit in range(len(self.targetisos[eleit])):
                            iso = self.targetisos[eleit][isoit]
                            pathtoiso = pathtoele + str(iso) + '/'
                            # the filename
                            finname = ele + str(iso) + self.projectile + 'X' + cos + '_' + str(ene) + 'MeV.sig'
                            csin = open(pathtoiso + finname, 'r')
                            csreadin = csin.read().split()
                            try:
                                cs = float(csreadin[0])
                                csunc = float(csreadin[1])
                                csisos.append((cs, csunc))
                            except ValueError:
                                print 'Problem moving to float for cross section or uncertainty in file ' + finname
                        # now i have all the isos for a given target element. Now normalize them according to
                        # the natural abundance, which is already 1 in the file
                        csnorm = 0.
                        csnormunc = 0.
                        for it in range(len(csisos)):
                            csnorm += csisos[it][0] * self.targetisoabu[eleit][it]
                            csnormunc += csisos[it][1] * self.targetisoabu[eleit][it]
                        # now write this energy to the csout file in the necessary format
                        if csnorm > 0:
                            # transform numbers to correct format:
                            enwrite = format(ene, '.3E')
                            cswrite = format(csnorm, '.3E')
                            csuncwrite = format(csnormunc, '.3E')
                            # write
                            csout.writelines('  ' + enwrite + '  0.000E+00  ' + cswrite + '  ' + csuncwrite +
                                             ' INCL4.6.6\n')

                    # now write the last 0 in there for ending this one
                    csout.writelines('  0.000E+00\n')

                except ValueError:
                    # so the cosmogenic nuclide is not in the list for this target, oh well, moving on to next target
                    continue

            # close the csout file
            csout.flush()
            csout.close()


# Additional routines, little helpers
def engrid_setup():
    """
    Returns the energy grid vector
    :return:
    """
    energygrid = range(1, 51, 1)
    i = 50
    while i < 100:
        i += 5
        energygrid.append(i)
    i = 100
    while i < 400:
        i += 20
        energygrid.append(i)
    i = 400
    while i < 10000:
        i += 200
        energygrid.append(i)

    # make to np array
    return np.array(energygrid)


def targetdict():
    """
    Returns a dictionary for cosmogenic nuclide versus target element
    :return:
    """
    targ = [['H',  ['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'Ca', 'Fe', 'Ni']],
            ['He', ['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'Ca', 'Fe', 'Ni']],
            ['Li', ['C', 'Si']],
            ['Be', ['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'Ca', 'Fe', 'Ni']],
            ['C',  ['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'Ca', 'Fe', 'Ni']],
            ['Ne', ['Na', 'Mg', 'Al', 'Si', 'Ca', 'Fe', 'Ni']],
            ['Al', ['Mg', 'Al', 'Si', 'Ca', 'Fe', 'Ni']],
            ['Cl', ['Cl', 'K', 'Ca', 'Ti', 'Fe', 'Ni']],
            ['Ar', ['K', 'Ca', 'Fe', 'Ni']],
            ['Ca', ['K', 'Ca', 'Ti', 'Fe', 'Ni']],
            ['Ti', ['Ti', 'Fe', 'Ni']],
            ['Mn', ['Fe', 'Co', 'Ni']],
            ['Fe', ['Ni']],
            ['Ni', ['Ni']],
            ['Kr', ['Sr', 'Rb', 'Y', 'Zr', 'Nb']],
            ['I',  ['Te', 'Ba', 'La']],
            ['Xe', ['Ba', 'La']]
            ]

    # make the dictionary
    targele = []
    for it in range(len(targ)):
        targele.append(targ[it][0].lower())

    targdict = dict(zip(targele, targ))
    return targdict


def isobardict():
    """
    Return an isobar and the assigned stuff
    :return:
    isodict: Dictionary for isobars
    """
    isob = [['he3', 3, [2, 1]],
            ['he4', 4, [2, 1]],
            ['li6', 6, [3, 2, 1, 4, 5]],
            ['li7', 7, [3, 2, 1, 4, 5]],  # including Be7
            ['be10', 10, [4, 3]],
            ['c14', 14, [6, 5, 4]],
            ['ne20', 20, [10, 9, 8, 7, 6, 11, 12]],
            ['ne21', 21, [10, 9, 8, 7, 6, 11, 12]],
            ['ne22', 22, [10, 9, 8, 7, 6, 11, 12, 13, 14]],  # contains also na22
            ['al26', 26, [13]],
            ['ar36', 36, [18, 17, 19, 20]],  # contains also cl36
            ['ar38', 38, [18, 17, 16, 15, 14, 19, 20]],
            ['ar40', 40, [18, 17, 16, 15, 14]],
            ['ca41', 41, [20, 21, 22]],
            ['ti44', 44, [22, 23, 24]],
            ['mn53', 53, [25, 26, 27, 28]],
            ['fe60', 60, [26, 25, 24, 23, 22]],
            ['ni59', 59, [28, 29, 30]],
            ['kr78', 78, [36, 35, 37, 38, 39, 40]],
            ['kr80', 80, [36, 35, 37, 38, 39, 40]],
            ['kr82', 82, [36, 35, 37, 38, 39, 40]],
            ['kr83', 83, [36, 35, 34, 33, 32, 37, 38, 39, 40]],
            ['kr84', 84, [36, 35, 34, 33, 32, 37]],
            ['kr86', 86, [36, 35, 34, 33]],
            ['xe124', 124, [54, 53, 55, 56, 57, 58]],
            ['xe126', 126, [54, 53, 55, 56, 57, 58]],
            ['xe128', 128, [54, 53, 55, 56, 57, 58]],
            ['xe129', 129, [54, 53, 52, 51, 50, 55, 56, 57, 58]],  # contains also i129
            ['xe130', 130, [54, 53, 55]],
            ['xe131', 131, [54, 53, 52, 51, 50, 55, 56, 57, 58]],
            ['xe132', 132, [54, 53, 52, 51, 50, 55]],
            ['xe134', 134, [54, 53, 52, 51, 50]],
            ['xe136', 136, [54, 53, 52, 51, 50, 55]]]
    # make a list with the isotope
    isotope = []
    for it in range(len(isob)):
        isotope.append(isob[it][0].lower())
    isodict = dict(zip(isotope, isob))
    return isodict


def isotopesdict():
    """
    Returns a dictionary for elements with all their isotopes and natural abundnaces
    Reads file isotopes.txt
    http://www.nndc.bnl.gov/nudat2/indx_sigma.jsp

    currently: split decays are not incorporated, but table would be ready for it

    :return:
    """
    # define what is in what column
    col_a = 0
    col_name = 1
    col_z = 2
    col_abu = 10

    # fead in file
    f = open('gdata/isotopes.txt', 'r')
    datain = []
    for line in f:
        datain.append(line.split('\t'))

    # now create a useful list with entries as following (per line): [elename, z, [a1, a2, a3], [abu1, abu2, abu3]]
    data = []
    zcomp = -1
    ele = None
    atmp = None
    aadd = []
    abuadd = []
    for it in range(1, len(datain)):   # skip header
        # if it's the same as last time
        if zcomp == int(datain[it][col_z]) and atmp != int(datain[it][col_a]):
            atmp = int(datain[it][col_a])
            aadd.append(atmp)
            abuadd.append(float(datain[it][col_abu].replace('%', '')) / 100.)
        elif zcomp != int(datain[it][col_z]):
            # add to list if not the first
            if zcomp != -1:   # take away the first one
                data.append([ele, zcomp, aadd, abuadd])
                aadd = []
                abuadd = []
            # make the settings for the first
            zcomp = int(datain[it][col_z])
            ele = datain[it][col_name].replace(' ', '')
            atmp = int(datain[it][col_a])
            aadd.append(atmp)
            abuadd.append(float(datain[it][col_abu].replace('%', '')) / 100.)

    # now get an element list and check that all the abundances sum up to 1 for each element
    elelist = []
    for it in range(len(data)):
        elelist.append(data[it][0].lower())
        # abutest
        if np.abs(np.sum(data[it][3]) - 1.) > 0.00001:
            print 'Problem in isotope abundance for element ' + data[it][0] + ', array: ' + str(data[it][3])

    # now make the dictionary
    dictele = dict(zip(elelist, data))

    return dictele


def cosnucsplitter(cnuc):
    """
    Splits the cosmogenic nuclide into a string and an integer, i.e., 'ne21' into 'ne', 21 and returns it
    :param cnuc: your cosmogenic nuclide as a string
    :return:
    """
    name = ''
    mass = ''
    for it in range(len(cnuc)):
        try:
            mass += str(int(cnuc[it]))
        except ValueError:
            name += cnuc[it]
    return name, int(mass)


def mylooper(rootfname, cfilefname, targele, targiso, ene):
    """
    Creates myLooper.C file from given stuff
    :return:
    """
    fname = 'myLooper_' + targele + str(int(targiso)) + '_' + str(ene) + 'MeV.C'

    towr = '{\n\n' \
           '// Read in normalization from global tree\n' \
           'TFile * file = new TFile("' + rootfname + '");\n' \
           'TTree * tree = (TTree *)file->Get("gt");\n' \
           'gt->GetEntry(gt->GetEntries()-1);\n' \
           'Float_t normalize = gt->GetLeaf("geometricCrossSection")->GetValue()\n' \
           '                    / gt->GetLeaf("nShots")->GetValue();\n\n' \
           '// write normalization and energy string to conftmp\n' \
           'ofstream fout;\n' \
           'string fname = "conftmp";\n' \
           'fout.open(fname.c_str(), ios_base::out);\n' \
           'if(!fout){\n' \
           '    cout << "Problem opening conftmp for writing" << endl;\n' \
           '}\n\n' \
           'fout << normalize;\n' \
           'fout.flush();\n' \
           'fout.close();\n\n' \
           '// Now run the C class\n' \
           '// compile the script\n' \
           'gROOT->ProcessLine(".L ' + cfilefname + '");\n\n' \
           '// look for the cross section stuff\n' \
           'gROOT->ProcessLine("INCLet k2");\n' \
           'k2.Loop();\n\n' \
           '}'

    # write the output file
    fout = open(fname, 'w')
    fout.write(towr)
    fout.flush()
    fout.close()

    return fname


def headerwriter(rootf, targele, targiso, ene, pts):
    """
    Writes the incl h file and gives back its filename
    :param rootf:   file name of the root file
    :param targele: target element ('fe')
    :param targiso: target isotope ('54')
    :param ene:     energy [MeV]
    :param pts:     The path to the script folder, since outside the class
    :return:        file name of the ICNLet header file
    """
    # file name for the header file
    fname = 'INCLet_' + targele + str(int(targiso)) + '_' + str(ene) + 'MeV.h'

    # read in part1 and 2 for the header file
    f = open(pts + 'gdata/INCLet_h_p1.txt', 'r')
    part1 = f.read()
    f.close()
    f = open(pts + 'gdata/INCLet_h_p2.txt', 'r')
    part2 = f.read()
    f.close()

    # now write the filling:
    filling = '      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject('\
              '\"' + rootf + '\");\n' \
                             '      if (!f || !f->IsOpen()) {\n' \
                             '         f = new TFile(\"' + rootf + '\"'

    # write the whole thing
    fout = open(fname, 'w')
    fout.write(part1)
    fout.write(filling)
    fout.write(part2)
    fout.flush()
    fout.close()

    return fname
