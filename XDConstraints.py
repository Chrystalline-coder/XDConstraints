# -----------------------------------------------------------------------------
#  "THE BEER-WARE LICENSE" (Revision 43):
#  <lkrause@gwdg.de> wrote this file. As long as you retain this notice you
#  can do whatever you want with this stuff. If we meet some day, and you think
#  this stuff is worth it, you can buy me a beer in return.
#  last change: 2014-07-29
#
#  README:
#  this script generates reset bond, Uij and riding hydrogen constraints
#  input: xd.mas
#  ----------------------------------------------------------------------------

# ##################################################################### #
#                                                                       #
#                         change values here!                           #
#                                                                       #
# ##################################################################### #

# round Uij values to X digits:
digits = 4

# electron distribution for all atoms the hydrogen constraints are generated for:
# SYNTAX: 'NAME':'electron distribution'
# 'NAME' doesn't need to match the scatter entry name! e.g. scatter entry is 'Cv' but 'NAME' can be 'Csp3'
# 'electron distribution' -> '2 -2 0 0 -2 (...)' would be: '2\s+-2\s+0\s+0\s+-2\s+(...)'
# the '\s+' term allows for any number of whitespaces to avoid parsing errors!
electron_distribution = {
   'H':'-1\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0',
'Csp2':'2\s+-2\s+0\s+0\s+-2\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0',
'Csp3':'2\s+-1\s+0\s+0\s+-3\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0',
   'N':'2\s+-2\s+0\s+0\s+-3\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0',
   'O':'2\s+-2\s+0\s+0\s+-4\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0'
                        }

# assign trivial names for fragments to allow for individual H-X distances and vibrational constraints
# SYNTAX: ('NAME','number of Hydrogen Atoms'):'Fragment Name'
# 'NAME' must match 'NAME' in electron_distribution!
# EXAMPLE CH3 group: Csp3 carbon atom with 3 hydrogen atoms attached -> ('Csp3',3):'CH3'
trivial_names = {
('Csp2',1):'Csp2-H',
('Csp2',2):'Csp2-H2 (check!)',
('Csp2',3):'Csp2-H3 (wrong!)',
('Csp3',1):'Csp3-H',
('Csp3',2):'Csp3-H2',
('Csp3',3):'Csp3-H3',
('Csp3',6):'Csp3-H3 (disordered!)',
('Csp3',9):'Csp3-H3 (disordered!)', #Added by CS 2015-06-30
   ('O',2):'OH2',
   ('O',1):'O-H',
   ('N',1):'N-H',
   ('N',2):'NH2',
   ('N',3):'NH3'
                }

# reset bond distances (H - Fragment)
# SYNTAX: 'Fragment Name':distance
# 'Fragment Name' must match 'Fragment Name' in trivial_names!
# REF: Allen and Bruno, Bond lengths revisited Acta Cryst. (2010). B66, 380-386
XH_distances = {
               'Csp2-H':1.083,
     'Csp2-H2 (check!)':1.082,
               'Csp3-H':1.099,
              'Csp3-H2':1.092,
              'Csp3-H3':1.077,
'Csp3-H3 (disordered!)':1.077,
                  'OH2':0.983,
                  'O-H':0.980,
                  'N-H':1.027,
                  'NH2':1.013,
                  'NH3':1.032
                }

# vibrational constraint for hydrogen atoms attached to 'Fragment Name', X times the vibration of 'Fragment Name':
# DEFAULT: 1.2000 (assumed if no value is specified in XH_Uij_values!)
# SYNTAX: 'Fragment Name':vibration
# comma seperated entries, 'Fragment Name' must match 'Fragment Name' in trivial_names!
XH_Uij_values = {
                'Csp3-H3':1.500,
  'Csp3-H3 (disordered!)':1.500
                }

# ##################################################################### #
#                                                                       #
#                     don't change values here!                         #
#                                                                       #
# ##################################################################### #

import re, math, collections
# ################################################ #
#                     classes                      #
# ################################################ #
class atom():
  def __init__(self):
    self.count    = 0
    self.number   = None # number in atom table
    self.name     = None # name
    self.bname    = None # name of atom the first principle axis points to; assumed to be the partner atom in case of hydrogen atoms
    self.bnumber  = None # number of partner atom
    self.scatname = None # scattering factor entry name
    self.scatnumb = None # scattering factor entry number
    self.dist     = 0.00 # Reset Bond distance for hydrogen atoms attached to this atom
    self.Uij      = 1.20 # Uij factor for hydrogen atoms attached to this atom
    self.element  = None # 'NAME' (assigned in electron_distribution) to avoid NameErrors concerning different scattering factor entry namings (e.g. Cv or CV or just C for sp3 carbon atoms)

  def assign(self,number,name,bname,bnumber,scatname,scatnumb,element):
    self.number   = number
    self.name     = name
    self.bname    = bname
    self.scatname = scatname
    self.scatnumb = scatnumb
    self.element  = element

  def update(self):
    try:
      self.bnumber = atoms[self.bname].number# get atom number of partner
    except KeyError:
      self.bnumber = 0# "bound" to dummy atom!
    try:
      self.type = [trivial_names[x] for x in trivial_names.keys() if x == (self.scatname,self.count)][0]# try to assign a trivial name for the fragment
    except IndexError:
      self.type = None# no special treatment!
    if self.element == 'H':
      try:
        self.type = atoms[self.bname].type
        self.dist = [XH_distances[x] for x in XH_distances.keys() if x == (atoms[self.bname].type)][0]
        self.Uij  = [XH_Uij_values[x] for x in XH_Uij_values.keys() if x == (atoms[self.bname].type)][0]
      except IndexError:
        pass  # no distance or Uij entry found!
      except AttributeError:
        pass  # 1st principle axis points on unknown atom!

# ################################################ #
#                    functions                     #
# ################################################ #
def Uij(xx,name,bnumber):
  xx = round(xx, digits)# avoid '1.2564332e-17' entries!
  if xx == 0:
    Uxx = ''
  else:
    Uxx = ' {:.{width}f} {}/{:3s}'.format(xx,name,str(bnumber),width=digits)
  return Uxx

def assign_atoms(electron_distribution):
  scatentry = []
  scatnumb = 0
  for i in scat_list:
    assigned = False
    scatnumb += 1
    for j in electron_distribution:
      if re.search(electron_distribution[j], i):
        assigned = True
        scatentry.append([str(scatnumb)]+[i.split()[0]]+[j])
    if assigned == False:
      scatentry.append([str(scatnumb)]+[i.split()[0]]+['X'])
  return scatentry

# ################################################ #
#                    header                        #
# ################################################ #
print '\n #######################################\n'
print '  THIS SCRIPT GENERATES RESET BOND, Uij\n'
print '     AND RIDING HYDROGEN CONSTRAINTS\n'
print ' #######################################\n'

# ################################################ #
#                  file handler                    #
# ################################################ #
master_name = raw_input(' enter filename [xd.mas]: ') or 'xd.mas'
if not re.search('\.mas', master_name):
  master_name = master_name + '.mas'
try:
  with open(master_name, 'r') as master_file:
    master = master_file.readlines()
except IOError:
  print ' error: {} not found!'.format(master_name)
  raise SystemExit

# ################################################ #
#                  project info                    #
# ################################################ #
scat_list = []
atom_list = []
cell_exp  = re.compile('^CELL ')
title_exp = re.compile('^TITLE ')
scat_exp  = re.compile('[a-zA-Z\+\-]+\s+CHFW\s+CHFW\s+CSZD\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+\s+-*\d+')
atom_exp  = re.compile('[a-zA-Z\(\)0-9\']+\s+[a-zA-Z\(\)0-9\']+\s+[XYZxyz]\s+[a-zA-Z\(\)0-9\']+\s+[a-zA-Z\(\)0-9\']+\s+[XYZxyz]\s+[RLrl]\s+\d\s+\d+\s+\d+\s+\d\s+[0-9a-zA-Z_]*[a-zA-Z\(\)0-9\']*')
for line in master:
  if scat_exp.search(line):
    scat_list.append(line)
  elif atom_exp.search(line):
    atom_list.append(line)
  elif title_exp.search(line):
    title = line.split()[1]
  elif cell_exp.search(line):
    a,b,c,alpha,beta,gamma = [float(j) for j in line.split()[1:]]

# ################################################ #
#               cell calculations                  #
# ################################################ #
volume = a*b*c*math.sqrt(1+2*math.cos(alpha*math.pi/180)*math.cos(beta*math.pi/180)*math.cos(gamma*math.pi/180)-math.cos(alpha*math.pi/180)*math.cos(alpha*math.pi/180)-math.cos(beta*math.pi/180)*math.cos(beta*math.pi/180)-math.cos(gamma*math.pi/180)*math.cos(gamma*math.pi/180))
a_ = b*c*math.sin(gamma*math.pi/180)/volume
b_ = a*c*math.sin(beta*math.pi/180)/volume
c_ = a*b*math.sin(alpha*math.pi/180)/volume

# ################################################ #
#                  introduction                    #
# ################################################ #
print '\n TITLE:  {}'.format(title)
print ' CELL:   {:.4f} {:.4f} {:.4f} {:.3f} {:.3f} {:.3f}'.format(a,b,c,alpha,beta,gamma)
print ' VOLUME: {:.2f} A^3'.format(volume)
Header_Out = '!TITLE: {}\n!CELL:  {:.4f} {:.4f} {:.4f} {:.3f} {:.3f} {:.3f}\n!set Kappa = 1.1 and Kappa\' = 1.18 for all Hydrogen Atoms in the corresponding *.inp file\n'.format(title,a,b,c,alpha,beta,gamma)

# ################################################ #
#              atom classification                 #
# ################################################ #
# link scattering factors from master file to atoms determined in electron_distribution 
scatentry = assign_atoms(electron_distribution)

# assign atoms to class atom() and fill the dictionnary atoms{}
atoms = collections.OrderedDict()
atom_number = 0
for atom_entry in atom_list:
  i = atom_entry.split()
  if len(i) > 4:# skip dummys
    atom_number += 1
    for j in scatentry:
      if i[8] == j[0]:# i[8] is the scatter entry number
        try:
          atoms[i[0]].assign(atom_number,i[0],i[1],None,j[2],j[0],j[2])
        except KeyError:# atom is not yet initiated as atom()!
          atoms[i[0]] = atom()
          atoms[i[0]].assign(atom_number,i[0],i[1],None,j[2],j[0],j[2])
        if atoms[i[0]].element == 'H':# Only the Hydrogen atoms are of interest
          try:
            atoms[i[1]].count += 1# count same atoms -> CH, CH2, CH3 etc
          except KeyError:
            atoms[i[1]] = atom()
            atoms[i[1]].count += 1

# ################################################ #
#                header and list                   #
# ################################################ #
print '\n {:<7} {:<3} {:<7} {:<4} {:<7} {:<7}'.format('ATOM','#','ATOM0','#','DIST','TYPE')

warnings = []
for i in atoms:
  atoms[i].update()
  ''' DEBUG 
  print ' {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(atoms[i].number,atoms[i].name,atoms[i].bname,atoms[i].bnumber,atoms[i].scatname,atoms[i].scatnumb,atoms[i].dist,atoms[i].Uij,atoms[i].element)
  DEBUG '''
  if atoms[i].element == 'H':
    if atoms[i].type == None:# unknown partner!
      warnings.append([atoms[i].name] + [atoms[i].bname])
    print ' {:<7} {:<3} {:<7} {:<4} {:<7.3f} {:<}'.format(atoms[i].name,atoms[i].number,atoms[i].bname,atoms[i].bnumber,atoms[i].dist,atoms[i].type)
print ' -\n'

# ################################################ #
#                     alerts                       #
# ################################################ #
Alerts_Out = ''
if warnings:
  Alerts_Out = '\n!Warning! {} alert(s) pending:\n'.format(str(len(warnings)))
  for i in warnings:
    Alerts_Out += '! - {:6} first principle axis points to unknown atom: {:6}\n'.format(str(i[0]),str(i[1]))
  Alerts_Out += '! open xd_constraints.py in an editor for help!\n'
  print Alerts_Out

# ################################################ #
#             main loop starts here                #
# ################################################ #
Uiso_Out     = '\n!Uiso CONSTRAINTS\n'
Riding_Out   = '\n!RIDING HYDROGEN CONSTRAINTS\n'
Reset_Out    = '\n!RESET BONDS\n'
Riding_Alert = ''

for i in atoms:
  if atoms[i].element == 'H':
  # ################################################ #
  #            hydrogen Uiso constraints             #
  # ################################################ #
    fac = atoms[i].Uij
    aa = fac*a*a*a_*a_/3
    ab = fac*2*a*a_*b_*b*math.cos(gamma*math.pi/180)/3
    ac = fac*2*a*a_*c*c_*math.cos(beta*math.pi/180)/3
    bb = fac*b*b*b_*b_/3
    bc = fac*b*b_*c*c_*math.cos(alpha*math.pi/180)/3
    cc = fac*c*c*c_*c_/3
    U11 = Uij(aa,'U11',atoms[i].bnumber)
    U12 = Uij(ab,'U12',atoms[i].bnumber)
    U13 = Uij(ac,'U13',atoms[i].bnumber)
    U22 = Uij(bb,'U22',atoms[i].bnumber)
    U23 = Uij(bc,'U23',atoms[i].bnumber)
    U33 = Uij(cc,'U33',atoms[i].bnumber)
    if [True for k in warnings if atoms[i].name == k[0]]:
      Uiso_Out += '!'
    Uiso_Out += 'CON{}{}{}{}{}{} -1 U11/{} = 0 !Atom {:<3} is {:6s}, atom {:<3} is {:6s}\n'.format(U11,U12,U13,U22,U23,U33,atoms[i].number,atoms[i].bnumber,atoms[i].bname,atoms[i].number,atoms[i].name)

  # ################################################ #
  #           riding hydrogen constraints            #
  # ################################################ #
    if [True for k in warnings if atoms[i].name == k[0]]:
      Riding_Alert = '!'
    else:
      Riding_Alert = ''
    X = '{}CON 1.0 X/{:<3} -1.0 X/{:<3} = 0'.format(Riding_Alert,atoms[i].bnumber,atoms[i].number)
    Y = '{}CON 1.0 Y/{:<3} -1.0 Y/{:<3} = 0'.format(Riding_Alert,atoms[i].bnumber,atoms[i].number)
    Z = '{}CON 1.0 Z/{:<3} -1.0 Z/{:<3} = 0'.format(Riding_Alert,atoms[i].bnumber,atoms[i].number)
    Riding_Out += '{} !Atom {} is {}, atom {} is {}\n{}\n{}\n'.format(X,atoms[i].bnumber,atoms[i].bname,atoms[i].number,atoms[i].name,Y,Z)

  # ################################################ #
  #               reset bond commands                #
  # ################################################ #
    if [True for k in warnings if atoms[i].name == k[0]]:
      Reset_Out += '!'
    Reset_Out += 'RESET BOND {:6} {:6} {:6.3f}\n'.format(atoms[i].bname,atoms[i].name,atoms[i].dist)

# ################################################ #
#               write output file                  #
# ################################################ #
with open('xd.const', 'w') as Out:
  Out.write(Header_Out + Alerts_Out + Uiso_Out + Riding_Out + Reset_Out)
  print ' {} written!'.format(Out.name)


#print ' __  __ ____     ____   ____   _   _  ____  _____  ____      _     _  _   _  _____  ____  '
#print ' \ \/ /|  _ \   / ___| / __ \ | \ | |/ ___||_   _||  _ \    / \   | || \ | ||_   _|/ ___| '
#print '  \  / | | | | | |    | |  | ||  \| |\___ \  | |  | |_) |  / _ \  | ||  \| |  | |  \___ \ '
#print '  /  \ | |_| | | |___ | |__| || |\  | ___) | | |  |  _ <  / ___ \ | || |\  |  | |   ___) |'
#print ' /_/\_\|____/   \____| \____/ |_| \_||____/  |_|  |_| \_\/_/   \_\|_||_| \_|  |_|  |____/ '