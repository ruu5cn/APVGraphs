#!/usr/bin/env python3

#Creates (Q2, x, Asymmetry) data points and writes them into a data file

#Used on 7/19
#7/13 update: scans (Q2, x) points in a specific triangular region instead of in a square
#7/19 update: replaces fixed-target formulas for s with dual beam ones. Fixed g1_g LH bug 
#7/24 update: tested out some EIC graphs
#7/25 update: Make sure to iron out stuff with M

#The most important function is Sigma(), which is short for d(sigma)/dxdy
#all units are in GeV, m

############################################# 
#important settings!
FIXED_TARGET = False
POSITRON_BEAM = False
############################################# 

import numpy as np
import lhapdf

#Universal constants
SIN2THETAW = 0.23 #1, sin squared of weak mixing angle. Slight changes have a surprisingly large effect   
GF = 1.16637e-5 #GeV^-2
HBARC = 1.97e-16 #GeV*m
ALPHA = 7.299e-3 #1
AHC = HBARC*ALPHA #GeV*m
MP = 0.93828 #GeV
MN = 0.93957 #GeV
MZ = 91.188 #GeV
PI = np.pi #1

#in-code constants
M = MP #To allow easier switching between MP and MN
QE = 7 
PRINT = False #Prints out structure functions values and each term when True

#PDF grid selection and creation
#unpol_pdf_name = "CT10"
unpol_pdf_name = "CT18NLO"
#unpol_pdf_name = "CT14nlo"
UNPOL_PDF = lhapdf.mkPDF(unpol_pdf_name,0)
pol_pdf_name = "NNPDFpol10_100"
POL_PDF = lhapdf.mkPDF(pol_pdf_name,0)


############################################# 
#coupling constants
def gv(q):
  #neutral vector coupling constant of lepton/quark to Z boson
  #uses common quark numbering convention (https://particle.wiki/wiki/PDG_particle_numbering_scheme) 
  if (((q==1) or (q==3)) or (q == 5)): #dsb
    return -0.5+ 2.0/3.0*SIN2THETAW
  if ((q==2) or (q==4) or (q==6)): #uct
    return 0.5 - (4.0/3.0)*SIN2THETAW
  if (q== QE): #electron
    return -0.5 + 2*SIN2THETAW
  return 0.0

def ga(q):
  #neutral axial vector coupling constant of lepton/quark to Z boson
  #uses common quark numbering convention (https://particle.wiki/wiki/PDG_particle_numbering_scheme) 
  if (((q==1) or (q==3)) or (q == 5)): #dsb
    return -0.5;
  if ((q==2) or (q==4) or (q==6)): #uct
    return 0.5
  if (q==7): #electron
    return -0.5;
  if (q == -7):
    return 0.5
  return 0.0;

def eq(q):
  #returns electric charge given quark number
  #uses common quark numbering convention (https://particle.wiki/wiki/PDG_particle_numbering_scheme)
  if (((q==1) or (q==3)) or (q == 5)): #dsb
      return (-1.0/3.0);
  if ((q==2) or (q==4) or (q==6)): #uct
      return (2.0/3.0);
  return 0.0

############################################# 
#various helpful variables
def UPSILON(x, Q2, y, qe, M):
  return 1 - y - (M*M*x*x*y*y)/Q2
def ZETA(x, Q2, y, qe, M):
  return 2 - y - (2*M*M*x*x*y*y)/Q2
def XI(x, Q2, y, qe, M):
  return 1 - y - (M*M*x*x*y)/Q2
def BIGXI(x, Q2, y, qe, M):
  return 1 + 2*M*M*x*x*y/Q2


#############################################
#Structure Functions
def F2(x, Q2, spr, qe):
  val = 0.0 #value to be returned
  for q0 in range(6):
    q = q0 + 1
    quarkVal = (UNPOL_PDF.xfxQ2(q, x, Q2) + UNPOL_PDF.xfxQ2(-q, x, Q2)) #don't forget the implied *x
    if (spr == "g"):
      val += eq(q)**2 * quarkVal
    elif (spr == "gz"):
      val += 2*eq(q)*gv(q)*quarkVal
    elif (spr == "z"):
      val += (gv(q)**2 + ga(q)**2)*quarkVal
    else:
      raise Exception("spr must be g, gz, or z")
  return val

def F1(x, Q2, spr, qe):
  return F2(x, Q2, spr, qe)/(2*x)

def F3(x, Q2, spr, qe):
  val =	0.0 #value to be returned                                                                                   
  for q0 in range(6):
    q =	q0 + 1
    quarkVal = (UNPOL_PDF.xfxQ2(q, x, Q2) - UNPOL_PDF.xfxQ2(-q, x, Q2))/x
    if (spr == "g"):
      val += 0 * quarkVal
    elif (spr == "gz"): 
      val += 2*eq(q)*ga(q) * quarkVal
    elif (spr == "z"):
      val += 2*gv(q)*ga(q) * quarkVal
    else:
      raise Exception("spr must be g, gz, or z")
  return val

def g1(x, Q2, spr, qe):
  val = 0.0 #value to be returned
  for q0 in range(6):
    q = q0 + 1
    quarkVal = (POL_PDF.xfxQ2(q, x, Q2) + POL_PDF.xfxQ2(-q, x, Q2))/x
    if (spr == "g"):
      val += eq(q)**2 * quarkVal
      #print(f"q: {q}; qval: {quarkVal}")
      #print(f"q: {q}; +val: {0.5*eq(q)**2 * quarkVal}")
    elif (spr == "gz"):
      val += 2*eq(q)*gv(q) * quarkVal
    elif (spr == "z"):
      val += (gv(q)**2 + ga(q)**2)*quarkVal
    else:
      raise Exception("spr must be g, gz, or z")
  val *= 0.5
  return val

def g2(x, Q2, spr, qe):
  #dubious but it might make things work
  val = 0.0 #value to be returned
  for q0 in range(6):
    q = q0 + 1
    quarkVal = (POL_PDF.xfxQ2(q, x, Q2) + POL_PDF.xfxQ2(-q, x, Q2))/x
    if (spr == "g"):
      val += 0 * quarkVal
    elif (spr == "gz"):
      val += 0 * quarkVal
    elif (spr == "z"):
      val += 0.5*(ga(q)**2) *quarkVal
    else:
      raise Exception("spr must be g, gz, or z")
  return val

def g5(x, Q2, spr, qe):
  val = 0.0 #value to be returned
  for q0 in range(6):
    q = q0 + 1
    quarkVal = (POL_PDF.xfxQ2(q, x, Q2) - POL_PDF.xfxQ2(-q, x, Q2))/x
    if (spr == "g"):
      val += 0 * quarkVal
    elif (spr == "gz"):
      val += eq(q)*ga(q) * quarkVal
    elif (spr == "z"):
      val += gv(q)*ga(q)*quarkVal
    else:
      raise Exception("spr must be g, gz, or z")
  return val

def g4(x, Q2, spr, qe):
  return g5(x, Q2, spr, qe)*2*x


####################################
#Function for Cross-section
def Sigma(x, Q2, s, qe, pol):
  global PRINT, M
  
  F1_g = F1(x, Q2, "g", qe)
  F1_gz = F1(x, Q2, "gz", qe)
  F1_z = F1(x, Q2, "z", qe)
  F2_g = F2(x, Q2, "g", qe)
  F2_gz = F2(x, Q2, "gz", qe)
  F2_z = F2(x, Q2, "z", qe)
  F3_g = F3(x, Q2, "g", qe)
  F3_gz = F3(x, Q2, "gz", qe)
  F3_z = F3(x, Q2, "z", qe)
  g1_g = g1(x, Q2, "g", qe)
  g1_gz = g1(x, Q2, "gz", qe)
  g1_z = g1(x, Q2, "z", qe)
  g4_g = g4(x, Q2, "g", qe)
  g4_gz = g4(x, Q2, "gz", qe)
  g4_z = g4(x, Q2, "z", qe)
  g5_g = g5(x, Q2, "g", qe)
  g5_gz = g5(x, Q2, "gz", qe)
  g5_z = g5(x, Q2, "z", qe)

  if PRINT:
    print(f"""
  F1_g: {F1_g}
  F1_gz: {F1_gz}
  F1_z: {F1_z}
  F2_g: {F2_g}
  F2_gz: {F2_gz}
  F2_z: {F2_z}
  F3_g: {F3_g}
  F3_gz: {F3_gz}
  F3_z: {F3_z}
  g1_g: {g1_g}
  g1_gz: {g1_gz}
  g1_z: {g1_z}
  g4_g: {g4_g}
  g4_gz: {g4_gz}
  g4_z: {g4_z}
  g5_g: {g5_g}
  g5_gz: {g5_gz}
  g5_z: {g5_z}
  """)
    
  eta_g = 1.0
  eta_gz = ((GF*MZ**2)/(2**1.5*PI*ALPHA)*(Q2)/(Q2 + MZ**2))
  eta_z = eta_gz**2

  y = Q2/((s-M*M)*x)
  upsilon  = UPSILON(x, Q2, y, qe, M)
  zeta = ZETA(x, Q2, y, qe, M)
  xi = XI(x, Q2, y, qe, M)
  bigxi = BIGXI(x, Q2, y, qe, M)
  
  t1 = 0.0
  t2 = 0.0
  t3 = 0.0
  t4 = 0.0
  t5 = 0.0

  if (pol=="0"):
    t1 = (4 * PI * y * (AHC**2))/Q2 * (F1_g - eta_gz * gv(QE) * F1_gz + eta_z * (gv(QE)**2 + ga(QE)**2 ) * F1_z)
    t2 = (4*PI*AHC**2)/(x*y*Q2)*upsilon * (F2_g - eta_gz*gv(QE)*F2_gz + eta_z*(gv(QE)**2 + ga(QE)**2)*F2_z)
    t3 = -(2*PI*AHC**2)*(2-y)/Q2*(eta_gz*ga(QE)*F3_gz - 2*eta_z*gv(QE)*ga(QE)*F3_z)
    if PRINT:
      print(f"bottom t1: {t1}")
      print(f"bottom t2: {t2}")
      print(f"bottom t3: {t3}")
  elif (pol=="L"):
    t1 = (4*PI*y*AHC**2)/Q2*(eta_gz*ga(0)*F1_gz - 2*eta_z*gv(0)*ga(0)*F1_z);
    t2 = (4*PI*AHC**2)/(x*y*Q2)*upsilon*(eta_gz*ga(0)*F2_gz-2*eta_z*gv(0)*ga(0)*F2_z);
    t3 = -(2*PI*AHC**2)*(2-y)/Q2*(-eta_gz*gv(0)*F3_gz + eta_z*(gv(0)*gv(0)+ga(0)*ga(0))*F3_z);
    if PRINT:
      print(f"top t1: {t1}")
      print(f"top t2: {t2}")
      print(f"top t3: {t3}")
  elif (pol=="H"):
    t1 = (4*PI*AHC**2)/Q2*zeta * (eta_gz*ga(QE)*g1_gz - 2*eta_z*gv(QE)*ga(QE)*g1_z)
    #t2, t3 = 0 since g2 = g3 = 0 under parton model
    t4 = -(4*PI*AHC**2)/(x*y*Q2)*bigxi*upsilon * (-eta_gz*gv(QE)*g4_gz + eta_z*(gv(QE)**2 + ga(QE)**2)*g4_z)
    t5 = -(4*PI*y*AHC**2)/Q2*bigxi * (-eta_gz*gv(QE)*g5_gz + eta_z*(gv(QE)**2 + ga(QE)**2)*g5_z)
    if PRINT:
      print(f"top t1: {t1}")
      print(f"top t4: {t4}")
      print(f"top t5: {t5}")
  elif (pol=="LH"):
    t1 = (4*PI*AHC**2)/Q2*zeta * (g1_g-eta_gz*gv(QE)*g1_gz + eta_z*(gv(QE)**2 + ga(QE)**2)*g1_z)
    #t2, t3 = 0 since g2 = g3 = 0 under parton model
    t4 = -(4*PI*AHC**2)/(x*y*Q2)*bigxi*upsilon * (eta_gz*ga(QE)*g4_gz - 2*eta_z*gv(QE)*ga(QE)*g4_z)
    t5 = -(4*PI*y*AHC**2)/Q2*bigxi * (eta_gz*ga(QE)*g5_gz - 2*eta_z*gv(QE)*ga(QE)*g5_z)
    if PRINT:
      print(f"top t1: {t1}")
      print(f"top t4: {t4}")
      print(f"top t5: {t5}")
  else:
    raise Exception("pol must be 0 (no polarization), L (lepton polarization), H (hadron polarization), or LH (both lepton and handron polarization)")

  """
  term_l = [t1, t2, t3, t4, t5]
  for t in range(len(term_l)):
    print(f't{t+1}: {term_l[t]}')
  """
  return t1 + t2 + t3 + t4 + t5


############################################# 
def main():
  global SIN2THETAW, M

  pol = "L"
  if FIXED_TARGET:
    E = 11.0
    s = 2*M*E + M*M
  else:
    Ek = 18.0
    Ep = 275.0
    s = 4*Ek*Ep
  
  MIN_Q2 = 1.0
  LG_ORD = 4 #largest order shown
  Q2_DENS = 4 #number of Q2 lines per order of 10
  N_Q2 = LG_ORD*Q2_DENS
  N_X = 50 #number of x-points on a graph
  
  Q2 = float(MIN_Q2)
  
  #data file that I write to is typically the day that I've run the program, in the form DDMon
  data = open("01Aug00.txt", "w")
  data.write(f"{N_X}\n")
  
  for i_q2 in range(N_Q2):
    min_x = Q2/(s-M**2)
    max_x = Q2/(4-M**2+Q2)
    x = max_x
    if max_x < min_x:
      break #stops data collection entirely at this point because something's clearly gone wrong
    for i_x in range(N_X):
      sigma_b = Sigma(x, Q2, s, -1, "0")
      sigma_t = Sigma(x, Q2, s, -1, pol)
      apv = sigma_t/sigma_b
      
      data.write(f"{Q2},{x},{apv}\n") #change this line if you don't want the data read to a file
      
      x *= (min_x/max_x)**(1.0/(N_X-1))
    #runs after all x points are collected
    Q2 *= 10**(1/Q2_DENS)
    
  data.close()
  print("APV5, Pol " + pol)

main()
