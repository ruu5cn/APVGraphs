#!/usr/bin/env python3

#Just the structure functions

#Used on 7/18
#7/18 update: This is the transverse lepton/unpolarized hadron case.
#stil working on what exactly k'*p means

#all units are in GeV, m,

import numpy as np
import lhapdf

GF = 1.16637e-5 #GeV^-2
HBARC = 1.97e-16 #GeV*m
ALPHA = 7.299e-3 #1
AHC = ALPHA#*HBARC #GeV*m
MP = 0.93828 #GeV
MN = 0.93957 #GeV
ME = 0.000511 #GeV #electron mass
MEC = ME*3.0e8 #Mass of electron times speed of light
MZ = 91.188 #GeV
PI = np.pi #1

QE = 0

GZ_ON = 1
Z_ON = 1

#PDF grid selection and creation
#unpol_pdf_name = "CT10"
unpol_pdf_name = "CT18NLO"
#unpol_pdf_name = "CT14nlo"
UNPOL_PDF = lhapdf.mkPDF(unpol_pdf_name,0)
pol_pdf_name = "NNPDFpol10_100"
POL_PDF = lhapdf.mkPDF(pol_pdf_name,0)

#sin squared of weak mixing angle
SIN2THETAW = 0.2232

PRINT = False

def gv(q):
  #neutral vector coupling constant of lepton/quark to Z boson

  #sin2thetaw =0.2232 # for initial cross-checks
  #sin2thetaw = 0.23
  #sin2thetaw = 0.232 #https://userweb.jlab.org/~xiaochao/memo/memo39-g3/memo39-polpv_20151110.pdf, page 3
  #sin2thetaw = 0.25 #Also above paper, but in page 4 
  
  if (((q==1) or (q==3)) or (q == 5)): #dsb
    return -0.5+ 2.0/3.0*SIN2THETAW
  if ((q==2) or (q==4) or (q==6)): #uct
    return 0.5 - (4.0/3.0)*SIN2THETAW
  if (q== QE): #electron
    return -0.5 + 2*SIN2THETAW
  return 0.0


def ga(q):
  #neutral axial vector coupling constant of lepton/quark to Z boson
  if (((q==1) or (q==3)) or (q == 5)): #dsb
    return -0.5;
  if ((q==2) or (q==4) or (q==6)): #uct
    return 0.5
  if (q==QE): #electron
    return -0.5;
  return 0.0;

        
def eq(q):
  #returns electric charge given quark number
  #uses common quark numbering convention (https://particle.wiki/wiki/PDG_particle_numbering_scheme)
  if (((q==1) or (q==3)) or (q == 5)): #dsb
      return (-1.0/3.0);
  if ((q==2) or (q==4) or (q==6)): #uct
      return (2.0/3.0);
  return 0.0

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
    quarkVal = (UNPOL_PDF.xfxQ2(q, x, Q2) + UNPOL_PDF.xfxQ2(-q, x, Q2))/x
    if (spr == "g"):
      val += eq(q)**2 * quarkVal
    elif (spr == "gz"):
      val += 2*eq(q)*gv(q)*quarkVal
    elif (spr == "z"):
      val += (gv(q)**2 + ga(q)**2)*quarkVal
    else:
      raise Exception("spr must be g, gz, or z")
  val *= x
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

def Sigma(x, Q2, s, qe, pol):
  global PRINT
  
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
  if PRINT and ((F1_gz > 10*F1_g) or (F2_gz > 10*F1_g) or (F3_gz > 10*F1_g)):
    print(f"""At x={x}, Q2={Q2};
    F1_g: {F1_g}
    F1_gz: {F1_g}
    F2_gz: {F1_gz}
    F3_gz: {F3_gz}
    """)
  if ((F1_gz > 100*F2_g) or (F2_gz > 100*F2_g) or (F3_gz > 100*F2_g)) and PRINT:
    print(f"""At x={x}, Q2={Q2};
    F2_g: {F2_g}
    F1_gz: {F1_g}
    F2_gz: {F1_gz}
    F3_gz: {F3_gz}
    """)
  if PRINT and ((F3_gz > 10*F1_gz)):
    print(f"""At x={x}, Q2={Q2};
    F1_gz: {F1_gz}
    F2_gz: {F2_gz}
    F3_gz: {F3_gz}
    """)
  if PRINT and False:
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
  """)
  eta_g = 1.0
  eta_gz = GZ_ON*((GF*MZ**2)/(2**1.5*PI*ALPHA)*(Q2)/(Q2 + MZ**2))
  eta_z = Z_ON*eta_gz**2

  M = MP
  m = ME
  mc = ME
  y = Q2/((s-M*M)*x)

  upsilon  = UPSILON(x, Q2, y, qe, M)
  zeta = ZETA(x, Q2, y, qe, M)
  xi = XI(x, Q2, y, qe, M)
  bigxi = BIGXI(x, Q2, y, qe, M)

  t1 = 0.0
  t2 = 0.0
  t3 = 0.0
  if (pol == "0"):
    t1 = -4*PI*y*AHC**2/(Q2**2)*((2*mc**2-Q2)*F1_g - 2*(2*mc**2-Q2)*gv(QE)*eta_gz*F1_gz  + ((gv(QE)**2 + ga(QE)**2)*(2*m**2 - Q2)-8*ga(QE)**2*m**2)*eta_z*F1_z)
    t2 = 4*PI*AHC**2/(x*y*Q2)*(upsilon*F2_g - 2*gv(QE)*eta_gz*upsilon*F2_gz + eta_z*((gv(QE)**2 + ga(QE)**2)*upsilon - 4*ga(QE)**2*m**2*M**2*x**2*y**2/(Q2**2))*F2_z)
    t3 = 4*PI*AHC**2*ga(QE)*(y-2)/Q2*(eta_gz*F3_gz - gv(QE)*eta_z*F3_z)
  if pol == "L":
    t1 = 16*PI*y*AHC**2*m*ga(QE)/(Q2**2)*(eta_gz*F1_gz) #- gv(QE)*eta_z*F1_z)
    t2 = -16*PI*m*M**2*x*y*AHC**2*ga(QE)/(Q2**3)*(eta_gz*F2_gz )#- gv(QE)*eta_z*F2_z)
    t3 = 4*PI*AHC**2*m/(Q2**2)*(-2*y*gv(QE)*eta_gz*F3_gz + (y*(gv(QE)**2 + ga(QE)**2) - 2*ga(QE)**2)*eta_z*F3_z)
    #kpDOTs = (s-M**2)/(2*M) #calling k' dot s = 1/2*E, which probably horribly wrong, but I think it's the best thing I could do right now
    kpDOTs = 1 #now we're trying this, just because
    #kpDOTs = x
    t1 *= kpDOTs #double-check later whether it's 1/2 or 1/4
    t2 *= kpDOTs
    t3 *= kpDOTs
  if PRINT:
    print(f"t1: {t1}")
    print(f"t2: {t2}")
    print(f"t3: {t3}")
    print(f"n_gz: {eta_gz}")
    print(f"y: {y}")
  return t1 + t2 + t3

def main():
  global SIN2THETAW
  
  E = 11.0
  s = 2*MP*E + MP**2 #This is unfortunately just fixed-target scattering, will need to \update for dual-beam collider
  MIN_Q2 = 1.0
  LG_ORD = 4 #largest order shown
  Q2_DENS = 8
  #MAX_Q2 = 12.0                                                                        
  N_Q2 = LG_ORD*Q2_DENS

  N_X = 50
  X_START = MIN_Q2/(s-MP**2)
  x = X_START
  Q2 = float(MIN_Q2)

  #data file that I write to is typically the day that I've run the program, in the for\mmDDMon                                                            
  data = open("10Aug00.txt", "w")
  data.write(f"{N_X}\n")
  #data.write(f"{N_Q2}\n")                                                             
  for i_q2 in range(N_Q2):
    min_x = Q2/(s-MP**2)
    max_x = Q2/(4-MP**2+Q2)
    x = max_x
    if max_x < min_x:
      break
    for i_x in range(N_X):
      sigma_b = Sigma(x, Q2, s, -1, "0")
      sigma_t = Sigma(x, Q2, s, -1, "L")
      apv = sigma_t/sigma_b
      data.write(f"{Q2},{x},{apv}\n")
      if N_X != 1: x *= (min_x/max_x)**(1.0/(N_X-1))
      #x += (MAX_X-MIN_X)/N_X                                                         
      #float(MIN_X) + (0.5/N_X)*(MAX_X - MIN_X)                                       
    Q2 *= 10**(1/Q2_DENS)
  print("Transversely polarized lepton, unpolarized hadron")
  data.close()

def singlePoint(x, Q2, s):
  s_t = Sigma(x, Q2, s, -1, "L")
  s_b = Sigma(x, Q2, s, -1, "0")

  APV = s_t/s_b
  print(f"For x={x}, Q2={Q2}, APV={APV}")
  
E = 11.0
s = 2*MP*E + MP*MP

x = 0.4
Q2 = 3.16

main()
#singlePoint(x,Q2,s)
