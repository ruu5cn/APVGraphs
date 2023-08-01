#Graph Drawer
#I would copy and paste this code into a Jupyer notebook to create the graph.

#Requires ROOT to be useable

from __future__ import print_function
from ROOT import TCanvas, TGraph
from ROOT import gROOT
import ROOT

from array import array
from decimal import Decimal

RF = -1.0e6 #factor that multiplies APV when it reads in

c1 = TCanvas( 'c1', 'Asymmetry Graph', 200, 10, 800, 800 )
c1.SetLeftMargin(0.22);
 
#c1.SetFillColor( 42 )
c1.SetGrid()
c1.SetLogx()
c1.SetLogy()

#Put the name of your data file here!
data = open("data/28Jul02.txt", "r")

#read in important data
N_X = int(data.readline())

graphList = []
xList = [] #not a list of x-values, each item in this list is a list of x-values
apvList = [] #ditto

Q2List = []  #This is actually a list of Q2 values
textList = []


line = data.readline()
p1 = line.find(",", 0, -1)
p2 = line.find(",", p1 + 1, -1)

maxX = 0.5
minX = 0.5 #helps create limits on graph

while line != "":
    Q2Val = float(line[0:p1])
    Q2List.append(float(line[0:p1]))
    xList.append(array('d'))
    apvList.append(array('d'))
    pnts = 0;
    for i in range(N_X):
        xVal = float(line[p1+1:p2])
        xList[-1].append(float(line[p1+1:p2]))
        #Q2 = float(line[p1 + 1: p2])
        apvList[-1].append(float(line[p2 + 1: -1])*RF)
        pnts += 1
        
        line = data.readline()
        p1 = line.find(",", 0, -1)
        p2 = line.find(",", p1 + 1, -1)
    #this runs after we've been through all the x's for a given Q2
    if max(xList[-1]) > maxX:
        maxX = max(xList[-1])
    if min(xList[-1]) < minX:
        minX = min(xList[-1])
    graphList.append(TGraph(pnts, xList[-1], apvList[-1]))

mg = ROOT.TMultiGraph()
mg.SetTitle("Hadron polarized (JLab);")
#mg.CenterTitle()
#mg.GetXaxis().SetTitle("x")

for i_g in range(len(graphList)):
    graphList[i_g].SetLineColor(4)
    graphList[i_g].SetLineWidth(3)
    #graphList[i_g].SetMarkerStyle(2)
    #add label
    #label = "Q2 = " + str(Q2List[i_g])
    mg.Add(graphList[i_g])
    textList.append(ROOT.TText(xList[i_g][0]*1.01, apvList[i_g][0], "Q2 = " + str(round(Q2List[i_g], 2))))
    textList[i_g].SetTextSize(0.02)
        #textList[i_g].Draw()
    
print(textList) 
mg.Draw("ACP")
mg.GetXaxis().SetTitle("x")
mg.GetYaxis().SetTitle("APV (-ppm)")
mg.GetXaxis().SetLimits(0.5*minX,2.0*maxX);

for i in range(len(graphList)):
    textList[i].Draw()

# TCanvas.Update() draws the frame, after which one can change it
c1.Update()
#c1.GetFrame().SetFillColor( 21 )
#c1.GetFrame().SetBorderSize( 12 )
c1.Modified()
c1.Update()
c1.Draw()

#Save your graph here!
c1.SaveAs("UnpolL-LongH-JLAB.pdf")
