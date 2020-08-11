## Particle mass and energy must be set to determine other decay properties.
## Add angular distribution for more information.

## Should be callable with lists of particles too?

import ROOT

# Units in GeV
mass = 200

energy = 500

# Units in ns
taus = [0.01,0.1,1,10,100]

# Store graph formatting
graph_formatting = {
  "total_dist" : {"xrange": (5e-5,20),
                  "yrange": (0,1),
                  "xlabel": "Total decay displacement [m]",
                  "ylabel": "Fraction of decays"},
  "decayed_by_dist" : {"xrange": (5e-5,20),
                  "yrange": (0,1),
                  "xlabel": "Total decay displacement [m]",
                  "ylabel": "Probability particle has decayed"},
  "transverse_dist" : {"xrange": (5e-5,20),
                  "yrange": (0,1),
                  "xlabel": "Transverse decay displacement [m]",
                  "ylabel": "Fraction of decays"},
  "dist" : {"xrange": (5e-5,20),
                  "yrange": (0,1),
                  "xlabel": "Decay displacement [m]",
                  "ylabel": "Fraction of decays"},                  
}

def make_graph(xvals, yvals) :

  graph = ROOT.TGraph()
  for x, y in zip(xvals,yvals) :
    graph.SetPoint(graph.GetN(),x,y)
  return graph

# Imports ATLAS style for plotting
import AtlasStyle
AtlasStyle.SetAtlasStyle()

def plot_graphs(graphs,legend_list,graphname,tag="") :

  # Make a canvas to put the plot on.
  # Format it to have log x axis.
  c = ROOT.TCanvas("canvas",'',0,0,600,600)
  c.SetLogx(True)
  c.SetLogy(False)
  c.SetGridx(0)
  c.SetGridy(0)

  # Retrieve info
  properties = graph_formatting[graphname]

  # Decide what colours to use.
  # These ones look decent, but obviously use
  # whatever you like best.
  goodColours = [ROOT.kRed-7,ROOT.kAzure+1,ROOT.kTeal+1,ROOT.kOrange-4,ROOT.kViolet-5]

  legend = ROOT.TLegend(0.2,0.2,0.52,0.4)
  # Make the text a nice fond, and big enough
  legend.SetTextFont(42)
  legend.SetTextSize(0.04)
  # A few more formatting things .....
  legend.SetBorderSize(0)
  legend.SetLineColor(0)
  legend.SetLineStyle(1)
  legend.SetLineWidth(1)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)

  for i,item in enumerate(graphs) :
    if type(item) is tuple :
      these = item
      these[1].SetLineStyle(2)
    else :
      these = (item,)
    for j,graph in enumerate(these) :
      graph.SetLineColor(goodColours[i])
      graph.SetLineWidth(3)
      graph.GetXaxis().SetLimits(properties["xrange"][0],properties["xrange"][1])
      graph.GetYaxis().SetRangeUser(properties["yrange"][0],properties["yrange"][1])    
      if i > 0 or j > 0 :
        graph.Draw("L")
      else :
        graph.GetXaxis().SetTitle(properties["xlabel"])
        graph.GetYaxis().SetTitle(properties["ylabel"])      
        graph.Draw("AL")
      if j == 0 :
        legend.AddEntry(graph,legend_list[i],"L")

  # Actually draw the legend
  legend.Draw()

  # This is one way to draw text on the plot
  myLatex = ROOT.TLatex()
  myLatex.SetTextColor(ROOT.kBlack)
  myLatex.SetNDC()
  myLatex.SetTextSize(0.04)
  myLatex.SetTextFont(42)
  myLatex.DrawLatex(0.2,0.46,"M: {0} GeV".format(mass))
  myLatex.DrawLatex(0.2,0.42,"E: {0} GeV".format(energy))

  # Update the canvas
  c.Update()

  # Save the output
  savename = graphname
  if tag : savename = savename + "_"+tag
  c.SaveAs(savename+".eps")

  # Save another without the log
  c.SetLogx(False)
  c.SaveAs(savename+"_linear.eps")

# And start
import numpy as np
from scipy import constants
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import scipy.special as sc

graphs = {}

for tau in taus :

  # Dict to hold my graphs
  tau_string = "{0}".format(tau)
  graphs[tau_string] = {}

  # tau is in ns; convert to seconds
  tau_s = tau * constants.nano
  print("Tau in seconds:",tau_s)

  ## Mean lifetime tau:
  # N(t) = N(0) exp(-t/tau)
  ## Mean distance xMean:
  # N(x) = N(0) exp(-x/xMean)
  # where xMean for relativistic particle = beta gamma c tau

  # Energy is total energy; mass is rest mass in GeV/c^2
  # Gamma = total energy/rest energy, where rest energy = rest mass * c^2
  gamma = float(energy)/float(mass)
  beta = np.sqrt(1.-(1./gamma)**2)
  xMean = beta * gamma * constants.c * tau_s
  print("Velocity is",beta,"of c.")
  print("Beta-gamma is",beta*gamma)
  print("ctau is",constants.c * tau_s,"m")
  print("xMean is", xMean,"m")

  # For defining graphs
  xvals_distance = np.logspace(-6,2,200,endpoint=True)

  ## Now create decay distribution.
  def N(x) :
    return np.exp(-x/xMean)

  # Save TGraph
  yvals_totaldist = N(xvals_distance)
  graphs[tau_string]["total_dist"] = make_graph(xvals_distance,yvals_totaldist)

  ## Probability particle has decayed at location x
  N_total = integrate.quad(N, 0, 200)[0]
  def N_cumulative(x) :
    integral = integrate.quad(N,0,x)[0]
    return integral/N_total

  # Save
  yvals_decayprob = [N_cumulative(x) for x in xvals_distance]
  graphs[tau_string]["decayed_by_dist"] = make_graph(xvals_distance,yvals_decayprob)

  ## Probability particle has decayed by *transverse* distance x
  def N_transverse(r) :
    integrand = lambda z : sc.expi(- np.sqrt(r**2+z**2)/xMean)
    # Appropriate z range for integral depends on lifetime.
    z_end = 5*tau
    val = integrate.quad(integrand,-z_end,z_end)    
    return -val[0]

  # Save
  norm_val = N_transverse(xvals_distance[0])
  print "Normalisation value is",norm_val
  yvals_transversedist = [N_transverse(x)/norm_val for x in xvals_distance]
  print yvals_transversedist[0:20]
  graphs[tau_string]["transverse_dist"] = make_graph(xvals_distance,yvals_transversedist)
  

# Now make plots.
legend_list = ["#tau = {0} ns".format(i) for i in sorted(graphs.keys())]
for graph_name in graphs["1"].keys() :

  graph_list = [graphs[i][graph_name] for i in sorted(graphs.keys())]
  plot_graphs(graph_list,legend_list,graph_name)

# Comparison graph.
graphs_total = [graphs[i]["total_dist"] for i in sorted(graphs.keys())]
graphs_transverse = [graphs[i]["transverse_dist"] for i in sorted(graphs.keys())]
plot_graphs(zip(graphs_total,graphs_transverse),legend_list,"dist",tag="compare")



# Matplotlib example
## Plot 1: Number of particles versus total decay distance
#xvals = np.logspace(-6,2,200,endpoint=True)
#yvals_cumulative = [N_cumulative(x) for x in xvals]
#plt.figure()
#plt.plot(xvals, N(xvals), 'b', xvals, yvals_cumulative, 'r')
#plt.xlim(6e-6, 200)
#plt.xscale('log')
#plt.show()



