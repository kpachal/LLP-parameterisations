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
                  "ylabel": "Relative decay probability"},
  "decayed_by_dist" : {"xrange": (5e-5,20),
                  "yrange": (0,1),
                  "xlabel": "Total decay displacement [m]",
                  "ylabel": "Probability particle has decayed"},
  "transverse_dist" : {"xrange": (5e-5,20),
                  "yrange": (0,1),
                  "xlabel": "Transverse decay displacement [m]",
                  "ylabel": "Relative decay probability"},
  "dist" : {"xrange": (5e-5,20),
                  "yrange": (0,1),
                  "xlabel": "Decay displacement [m]",
                  "ylabel": "Relative decay probability"},
}

# Everything in m
detector_info = {
  "Inner detector" : {"abs_z1" : 0.3315, 
                      "abs_z2" : 2.72,
                      "r1" : 0.033,
                      "r2" : 1.066  },
  "Calorimeters" : {
                    # Tile cal: r = 2.3 to 3.9
                    # total z length: 12m, so abs_z1 = abs_z2 = 6 m
                    # LAr cal: r = 1.15 to 2.25
                    # Approximate whole as LAr inner radius to Tile outer radius,
                    # and 6 m long to include endcaps? -> yes
                    "abs_z1" : 6.0,
                    'abs_z2' : 6.0,
                    "r1" : 1.15,
                    "r2" : 3.9 },
  "Muon spectrometer" : {
                    "abs_z1" : 6.,
                    "abs_z2" : 13.,
                    "r1" : 5.,
                    "r2": 10.,
  }
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
detector_data = {}

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

  # For transverse distances
  def integrand_z(z, r, dummy1, dummy2, dummy3, dummy4) :
    return sc.expi(- np.sqrt(r**2+z**2)/xMean)

  ## Probability particle has decayed by *transverse* distance x
  def N_transverse(r) :
    # Appropriate z range for integral depends on lifetime.
    z_end = 5*tau
    val = integrate.quad(integrand_z,-z_end,z_end,args=(r, 0,0,0,0))
    return -val[0]

  # Save
  norm_val = N_transverse(xvals_distance[0])
  yvals_transversedist = [N_transverse(x)/norm_val for x in xvals_distance]
  graphs[tau_string]["transverse_dist"] = make_graph(xvals_distance,yvals_transversedist)

  def lim_z(r, abs_z1, abs_z2, r1, r2) :
    # z limits for integration come from a straight line between z1, r1, and z2, r2 at r.
    # Use python's linear interpolation.
    z = np.interp(r,[r1,r2],[abs_z1,abs_z2])
    return [-z,z]


  def N_subdetector(abs_z1, abs_z2, r1, r2) :
    args = (abs_z1, abs_z2, r1, r2)
    val = integrate.nquad(integrand_z,[lim_z,[r1,r2]], args=args)
    return -val[0]

  # Total in 2D is integral of N_transverse out to appropriate r.
  # Let's start from value of first encountering detector.
  N_transverse_total = N_subdetector(5*tau,5*tau,0,5*tau)

  test = N_subdetector(1,1,1,2)
  print "Frac in subdector from -1 to 1 and 1 to 2:",test/N_transverse_total
  print "Total:",N_transverse_total
  print "N in here:",test

  print "Test limit. Interpolate from z=1 at r=1 to z=3 at r=2, evaluate for r=1.5"
  print "z_lims are",lim_z(1.5,1,3,1,2)

  # Calculation per subdetector
  for detector in detector_info.keys() :

    dimensions = detector_info[detector]
    #n_particles = N_subdetector()



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



