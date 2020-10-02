import ROOT

# Create colors
#color1 = ROOT.TColor().GetColor("#5fa55a")
color1 = ROOT.TColor().GetColor("#4C8448")
color2 = ROOT.TColor().GetColor("#01b4bc")
color3 = ROOT.TColor().GetColor("#f6d51f")
color4 = ROOT.TColor().GetColor("#fa8925")
#color5 = ROOT.TColor().GetColor("#fa5457")
color5 = ROOT.TColor().GetColor("#E14B4E")

def plot_graphs(graphs,legend_list,graphname,properties,chief_line_1="",chief_line_2="",extraline="",tag="",legend_low=True) :

  # Make a canvas to put the plot on.
  # Format it to have log x axis.
  c = ROOT.TCanvas("canvas",'',0,0,600,600)
  c.SetLogx(True)
  c.SetLogy(False)
  c.SetGridx(0)
  c.SetGridy(0)

  # Decide what colours to use.
  # These ones look decent, but obviously use
  # whatever you like best.
  #goodColours = [ROOT.kRed-7,ROOT.kAzure+1,ROOT.kTeal+1,ROOT.kOrange-4,ROOT.kViolet-5]
  goodColours = [color1,color2,color3,color4,color5]
  if len(graphs) < 4 :
    #goodColours = [ROOT.kRed-7,ROOT.kAzure+1,ROOT.kOrange-4]
    goodColours = [color2,color3,color5]

  leg_height = 0.04*len(graphs)
  leg_top = 0.2+leg_height if legend_low else 0.8
  leg_bottom = 0.2 if legend_low else 0.8-leg_height
  legend = ROOT.TLegend(0.2,leg_bottom,0.52,leg_top)
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
      graph.SetLineWidth(4)
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
  if extraline :
    myLatex.DrawLatex(0.2,leg_top+0.02,extraline)
    leg_top = leg_top+0.04
  myLatex.DrawLatex(0.2,leg_top+0.06,chief_line_1)
  myLatex.DrawLatex(0.2,leg_top+0.02,chief_line_2)
  #myLatex.DrawLatex(0.2,leg_top+0.06,"M: {0} GeV".format(mass))
  #myLatex.DrawLatex(0.2,leg_top+0.02,"E: {0} GeV".format(energy))  

  # Update the canvas
  c.Update()

  # Save the output
  savename = graphname
  if tag : savename = savename + "_"+tag
  c.SaveAs("plots/"+savename+".eps")

  # Save another without the log
  c.SetLogx(False)
  c.SaveAs("plots/"+savename+"_linear.eps")
