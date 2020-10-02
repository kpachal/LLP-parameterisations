import ROOT
import math

from plotter import plot_graphs

# Ref: https://www.pp.rhul.ac.uk/~cowan/stat/notes/medsigNote.pdf
def estimate_sensitivity(s, b, db) :

    num_term_1 = (s + b)*(b + db**2)
    denom_term_1 = b**2 + (s+b)*db**2
    term_1 = (s+b)*math.log(num_term_1/denom_term_1)

    arg_term2 = 1 + s*db**2/(b*(b + db**2))
    term_2 = b**2/db**2 * math.log(arg_term2)

    ZA = math.sqrt(2*(term_1 - term_2))
    return ZA

# Ref: https://root.cern.ch/doc/v608/rs__numbercountingutils_8C.html
def estimate_sensitivity_root(s, b, rel_db) :
  pval = ROOT.RooStats.NumberCountingUtils.BinomialExpP(s,b,rel_db)
  zval = ROOT.RooStats.PValueToSignificance(pval/2.)
  return zval

def observed_significance_root(obs,b,rel_db) :
    pval = ROOT.RooStats.NumberCountingUtils.BinomialObsP(obs,b,rel_db)
    zval = ROOT.RooStats.PValueToSignificance(pval/2.)
    return zval

# O signal events at 140 ifb:
print("0,0.1,0.1:", estimate_sensitivity(0,0.1,0.1))
print("using root estimate:", estimate_sensitivity_root(0,0.1,1))

## 1 signal event at 140 ifb:
print("1,0.1,0.1:", estimate_sensitivity(1,0.1,0.1))
print("using root estimate:", estimate_sensitivity_root(1,0.1,1))

## 2 signal events at 140 ifb:
print("2,0.1,0.1:", estimate_sensitivity(2,0.1,0.1))
print("using root estimate:", estimate_sensitivity_root(2,0.1,1))

## Let's assume a 10 percent impact on signal efficiency combined with
## s proportional to L otherwise.
## With 300 ifb:
s = 300./140. * 0.9 * 2
print("Correct s is",s)
print("3,0.1,0.1:", estimate_sensitivity(3.8,0.1,0.1))


L_start = 36
L_end = 300
s_start = 1

b_vals = [0.01,0.1,1]
rel_db_vals = [0.1,1,2]

# Store graph formatting
graph_formatting = {
  "sensitivity" : {"xrange": (30,320),
                  "yrange": (0,5),
                  "xlabel": "Luminosity [fb-1]",
                  "ylabel": "Significance [#sigma]"},
}

# By what relative efficiency should we assume signal efficiency decreases
# as we tighten cuts?
# This will be: each time we double luminosity this is the amount lost.
efficiency_hit = 0.1

graphs_dict = {}

for b in b_vals :
  graphs_dict[b] = {}
  for rel_db in rel_db_vals :
    graphs_dict[b][rel_db] = {}
    discovery_sig = ROOT.TGraph()
    exclusion_sig = ROOT.TGraph()
    for luminosity in range(L_start, L_end+1) :

      # TODO add efficiency degradation
      s = s_start * luminosity/L_start

      d_sig = estimate_sensitivity_root(s,b,rel_db)
      discovery_sig.SetPoint(discovery_sig.GetN(),luminosity,d_sig)

      d_excl = observed_significance_root(s+b,b,rel_db)
      exclusion_sig.SetPoint(exclusion_sig.GetN(),luminosity,d_sig)
    graphs_dict[b][rel_db]["discovery"] = discovery_sig
    graphs_dict[b][rel_db]["exclusion"] = exclusion_sig

  # Make plot comparing background effects
  graphs_disc = [graphs_dict[b][i]["discovery"] for i in rel_db_vals]
  plot_graphs(graphs_disc,["{0}".format(i) for i in rel_db_vals],"sensitivity",graph_formatting["sensitivity"],tag="b{0}".format(b),legend_low=False)


# Now do something to compare to
graphs_dict_standard = {}

# For 100 background events with 10% relative uncertainty
# a signal value of 16.7 events
# gives the same significance we saw for 1,0.1,1
b_start = 100
s_start = 16.7
# And this if we want to match 1,0.1,2
s_start = 8
#rel_db = 0.1
comparison_graph = ROOT.TGraph()
for luminosity in range(L_start, L_end+1) :

  # any changes? Unsure
  s = s_start * luminosity/L_start
  b = b_start * luminosity/L_start
  rel_db = math.sqrt(b)/b

  d_sig = estimate_sensitivity_root(s,b,rel_db)
  comparison_graph.SetPoint(comparison_graph.GetN(),luminosity,d_sig)

compare_list = [graphs_dict[0.1][2]["discovery"],comparison_graph]
plot_graphs(compare_list,["Zero-background analysis","High-background analysis"],"comparison",graph_formatting["sensitivity"],legend_low=False)







