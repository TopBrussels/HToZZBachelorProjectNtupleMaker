# HToZZBachelorProjectNtupleMaker

Very simple ntuple maker class.
Should compile out of the box using the "source compile.sh" command.
Ntupler takes TopTree xml files as single argument and produces a very simple root tuple
maked0plot.py is a very simple pyRoot analyser that books histograms and loops over the output tree.

use:
> source compile.sh
> ./Ntupler <myxmlfile.xml>
> python -i maked0plot.py #
