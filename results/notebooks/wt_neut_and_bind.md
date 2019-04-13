
<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Serum-neutralization-and-binding-to-wildtype-HA" data-toc-modified-id="Serum-neutralization-and-binding-to-wildtype-HA-1">Serum neutralization and binding to wildtype HA</a></span><ul class="toc-item"><li><span><a href="#Overview" data-toc-modified-id="Overview-1.1">Overview</a></span></li><li><span><a href="#Import-Python-packages" data-toc-modified-id="Import-Python-packages-1.2">Import Python packages</a></span></li><li><span><a href="#Configuration-and-setup" data-toc-modified-id="Configuration-and-setup-1.3">Configuration and setup</a></span></li><li><span><a href="#Read-neutralization-data" data-toc-modified-id="Read-neutralization-data-1.4">Read neutralization data</a></span></li><li><span><a href="#Fit-and-plot-neutralization-curves" data-toc-modified-id="Fit-and-plot-neutralization-curves-1.5">Fit and plot neutralization curves</a></span></li><li><span><a href="#Compare-neutralization-to-binding-data" data-toc-modified-id="Compare-neutralization-to-binding-data-1.6">Compare neutralization to binding data</a></span></li><li><span><a href="#Examine-IC50-/-IC95-for-mutational-antigenic-profiling" data-toc-modified-id="Examine-IC50-/-IC95-for-mutational-antigenic-profiling-1.7">Examine IC50 / IC95 for mutational antigenic profiling</a></span></li></ul></li></ul></div>

# Serum neutralization and binding to wildtype HA

## Overview
Here we analyze neutralization assays of the serum against virus with the wildtype Perth/2009 HA.
We then compare the IC50s to the binding titers measured as part of the [Nhat et al (2017)](https://www.nature.com/articles/s41598-017-06177-0) study that characterized these IC50s.

The neutralization assays where set up using the "Rachel-style 2019" format on the Bloom lab plate reader as [described here](https://jbloomlab.github.io/neutcurve/rachelstyle2019_example.html), so we can analyze the raw Excel data off the plate reader.

## Import Python packages
We use the Bloom lab [neutcurve](https://jbloomlab.github.io/neutcurve) package for fitting neutralization curves, and [plotnine](https://plotnine.readthedocs.io) for ggplot2-like plotting:


```python
import os
import warnings

from IPython.display import display, HTML

import numpy

import pandas as pd

from plotnine import *

import yaml

import neutcurve
from neutcurve.colorschemes import CBPALETTE
import neutcurve.parse_excel

print(f"Using neutcurve version {neutcurve.__version__}")
```

    Using neutcurve version 0.1.dev0


Set output format of pandas Data Frames:


```python
pd.set_option('display.float_format', '{:.3g}'.format)
```

Hide warnings that can clutter output:


```python
warnings.simplefilter('ignore')
```

## Configuration and setup
Read general configuration from [config.yaml](config.yaml):


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Read the neutralization assay configuration from the specified file:


```python
print(f"Reading neutralization assay setup from {config['wt_neut_config']}")

with open(config['wt_neut_config']) as f:
    neut_config = yaml.safe_load(f)
```

    Reading neutralization assay setup from data/wt_neut_config.yaml


Get the output directory:


```python
outdir = config['wt_neut_and_bind_outdir']
os.makedirs(outdir, exist_ok=True)
print(f"Output will be written to {outdir}")
```

    Output will be written to results/wt_neut_and_bind


## Read neutralization data
Next, for each dict in *neut_config*, we use
[neutcurve.parse_excel.parseRachelStyle2019](https://jbloomlab.github.io/neutcurve/neutcurve.parse_excel.html#neutcurve.parse_excel.parseRachelStyle2019) to create a tidy
data frame appropriate for passing to
[neutcurve.CurveFits](https://jbloomlab.github.io/neutcurve/neutcurve.curvefits.html#neutcurve.curvefits.CurveFits).
We then concatenate all the
tidy data frames to get our neutralization data:


```python
neutdata = []  # store all data frame, then concatenate at end

for sampledict in neut_config:
    assert len(sampledict) == 1
    sampleset, kwargs = list(sampledict.items())[0]
    print(f"Parsing data for {sampleset}...")
    neutdata.append(neutcurve.parse_excel.parseRachelStyle2019(**kwargs))

neutdata = pd.concat(neutdata)
print(f"Read data for {len(neutdata.groupby('serum'))} sera and "
      f"{len(neutdata.groupby(['serum', 'virus']))} serum / virus pairs.")
      
assert neutdata['virus'].unique() == ['wt'], 'expected data only for wildtype virus'
```

    Parsing data for 19Feb19...
    Parsing data for 22Feb19...
    Parsing data for 25Feb19...
    Parsing data for 26Feb19...
    Parsing data for 27Feb19...
    Parsing data for 04Mar19...
    Read data for 48 sera and 48 serum / virus pairs.


We write the neutralization data to a CSV file in our output directory:


```python
neutdatafile = os.path.join(outdir, 'neutdata.csv')
neutdata.to_csv(neutdatafile, index=False)
print(f"Wrote neutralization data to {neutdatafile}")
```

    Wrote neutralization data to results/wt_neut_and_bind/neutdata.csv


## Fit and plot neutralization curves

Now we fit the neutralization curves with a [neutcurve.CurveFits](https://jbloomlab.github.io/neutcurve/neutcurve.curvefits.html#neutcurve.curvefits.CurveFits):


```python
fits = neutcurve.CurveFits(neutdata)
```

Make plots that show the curves for all replicates for each serum:


```python
fig_reps, _ = fits.plotReplicates(legendtitle='replicate',
                                  xlabel='serum dilution')
```


![png](wt_neut_and_bind_files/wt_neut_and_bind_22_0.png)


Make plots that just show the across-replicate averages for each serum /
virus combination:


```python
fig_avgs, _ = fits.plotAverages(xlabel='serum dilution')
```


![png](wt_neut_and_bind_files/wt_neut_and_bind_24_0.png)


Save the above neutralization curve plots to PDF files:


```python
for plotname, fig in [('replicates_plot', fig_reps),
                      ('average_plot', fig_avgs)]:
    plotfile = os.path.join(outdir, plotname + '.pdf')
    print(f"Creating plot {plotfile}")
    fig.savefig(plotfile)
```

    Creating plot results/wt_neut_and_bind/replicates_plot.pdf
    Creating plot results/wt_neut_and_bind/average_plot.pdf


## Compare neutralization to binding data
Now we want to compare the IC50s from the neutralization curves to the binding titers.
First, we read in the serum information that has the binding titers:


```python
print(f"Reading serum information from {config['serum_info']}")
sera = pd.read_csv(config['serum_info'])
print(f"Read information for {len(sera)} sera.")
```

    Reading serum information from data/serum_info.csv
    Read information for 48 sera.


Now we get the neutralization curve fit parameters (including the IC50 as well as the IC95), make sure we have serum information for all sera for which we have curves, and then merge them into a single data frame:


```python
fitparams = fits.fitParams(ics=[50, 95])

assert set(fitparams['serum']) <= set(sera['serum'])

neut_and_bind = sera.merge(fitparams, on='serum')
```

Now we plot the correlation between the neutralization titer (reciprocal IC50) and each of the binding titers.
We color points differently depending on whether they are child (<10 years) or adult ($\ge$10 years), and set their shapes based on whether the neutralization titer (IC50) is interpolated (in the range of the data) or simply an upper / lower bound:


```python
ic50_vs_bind_plot = (
    neut_and_bind
    .melt(id_vars=['age', 'ic50', 'ic50_bound'],
          value_vars=['H3_2011_binding', 'H1_2009_binding'],
          var_name='bind_type',
          value_name='binding titer'
          )
    .assign(bind_type=lambda x: x['bind_type'].str.replace('_', ' '),
            neut_titer=lambda x: 1 / x['ic50'],
            neut_titer_type=lambda x: x['ic50_bound']
                                      .map({'interpolated': 'interpolated',
                                            'upper': 'lower bound',
                                            'lower': 'upper bound'}),
            age=lambda x: numpy.where(x['age'] < 10, 'child', 'adult')
            )
    .pipe(ggplot, aes('neut_titer', 'binding titer',
                      shape='neut_titer_type')) +
          geom_point(aes(color='age'), size=3, alpha=0.6, stroke=0) +
          geom_smooth(method='lm', size=0.3) +
          scale_y_log10() +
          scale_x_log10(name='neutralization titer') +
          scale_color_manual(values=CBPALETTE[1: ]) +
          facet_wrap('~ bind_type') +
          theme_bw() +
          theme(figure_size=(4, 2),
                legend_key=element_blank()) +
          guides(shape=guide_legend(title='neutralization titer'))
          
    )

_ = ic50_vs_bind_plot.draw()

ic50_vs_bind_plotfile = os.path.join(outdir, 'ic50_vs_bind.pdf')
print(f"Saving plot to {ic50_vs_bind_plotfile}")
ic50_vs_bind_plot.save(ic50_vs_bind_plotfile)
```

    Saving plot to results/wt_neut_and_bind/ic50_vs_bind.pdf



![png](wt_neut_and_bind_files/wt_neut_and_bind_32_1.png)


Overall, the results in the above plot make sense:
  - There is little relationship between neutralization titer on the H3 Perth/2009 HA and the H1 binding.
  - There is a correlation between H3 binding and H3 neutralization.
  - All the sera that neutralizes binds well, as expected.
  - Some sera bind but do not neutralize well; these sera may bind non-neutralizing epitopes.
  - It appears that there is a stronger association of high binding with high neutralization for children. This
    This is consistent with the basic idea suggested in [Ranjeva et al (2019)](https://www.nature.com/articles/s41467-019-09652-6) that sera from children tends to target HI / neutralizing epitopes, whereas for adults this isn't always the case.

## Examine IC50 / IC95 for mutational antigenic profiling
Below we tabulate the serum by IC95 to see which ones might have appropriate concentrations for mutational antigenic profiling:


```python
display(HTML(neut_and_bind
             .sort_values(['ic95', 'ic50'])
             [['serum', 'ic95_str', 'ic50_str', 'age', 'H3_2011_binding',
               'H1_2009_binding']]
             .to_html(index=False)
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>ic95_str</th>
      <th>ic50_str</th>
      <th>age</th>
      <th>H3_2011_binding</th>
      <th>H1_2009_binding</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>HC070097</td>
      <td>0.000401</td>
      <td>7.14e-05</td>
      <td>31.5</td>
      <td>1.81e+03</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC070021</td>
      <td>0.000478</td>
      <td>0.000109</td>
      <td>2.4</td>
      <td>1.81e+03</td>
      <td>22.5</td>
    </tr>
    <tr>
      <td>HC070041</td>
      <td>0.000849</td>
      <td>0.000396</td>
      <td>2.1</td>
      <td>1.81e+03</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC150036</td>
      <td>0.0012</td>
      <td>0.000388</td>
      <td>2.2</td>
      <td>1.81e+03</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC150044</td>
      <td>0.00121</td>
      <td>0.000341</td>
      <td>3.4</td>
      <td>1.81e+03</td>
      <td>21.3</td>
    </tr>
    <tr>
      <td>HC060002</td>
      <td>0.00178</td>
      <td>0.000562</td>
      <td>3.5</td>
      <td>773</td>
      <td>10</td>
    </tr>
    <tr>
      <td>HC140028</td>
      <td>0.00188</td>
      <td>0.00105</td>
      <td>2.5</td>
      <td>1.18e+03</td>
      <td>605</td>
    </tr>
    <tr>
      <td>HC120043</td>
      <td>0.00256</td>
      <td>0.000573</td>
      <td>2.5</td>
      <td>1.81e+03</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC080054</td>
      <td>0.00324</td>
      <td>0.000503</td>
      <td>33.5</td>
      <td>1.81e+03</td>
      <td>528</td>
    </tr>
    <tr>
      <td>HC150108</td>
      <td>0.00325</td>
      <td>0.000625</td>
      <td>30.5</td>
      <td>1.22e+03</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC080048</td>
      <td>0.00473</td>
      <td>0.0016</td>
      <td>6.3</td>
      <td>514</td>
      <td>10</td>
    </tr>
    <tr>
      <td>HC110074</td>
      <td>0.00589</td>
      <td>0.00176</td>
      <td>38.8</td>
      <td>495</td>
      <td>90.6</td>
    </tr>
    <tr>
      <td>HC080043</td>
      <td>0.00605</td>
      <td>0.000684</td>
      <td>3.3</td>
      <td>982</td>
      <td>57.2</td>
    </tr>
    <tr>
      <td>HC150099</td>
      <td>0.00614</td>
      <td>0.0011</td>
      <td>30.5</td>
      <td>486</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC070072</td>
      <td>0.00617</td>
      <td>0.00118</td>
      <td>3.3</td>
      <td>1.81e+03</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC140010</td>
      <td>0.00629</td>
      <td>0.00262</td>
      <td>2.1</td>
      <td>830</td>
      <td>813</td>
    </tr>
    <tr>
      <td>HC080004</td>
      <td>0.00731</td>
      <td>0.00129</td>
      <td>3.5</td>
      <td>1.06e+03</td>
      <td>10</td>
    </tr>
    <tr>
      <td>HC050031</td>
      <td>0.00736</td>
      <td>0.00223</td>
      <td>5.5</td>
      <td>1.1e+03</td>
      <td>703</td>
    </tr>
    <tr>
      <td>HC150111</td>
      <td>0.00841</td>
      <td>0.00224</td>
      <td>34.5</td>
      <td>821</td>
      <td>92.9</td>
    </tr>
    <tr>
      <td>HC060077</td>
      <td>0.00864</td>
      <td>0.00174</td>
      <td>4.5</td>
      <td>1.81e+03</td>
      <td>50.4</td>
    </tr>
    <tr>
      <td>HC080059</td>
      <td>&gt;0.00926</td>
      <td>0.0019</td>
      <td>38.5</td>
      <td>1.81e+03</td>
      <td>511</td>
    </tr>
    <tr>
      <td>HC150140</td>
      <td>&gt;0.00926</td>
      <td>0.00217</td>
      <td>37.5</td>
      <td>1.81e+03</td>
      <td>805</td>
    </tr>
    <tr>
      <td>HC110090</td>
      <td>&gt;0.00926</td>
      <td>0.00219</td>
      <td>34.5</td>
      <td>1.13e+03</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC090064</td>
      <td>&gt;0.00926</td>
      <td>0.00274</td>
      <td>39.5</td>
      <td>1.81e+03</td>
      <td>793</td>
    </tr>
    <tr>
      <td>HC060106</td>
      <td>&gt;0.00926</td>
      <td>0.00395</td>
      <td>3.5</td>
      <td>1.81e+03</td>
      <td>577</td>
    </tr>
    <tr>
      <td>HC090070</td>
      <td>&gt;0.00926</td>
      <td>0.00431</td>
      <td>38.5</td>
      <td>1.81e+03</td>
      <td>1.2e+03</td>
    </tr>
    <tr>
      <td>HC100072</td>
      <td>&gt;0.00926</td>
      <td>0.00453</td>
      <td>26.6</td>
      <td>737</td>
      <td>1.31e+03</td>
    </tr>
    <tr>
      <td>HC040090</td>
      <td>&gt;0.00926</td>
      <td>0.00626</td>
      <td>26.5</td>
      <td>658</td>
      <td>831</td>
    </tr>
    <tr>
      <td>HC130149</td>
      <td>&gt;0.00926</td>
      <td>0.00682</td>
      <td>22.5</td>
      <td>648</td>
      <td>609</td>
    </tr>
    <tr>
      <td>HC100080</td>
      <td>&gt;0.00926</td>
      <td>0.00825</td>
      <td>38.5</td>
      <td>1.11e+03</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC150124</td>
      <td>&gt;0.00926</td>
      <td>0.00925</td>
      <td>20.5</td>
      <td>579</td>
      <td>1.37e+03</td>
    </tr>
    <tr>
      <td>HC070063</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>2.2</td>
      <td>10</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC070009</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>2.3</td>
      <td>21.4</td>
      <td>10</td>
    </tr>
    <tr>
      <td>HC150021</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>2.3</td>
      <td>10</td>
      <td>812</td>
    </tr>
    <tr>
      <td>HC060062</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>2.5</td>
      <td>10</td>
      <td>824</td>
    </tr>
    <tr>
      <td>HC150047</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>2.5</td>
      <td>10</td>
      <td>10</td>
    </tr>
    <tr>
      <td>HC070005</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>2.5</td>
      <td>44.5</td>
      <td>1.01e+03</td>
    </tr>
    <tr>
      <td>HC150080</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>2.5</td>
      <td>10</td>
      <td>1.22e+03</td>
    </tr>
    <tr>
      <td>HC070047</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>2.9</td>
      <td>10</td>
      <td>10</td>
    </tr>
    <tr>
      <td>HC110020</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>3.1</td>
      <td>10</td>
      <td>1.81e+03</td>
    </tr>
    <tr>
      <td>HC060087</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>3.5</td>
      <td>10</td>
      <td>702</td>
    </tr>
    <tr>
      <td>HC130135</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>21.2</td>
      <td>1.81e+03</td>
      <td>1.21e+03</td>
    </tr>
    <tr>
      <td>HC040121</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>26.5</td>
      <td>605</td>
      <td>531</td>
    </tr>
    <tr>
      <td>HC080077</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>30.5</td>
      <td>1.22e+03</td>
      <td>1.05e+03</td>
    </tr>
    <tr>
      <td>HC130143</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>30.5</td>
      <td>1.1e+03</td>
      <td>901</td>
    </tr>
    <tr>
      <td>HC080061</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>30.5</td>
      <td>1.81e+03</td>
      <td>509</td>
    </tr>
    <tr>
      <td>HC120108</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>31</td>
      <td>1.18e+03</td>
      <td>10</td>
    </tr>
    <tr>
      <td>HC100062</td>
      <td>&gt;0.00926</td>
      <td>&gt;0.00926</td>
      <td>32.5</td>
      <td>656</td>
      <td>548</td>
    </tr>
  </tbody>
</table>


We also write all of the data on the curve fits and binding to a file:


```python
neut_and_bind_file = os.path.join(outdir, 'neut_and_bind.csv')
print(f"Writing neutralization and binding data to {neut_and_bind_file}")
neut_and_bind.to_csv(neut_and_bind_file)
```

    Writing neutralization and binding data to results/wt_neut_and_bind/neut_and_bind.csv



```python

```
