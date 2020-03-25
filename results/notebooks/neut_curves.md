# Neutralization assays on HA mutants

## Overview
Here we analyze neutralization assays of the serum against virus with the wildtype Perth/2009 HA and strongly selected mutants identified in the mutational antigenic profiling.

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

    Using neutcurve version 0.3.0


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
print(f"Reading neutralization assay setup from {config['neut_config']}")

with open(config['neut_config']) as f:
    neut_config = yaml.safe_load(f)
```

    Reading neutralization assay setup from data/neut_config.yaml


Get the output directory:


```python
outdir = config['neut_curves_outdir']
os.makedirs(outdir, exist_ok=True)
print(f"Output will be written to {outdir}")
```

    Output will be written to results/neut_curves


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

display(HTML(neutdata.head().to_html(index=False)))
```

    Parsing data for VIDD4...
    Parsing data for HC070021...
    Parsing data for HC150044...
    Parsing data for 2HC080043...
    Parsing data for HC070041...
    Parsing data for HC150036...
    Parsing data for HC080054...
    Parsing data for HC140028...
    Parsing data for HC060002...
    Parsing data for HC120043...
    Parsing data for ferret-Pitt2...
    Read data for 11 sera and 44 serum / virus pairs.



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>VIDD4</td>
      <td>wt</td>
      <td>1</td>
      <td>6.78e-06</td>
      <td>1</td>
    </tr>
    <tr>
      <td>VIDD4</td>
      <td>wt</td>
      <td>1</td>
      <td>1.36e-05</td>
      <td>1.01</td>
    </tr>
    <tr>
      <td>VIDD4</td>
      <td>wt</td>
      <td>1</td>
      <td>2.71e-05</td>
      <td>0.955</td>
    </tr>
    <tr>
      <td>VIDD4</td>
      <td>wt</td>
      <td>1</td>
      <td>5.43e-05</td>
      <td>0.919</td>
    </tr>
    <tr>
      <td>VIDD4</td>
      <td>wt</td>
      <td>1</td>
      <td>0.000109</td>
      <td>0.858</td>
    </tr>
  </tbody>
</table>


Now read in the serum info which includes names and groups for the sera, and add this to `neutdata` data frame:


```python
with open(config['serum_info']) as f:
    serum_info = (pd.DataFrame(yaml.safe_load(f))
                  .transpose()
                  .rename_axis('serum')
                  .reset_index()
                  )
    
neutdata = neutdata.merge(serum_info,
                          how='left',
                          on='serum',
                          validate='many_to_one')

display(HTML(neutdata.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>name</th>
      <th>description</th>
      <th>group</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>VIDD4</td>
      <td>wt</td>
      <td>1</td>
      <td>6.78e-06</td>
      <td>1</td>
      <td>age-64-Hutch</td>
      <td>collected at Hutch in 11/2008 from person born in 1945</td>
      <td>VIDD_sera</td>
    </tr>
    <tr>
      <td>VIDD4</td>
      <td>wt</td>
      <td>1</td>
      <td>1.36e-05</td>
      <td>1.01</td>
      <td>age-64-Hutch</td>
      <td>collected at Hutch in 11/2008 from person born in 1945</td>
      <td>VIDD_sera</td>
    </tr>
    <tr>
      <td>VIDD4</td>
      <td>wt</td>
      <td>1</td>
      <td>2.71e-05</td>
      <td>0.955</td>
      <td>age-64-Hutch</td>
      <td>collected at Hutch in 11/2008 from person born in 1945</td>
      <td>VIDD_sera</td>
    </tr>
    <tr>
      <td>VIDD4</td>
      <td>wt</td>
      <td>1</td>
      <td>5.43e-05</td>
      <td>0.919</td>
      <td>age-64-Hutch</td>
      <td>collected at Hutch in 11/2008 from person born in 1945</td>
      <td>VIDD_sera</td>
    </tr>
    <tr>
      <td>VIDD4</td>
      <td>wt</td>
      <td>1</td>
      <td>0.000109</td>
      <td>0.858</td>
      <td>age-64-Hutch</td>
      <td>collected at Hutch in 11/2008 from person born in 1945</td>
      <td>VIDD_sera</td>
    </tr>
  </tbody>
</table>


We write the neutralization data to a CSV file in our output directory:


```python
neutdatafile = os.path.join(outdir, 'neutdata.csv')
neutdata.to_csv(neutdatafile, index=False)
print(f"Wrote neutralization data to {neutdatafile}")
```

    Wrote neutralization data to results/neut_curves/neutdata.csv


## Fit and plot neutralization curves

Now we fit the neutralization curves with a [neutcurve.CurveFits](https://jbloomlab.github.io/neutcurve/neutcurve.curvefits.html#neutcurve.curvefits.CurveFits):


```python
fits = neutcurve.CurveFits(neutdata,
                           serum_col='name')
```

Make plots that show the curves for all replicates for each serum:


```python
fig_reps, _ = fits.plotReplicates(legendtitle='replicate',
                                  xlabel='serum dilution')
```


![png](neut_curves_files/neut_curves_23_0.png)


Now make nice plots that show all curves for each mutant:


```python
fig_sera, _ = fits.plotSera(xlabel='serum dilution')
```


![png](neut_curves_files/neut_curves_25_0.png)


Save a PDF of the above plot:


```python
plotfile = os.path.join(outdir, 'all_curves_by_sera.pdf')
print(f"Creating plot {plotfile}")
fig_sera.savefig(plotfile)
```

    Creating plot results/neut_curves/all_curves_by_sera.pdf


Display fit parameters:


```python
display(HTML(fits.fitParams().to_html()))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>age-64-Hutch</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>0.000236</td>
      <td>interpolated</td>
      <td>0.000236</td>
      <td>0.000236</td>
      <td>1.92</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>age-64-Hutch</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>0.000244</td>
      <td>interpolated</td>
      <td>0.000244</td>
      <td>0.000244</td>
      <td>1.96</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>age-64-Hutch</td>
      <td>F159G</td>
      <td>average</td>
      <td>3</td>
      <td>0.00602</td>
      <td>interpolated</td>
      <td>0.00602</td>
      <td>0.00602</td>
      <td>1.47</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>age-64-Hutch</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.000301</td>
      <td>interpolated</td>
      <td>0.000301</td>
      <td>0.000301</td>
      <td>1.87</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>age-2.4</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>6.69e-05</td>
      <td>interpolated</td>
      <td>6.69e-05</td>
      <td>6.69e-05</td>
      <td>1.56</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>age-2.4</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>7.16e-05</td>
      <td>interpolated</td>
      <td>7.16e-05</td>
      <td>7.16e-05</td>
      <td>1.79</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>age-2.4</td>
      <td>F159G</td>
      <td>average</td>
      <td>3</td>
      <td>0.00278</td>
      <td>lower</td>
      <td>&gt;0.00278</td>
      <td>0.00308</td>
      <td>1.78</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>age-2.4</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.000748</td>
      <td>interpolated</td>
      <td>0.000748</td>
      <td>0.000748</td>
      <td>1.86</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>age-3.4</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>0.000172</td>
      <td>interpolated</td>
      <td>0.000172</td>
      <td>0.000172</td>
      <td>1.77</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>9</th>
      <td>age-3.4</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>0.000149</td>
      <td>interpolated</td>
      <td>0.000149</td>
      <td>0.000149</td>
      <td>1.81</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>10</th>
      <td>age-3.4</td>
      <td>F159G</td>
      <td>average</td>
      <td>3</td>
      <td>0.00556</td>
      <td>lower</td>
      <td>&gt;0.00556</td>
      <td>0.00667</td>
      <td>1.47</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>11</th>
      <td>age-3.4</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.00113</td>
      <td>interpolated</td>
      <td>0.00113</td>
      <td>0.00113</td>
      <td>1.72</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>12</th>
      <td>age-3.3</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>0.00014</td>
      <td>interpolated</td>
      <td>0.00014</td>
      <td>0.00014</td>
      <td>1.1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>13</th>
      <td>age-3.3</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>0.000195</td>
      <td>interpolated</td>
      <td>0.000195</td>
      <td>0.000195</td>
      <td>1.07</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>14</th>
      <td>age-3.3</td>
      <td>F159G</td>
      <td>average</td>
      <td>3</td>
      <td>0.0121</td>
      <td>interpolated</td>
      <td>0.0121</td>
      <td>0.0121</td>
      <td>0.555</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>15</th>
      <td>age-3.3</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.000635</td>
      <td>interpolated</td>
      <td>0.000635</td>
      <td>0.000635</td>
      <td>0.458</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>16</th>
      <td>age-2.1</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>0.000193</td>
      <td>interpolated</td>
      <td>0.000193</td>
      <td>0.000193</td>
      <td>2.9</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>17</th>
      <td>age-2.1</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>0.00018</td>
      <td>interpolated</td>
      <td>0.00018</td>
      <td>0.00018</td>
      <td>2.76</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>18</th>
      <td>age-2.1</td>
      <td>K189D</td>
      <td>average</td>
      <td>3</td>
      <td>0.000359</td>
      <td>interpolated</td>
      <td>0.000359</td>
      <td>0.000359</td>
      <td>2.41</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>19</th>
      <td>age-2.1</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.000294</td>
      <td>interpolated</td>
      <td>0.000294</td>
      <td>0.000294</td>
      <td>2.2</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>20</th>
      <td>age-2.2</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>0.000237</td>
      <td>interpolated</td>
      <td>0.000237</td>
      <td>0.000237</td>
      <td>2.33</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>21</th>
      <td>age-2.2</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>0.00023</td>
      <td>interpolated</td>
      <td>0.00023</td>
      <td>0.00023</td>
      <td>2.3</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>22</th>
      <td>age-2.2</td>
      <td>K189D</td>
      <td>average</td>
      <td>3</td>
      <td>0.000609</td>
      <td>interpolated</td>
      <td>0.000609</td>
      <td>0.000609</td>
      <td>2.01</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>23</th>
      <td>age-2.2</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.000574</td>
      <td>interpolated</td>
      <td>0.000574</td>
      <td>0.000574</td>
      <td>1.9</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>24</th>
      <td>age-33.5</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>0.000806</td>
      <td>interpolated</td>
      <td>0.000806</td>
      <td>0.000806</td>
      <td>1.76</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>25</th>
      <td>age-33.5</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>0.000716</td>
      <td>interpolated</td>
      <td>0.000716</td>
      <td>0.000716</td>
      <td>1.86</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>26</th>
      <td>age-33.5</td>
      <td>K189D</td>
      <td>average</td>
      <td>3</td>
      <td>0.00193</td>
      <td>interpolated</td>
      <td>0.00193</td>
      <td>0.00193</td>
      <td>1.47</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>27</th>
      <td>age-33.5</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.0012</td>
      <td>interpolated</td>
      <td>0.0012</td>
      <td>0.0012</td>
      <td>1.73</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>28</th>
      <td>age-2.5-b</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>0.000237</td>
      <td>interpolated</td>
      <td>0.000237</td>
      <td>0.000237</td>
      <td>1.03</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>29</th>
      <td>age-2.5-b</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>0.000508</td>
      <td>interpolated</td>
      <td>0.000508</td>
      <td>0.000508</td>
      <td>3.84</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>30</th>
      <td>age-2.5-b</td>
      <td>K189D</td>
      <td>average</td>
      <td>3</td>
      <td>0.00115</td>
      <td>interpolated</td>
      <td>0.00115</td>
      <td>0.00115</td>
      <td>0.688</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>31</th>
      <td>age-2.5-b</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.00103</td>
      <td>interpolated</td>
      <td>0.00103</td>
      <td>0.00103</td>
      <td>2.07</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>32</th>
      <td>age-3.5</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>0.000321</td>
      <td>interpolated</td>
      <td>0.000321</td>
      <td>0.000321</td>
      <td>2.32</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>33</th>
      <td>age-3.5</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>0.000257</td>
      <td>interpolated</td>
      <td>0.000257</td>
      <td>0.000257</td>
      <td>2.12</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>34</th>
      <td>age-3.5</td>
      <td>K189D</td>
      <td>average</td>
      <td>3</td>
      <td>0.00108</td>
      <td>interpolated</td>
      <td>0.00108</td>
      <td>0.00108</td>
      <td>1.79</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>35</th>
      <td>age-3.5</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.00113</td>
      <td>interpolated</td>
      <td>0.00113</td>
      <td>0.00113</td>
      <td>1.6</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>36</th>
      <td>age-2.5</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>0.00052</td>
      <td>interpolated</td>
      <td>0.00052</td>
      <td>0.00052</td>
      <td>1.57</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>37</th>
      <td>age-2.5</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>0.000387</td>
      <td>interpolated</td>
      <td>0.000387</td>
      <td>0.000387</td>
      <td>1.36</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>38</th>
      <td>age-2.5</td>
      <td>K189D</td>
      <td>average</td>
      <td>3</td>
      <td>0.00463</td>
      <td>lower</td>
      <td>&gt;0.00463</td>
      <td>0.00719</td>
      <td>1.88</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>39</th>
      <td>age-2.5</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.00353</td>
      <td>interpolated</td>
      <td>0.00353</td>
      <td>0.00353</td>
      <td>2.11</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>40</th>
      <td>ferret-Pitt2</td>
      <td>wt</td>
      <td>average</td>
      <td>3</td>
      <td>0.000178</td>
      <td>interpolated</td>
      <td>0.000178</td>
      <td>0.000178</td>
      <td>1.64</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>41</th>
      <td>ferret-Pitt2</td>
      <td>F193F</td>
      <td>average</td>
      <td>3</td>
      <td>0.000129</td>
      <td>interpolated</td>
      <td>0.000129</td>
      <td>0.000129</td>
      <td>1.29</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>42</th>
      <td>ferret-Pitt2</td>
      <td>K189D</td>
      <td>average</td>
      <td>3</td>
      <td>0.000508</td>
      <td>interpolated</td>
      <td>0.000508</td>
      <td>0.000508</td>
      <td>1.51</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>43</th>
      <td>ferret-Pitt2</td>
      <td>F193D</td>
      <td>average</td>
      <td>3</td>
      <td>0.000692</td>
      <td>interpolated</td>
      <td>0.000692</td>
      <td>0.000692</td>
      <td>1.73</td>
      <td>1</td>
      <td>0</td>
    </tr>
  </tbody>
</table>

