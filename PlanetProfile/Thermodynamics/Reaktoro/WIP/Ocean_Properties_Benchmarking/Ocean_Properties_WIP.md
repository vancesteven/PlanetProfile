# Ocean Properties WIP

# Properties Benchmark

We benchmark the Reaktoro-simulated ocean properties against the existing laboratory-measured ocean properties currently integrated within Reaktoro.

Importantly, we apply a Seafreeze $H_2O$ correction to Reaktoro's simulated properties to align Reaktoro's and Seafreeze's Pure $H_2O$ data, according to the equation S1:

$S1: Property_{corrected(P,T)}=Property_{solution(P,T)}-Property_{ReaktoroH_2O(P,T)}+Property_{SeafreezeH_2O(P,T)}$

For more details, see the manuscript and [Seafreeze Pure $H_2O$ section](#seafreeze-pure-h_2o)
## Seafreeze Pure $H_2O$
The Seafreeze Pure $H_2O$ is a Gibbs parameterization of pure $H_2O$ and ice polymorphs that computes the thermodynamic and elastic properties up to 100 GPa and 10,000K.
The SeaFreeze package is available at [Github link](https://github.com/Bjournaux/SeaFreeze).

Below, we benchmark Reaktoro's Pure $H_2O$ calculations to Seafreeze's Pure $H_2O$ data. Subsequently, we apply a correction to achieve agreement between Reaktoro's and Seafreeze's Pure $H_2O$ properties, which we apply to Reaktoro's aqueous solution calculations as outlined by equation S1.
### Seafreeze Pure $H_2O$ benchmark without Seafreeze correction:


### Seafreeze Pure $H_2O$ with Seafreeze correction:


## Seafreeze NaCl
## Gibb's Seawater
## Vance et al.'s $MgSO_4$
