### In Silico Experiment on Molecular Communication Systems

#### Overview
This README documents an in-silico experiment conducted on a molecular communication system, inspired by [Basu05]'s molecular communication framework. This experiment, part of an academic project under the Department of Basic Informatics at the Graduate School of Information Science, aims to computationally replicate and understand complex biological interactions and communications at a molecular level.

#### Experiment Context
Molecular communication involves the use of molecules to relay information in biological systems. This concept is critical in understanding natural phenomena like Quorum Sensing, where organisms communicate based on molecular concentration, influencing gene expression and behavior in a population.

The experiment focuses on AHL (N-acyl homoserine lactone), a molecule involved in bacterial communication. By varying AHL concentrations, the experiment observes the response in gene expression, particularly the production of GFP (Green Fluorescent Protein), a commonly used marker in molecular biology.

#### Mathematical Models
- **Hill Function**: To model the gene expression dynamics, the Hill function is utilized, reflecting how protein production rates depend on the concentration of certain molecules.
- **Decay Model**: Accounts for the natural decay of molecules over time, using parameters like initial concentration and decay rate.
- **Steady State Model**: Calculates the steady-state concentrations of various molecular species in the system under given conditions.

#### Simulation Components
- `GFPExpressionCalculator`: Calculates GFP production based on AHL concentration.
- `AHLConcentrationCalculator`: Estimates AHL concentration at different locations from the source, assuming exponential decay with distance.
- `GFPHighConcentrationLocator`: Identifies regions with high GFP concentration, which can be critical for understanding pattern formation in molecular communication.

#### Experimental Setup
- Different patterns of AHL sources are set up to observe the resultant GFP expression patterns.
- The concentration of AHL and GFP at various points in the system is calculated and plotted, offering insights into how molecular communication manifests spatially.

#### Implementation
- The experiment is implemented using Python, utilizing libraries like NumPy, SciPy, and Matplotlib for mathematical computations and visualizations.
- The code structure involves defining classes for each model and a simulation setup that integrates these models to study the system's behavior under various scenarios.

#### Results and Analysis
- The results are compared with in-vitro experiments to validate the in-silico model.
- Patterns of GFP expression are analyzed to understand how molecular signals propagate and influence cellular behavior in a synthetic multicellular system.

<img src="https://github.com/KeishiNishio/In_Silico_experiment/blob/main/patternA.png" width="200"><img src="https://github.com/KeishiNishio/In_Silico_experiment/blob/main/patternB.png" width="200"><img src="https://github.com/KeishiNishio/In_Silico_experiment/blob/main/patternC.png" width="200"><img src="https://github.com/KeishiNishio/In_Silico_experiment/blob/main/patternD.png" width="200">

#### Conclusion
This in-silico experiment provides valuable insights into molecular communication, offering a computational approach to study complex biological systems. It demonstrates the potential of computational biology in augmenting traditional experimental methods, allowing for more nuanced exploration of molecular dynamics in biological systems.

---

#### References
- S. Basu et al., “A synthetic multicellular system for programmed pattern formation,” Nature (2005).
