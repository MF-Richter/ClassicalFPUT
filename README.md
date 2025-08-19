# ClassicalFPUT

~Brief overview of scientific purpose (1–2 sentences).~

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ClassicalFPUT



## Author
- Moritz F. Richter



## Project Structure

- `src/` — core code modules
- `scripts/` — main analysis scripts
- `data/` — raw and processed data (see notes below)
- `results/` — generated outputs


## Installation & Setup

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "ClassicalFPUT"
```
which auto-activate the project and enable local path handling from DrWatson.


## How to Reproduce Results

- List main scripts to run and what they produce (e.g., “Run `scripts/plot_results.jl` for main figures”)


## Data Availability

- “Raw data not included. Please request from author.”


## Citation

- “If you use this code, please cite: [citation]”


## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
