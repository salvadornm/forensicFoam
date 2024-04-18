# Documentation

AMReX-based Navier-Stokes Reactive compressible LES code 

--8<-- "README.md"


## Installation

Download the latest version on the master branch using GitHub CLI

```bash
$ gh repo clone salvadornm/forensicFoam
```


or Download the latest release from [GitHub](https://github.com/salvadornm/forensicFoam/releases)



## Install OpenFOAM


### Visualization
[Paraview](https://www.paraview.org), 

[VisIt](https://visit-dav.github.io/visit-website/), 


## Documentation Editing

For help editing the documentation visit [mkdocs.org](https://www.mkdocs.org). To generate the docs locally `mkdocs serve`
and point the browser to [127.0.0.1.8000](http://127.0.0.1:8000)

You will need to install the `python-markdown-math` extension for rendering equations and the `markdown-callouts` extension for correctly displaying the warning and note blocks. All requirements can be installed automatically using

```bash
$ pip install -r docs/requirements.txt
```

You may need to install

```bash
$ pip install pip-tools
```

if you add new markdown extensions, edit the requirements.in

```bash
$ pip-compile requirements.in
```

