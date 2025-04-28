# CellModeller

[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-green)](./licence.txt)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Docs](https://img.shields.io/badge/docs-doxygen-lightgrey.svg)](#-documentation)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](Community/code_of_conduct.md)

**CellModeller** is a Python-based framework for modeling large-scale multicellular
systems such as biofilms, plant and animal tissues. It enables simulation of
cellular growth, division, biophysics, gene regulation, and intercellular
signalling through customizable models written in Python.

---

## 📚 Table of Contents

- [Features](#-features)
- [Installation](#-installation)
- [Running the GUI](#-running-the-gui)
- [Writing Your Own Models](#-writing-your-own-models)
- [Documentation](#-documentation)
- [Troubleshooting](#-troubleshooting)
- [Citation](#-citation)
- [Contributing](#-contributing)
- [To-Do](#-todo)

---

## 🚀 Features

- Modular, Python-based simulation framework
- Support for intracellular processes, gene networks, and signaling
- Graphical User Interface (GUI) for interactive simulation
- Example models and tutorials to help you get started
- Fully scriptable for command-line use and automation
- Doxygen-compatible source code documentation

---

## 📦 Installation

We recommend using **Conda** to manage your environment:

```bash
conda create -n cellmodeller python=3.10
conda activate cellmodeller
```

Clone and install the project:

```bash
git clone https://github.com/ingallslab/CellModeller-ingallslab.git
cd CellModeller
pip install -e .
```

---

## 💻 Running the GUI

Launch the GUI:

```bash
python CellModeller/Scripts/CellModellerGUI.py
```

### GUI Buttons
- **Load Model** – Choose a model script to simulate
- **Load Pickle** – Load saved simulation states
- **Reset Simulation** – Restart the simulation
- **Run** – Start/resume the simulation
- **Save Pickle** – Toggle saving simulation data

### Keyboard Controls
- Move: `WASD`, `QE`, arrow keys
- Rotate: `IJKL`, `UO`
- Zoom: Scroll wheel
- Mouse drag: look around
- Right-click: move camera

### Running Models from CLI

You can also run models directly:

```bash
python Scripts/CellModellerGUI.py Examples/KyleP/Tutorial_1a.py
```

---

## 🧬 Writing Your Own Models

A CellModeller *model* is a Python script that defines:
- Simulation parameters
- Initial cell states and types
- Growth and interaction rules
- Modules (e.g. biophysics, signaling)

Explore the [`Examples/`](./Examples/KyleP) folder for real model templates.

> ⚠️ Our tutorials are a work in progress and may be slightly outdated — updates are on the way!

---

## 🧰 Documentation

You can generate local documentation using [Doxygen](https://www.doxygen.nl/):

```bash
doxygen ./doxygenconfig
```

Then open the docs at:

```bash
html/index.html
```

---

## 🔧 Troubleshooting

<details>
<summary>
No UserModule found in ...
</summary>
Note that a lot of the old examples use deprecated parts of the application, so
will not work until rewritten using the new API.

If you are using the new API:
1. Check if there is a subclass of UserModule in the given path.
2. Check if it is a concrete implementation or an abstract class by
instantiating it either at the top level of the file or in an interactive shell.
> If not, remember to implement the abstract methods sim_step() and
on_sim_ready(), even if they are empty functions!
</details>

<details>
<summary>
pyopencl has no attribute SVM
</summary>
This error has appeared on some Apple silicon Mac computers. Currently no
workarounds exist, rollback to v4.3.1. We are working on a fix.
</details>

Think your problem runs deeper than that?
Open an [issue](https://github.com/ingallslab/CellModeller-ingallslab/issues) on GitHub,
mentioning what you've tried and what you expect to happen.

---

## 📚 Citation

If you use CellModeller in academic work, please cite:

- **Rudge et al. (2012)** – *Computational modeling of synthetic microbial biofilms*, ACS Synthetic Biology 1 (8), 345–352
- **Rudge et al. (2013)** – *Cell polarity-driven instability generates self-organized, fractal patterning of cell layers*, ACS Synthetic Biology 2 (12), 705–714

---

## 🤝 Contributing

We welcome contributions of all kinds!
Check out our [CONTRIBUTING.md](./Community/contributing.md) for guidelines.

---

## ✅ TODO

See [TODO.md](./TODO.md) for plans and upcoming improvements.
Also check out the [CHANGELOG.md](./CHANGELOG.md) for recent updates.

[^1]

[^1]: This document was made with GenAI: ChatGPT using the GPT-4-turbo model,
April 2025 release.
