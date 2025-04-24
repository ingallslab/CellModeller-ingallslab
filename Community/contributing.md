# How to Contribute to the CellModeller Project

Thank you for your interest in contributing to CellModeller! Whether you’ve
found a bug, have an idea for a new feature, or want to help improve the
codebase, we welcome your contributions.

---

## Reporting Issues

One of the easiest ways to contribute is by using the software and reporting
any problems or ideas for improvement.

 - Found a bug? Have a feature request? Start by checking
 [existing issues](https://github.com/cellmodeller/CellModeller/issues).
 - Use the provided templates to file a new
 [bug report or feature request](https://github.com/cellmodeller/CellModeller/issues/new/choose).

Clear, detailed reports help us understand and address issues quickly.

---

## Making Code Changes

### 1. **Fork the repository**

 - Click the “Fork” button at the top right of the repo.
 - Clone your fork locally.

```bash
git clone https://github.com/your-username/CellModeller.git
cd CellModeller
```

### 2. **Create a new branch**

 - Name your branch descriptively (e.g., `fix/gui-launch-error` or
 `feature/move-to-pyproject`).

```bash
git checkout -b your-branch-name
```

### 3. **Make your changes**

When editing or adding code, please follow these project guidelines:

 - **Code formatting** is enforced via [**pre-commit**](https://pre-commit.com/).
 We use a `.pre-commit-config.yaml` to manage hooks and apply them only to
 relevant files/directories.

Before committing, run:
```bash
pre-commit install  # Run once to install the hook
pre-commit run --all-files
```

 - Write **meaningful commit messages** using the
 [**Conventional Commits**](https://www.conventionalcommits.org/) spec.[^2]

Example:
```
feat(gui): add launch button for CellModellerGUI
fix(simulation): correct time step rounding error
```

[^2]: ### Commit Message Guidelines
    The commit message should be structured as follows:

    ```
    <type>[optional scope]: <description>

    [optional body]

    [optional footer(s)]
    ```

    Type should be one of the following:
    type    | Use case
    --------|---------
    project | Change workflows, build systems, or update dependencies
    docs    | Only affects documentation
    feat    | New features/functionality
    fix     | Bug fixes
    refactor| Not a feature or a fix, retains logic
    style   | Does not change meaning of code
    test    | Only affects tests
    If a commit falls under multiple categories, it can most likely be split
    into several smaller commits. If not, pick the most appropriate.

    Valid scopes:
     - scripts
     - examples
     - simulation (files in CellModeller and not in a more specific scope)
     - gui
     - biophysics
     - utils

    Additional scopes may be added here as new base modules are implemented.

    Description should start with a lowercase letter, use imperative, present
    tense (i.e. "change" not "changed" nor "changes"), and not end with a dot (.)

    The body should also use imperative, present tense and include motivation for
    changes as well as contrast with previous behaviour if applicable.

    The footer contains references to issues with "Closes" and any information
    about a "BREAKING CHANGE:", following which the rest of the commit message
    is used for details.

 - Document **public API functions and classes** with **Doxygen-style docstrings**,
 which are used to generate reference documentation.

Example:
```python
def simulate_growth(cells, time_step):
  """
  @brief Simulate one time step of cell growth.
  @param cells List of Cell objects to update.
  @param time_step Duration of simulation step.
  @return Updated list of Cell objects.
  """
  pass
```

 - Keep `TODO.md` up to date with a note on what you're working on in your branch.
 - Test your changes locally and add unit tests if applicable.

### 4. **Submit a pull request**

 - Push your branch and open a PR to the `master` branch.
 - Add a bullet to the `CHANGELOG.md` under the `Unreleased` section summarizing
 your change.
 - Use a clear, concise PR title and description.

If your PR is related to an existing issue, mention it in the description (e.g.,
“Closes #123”).

---

## Release & Versioning Practices

We're committed to maintaining a healthy and transparent release process.
Here's how we manage that:

 - We follow **[Semantic Versioning](https://semver.org/)**
(e.g., `v1.2.0` -> `v2.0.0` for breaking changes).
 - Our `CHANGELOG.md` follows the
[**Keep a Changelog**](https://keepachangelog.com/en/1.0.0/) format.
 - PRs should include changelog entries summarizing user-facing changes.
 - All new releases should be **tested thoroughly** before final PR submission.

---

## For Frequent Contributors

If you contribute regularly, we’ll invite you to join the organization. This
grants you direct access to create branches and collaborate more closely with
the team.

Even if you have push access, please continue to:
 - Use feature branches and open pull requests for any changes.
 - Avoid committing directly to `master`.

---

## More Resources

 - [How to Contribute to Open Source (GitHub Guide)](https://opensource.guide/how-to-contribute/)
 - [Doxygen Documentation Generator](https://www.doxygen.nl/manual/docblocks.html)
 - [Black Code Formatter](https://black.readthedocs.io/)
 - [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/)
 - [Semantic Versioning](https://semver.org/)
 - [Keep a Changelog](https://keepachangelog.com/)

---

Thanks again for helping improve CellModeller![^1]

[^1]: This document was generated by ChatGPT using the GPT-4-turbo model, April
2025 release.
