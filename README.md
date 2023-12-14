# AiiDA-WorkTree
[![PyPI version](https://badge.fury.io/py/aiida-worktree.svg)](https://badge.fury.io/py/aiida-worktree)
[![Unit test](https://github.com/superstar54/aiida-worktree/actions/workflows/ci.yaml/badge.svg)](https://github.com/superstar54/aiida-worktree/actions/workflows/ci.yaml)
[![Docs status](https://readthedocs.org/projects/aiida-worktree/badge)](http://aiida-worktree.readthedocs.io/)

Provides the third workflow component: `WorkTree`, to design flexible node-based workflows using AiiDA.

In AiiDA, there are two workflow components: `workfunction` and `WorkChain`. Workfunction is easy to implement but it does not support automatic checkpointing, which is important for long-running calculations. Workchain supports automatic checkpointing but it is difficult to implement and also not as flexible as the `workfunction`. AiiDA-WorkTree provides the third component: `WorkTree`. It is easy to implement and supports automatic checkpointing. It is also flexible and can be used to design complex workflows.


Here is a detailed comparison between the ``WorkTree`` with two AiiDA built-in workflow components.


| Aspect                   | WorkFunction           | WorkChain                     | WorkTree               |
| ------------------------ | ---------------------- | ----------------------------- | ---------------------- |
| Use Case                 | Short-running jobs     | Long-running jobs             | Long-running jobs      |
| Checkpointing            | ``No``                 | Yes                           | Yes                    |
| Execution order          | ``Sequential``         | ``Hybrid Sequential-Parallel``| Directed Acyclic Graph |
| Non-blocking             | ``No``                 | Yes                           | Yes                    |
| Implementation           | Easy                   | ``Difficult``                 | Easy                   |
| Dynamic                  | ``No``                 | ``No``                        | Yes                    |
| Ready to Use             | Yes                    | ``Need PYTHONPATH``           | Yes                    |
| Subprocesses Handling    | ``No``                 | Launches & waits              | Launches & waits       |
| Flow Control             | All                    | `if`, `while`                 | `if`, `while`, `match` |
| Termination              | ``Hard exit``          | ExitCode                      | ExitCode               |
| Data Passing             | Direct passing         | Context                       | Link & Context         |
| Output Recording         | Limited support        | Out & validates               | Out                    |
| Port Exposing            | Limited support        | Manual & automatic            | Manual                 |



## Installation

```console
    pip install aiida-worktree
```


## Documentation
Check the [docs](https://aiida-worktree.readthedocs.io/en/latest/) and learn about the features.

## Examples

Create calcfunction nodes:

```python
from aiida_worktree import node

# define add calcfunction node
@node.calcfunction()
def add(x, y):
    return x + y

# define multiply calcfunction node
@node.calcfunction()
def multiply(x, y):
    return x*y

```

Create a worktree to link the nodes.

```python
from aiida_worktree import WorkTree
from aiida import load_profile
from aiida.orm import Int
load_profile()


wt = WorkTree("test_add_multiply")
wt.nodes.new(add, name="add1", x=Int(2.0), y=Int(3.0))
wt.nodes.new(multiply, name="multiply1", y=Int(4.0))
wt.links.new(wt.nodes["add1"].outputs[0], wt.nodes["multiply1"].inputs["x"])
wt.submit(wait=True)
```

Start the web app, open a terminal and run:
```console
worktree web start
```

Then visit the page http://127.0.0.1:8000/worktree, you should find a `first_workflow` Worktree, click the pk and view the WorkTree.

<img src="docs/source/_static/images/first-workflow.png" />


One can also generate the node graph from the process:
```console
verdi node generate pk
```

<img src="docs/source/_static/images/add_multiply.png"/>


## Development

### Pre-commit and Tests
To contribute to this repository, please enable pre-commit so the code in commits are conform to the standards.
```console
pip install -e .[tests, pre-commit]
pre-commit install
```

### Web app
See the [README.md](https://github.com/superstar54/aiida-worktree/blob/main/aiida_worktree/web/README.md)

### Build and publish
Build package:
```console
pip install build
python -m build
```
Upload to PyPI:
```console
pip install twine
twine upload dist/*
```

## License
[MIT](http://opensource.org/licenses/MIT)
