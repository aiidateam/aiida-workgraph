# AiiDA-WorkTree
[![PyPI version](https://badge.fury.io/py/aiida-worktree.svg)](https://badge.fury.io/py/aiida-worktree)
[![Unit test](https://github.com/superstar54/aiida-worktree/actions/workflows/ci.yaml/badge.svg)](https://github.com/superstar54/aiida-worktree/actions/workflows/ci.yaml)
[![Docs status](https://readthedocs.org/projects/aiida-worktree/badge)](http://aiida-worktree.readthedocs.io/)

Provides the third workflow component: `WorkTree`, to design flexible node-based workflows using AiiDA.

In AiiDA, there are two workflow components: `workfunction` and `WorkChain`. Workfunction is easy to implement but it does not support automatic checkpointing, which is important for long-running calculations. Workchain supports automatic checkpointing but it is difficult to implement and also not as flexible as the `workfunction`. AiiDA-Nodetree provides the third component: `WorkTree`. It is easy to implement and supports automatic checkpointing. It is also flexible and can be used to design complex workflows.


Here is a detailed comparison between the ``WorkTree`` with two AiiDA built-in workflow components.


| Aspect                   | WorkFunction           | WorkChain              | WorkTree               |
| ------------------------ | ---------------------- | ---------------------- | ---------------------- |
| Use Case                 | Short-running jobs     | Long-running jobs      | Long-running jobs      |
| Checkpointing            | ``No``                 | Yes                    | Yes                    |
| Non-blocking             | ``No``                 | Yes                    | Yes                    |
| Implementation           | Easy                   | ``Difficult``          | Easy                   |
| Dynamic                  | ``No``                 | ``No``                 | Yes                    |
| Ready to Use             | Yes                    | ``No``,Need PYTHONPATH | Yes                    |
| Subprocesses Handling    | ``No``                 | Launches & waits       | Launches & waits       |
| Flow Control             | All                    | `if`, `while`          | `if`, `while`          |
| Termination              | ``Hard exit``          | ExitCode               | ExitCode               |
| Capabilities             | Calls calcs and works  | Calls any process      | Calls any process      |
| Data Passing             | Direct passing         | Context dictionary     | Link                   |
| Output Recording         | Limited support        | out & validates        | out                    |
| Port Exposing            | Limited support        | Supports automatic     | Limited support        |



## Installation

```console
    pip install --upgrade aiida-worktree
```


## Documentation
Check the [docs](https://aiida-worktree.readthedocs.io/en/latest/) and learn about the features.

## Examples

Create nodes from `calcfunction`.

```python
from aiida_worktree import node
from aiida.engine import calcfunction

# define add calcfunction node
@node()
@calcfunction
def add(x, y):
    return x + y

# define multiply calcfunction node
@node()
@calcfunction
def multiply(x, y):
    return x*y

```

Create a worktree to link the nodes.

```python
from aiida_worktree import WorkTree
from aiida import load_profile
from aiida.orm import Int
load_profile()

x = Int(2.0)
y = Int(3.0)
z = Int(4.0)

nt = WorkTree("test_add_multiply")
nt.nodes.new(add, name="add1", x=x, y=y)
nt.nodes.new(multiply, name="multiply1", y=z)
nt.links.new(nt.nodes["add1"].outputs[0], nt.nodes["multiply1"].inputs["x"])
nt.submit(wait=True)
```

The node graph from the worktree process:

<img src="docs/source/_static/images/add_multiply.png"/>


## TODO
- For the moment, I did not create a `WorkTreeNode` for the `WorkTree` process. I used the `WorkChainNode`, because AiiDA hard codes the `WorkChainNode` for the command (report), graph etc.


## Bugs
- the `report` does not work.

## License
[MIT](http://opensource.org/licenses/MIT)
